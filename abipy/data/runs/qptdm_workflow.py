#!/usr/bin/env python
from __future__ import division, print_function

import os
import numpy as np
import yaml
import collections
from collections import OrderedDict
import cPickle as pickle
import abipy.abilab as abilab
Workflow = abilab.Workflow
from abipy.tools import AttrDict

from pydispatch import dispatcher
from pymatgen.io.abinitio.task import Task, Dependency
from pymatgen.io.abinitio.workflow import BandStructureWorkflow, IterativeWorkflow

import logging
logger = logging.getLogger(__name__)


class OrderedDefaultdict(collections.OrderedDict):
    """
    Taken from:
    http://stackoverflow.com/questions/4126348/how-do-i-rewrite-this-function-to-implement-ordereddict/4127426#4127426
    """
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultdict, self).__init__(*args, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

    def __reduce__(self):  # optional, for pickle support
        args = self.default_factory if self.default_factory else tuple()
        return type(self), args, None, None, self.items()


class DefaultOrderedDict(OrderedDict):
    """
    Taken from
    http://stackoverflow.com/questions/6190331/can-i-do-an-ordered-default-dict-in-python
    """
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and not callable(default_factory)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))
    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                        OrderedDict.__repr__(self))



def hello(signal, sender):
    print("on_hello with sender %s, signal %s" % (sender, signal))

class Callback(object):

    def __init__(self, func, work, deps, cbk_data):
        self.func  = func
        self.work = work
        self.deps = deps
        self.cbk_data = cbk_data or {}
        self._disabled = False

    def __call__(self, flow):
        if self.can_execute():
            print("in callback")
            #print("in callback with sender %s, signal %s" % (sender, signal))
            cbk_data = self.cbk_data.copy()
            #cbk_data["_w_idx"] = self.w_idx
            return self.func(flow=flow, work=self.work, cbk_data=cbk_data)
        else:
            raise Exception("Cannot execute")

    def can_execute(self):
        return not self._disabled and [dep.status == Task.S_OK  for dep in self.deps]

    def disable(self):
        self._disabled = True

    def handle_sender(self, sender):
        return sender in [d.node for d in self.deps]


def bandstructure_flow(workdir, manager, scf_input, nscf_input):
    # Create the container that will manage the different workflows.
    flow = AbinitFlow(workdir, manager)
    flow.register_work(BandStructureWorkflow(scf_input, nscf_input))
    return flow.allocate()


def qptdm_workflow(wfk_file, scr_input, workdir, manager):
    """
    Factory function that builds a `QptdmWorkflow`.

    Args:
        scr_input:
            Input for the screening calculation.
        wfk_file:
            Path to the ABINIT WFK file to use for the computation of the screening.
        workdir:
            String defining the working directory. 
        manager:
            `TaskManager` object.

    Return
        `QptdmWorflow` object.
    """
    manager = work.manager

    # Build a temporary workflow with a shell manager just 
    # to run ABINIT to get the list of q-points for the screening.
    shell_manager = manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=0))

    fake_input = scr_input.deepcopy()
    tmpdir = os.path.join(workdir, "_qptdm_run")

    w = Workflow(workdir=tmpdir, manager=shell_manager)
    w.register(fake_input)
    w.build()

    # Create the symbolic link and add the magic value 
    # nqpdm = -1 to get the list of q-points
    fake_task = w[0]
    fake_task.inlink_file(wfk_file)
    fake_task.strategy.add_extra_abivars({"nqptdm": -1})

    w.start()

    # Parse the section with the q-points
    qpoints = yaml_kpoints(fake_task.log_file.path, tag="<QPTDM>")
    #print(qpoints)
    #w.remove()

    # Now we can build the final workflow.
    work = QptdmWorkflow(wfk_file, workdir, manager)

    for qpoint in qpoints:
        qptdm_input = scr_input.deepcopy()
        qptdm_input.set_variables(
            nqptdm=1,
            qptdm=qpoint
        )
        work.register(qptdm_input)

    return work


class QptdmWorkflow(Workflow):
    """
    This workflow parallelizes the calculation of the q-points of the screening. 
    It also provides the callback on_all_ok that calls mrgscr to merge 
    all the partial screening files produced.
    """
    #def __init__(self, wfk_file, workdir, manager):
    #    """
    #    Args:
    #        wfk_file:
    #            The path of the WFK file.
    #        workdir:
    #            String defining the working directory.
    #        manager:
    #            `TaskManager` object.
    #    """
    #    super(QptdmWorkflow, self).__init__(workdir, manager)

    #    self.wfk_file = os.path.abspath(wfk_file)

    def build(self):
        """Build files and directories."""
        super(QptdmWorkflow, self).build()

        # Create symbolic links to the WFK file.
        for task in self:
            task.inlink_file(self.wfk_file)

    def fill(self, wfk_file, scr_input):
        """
        Factory function that builds a `QptdmWorkflow`.

        Args:
            scr_input:
                Input for the screening calculation.
            wfk_file:
                Path to the ABINIT WFK file to use for the computation of the screening.

        Return
            `QptdmWorflow` object.
        """
        wfk_file = self.wfk_file = os.path.abspath(wfk_file)
        # Build a temporary workflow with a shell manager just 
        # to run ABINIT to get the list of q-points for the screening.
        shell_manager = self.manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=0))

        fake_input = scr_input.deepcopy()
        tmpdir = os.path.join(self.workdir, "_qptdm_run")

        w = Workflow(workdir=tmpdir, manager=shell_manager)
        w.register(fake_input)
        w.build()

        # Create the symbolic link and add the magic value 
        # nqpdm = -1 to get the list of q-points
        fake_task = w[0]
        fake_task.inlink_file(wfk_file)
        fake_task.strategy.add_extra_abivars({"nqptdm": -1})

        w.start()

        # Parse the section with the q-points
        qpoints = yaml_kpoints(fake_task.log_file.path, tag="<QPTDM>")
        #print(qpoints)
        #w.remove()

        # Now we can register the task for the different q-points 
        for qpoint in qpoints:
            qptdm_input = scr_input.deepcopy()
            qptdm_input.set_variables(
                nqptdm=1,
                qptdm=qpoint
            )
            self.register(qptdm_input)

        return self

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        In run `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Workflow`.
        """
        logger.info("about to call mrgscr in on_all_ok")

        scr_files = filter(None, [task.outdir.has_abiext("SCR") for task in self])
        logger.debug("scr_files:\n%s" % str(scr_files))
        assert len(scr_files) == len(self)

        mrgscr = abilab.Mrgscr(verbose=1)

        final_scr = mrgscr.merge_qpoints(scr_files, out_prefix="out", cwd=self.outdir.path)

        results = dict(
            returncode=0,
            message="mrgscr done",
            final_scr=final_scr,
        )

        return results


def cbk_qptdm_workflow(flow, work, cbk_data):
    #w_idx = cbk_data["_w_idx"]
    #workdir = flow[w_idx].workdir
    scr_input = cbk_data["input"]
    #wfk_file = cbk_data["wfk_file"]
    wfk_file = "/Users/gmatteo/Coding/abipy/abipy/data/runs/WORKS/work_0/task_1/outdata/out_WFK"

    work.fill(wfk_file, scr_input)
    work.allocate()
    #work.add_deps(cbk.deps)
    work.connect_signals()
    work.build()

    return work

    #manager = flow.manager
    ## TODO
    ## Index of the workflow and of the task that 
    ## has produced the WFK needed to generate the `QptdmWorkflow`.
    #wi, ti = (0,1)
    #task = flow[wi][ti]
    #wfk_file = task.outdir.has_abiext("WFK")
    #
    #try:
    #    return qptdm_workflow(wfk_file, scr_input, workdir, manager)

    #except Exception as exc:
    #    raise 


def yaml_kpoints(filename, tag="<KPOINTS>"):
    end_tag = tag.replace("<", "</")

    with open(filename, "r") as fh:
        lines = fh.readlines()

    start, end = None, None
    for i, line in enumerate(lines):
        if tag in line:
            start = i
        elif end_tag in line:
            end = i
            break
                                                                                                             
    if start is None or end is None:
        raise ValueError("%s\n does not contain any valid %s section" % (filename, tag))
                                                                                                             
    if start == end:
        # Empy section ==> User didn't enable Yaml in ABINIT.
        raise ValueError("%s\n contains an empty RUN_HINTS section. Enable Yaml support in ABINIT" % filename)

    s = "".join(l for l in lines[start+1:end])

    try:
        d = yaml.load(s)
    except Exception as exc:
        raise ValueError("Malformatted Yaml section in file %s:\n %s" % (filename, str(exc)))

    return np.array(d["reduced_coordinates_of_qpoints"])
    #return KpointList(reciprocal_lattice, frac_coords, weights=None, names=None)


def yaml_irred_perts(filename, tag="<KPOINTS>"):
    end_tag = tag.replace("<", "</")

    with open(filename, "r") as fh:
        lines = fh.readlines()

    start, end = None, None
    for i, line in enumerate(lines):
        if tag in line:
            start = i
        elif end_tag in line:
            end = i
            break

    if start is None or end is None:
        raise ValueError("%s\n does not contain any valid %s section" % (filename, tag))
                                                                                                             
    if start == end:
        # Empy section ==> User didn't enable Yaml in ABINIT.
        raise ValueError("%s\n contains an empty RUN_HINTS section. Enable Yaml support in ABINIT" % filename)

    s = "".join(l for l in lines[start+1:end])

    try:
        d = yaml.load(s)

    except Exception as exc:
        raise ValueError("Malformatted Yaml section in file %s:\n %s" % (filename, str(exc)))

    return np.array(d["reduced_coordinates_of_qpoints"])
    #return KpointList(reciprocal_lattice, frac_coords, weights=None, names=None)


class PhononWorkflow(Workflow):
    """
    This workflow is used for parallelizing the calculation of the 
    q-points of the screening. It also provides a on_all_ok method 
    that calls mrgddb to merge the partial DDB files.
    """
    def __init__(self, wfk_file, workdir, manager):
        super(PhononWorkflow, self).__init__(workdir, manager)

        self.wfk_file = os.path.abspath(wfk_file)

    def build(self):
        super(PhononWorkflow, self).build()

        # Create symbolic links to the WFK file.
        for task in self:
            task.inlink_file(self.wfk_file)

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        In run `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Workflow`.
        """
        logger.info("about to call anaddb.")

        ddb_files = filter(None, [task.outdir.has_abiext("DDB") for task in self])
        assert len(ddb_files) == len(self)
        logger.debug("ddb_files:\n%s" % str(ddb_files))

        mrgddb = abilab.Mrgddb(verbose=1)

        mrgddb.merge_qpoints(ddb_files, out_prefix="out", cwd=self.outdir.path)

        results = dict(
            returncode=0,
            message="merged done",
        )

        return results

#import collections
#AbinitPertubation = collections.namedtuple("IrredPert", "qpoint idir ipert")

def phonon_work(workdir, manager, ph_input, wfk_file): #, with_loto=):
    """
    Factory function that builds a `list of PhononWorkflow` objects.

    Args:
        workdir:
            Working directory. 
        manager:
            `TaskManager` object.
        ph_input:
            Input for the phonon calculation.
        wfk_file:
            Path to the ABINIT wavefunction file to use in the DFPT runs.

    Return:
        List of `PhononWorflow` objects. Each workflow computes all the irreducible
        atomic perturbations for s given q.
    """
    # Build a temporary workflow with a shell manager just to run 
    # ABINIT to get the list of q-points for the phonons.
    shell_manager = manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=0))

    #fake_input = ph_input.deepcopy()
    #w = Workflow(workdir=os.path.join(workdir, "phonon_run"), manager=shell_manager)
    #w.register(fake_input)
    #w.build()

    ## Create the symbolic link and add the magic value nqpdm = -1 to 
    ## get the list of q-points
    #fake_task = w[0]
    #fake_task.inlink_file(wfk_file)
    #fake_task.strategy.add_extra_abivars({"nqptdm": -1})

    #w.start()

    ## Parse the section with the q-points in the irreducible zone.
    #qpoints = yaml_kpoints(fake_task.log_file.path, tag="<QPTDM>")
    #print(qpoints)

    # Now we can build the final list of workflows:
    # One workflow per q-point, each workflow computes all 
    # the irreducible perturbations for a singe q-point.

    qpoints = np.reshape([
        0.0, 0.0, 0.0,
        ], (-1,3))

    works = []

    for iq, qpoint in enumerate(qpoints):

        # Get the irreducible perturbations for this q-point.
        fake_input = ph_input.deepcopy()
                                                                                                
        w = Workflow(workdir=os.path.join(workdir, "irred_pert_run"), manager=shell_manager)
        w.register(fake_input)
        w.build()
                                                                                                
        # Create the symbolic link and add the magic value 
        # nqpdm = -1 to get the list of q-points.
        fake_task = w[0]
        fake_task.inlink_file(wfk_file)
        vars = dict(
            qpt=qpoint
            #??
            )
        fake_task.strategy.add_extra_abivars(vars)
                                                                                                
        w.start()

        # Parse the file to get the perturbations.
        irred_perts = yaml_irred_perts(fake_task.log_file.path, tag="<QPTDM>")

        # Compute all the irreducible perturbations.
        dirpath = os.path.join(workdir, "QPT_" + str(i))
        work_q = PhononWorkflow(dirpath, manager, wfk_file)

        for irred_pert in irred_perts:
            new_input = ph_input.deepcopy()
            new_input.set_variables(
                #rfpert=1,
                qpt=irred_pert.qpoint,
                idir=irred_pert.idir,
                ipert=irred_pert.ipert,
            )
            work_q.register(new_input)

        works.append(work_q)

    return works


def gw_workflow(workdir, manager, scf_input, nscf_input, scr_input, sigma_input):
    # NEW version

    # Create the container that will manage the different workflows.
    flow = AbinitFlow(workdir, manager)

    # Register the first workflow (GS + NSCF calculation)
    bands_work = flow.register_work(BandStructureWorkflow(scf_input, nscf_input))

    assert not bands_work.scf_task.depends_on(bands_work.scf_task)
    assert bands_work.nscf_task.depends_on(bands_work.scf_task)

    # Register the callback that will be executed to build another workflow
    scr_work = flow.register_cbk(cbk=cbk_qptdm_workflow, cbk_data={"input": scr_input},
                                 deps={bands_work.nscf_task: "WFK"}, work_class=QptdmWorkflow
                                )
                             

    #assert scr_work.depends_on(bands_work.nscf_task)
    #assert not scr_work.depends_on(bands_work.scf_task)

    #sigma_work = Workflow(sigma_input, deps={bands_work.nscf_task: "WFK", scr_work: "SCR"})
    #flow.register(sigma_work)

    # The last workflow is a SIGMA run that will use 
    # the data produced in the previous two workflows.
    sigma_work = Workflow()
    #sigma_task = sigma_work.register(sigma_input, deps={bands_work.nscf_task: "WFK"})
    sigma_task = sigma_work.register(sigma_input, deps={bands_work.nscf_task: "WFK", scr_work: "SCR"})
    flow.register_work(sigma_work)

    flow.allocate()
    flow.show_dependencies()

    assert sigma_task.depends_on(bands_work.nscf_task)
    assert not sigma_task.depends_on(bands_work.scf_task)
    print("sigma_work.deps", sigma_work.deps)
    print("sigma_task.deps", sigma_task.deps)
    assert sigma_task.depends_on(scr_work)
    #assert sigma_work.depends_on(scr_work)

    return flow


class GW0_Workflow(IterativeWorkflow):
    def __init__(self, wfk_file, sigma_input, workdir, manager):
        self.wfk_file = os.path.abspath(wfk_file)
        max_niter = 25

        inputs = []
        for i in range(max_niter):
            new = sigma_input.copy()
            inputs.append(new)

        IterativeWorkflow.__init__(sigma_generator, max_niter=max_niter, workdir=workdir, manager=manager)

    def build(self):
        """Build files and directories."""
        super(GW0_Workflow, self).build()

        # Create symbolic links to the WFK file.
        for task in self:
            task.inlink_file(self.wfk_file)

    def exit_iteration(self, *args, **kwargs):
        """
        Return a dictionary with the results produced at the given iteration.
        The dictionary must contains an entry "converged" that evaluates to
        True if the iteration should be stopped.
        """

    #def on_all_ok(self):
        # Create the container that will manage the different workflows.
        #work = Workflow()

        ## Register the first workflow (GS + NSCF calculation)
        #bands_work = flow.register(BandStructureWorkflow(scf_input, nscf_input))

        ## FIXME there's a problem with the ids if I want to add a task at runtime
        ## e.g if I want to iterate the SCHF task.
        ##sigma_work = Workflow(hf_input, deps={bands_work.nscf_task: "WFK"})
        ##flow.register(sigma_work)

        #flow.allocate()
        #flow.show_dependencies()
        #return flow


class AbinitFlow(collections.Iterable):
    """
    This object is a container of workflows. Its main task is managing the 
    possible inter-depencies among the workflows and the generation of
    dynamic worflows that are generates by callabacks registered by the user.
    """
    #PICKLE_FNAME = "__AbinitFlow__.pickle"
    PICKLE_FNAME = "__workflow__.pickle"

    def __init__(self, workdir, manager):
        self.workdir = os.path.abspath(workdir)

        self.manager = manager.deepcopy()

        # List of workflows handled by self.
        self._works = []

        # List of callbacks that must be executed when the dependencies reach S_OK
        self._callbacks = []

        # Dictionary [work] --> list of deps with the dependencies of each workflow.
        #self._deps_dict = collections.defaultdict(list)

    def __repr__(self):
        return "<%s at %s, workdir=%s>" % (self.__class__.__name__, id(self), self.workdir)

    def __str__(self):
        return "<%s, workdir=%s>" % (self.__class__.__name__, self.workdir)

    def __len__(self):
        return len(self._works)

    def __iter__(self):
        return self._works.__iter__()

    def __getitem__(self, slice):
        return self._works[slice]

    def used_ids(self):
        """
        Returns a set with all the ids used so far to identify `Task` and `Workflow`.
        """
        ids = []
        for work in self:
            ids.append(work.node_id)
            for task in work:
                ids.append(task.node_id)

        used_ids = set(ids)
        assert len(ids_set) == len(ids)

        return used_ids

    def generate_new_nodeid(self):
        """Returns an unused node identifier."""
        used_ids = self.used_ids()

        import itertools
        for nid in intertools.count():
            if nid not in used_ids:
                return nid

    def check_status(self):
        """Check the status of the workflows in self."""
        for work in self:
            work.check_status()

    def build(self, *args, **kwargs):
        for work in self:
            work.build(*args, **kwargs)

    def build_and_pickle_dump(self, protocol=-1):
        self.build()
        self.pickle_dump(protocol=protocol)

    def pickle_dump(self, protocol=-1):
        """Save the status of the object in pickle format."""
        filepath = os.path.join(self.workdir, self.PICKLE_FNAME)

        with open(filepath, mode="w" if protocol == 0 else "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)

        # Atomic transaction.
        #import shutil
        #filepath_new = filepath + ".new"
        #filepath_save = filepath + ".save"
        #shutil.copyfile(filepath, filepath_save)

        #try:
        #    with open(filepath_new, mode="w" if protocol == 0 else "wb") as fh:
        #        pickle.dump(self, fh, protocol=protocol)

        #    os.rename(filepath_new, filepath)
        #except:
        #    os.rename(filepath_save, filepath)
        #finally:
        #    try
        #        os.remove(filepath_save)
        #    except:
        #        pass
        #    try
        #        os.remove(filepath_new)
        #    except:
        #        pass

    @staticmethod
    def pickle_load(filepath):
        with open(filepath, "rb") as fh:
            flow = pickle.load(fh)

        flow.connect_signals()
        return flow

    def register_work(self, work, deps=None, manager=None):
        """
        Register a new `Workflow` and add it to the internal list, 
        taking into account possible dependencies.

        Args:
            work:
                `Workflow` object.
            deps:
                List of `Dependency` objects specifying the dependency of this node.
                An empy list of deps implies that this node has no dependencies.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the workflow.

        Returns:   
            The workflow.
        """
        # Directory of the workflow.
        work_workdir = os.path.join(self.workdir, "work_" + str(len(self)))

        # Make a deepcopy since manager is mutable and we might change it at run-time.
        manager = self.manager.deepcopy() if manager is None else manager.deepcopy()

        work.set_workdir(work_workdir)
        work.set_manager(manager)

        self._works.append(work)

        if deps:
            deps = [Dependency(node, exts) for node, exts in deps.items()]
            work.add_deps(deps)
            #self._deps_dict[work].extend(deps)

        return work

    def register_cbk(self, cbk, cbk_data, deps, work_class, manager=None):
        """
        Registers a new workflow and add it to the internal list, 
        taking into account possible dependencies.
                                                                                                            
        Args:
            cbk:
                Callback function.
            cbk_data
            deps:
                List of `Dependency` objects specifying the dependency of this node.
            work_class:
                `Workflow` class to instantiate.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the `Flow`.
                                                                                                            
        Returns:   
            The `Workflow` that will be finalized by the callback.
        """
        # Directory of the workflow.
        work_workdir = os.path.join(self.workdir, "work_" + str(len(self)))
                                                                                                            
        # Make a deepcopy since manager is mutable and we might change it at run-time.
        manager = self.manager.deepcopy() if manager is None else manager.deepcopy()
                                                                                                            
        # create an empty workflow and register the callback
        work = work_class(workdir=work_workdir, manager=manager)
        
        self._works.append(work)
                                                                                                            
        deps = [Dependency(node, exts) for node, exts in deps.items()]
        if not deps:
            raise ValueError("A callback must have deps!")

        work.add_deps(deps)

        # Wrap the callable in a Callback object and save 
        # useful info such as the index of the workflow and the callback data.
        cbk = Callback(cbk, work, deps=deps, cbk_data=cbk_data)
                                                                                                            
        self._callbacks.append(cbk)
                                                                                                            
        return work

    def allocate(self):
        for work in self:
            work.allocate()

        return self

    def show_dependencies(self):
        #for node, deps in self._deps_dict.items():
        #    print("Node %s had dependencies:" % str(node))
        #    for i, dep in enumerate(deps):
        #        print("%d) %s, status=%s" % (i, str(dep.node), str(dep.status)))

        for work in self:
            work.show_intrawork_deps()

    def on_dep_ok(self, signal, sender):
        print("on_dep_ok with sender %s, signal %s" % (str(sender), signal))

        for i, cbk in enumerate(self._callbacks):

            if not cbk.handle_sender(sender):
                print("Do not handle")
                continue

            if not cbk.can_execute():
                print("cannot execute")
                continue 

            # Execute the callback to generate the workflow.
            print("about to build new workflow")
            #empty_work = self._works[cbk.w_idx]

            # TODO better treatment of ids
            # Make sure the new workflow has the same id as the previous one.
            #new_work_idx = cbk.w_idx
            work = cbk(flow=self)
            #new_work.set_node_id(empty_work.node_id)
            ##new_work.set_task_nodeids(self, id_generator)
            #new_work.allocate()
            #new_work.add_deps(cbk.deps)
            #new_work.connect_signals()
            #new_work.build()
            work.add_deps(cbk.deps)

            #print("new: ", str(w))
            # Replace the empty workflow with the new one and 
            # register the callback for removal
            #self._works[cbk.w_idx] = new_work

            #dispatcher.send(signal=Task.S_OK, sender=work)
            #except Exception as exc:

            # Disable the callback.
            cbk.disable()

        # Update the database.
        self.pickle_dump()

    def connect_signals(self):
        """
        Connect the signals within the workflow.
        self is responsible for catching the important signals raised from 
        its task and raise new signals when some particular condition occurs.
        """
        # Connect the signals inside each Workflow.
        for work in self:
            work.connect_signals()

        # Observe the nodes that must reach S_OK in order to call the callbacks.
        for cbk in self._callbacks:
            for dep in cbk.deps:
                print("connecting %s \nwith sender %s, signal %s" % (str(cbk), dep.node, dep.node.S_OK))
                dispatcher.connect(self.on_dep_ok, signal=dep.node.S_OK, sender=dep.node, weak=False)

        self.show_receivers()

    def show_receivers(self, sender=dispatcher.Any, signal=dispatcher.Any):
        print("*** live receivers ***")
        for rec in dispatcher.liveReceivers(dispatcher.getReceivers(sender, signal)):
            print("receiver -->", rec)
        print("*** end live receivers ***")

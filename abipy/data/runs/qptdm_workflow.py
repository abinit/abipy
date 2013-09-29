#!/usr/bin/env python
from __future__ import division, print_function

import os
import numpy as np
import yaml
import collections
import cPickle as pickle
import abipy.abilab as abilab
from abipy.tools import AttrDict

from pydispatch import dispatcher
from pymatgen.io.abinitio.task import Task 
from pymatgen.io.abinitio.workflow import BandStructure, Dependency

import logging
logger = logging.getLogger(__name__)


def show_receivers(sender=dispatcher.Any, signal=dispatcher.Any):

    print("*** live receivers ***")
    for rec in dispatcher.liveReceivers(dispatcher.getReceivers(sender, signal)):
        print("receiver -->", rec)
    print("*** end live receivers ***")


def hello(signal, sender):
    print("on_hello with sender %s, signal %s" % (sender, signal))

class Callback(object):

    def __init__(self, func, w_idx, deps, cb_data):
        self.func  = func
        self.w_idx = w_idx
        self.deps = deps
        self.cb_data = cb_data or {}

    def __call__(self, works, **kwargs):
        print("in callback")
        #sender = kwargs["sender"]
        #signal = kwargs["sender"]
        #print("in callback with sender %s, signal %s" % (sender, signal))
        if self.can_execute():
            cb_data = self.cb_data.copy()
            cb_data["_w_idx"] = self.w_idx
            return self.func(works, cb_data=cb_data, **kwargs)
        else:
            raise Exception("Cannot execute")

    def can_execute(self):
        return [dep.status == Task.S_OK  for dep in self.deps]

    def handle_sender(self, sender):
        return sender in self.deps


class QptdmWorkflow(abilab.Workflow):
    """
    This workflow parallelizes the calculation of the q-points of the screening. 
    It also provides the callback on_all_ok that calls mrgscr to merge 
    all the partial screening files produced.
    """
    def __init__(self, workdir, manager, wfk_file):
        """
        Args:
            workdir:
                String defining the working directory.
            manager:
                `TaskManager` object.
            wfk_file:
                The path of the WFK file.
        """
        super(QptdmWorkflow, self).__init__(workdir, manager)

        self.wfk_file = os.path.abspath(wfk_file)

    def build(self):
        """Buil files and directories."""
        super(QptdmWorkflow, self).build()

        # Create symbolic links to the WFK file.
        for task in self:
            task.inlink_file(self.wfk_file)

    def on_all_ok(self):
        logger.info("about to call mrgscr in on_all_ok")

        scr_files = []
        for task in self:
            scr = task.outdir.has_abiext("SCR")
            assert scr 
            scr_files.append(scr)

        logger.debug("scr_files:\n%s" % str(scr_files))
        mrgscr = abilab.Mrgscr(verbose=1)

        out_prefix = "out_hello"
        final_scr = mrgscr.merge_qpoints(scr_files, out_prefix=out_prefix, cwd=self.outdir.path)

        results = dict(
            returncode=0,
            message="mrgscr done",
            final_scr=final_scr,
        )

        return results


def cb_build_qptdm_workflow(works, cb_data):
    w_idx = cb_data["_w_idx"]
    workdir = works[w_idx].workdir
    scr_input = cb_data["scr_input"]

    manager = works.manager
    # TODO
    # Index of the workflow and of the task that 
    # has produced the WFK needed to generate the `QptdmWorkflow`.
    wi, ti = (0,1)
    task = works[wi][ti]
    wfk_file = task.outdir.has_abiext("WFK")
    
    try:
        work = build_qptdm_workflow(workdir, manager, scr_input, wfk_file)
        return work

    except Exception as exc:
        raise 


def build_qptdm_workflow(workdir, manager, scr_input, wfk_file):
    """
    Factory function that builds a `QptdmWorkflow`.

    Args:
        workdir:
            String defining the working directory. 
        manager:
            `TaskManager` object.
        scr_input:
            Input for the screening calculation.
        wfk_file:
            Path to the ABINIT WFK file to use for the computation of the screening.

    Return
        `QptdmWorflow` object.
    """
    # Build a temporary workflow with a shell manager just 
    # to run ABINIT to get the list of q-points for the screening.
    shell_manager = manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=0))

    fake_input = scr_input.deepcopy()

    w = abilab.Workflow(workdir=os.path.join(workdir, "_qptdm_run"), manager=shell_manager)
    w.register(fake_input)
    w.build()

    # Create the symbolic link and add the magic value 
    # nqpdm = -1 to get the list of q-points
    fake_task = w[0]
    fake_task.inlink_file(wfk_file)
    fake_task.strategy.add_extra_abivars({"nqptdm": -1})

    #print(scr_input)
    #print(fake_task.strategy.make_input())
    w.start()

    # Parse the section with the q-points
    qpoints = yaml_kpoints(fake_task.log_file.path, tag="<QPTDM>")
    #print(qpoints)
    #w.remove()

    # Now we can build the final workflow.
    work = QptdmWorkflow(workdir, manager, wfk_file)

    for qpoint in qpoints:
        qptdm_input = scr_input.deepcopy()
        qptdm_input.set_variables(
            nqptdm=1,
            qptdm=qpoint
            #qptdm=qpoint.frac_coords,
        )
        work.register(qptdm_input)

    return work


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


class PhononWorkflow(abilab.Workflow):
    """
    This workflow is used for parallelizing the calculation of the 
    q-points of the screening. It also provides a on_all_ok method 
    that calls mrgddb to merge the partial DDB files.
    """
    def __init__(self, workdir, manager, wfk_file):
        super(PhononWorkflow, self).__init__(workdir, manager)

        self.wfk_file = os.path.abspath(wfk_file)

    def build(self):
        super(PhononWorkflow, self).build()

        # Create symbolic links to the WFK file.
        for task in self:
            task.inlink_file(self.wfk_file)

    def on_all_ok(self):
        logger.info("about to call anaddb.")

        ddb_files = []
        for task in self:
            ddb = task.outdir.has_abiext("DDB")
            assert ddb
            ddb_files.append(scr)

        logger.debug("ddb_files:\n%s" % str(ddb_files))
        mrgddb = abilab.Mrgddb(verbose=1)

        out_prefix = "out_hello"
        mrgddb.merge_qpoints(ddb_files, out_prefix=out_prefix, cwd=self.outdir.path)

        results = dict(
            returncode=0,
            message="merged done",
        )

        return results

#import collections
#AbinitPertubation = collections.namedtuple("IrredPert", "qpoint idir ipert")

def build_phonon_workflow(workdir, manager, ph_input, wfk_file): #, with_loto=):
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
    #w = abilab.Workflow(workdir=os.path.join(workdir, "phonon_run"), manager=shell_manager)
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
                                                                                                
        w = abilab.Workflow(workdir=os.path.join(workdir, "irred_pert_run"), manager=shell_manager)
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


def build_gw_workflow(workdir, manager, scf_input, nscf_input, scr_input, sigma_input):

    # Create the container that will manage the different workflows.
    works = AbinitWorks(workdir, manager)

    # Register the first workflow (GS + NSCF calculation)
    wfk_work = abilab.Workflow()
    works.register(wfk_work)

    # Pass the information needed to build the GS + NSCF workflow:
    #   1) Input files 
    #   2) Dependency between the NSCF and the SCF task 
    scf_task = wfk_work.register(scf_input)
    nscf_task = wfk_work.register(nscf_input, deps={scf_task: "DEN"})
    #wfk_work = BandStructure(workdir + "/WFK", manager, scf_input, nscf_input)

    # Register a function that will be executed to build another workflow
    # the callback will have access the all the workflows that have been 
    # registered so far. The workflow will be generated at runtime and will
    # depend on the previous workflows specified in deps.
    scr_work = works.register(cb_build_qptdm_workflow, deps={nscf_task: "WFK"}, 
        cb_data={"scr_input": scr_input}) # wi, ti = (0,1)

    # The last workflow is a SIGMA run that will use 
    # the data produced in the previous two workflows.
    sigma_work = abilab.Workflow()
    sigma_task = works.register(sigma_work, deps={nscf_task: "WFK", scr_work: "SCR"})
    sigma_work.register(sigma_input)

    return works

class AbinitWorks(collections.Iterable):
    """
    This object is a container of workflows. Its main task is managing the 
    possible inter-depencies among the workflows and the generation of
    dynamic worflows that are generates by callabacks registered by the user.
    """
    #PICKLE_FNAME = "__AbinitWorks__.pickle"
    PICKLE_FNAME = "__workflow__.pickle"

    def __init__(self, workdir, manager):
        self.workdir = os.path.abspath(workdir)

        self.manager = manager.deepcopy()

        # List of workflows handled by self.
        self._works = []

        # List of callbacks that must be executed when the dependencies reach S_OK
        self._callbacks = []

        # Dictionary [work] --> list of deps 
        # with the dependencies of each workflow.
        self._deps_dict = collections.defaultdict(list)

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
            works = pickle.load(fh)

        works.connect_signals()
        return works

    def register(self, obj, deps=None, manager=None, cb_data=None):
        """
        Registers a new workflow and add it to the internal list, 
        taking into account possible dependencies.

        Args:
            obj:
                `Strategy` object or `AbinitInput` instance.
                if obj is a `Strategy`, we create a new `AbinitTask` from 
                the input strategy and add it to the list.
            deps:
                List of `Dependency` objects specifying the dependency of this node.
                An empy list of deps implies that this node has no dependencies.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the workflow.

        Returns:   
            `Dependency` object
        """
        # Directory of the workflow.
        work_workdir = os.path.join(self.workdir, "work_" + str(len(self)))

        # Make a deepcopy since manager is mutable and we might change it at run-time.
        manager = self.manager.deepcopy() if manager is None else manager.deepcopy()

        if not callable(obj):
            # We have a workflow object. 
            assert cb_data is None
            work = obj
            work.set_workdir(work_workdir)
            work.set_manager(manager)
            callback = None

        else:
            # We received a callback --> create an empty workflow and register the callback
            work = abilab.Workflow(work_workdir, manager)
            callback = obj
        
        self._works.append(work)

        if deps:
            deps = [Dependency(node, exts) for node, exts in deps.items()]
            self._deps_dict[work].extend(deps)

            print(80*"*")
            print("%s needs:\n%s" % (work, "\n".join(str(l) for l in deps)))
            print(80*"*")

        if callback is not None:
            if not deps:
                raise ValueError("A callback must have deps!")
            print("Callback dependencies:", deps)

            # Wraps the callable in a Callback object a save 
            # useful info such as the index of the workflow in self.
            cb = Callback(callback, w_idx=len(self)-1, deps=deps, cb_data=cb_data)

            self._callbacks.append(cb)

        return work

    def on_dep_ok(self, signal, sender):
        print("on_dep_ok with sender %s, signal %s" % (str(sender), signal))

        for i, cbk in enumerate(self._callbacks):
            if not cbk.handle_sender(sender):
                continue

            if not cbk.can_execute():
                continue 

            print("about to build new workflow")
            w_idx = cbk.w_idx
            w = cbk(works=self)
            w.connect_signals()
            w.build()

            # Replace the empty workflow with the new one 
            # and remove the callback.
            self._works[w_idx] = w
            self._callbacks[i] = None

            # Update the database.
            self.pickle_dump()

        while True:
            try:
                self._callbacks.remove(None)
            except ValueError:
                break

    def connect_signals(self):
        print("Connecting AbinitWorks")
        # Connect the signals inside each Workflow.
        for work in self:
            work.connect_signals()

        # Handle possible callbacks.
        for c in self._callbacks:
            for dep in c.deps:
                print("connecting %s \nwith sender %s, signal %s" % (str(c), dep.node, dep.node.S_OK))
                dispatcher.connect(self.on_dep_ok, signal=dep.node.S_OK, sender=dep.node, weak=False)
                #dispatcher.connect(hello, signal=dep._node.S_OK)
                                                                                                        
        show_receivers()

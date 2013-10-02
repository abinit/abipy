#!/usr/bin/env python
from __future__ import division, print_function

import os
import numpy as np
import yaml
import collections
from collections import OrderedDict
import cPickle as pickle
from abipy.abilab import Workflow, AbinitFlow, Mrgscr, Mrgddb
from abipy.tools import AttrDict

from pydispatch import dispatcher
from pymatgen.io.abinitio.tasks import Task, Dependency
from pymatgen.io.abinitio.workflows import BandStructureWorkflow, IterativeWorkflow

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

def bandstructure_flow(workdir, manager, scf_input, nscf_input):
    flow = AbinitFlow(workdir, manager)
    flow.register_work(BandStructureWorkflow(scf_input, nscf_input))
    return flow.allocate()


def g0w0_flow(workdir, manager, scf_input, nscf_input, scr_input, sigma_input):
    flow = AbinitFlow(workdir, manager)
    flow.register_work(G0W0_Workflow(scf_input, scf_input, nscf_input, scr_input, sigma_input))
    return flow.allocate()


class QptdmWorkflow(Workflow):
    """
    This workflow parallelizes the calculation of the q-points of the screening. 
    It also provides the callback on_all_ok that calls mrgscr to merge 
    all the partial screening files produced.
    """

    def build(self):
        """Build files and directories."""
        super(QptdmWorkflow, self).build()

        # Create symbolic links to the WFK file.
        for task in self:
            task.inlink_file(self.wfk_file)

    def fill(self, wfk_file, scr_input):
        """
        Create the tasks and register them in self.

        Args:
            wfk_file:
                Path to the ABINIT WFK file to use for the computation of the screening.
            scr_input:
                Input for the screening calculation.

        Return:
            self
        """
        assert len(self) == 0
        wfk_file = self.wfk_file = os.path.abspath(wfk_file)

        # Build a temporary workflow in the tmpdir that will use a shell manager 
        # to run ABINIT in order to get the list of q-points for the screening.
        shell_manager = self.manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=0))

        w = Workflow(workdir=self.tmpdir.path_join("_qptdm_run"), manager=shell_manager)

        fake_input = scr_input.deepcopy()
        fake_task = w.register(fake_input)
        w.build()

        # Create the symbolic link and add the magic value 
        # nqpdm = -1 to the input to get the list of q-points.
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
            self.register(qptdm_input, manager=self.manager)

        return self.allocate()

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        In run `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Workflow`.
        """
        scr_files = filter(None, [task.outdir.has_abiext("SCR") for task in self])

        msg = "on_all_ok: will call mrgscr to merge %s:\n" % str(scr_files)
        logger.debug(msg)
        assert len(scr_files) == len(self)

        mrgscr = Mrgscr(verbose=1)

        final_scr = mrgscr.merge_qpoints(scr_files, out_prefix="out", cwd=self.outdir.path)

        results = dict(
            returncode=0,
            message="mrgscr done",
            final_scr=final_scr,
        )

        return results


def phonon_flow(workdir, manager, scf_input, ph_input):
    # Create the container that will manage the different workflows.
    flow = AbinitFlow(workdir, manager)

    # Register the first workflow (GS calculation)
    #scf_work = flow.register_work(Workflow().register(scf_input, task_class=ScfTask)

    # Register the callback that will be executed to build the other workflows:
    # one workflow for q-point
    #scr_work = flow.register_cbk(cbk=cbk_phonons_workflow, cbk_data={"input": ph_input},
    #                             deps={scf_task: "WFK"}, work_class=QptdmWorkflow
    #                            )
    
    return flow.allocate()


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

    s = "".join(lines[start+1:end])

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

    s = "".join(lines[start+1:end])

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
        ddb_files = filter(None, [task.outdir.has_abiext("DDB") for task in self])

        msg = "on_all_ok: will call mrgddb to merge %s:\n" % str(ddb_files)
        logger.debug(msg)
        assert len(ddb_files) == len(self)

        mrgddb = Mrgddb(verbose=1)

        mrgddb.merge_qpoints(ddb_files, out_prefix="out", cwd=self.outdir.path)

        results = dict(
            returncode=0,
            message="merged done",
        )

        return results

    def fill(self, wfk_file, ph_input): #, with_loto=):
        """
        Factory function that builds a `list of PhononWorkflow` objects.

        Args:
            wfk_file:
                Path to the ABINIT wavefunction file to use in the DFPT runs.
            ph_input:
                Input for the phonon calculation.

        Return:
            List of `PhononWorflow` objects. Each workflow computes all the irreducible
            atomic perturbations for s given q.
        """
        assert len(self) == 0
        wfk_file = self.wfk_file = os.path.abspath(wfk_file)

        # Build a temporary workflow with a shell manager just to run 
        # ABINIT to get the list of q-points for the phonons.
        shell_manager = self.manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=0))

        #fake_input = ph_input.deepcopy()
        #w = Workflow(workdir=self.tmpdir.path_join("_ph_run"),, manager=shell_manager)
        #w.register(fake_input)
        #w.build()

        # Create the symbolic link and add the magic value 
        # nqpdm = -1 to get the list of q-points
        #fake_task = w[0]
        #fake_task.inlink_file(wfk_file)
        #fake_task.strategy.add_extra_abivars({"nqptdm": -1})

        #w.start()

        ## Parse the section with the q-points in the irreducible zone.
        #qpoints = yaml_kpoints(fake_task.log_file.path, tag="<QPTDM>")
        #print(qpoints)

        # Now we can build the final list of tasks:
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

        #work.allocate()

        return works


def cbk_qptdm_workflow(flow, work, cbk_data):
    scr_input = cbk_data["input"]
    # Use the WFK file produced by the second 
    # Task in the first Workflow (NSCF step).
    nscf_task = flow[0][1]
    wfk_file = nscf_task.outdir.has_abiext("WFK")

    work.fill(wfk_file, scr_input)
    #work.add_deps(cbk.deps)
    work.connect_signals()
    work.build()

    return work


def g0w0_flow_with_qptdm(workdir, manager, scf_input, nscf_input, scr_input, sigma_input):
    # Create the container that will manage the different workflows.
    print(manager)
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



#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.abilab as abilab
import abipy.data as data  

from abipy.data.runs import decorate_main

from abipy.data.runs.qptdm_workflow import *

def make_base_input():
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    inp.set_structure(data.structure_from_ucell("si"))

    # Global variables
    global_vars = dict(ecut=12,
                       nsppol=1,
                       nband=20,
                       timopt=-1,
                       paral_kgb=0,
                       chksymbreak=0,
                       prtwf=0
                       #accesswff=3,
                    )

    inp.set_variables(**global_vars)

    # Simple GS run.
    inp.set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])
    #inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    #inp.set_kmesh(ngkpt=[12,12,12], shiftk=[0.1,0.2,0.3])
    inp.tolvrs = 1e-8

    return inp


def all_inputs():
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    pseudos = data.pseudos("14si.pspnc")

    ecut = ecutwfn = 6

    global_vars = dict(
        ecut=ecut,
        timopt=-1,
        istwfk = "*1",
    )

    gs = abilab.AbiInput(pseudos=pseudos)
    gs.set_structure(structure)

    # This grid is the most economical, but does not contain the Gamma point.
    gs_kmesh = dict(
        ngkpt=[2,2,2],
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap. 
    gw_kmesh = dict(
        ngkpt=[2,2,2],
        shiftk=[0.0, 0.0, 0.0,  
                0.0, 0.5, 0.5,  
                0.5, 0.0, 0.5,  
                0.5, 0.5, 0.0]
    )

    # Dataset 1 (GS run)
    gs.set_variables(**global_vars)
    gs.set_kmesh(**gs_kmesh)
    gs.set_variables(tolvrs=1e-6,
                     nband=4,
                    )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf = abilab.AbiInput(pseudos=pseudos)
    nscf.set_structure(structure)
    nscf.set_variables(**global_vars)
    nscf.set_kmesh(**gw_kmesh)

    nscf.set_variables(iscf=-2,
                       getden=1,
                       tolwfr=1e-12,
                       nband=35,
                       nbdbuf=5,
                       )

    # Dataset3: Calculation of the screening.
    scr = abilab.AbiInput(pseudos=pseudos)
    scr.set_structure(structure)
    scr.set_kmesh(**gw_kmesh)
    scr.set_variables(**global_vars)

    scr.set_variables(
        optdriver=3,   
        getkss=2,      
        nband=25,    
        ecutwfn=ecutwfn,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
        ppmfrq="16.7 eV",
    )

    return gs, nscf, scr

@decorate_main
def main():
    gs, nscf, scr_input = all_inputs()

    policy = dict(autoparal=0, max_ncpus=2)
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=policy)

    workdir = "QPTDM"

    # This is to produce the out_WFK file
    #wfk_work = abilab.Workflow(workdir, manager)
    #gs_link = wfk_work.register(gs)
    #nscf_link = wfk_work.register(nscf, links=gs_link.produces_exts("DEN"))
    #wfk_work.start()
    #return 

    wfk_file = os.path.join(os.getcwd(), "out_WFK")
    qptdm_work = build_qptdm_workflow(workdir, manager, scr_input, wfk_file)
    #qptdm_work.connect(nscf_link.produces_exts("WFK"))
    #qptdm_work.connect(nscf_link, exts="WFK")

    qptdm_work.build_and_pickle_dump()

    return 

    works = AbinitWorks(workdir, manager)
    # One can register a workflow object.
    wlink0 = works.register(wfk_work)

    # Register a function that will be executed to build another workflow
    # the call back will have access the all the workflows that have been 
    # registered so far. The workflow will be generated at runtime and will
    # depend on the previous workflows specified in links.
    #wlink1 = works.register_callback(qptdm_work, links=wlink0)

    #wlink1 = works.register(sigma_work, links=[nscf_link.produce("WFK"), wlink1.produce("SCR")])

    #works.build_and_pickle_dump()

import collections
class AbinitWorks(collections.Iterable):
    """
    This object is a container of workflows. Its main task is managing the 
    possible inter-depencies among the workflows and the generation of the 
    dynamic worflows whose creations is done via callabacks.
    """
    def __init__(self, workdir, manager):
        self.workdir = os.path.abspath(workdir)

        self.manager = manager.deepcopy()

        self._works = []

        # Dict with the dependencies of each task, indexed by task.id
        self._links_dict = collections.defaultdict(list)

    def __len__(self):
        return len(self._works)

    def __iter__(self):
        return len(self._works)

    def __getitem__(self, slice):
        return self._works[slice]

    def register(self, obj, links=(), manager=None)
        """
        Registers a new workflow and add it to the internal list, taking into account possible dependencies.

        Args:
            obj:
                `Strategy` object or `AbinitInput` instance.
                if Strategy object, we create a new `AbinitTask` from the input strategy and add it to the list.
            links:
                List of `Link` objects specifying the dependency of this node.
                An empy list of links implies that this node has no dependencies.
            manager:
                The `TaskManager` responsible for the submission of the task. If manager is None, we use 
                the `TaskManager` specified during the creation of the workflow.

        Returns:   
            `Link` object
        """
        if links and not isinstance(links, (list, tuple)):
            links = [links]

        work_id = len(self)
        work_workdir = os.path.join(self.workdir, "work_" + str(work_id))

        # Make a deepcopy since manager is mutable and we might change it at run-time.
        manager = self.manager.deepcopy() if manager is None else manager.deepcopy()

        work = obj
        #if isinstance(obj, Workflow):
        #    work = obj
        #    work.set_observer(self)
        #else:
        #    # Callback.
    
        self._works.append(work)

        if links:
            self._links_dict[work_id].extend(links)
            logger.debug("work_id %s needs\n %s" % (work_id, [str(l) for l in links]))

        return Link(task)

        def notify(self, event, workflow):
            #if work.status == 

            # Finalize the workflow.
            retcode, message = workflow.finalize()

            for i, w in enumerate(self):
                if not w.depends_on(workflow): 
                    continue

                if isinstance(w, Workflow):
                    w.set_status(w.READY)
                    #w.connect()

                else:
                    # Assume callback: build the workflow here, register it 
                    # and set its status to READY.
                    w = w(self)
                    #w.connect()
                    w.set_status(w.READY)
                    self.works[i] = w



if __name__ == "__main__":
    import sys
    sys.exit(main())

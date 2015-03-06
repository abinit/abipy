#!/usr/bin/env python
"""Integration tests for the scheduler."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import abipy.data as abidata  
import abipy.abilab as abilab

from pymatgen.io.abinitio.mocks import mock_task_start


def make_scf_nscf_inputs(paral_kgb=1):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos("14si.pspnc")

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    structure = inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    ecut = 4
    global_vars = dict(ecut=ecut,
                       nband=8,
                       nstep=15,
                       paral_kgb=paral_kgb,
                    )

    if inp.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    inp.set_vars(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[2,2,2], shiftk=[0,0,0])
    inp[1].set_vars(tolvrs=1e-4)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    inp[2].set_kpath(ndivsm=2, kptbounds=kptbounds)
    inp[2].set_vars(tolwfr=1e-4)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()
    return scf_input, nscf_input


def itest_flow_with_deadlocks(fwp):
    """
    Test the behaviour of the scheduler in the presence of a deadlock
    when we ignore errored tasks and we try to run all tasks in the flow.
    The scheduler should detect the deadlock and exit when no other task can be executed.
    """
    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs()

    # Build the flow.
    flow = abilab.Flow(fwp.workdir, manager=fwp.manager)
    work0 = abilab.BandStructureWork(scf_input, nscf_input, dos_inputs=nscf_input)
    flow.register_work(work0)
    scf_task, nscf_task, dos_task = work0[0], work0[1], work0[2]

    work1 = abilab.Work()
    work1.register_nscf_task(nscf_input, deps={scf_task: "DEN", dos_task: "WFK"})
    # This task will deadlock when nscf_task reaches S_ERROR.
    work1.register_nscf_task(nscf_input, deps={scf_task: "DEN", nscf_task: "WFK"})
    flow.register_work(work1)

    flow.allocate()

    # Mock an Errored nscf_task. This will cause a deadlock in the flow.
    nscf_task = mock_task_start(nscf_task)

    # Here we set max_num_abierrs to a very large number.
    sched = flow.make_scheduler()
    sched.max_num_abierrs = 10000
    sched.start()
    flow.check_status(show=True)

    assert not flow.all_ok
    assert all(task.status == task.S_OK for task in [scf_task, dos_task, work1[0]])
    assert all(task.status == task.S_ERROR for task in [nscf_task])
    deadlocked, runnables, running = flow.deadlocked_runnables_running()
    assert deadlocked and not runnables and not running
    assert work1[1] in deadlocked
    #assert 0


def itest_flow_without_runnable_tasks(fwp):
    """
    Test the behaviour of the scheduler when we ignore errrors and 
    all the task that can be executed have been submitted.
    The scheduler should detect this condition and exit.
    """
    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs()

    # Build the flow.
    flow = abilab.Flow(fwp.workdir, manager=fwp.manager)
    work0 = abilab.BandStructureWork(scf_input, nscf_input)
    flow.register_work(work0)
    scf_task, nscf_task = work0.scf_task, work0.nscf_task

    flow.allocate()

    # Mock an Errored nscf_task. This will cause a deadlock in the flow.
    nscf_task = mock_task_start(nscf_task)

    sched = flow.make_scheduler()
    sched.max_num_abierrs = 10000
    sched.start()
    flow.check_status(show=True)

    assert not flow.all_ok
    assert scf_task.status == scf_task.S_OK
    assert nscf_task.status == nscf_task.S_ERROR

    deadlocked, runnables, running = flow.deadlocked_runnables_running()
    assert not deadlocked and not runnables and not running
    #assert 0

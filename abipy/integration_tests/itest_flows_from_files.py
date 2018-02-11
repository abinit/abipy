"""
Integration tests for flows/works/tasks that rely on external files e.g. DEN --> NscfTask.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import pytest
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk
#from abipy.core.testing import has_abinit, has_matplotlib


def make_scf_nscf_inputs(paral_kgb=1):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)

    # Global variables
    ecut = 6
    global_vars = dict(
        ecut=ecut,
        nband=8,
        nstep=15,
        paral_kgb=paral_kgb,
    )

    if multi.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(ngkpt=[8, 8, 8], shiftk=[0, 0, 0])
    multi[0].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
    multi[1].set_vars(tolwfr=1e-12)

    # Generate two input files for the GS and the NSCF run
    scf_input, nscf_input = multi.split_datasets()

    return scf_input, nscf_input


def itest_nscf_from_denfile(fwp, tvars):
    """Testing NscTask from pre-existent DEN file."""
    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs(paral_kgb=1)

    # Build the flow.
    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

    # Create a Work, all tasks in work will start from the DEN file.
    # Note that the file must exist when the work is created
    # Use the standard approach based on tasks and works if
    # there's a node who needs a file produced in the future.
    #work = flowtk.Work()
    den_filepath = abidata.ref_file("si_DEN.nc")
    flow.register_nscf_task(nscf_input, deps={den_filepath: "DEN"}, append=True)

    flow.allocate(use_smartio=True)
    assert len(flow) == 1 and len(flow[0]) == 1
    # Will remove output files (WFK)
    #flow.set_garbage_collector()

    scheduler = flow.make_scheduler()
    assert scheduler.start() == 0
    assert not scheduler.exceptions
    assert scheduler.nlaunch == 1

    flow.check_status(show=True)
    assert flow.all_ok
    assert all(work.finalized for work in flow)

    # The WFK files should have been removed because we called set_garbage_collector
    # but the input DEN.nc should exist
    #for task in flow[0]:
    #    assert not task.outdir.has_abiext("WFK")
    assert os.path.isfile(den_filepath)

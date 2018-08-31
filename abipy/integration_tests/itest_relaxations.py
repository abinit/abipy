"""Integration tests for structural relaxations."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import numpy as np
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk

from abipy.core.testing import has_abinit, has_matplotlib


def ion_relaxation(tvars, ntime=50):
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    # Perturb the structure (random perturbation of 0.1 Angstrom)
    structure.perturb(distance=0.02)

    global_vars = dict(
        ecut=6,
        ngkpt=[2,2,2],
        shiftk=[0,0,0],
        nshiftk=1,
        chksymbreak=0,
        paral_kgb=tvars.paral_kgb,
    )

    inp = abilab.AbinitInput(structure, pseudos=abidata.pseudos("14si.pspnc"))

    # Global variables
    inp.set_vars(global_vars)

    # Dataset 1 (Atom Relaxation)
    #inp[1].set_vars(
    # FIXME here there's a bug
    inp.set_vars(
        optcell=0,
        ionmov=2,
        tolrff=0.02,
        tolmxf=5.0e-5,
        ntime=ntime,
        #dilatmx=1.05, # FIXME: abinit crashes if I don't use this
    )

    return inp


def itest_atomic_relaxation(fwp, tvars):
    """Test atomic relaxation with automatic restart."""
    # Build the flow
    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

    ion_input = ion_relaxation(tvars, ntime=2)
    work = flow.register_task(ion_input, task_class=flowtk.RelaxTask)

    flow.allocate()
    flow.build_and_pickle_dump(abivalidate=True)

    # Run t0, and check status
    t0 = work[0]
    t0.start_and_wait()
    assert t0.returncode == 0
    t0.check_status()
    assert t0.status == t0.S_UNCONVERGED

    assert t0.initial_structure == ion_input.structure
    unconverged_structure = t0.get_final_structure()
    assert unconverged_structure != t0.initial_structure

    # Use the default value ntime=50 and we can converge the calculation.
    t0.input.set_vars(ntime=50)

    t0.build()
    assert t0.restart()
    t0.wait()
    assert t0.num_restarts == 1

    # At this point, we should have reached S_OK.
    assert t0.status == t0.S_DONE
    t0.check_status()
    assert t0.status == t0.S_OK

    final_structure = t0.get_final_structure()
    assert final_structure != unconverged_structure

    flow.show_status()
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    # post-processing tools
    if has_matplotlib():
        assert t0.inspect(show=False) is not None

    with t0.open_hist() as hist:
        hist.to_string(verbose=2)
        # from_file accepts HIST files as well.
        assert np.all(hist.structures[-1].frac_coords == abilab.Structure.from_file(hist.filepath).frac_coords)

    with t0.open_gsr() as gsr:
        gsr.to_string(verbose=2)
        gsr.pressure == 1.8280

    t0.get_results()


def make_ion_ioncell_inputs(tvars, dilatmx, scalevol=1, ntime=50):
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    # Perturb the structure (random perturbation of 0.1 Angstrom)
    #structure.perturb(distance=0.01)

    # Compress the lattice so that ABINIT complains about dilatmx
    structure.scale_lattice(structure.volume * scalevol)

    global_vars = dict(
        ecut=6,
        ecutsm=0.5,
        ngkpt=[4,4,4],
        shiftk=[0,0,0],
        nshiftk=1,
        chksymbreak=0,
        paral_kgb=tvars.paral_kgb,
    )

    multi = abilab.MultiDataset(structure, pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)

    # Global variables
    multi.set_vars(global_vars)

    # Dataset 1 (Atom Relaxation)
    multi[0].set_vars(
        optcell=0,
        ionmov=2,
        tolrff=0.02,
        tolmxf=5.0e-5,
        ntime=ntime,
    )

    # Dataset 2 (Atom + Cell Relaxation)
    multi[1].set_vars(
        optcell=1,
        ionmov=2,
        dilatmx=dilatmx,
        tolrff=0.02,
        tolmxf=5.0e-5,
        strfact=100,
        ntime=ntime,
    )

    ion_inp, ioncell_inp = multi.split_datasets()
    return ion_inp, ioncell_inp


def itest_relaxation_with_restart_from_den(fwp, tvars):
    """Test structural relaxations with automatic restart from DEN files."""
    # Build the flow
    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

    # Use small value for ntime to trigger restart, then disable the output of the WFK file.
    ion_input, ioncell_input = make_ion_ioncell_inputs(tvars, dilatmx=1.1, ntime=3)
    ion_input.set_vars(prtwf=0)
    ioncell_input.set_vars(prtwf=0)

    relax_work = flowtk.RelaxWork(ion_input, ioncell_input)
    flow.register_work(relax_work)

    assert flow.make_scheduler().start() == 0
    flow.show_status()

    assert all(work.finalized for work in flow)
    assert flow.all_ok

    # we should have (0, 1) restarts and no WFK file in outdir.
    for i, task in enumerate(relax_work):
        assert task.status == task.S_OK
        assert task.num_restarts == i
        assert task.num_corrections == 0
        assert not task.outdir.has_abiext("WFK")

    if has_matplotlib:
        assert relax_work.plot_ion_relaxation(show=False) is not None
        assert relax_work.plot_ioncell_relaxation(show=False) is not None

    flow.rmtree()


def itest_dilatmx_error_handler(fwp, tvars):
     """
     Test cell relaxation with automatic restart in the presence of dilatmx error.
     """
     # Build the flow
     flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

     # Decrease the volume to trigger DilatmxError
     ion_input, ioncell_input = make_ion_ioncell_inputs(tvars, dilatmx=1.01, scalevol=0.8)

     work = flowtk.Work()
     work.register_relax_task(ioncell_input)

     flow.register_work(work)
     flow.allocate()
     assert flow.make_scheduler().start() == 0
     flow.show_status()

     assert all(work.finalized for work in flow)
     assert flow.all_ok

     # t0 should have reached S_OK, and we should have DilatmxError in the corrections.
     t0 = work[0]
     assert t0.status == t0.S_OK
     print(t0.corrections)
     assert t0.num_corrections > 0
     assert t0.corrections[0]["event"]["@class"] == "DilatmxError"


def itest_relaxation_with_target_dilatmx(fwp, tvars):
    """Test structural relaxations with automatic restart from DEN files."""
    # Build the flow
    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

    # Use small value for ntime to trigger restart, then disable the output of the WFK file.
    ion_input, ioncell_input = make_ion_ioncell_inputs(tvars, dilatmx=1.05)

    target_dilatmx = 1.03
    relax_work = flowtk.RelaxWork(ion_input, ioncell_input, target_dilatmx=target_dilatmx)
    flow.register_work(relax_work)

    assert flow.make_scheduler().start() == 0
    flow.show_status()

    assert all(work.finalized for work in flow)
    assert flow.all_ok
    #assert relax_work.last_dilatmx <= target_dilatmx

    # we should have (0, 1) restarts
    for i, task in enumerate(relax_work):
        assert task.status == task.S_OK
        assert task.num_restarts == i
        assert task.num_corrections == 0

    assert relax_work[1].input["dilatmx"] == 1.03

    flow.rmtree()

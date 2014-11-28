"""Integration tests for structural relaxations."""
from __future__ import print_function, division, unicode_literals

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.core.testing import has_abinit

# Tests in this module require abinit >= 7.9.0
#pytestmark = pytest.mark.skipif(not has_abinit("7.9.0"), reason="Requires abinit >= 7.9.0")


def ion_relaxation(tvars, ntime=50):
    pseudos = abidata.pseudos("14si.pspnc")
    cif_file = abidata.cif_file("si.cif")
    structure = abilab.Structure.from_file(cif_file)

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

    inp = abilab.AbiInput(pseudos=pseudos)
    inp.set_structure(structure)

    # Global variables
    inp.set_variables(**global_vars)

    # Dataset 1 (Atom Relaxation)
    inp[1].set_variables(
        optcell=0,
        ionmov=2,
        tolrff=0.02,
        tolmxf=5.0e-5,
        ntime=ntime,
        #dilatmx=1.05, # FIXME: abinit crashes if I don't use this
    )

    return inp


def make_ion_ioncell_inputs(tvars):
    cif_file = abidata.cif_file("si.cif")
    structure = abilab.Structure.from_file(cif_file)

    # Perturb the structure (random perturbation of 0.1 Angstrom)
    structure.perturb(distance=0.01)

    pseudos = abidata.pseudos("14si.pspnc")

    global_vars = dict(
        ecut=6,
        ngkpt=[4,4,4],
        shiftk=[0,0,0],
        nshiftk=1,
        chksymbreak=0,
        paral_kgb=tvars.paral_kgb,
    )

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    inp.set_structure(structure)

    # Global variables
    inp.set_variables(**global_vars)

    # Dataset 1 (Atom Relaxation)
    inp[1].set_variables(
        optcell=0,
        ionmov=2,
        tolrff=0.02,
        tolmxf=5.0e-5,
        ntime=50,
        #ntime=5, To test the restart
        dilatmx=1.05, # FIXME: abinit crashes if I don't use this
    )

    # Dataset 2 (Atom + Cell Relaxation)
    inp[2].set_variables(
        optcell=1,
        ionmov=2,
        ecutsm=0.5,
        dilatmx=1.05,
        tolrff=0.02,
        tolmxf=5.0e-5,
        strfact=100,
        ntime=50,
        #ntime=5, To test the restart
        )

    ion_inp, ioncell_inp = inp.split_datasets()
    return ion_inp, ioncell_inp


def itest_simple_atomic_relaxation(fwp, tvars):
    """
    Test ion relaxation with automatic restart.
    """
    # Build the flow
    flow = abilab.AbinitFlow(fwp.workdir, manager=fwp.manager)

    ion_input = ion_relaxation(tvars, ntime=2)
    work = flow.register_task(ion_input, task_class=abilab.RelaxTask)
    t0 = work[0]

    flow.allocate()
    flow.build_and_pickle_dump()

    # Run t0, and check status
    assert t0.returncode == 0
    t0.start_and_wait()
    assert t0.returncode == 0
    t0.check_status()
    assert t0.status == t0.S_UNCONVERGED

    # Remove ntime from the input so that the next run will
    # use the default value ntime=50 and we can converge the calculation.
    # This one does not work
    #t0.strategy.remove_extra_abivars(["ntime"])
    t0.strategy.add_extra_abivars(dict(ntime=50))
    #t0.strategy.abinit_input.pop("ntime")
    #t0.strategy.abinit_input.remove_variables("ntime")
    #t0.strategy.abinit_input.set_variables(ntime=50)
    print("new input:\n", t0.strategy.abinit_input)
    t0.build()
    assert t0.restart()
    t0.wait()
    assert t0.num_restarts == 1

    assert t0.status == t0.S_DONE
    t0.check_status()
    assert t0.status == t0.S_OK

    flow.show_status()
    assert all(work.finalized for work in flow)
    assert flow.all_ok
    #assert flow.validate_json_schema()

#!/usr/bin/env python
"""
Integration tests for flows (require pytest, ABINIT and a properly configured environment)
"""
from __future__ import division, print_function

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from pymatgen.core.design_patterns import AttrDict
from abipy.core.testing import has_abinit


has_abinit_ge790 = has_abinit("7.9.0")


def make_scf_nscf_inputs(paral_kgb, pp_paths):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos(pp_paths)

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    inp.set_structure_from_file(abidata.cif_file("si.cif"))

    # Global variables
    ecut = 4
    global_vars = dict(ecut=ecut,
                       nband=6,
                       paral_kgb=paral_kgb)

    if inp.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[1, 1, 1], shiftk=[0, 0, 0])
    inp[1].set_variables(tolvrs=1e-4)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0],  # L point
        [0.0, 0.0, 0.0],  # Gamma point
        [0.0, 0.5, 0.5],  # X point
    ]

    inp[2].set_kpath(ndivsm=2, kptbounds=kptbounds)
    inp[2].set_variables(tolwfr=1e-6)
    
    # Generate two input files for the GS and the NSCF run.
    scf_input, nscf_input = inp.split_datasets()

    return scf_input, nscf_input


@pytest.mark.skipif(not has_abinit_ge790, reason="Requires abinit >= 7.9.0")
@pytest.mark.parametrize("inp", [{"paral_kgb": 0}, {"paral_kgb": 1}])
def test_bandstructure_flow(fwp, inp):
    """
    Test the building of a bandstructure flow and autoparal.
    Simple flow with one dependency: SCF -> NSCF.
    """
    workdir = fwp.workdir
    inp = AttrDict(inp)
    print("Input:", inp)

    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs(paral_kgb=inp.paral_kgb, pp_paths="14si.pspnc")

    # Build the flow and create the database.
    flow = abilab.bandstructure_flow(workdir, fwp.manager, scf_input, nscf_input)

    flow.build_and_pickle_dump()

    t0 = flow[0][0]
    t1 = flow[0][1]

    # Test initialization status, task properties and links (t1 requires t0)
    assert t0.status == t0.S_INIT
    assert t1.status == t1.S_INIT
    assert t0.can_run
    assert not t1.can_run
    assert t1.depends_on(t0)
    assert not t1.depends_on(t1)
    assert t0.isnc
    assert not t0.ispaw

    # Flow properties and methods
    assert flow.num_tasks == 2
    assert not flow.all_ok
    assert flow.ncpus_reserved == 0
    assert flow.ncpus_allocated == 0
    assert flow.ncpus_inuse == 0
    flow.check_dependencies()
    flow.show_status()
    flow.show_receivers()

    # Testing the autoparal feature
    for work in flow:
        for task in work:
            print(task)
            #confs, optimal = task.autoparal_fake_run()
            #print(confs, optimal)
            #aequal(task.status, task.S_INIT)
            #afalse(task.is_completed)
            #self.assertTrue(task.can_run)

    # Run t0, and check status
    assert t0.returncode == 0
    t0.start_and_wait()
    assert t0.returncode == 0
    assert t0.status == t0.S_DONE
    t0.check_status()
    assert t0.status == t0.S_OK
    assert t1.can_run

    # Cannot start t0 twice...
    with pytest.raises(t0.Error):
        t0.start()

    # But we can restart it.
    fired = t0.restart()
    assert fired
    t0.wait()
    assert t0.num_restarts == 1
    assert t0.status == t0.S_DONE
    t0.check_status()
    assert t0.status == t0.S_OK
    assert not t0.can_run

    # Now we can run t1.
    t1.start_and_wait()
    assert t1.status == t0.S_DONE
    t1.check_status()
    assert t1.status == t1.S_OK
    assert not t1.can_run

    # This one does not work yet
    #fired = t1.restart()
    #atrue(fired)
    #t1.wait()
    #aequal(t1.num_restarts, 1)
    #aequal(t1.status, t1.S_DONE)
    #t1.check_status()
    #aequal(t1.status, t1.S_OK)
    #afalse(t1.can_run)

    assert flow.all_ok
    #assert 0

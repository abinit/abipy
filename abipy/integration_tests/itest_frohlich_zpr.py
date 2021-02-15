"""Integration tests for phonon flows."""

#import numpy as np
import abipy.data as abidata
import abipy.abilab as abilab
#import abipy.flowtk as flowtk

from abipy.dfpt.ddb import DdbFile


def make_scf_input(usepaw=0):
    """Returns the GS input file"""
    # Here we use parameters similar to https://docs.abinit.org/tests/v8/Input/t57.in
    pseudos = abidata.pseudos("Ca.psp8", "O.psp8")

    structure = dict(
        acell=3 * [9.136],
        xred=[
           0.0000000000, 0.0000000000, 0.0000000000,
           0.5000000000, 0.5000000000, 0.5000000000],
        rprim=[
           0  , 0.5, 0.5,
           0.5, 0  , 0.5,
           0.5, 0.5, 0],
        typat=[1, 2],
        natom=2,
        ntypat=2,
        znucl=[20, 8],
    )

    scf_input = abilab.AbinitInput(structure=structure, pseudos=pseudos)

    scf_input.set_vars(
        nband=10,
        nbdbuf=2,
        diemac=6,
        #ecut=30,               # Underconverged ecut.
        ecut=15,
        nstep=100,
        tolvrs=1e-6,
        kptrlatt=[-2,  2,  2,  # In cartesian coordinates, this grid is simple cubic
                   2, -2,  2,
                   2,  2, -2],
    )

    return scf_input


def itest_frohlich_zpr_flow(fwp, tvars):
    """
    """
    # Build the SCF input.
    scf_input = make_scf_input()

    # Build the flow.
    from abipy.flowtk.effmass_works import FrohlichZPRFlow
    flow = FrohlichZPRFlow.from_scf_input(fwp.workdir, scf_input, ddb_node=None, ndivsm=2, tolwfr=1e-10,
                                          manager=fwp.manager, metadata=dict(mp_id="mp-123"))

    scheduler = flow.make_scheduler()
    assert scheduler.start() == 0
    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()
    assert flow.on_all_ok_num_calls == 1

    # Reconstruct python objects from JSON file.
    data = abilab.mjson_load(flow.outdir.path_in("zprfrohl_results.json"))

    assert data["metadata"]["mp_id"] == "mp-123"
    assert data["structure"].formula == 'Ca1 O1'
    ebands_kpath = data["ebands_kpath"]
    assert ebands_kpath.nsppol == 1
    assert ebands_kpath.kpoints.is_path
    assert ebands_kpath.homos[0].kpoint == [0, 0, 0]
    assert ebands_kpath.lumos[0].kpoint == [0.5, 0, 0.5]
    with DdbFile.from_string(data["ddb_string"]) as ddb:
        assert ddb.structure.formula == 'Ca1 O1'
        assert ddb.has_bec_terms(select="at_least_one")
        assert ddb.has_epsinf_terms(select="at_least_one_diagoterm")
    assert data["epsinf_cart"].shape == data["eps0_cart"].shape

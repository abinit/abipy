"""Integration tests for phonon flows."""

import numpy as np
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk
from abipy.core.testing import has_matplotlib


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
        nband=12,
        nbdbuf=2,
        diemac=6,
        #ecut=30,               # Underconverged ecut.
        ecut=15,
        nstep=100,
        tolvrs=1e-8,
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
    flow = FrohlichZPRFlow.from_scf_input(fwp.workdir, scf_input, ddb_node=None, ndivsm=3, tolwfr=1e-10,
                                           manager=fwp.manager, metadata=dict(mpi_id=123))

    scheduler = flow.make_scheduler()
    assert scheduler.start() == 0
    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    #with flow.open_final_ddb() as ddb:
    #    ddb_path = ddb.filepath
    #    ddb.to_string(verbose=2)
    #    assert len(ddb.structure) == 2

    ## These output files should be produced in the task workdir.
    ## Actually they should be in outdir but anaddb uses different conventions.
    #assert len(atask.wdir.list_filepaths(wildcard="*PHBST.nc")) == 1
    #assert len(atask.wdir.list_filepaths(wildcard="*PHDOS.nc")) == 1

    data = abilab.mjson_load(flow.outdir.path_in("zprfrohl_results.json"))

    assert data["mp_id"] == 123
    assert data["structure"].formula == "CaO"
    assert data["ebands_kpath"].nsppol == 1
    assert data["ebands_kpath"].kpoints.is_kpath
    #assert data["ddb"].structure == 2
    #assert data["epsinf_cart"].shape == data["eps0_cart"].shape

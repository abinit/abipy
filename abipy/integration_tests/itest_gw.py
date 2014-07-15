"""Integration tests for GW flows."""
from __future__ import division, print_function

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from pymatgen.io.abinitio.calculations import bse_with_mdf
from abipy.core.testing import has_abinit


# Tests in this module require abinit >= 7.9.0
pytestmark = pytest.mark.skipif(not has_abinit("7.9.0"), reason="Requires abinit >= 7.9.0")


def make_inputs(ngkpt, tvars):
    """
    Calculation of the GW correction
    Dataset 1: ground state calculation
    Dataset 2: NSCF calculation
    Dataset 3: calculation of the screening
    Dataset 4-5-6: Self-Energy matrix elements (GW corrections) with different values of nband
    """
    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=4)
    inp.set_structure_from_file(abidata.cif_file("si.cif"))

    # This grid is the most economical, but does not contain the Gamma point.
    scf_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap.
    gw_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.0, 0.0, 0.0,
                0.0, 0.5, 0.5,
                0.5, 0.0, 0.5,
                0.5, 0.5, 0.0]
    )

    # Global variables. gw_kmesh is used in all datasets except DATASET 1.
    ecut = 4

    inp.set_variables(
        ecut=ecut,
        pawecutdg=ecut*2 if inp.pseudos.allpaw else None,
        istwfk="*1",
        paral_kgb=tvars.paral_kgb,
        gwpara=2,
    )
    inp.set_kmesh(**gw_kmesh)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(**scf_kmesh)
    inp[1].set_variables(
        tolvrs=1e-6,
        nband=4)

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    inp[2].set_variables(iscf=-2,
                         tolwfr=1e-10,
                         nband=10,
                         nbdbuf=2)

    # Dataset3: Calculation of the screening.
    inp[3].set_variables(
        optdriver=3,
        nband=8,
        ecutwfn=ecut,
        symchi=1,
        inclvkb=0,
        ecuteps=2.0,
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    kptgw = [
         -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
         -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
          5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
         -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
          5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
          0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
      ]

    inp[4].set_variables(
            optdriver=4,
            nband=10,
            ecutwfn=ecut,
            ecuteps=2.0,
            ecutsigx=2.0,
            symsigma=1,
            #nkptgw=0,
            #gw_qprange=0,
    )

    bdgw = [4, 5]
    inp[4].set_kptgw(kptgw, bdgw)

    return inp.split_datasets()


def itest_g0w0_flow(fwp, tvars):
    """Test flow for G0W0 calculations."""
    scf, nscf, scr, sig = make_inputs(ngkpt=[2, 2, 2], tvars=tvars)

    flow = abilab.g0w0_flow(fwp.workdir, fwp.manager, scf, nscf, scr, sig)
    flow.build_and_pickle_dump()

    for task in flow[0]:
        task.start_and_wait()

    flow.check_status()
    flow.show_status()
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    # The sigma task should produce a SIGRES file.
    assert len(flow[0][-1].outdir.list_filepaths(wildcard="*SIGRES.nc")) == 1


#def itest_bse_with_mdf(fwp, tvars):
#    pseudos = abidata.pseudos("14si.pspnc")
#    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
#
#    kppa = scf_kppa = 1
#    nscf_nband = 6
#    nscf_ngkpt = [4,4,4]
#    nscf_shiftk = [0.1, 0.2, 0.3]
#    bs_loband = 2
#    bs_nband = nscf_nband
#    soenergy = 0.7
#    mdf_epsinf = 12
#    max_ncpus = 1
#    ecuteps = 2
#
#    extra_abivars = dict(
#        ecut=12,
#        istwfk="*1",
#    )
#
#    flow = abilab.AbinitFlow(workdir=fwp.workdir, manager=fwp.manager)
#
#    # BSE calculation with model dielectric function.
#    work = bse_with_mdf(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
#                       ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf,
#                       accuracy="normal", spin_mode="unpolarized", smearing=None,
#                       charge=0.0, scf_solver=None, **extra_abivars)
#
#    flow.register_work(work)
#    flow.allocate()
#    flow.build_and_pickle_dump()
#
#    for task in flow:
#        task.start_and_wait()
#        assert task.status == task.S_DONE
#
#    flow.check_status()
#    flow.show_status()
#    assert all(work.finalized for work in flow)
#    assert flow.all_ok


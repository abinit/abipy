"""
Integration tests for flows (require pytest, ABINIT and a properly configured environment)
"""
import os
import unittest
#import numpy.testing.utils as nptu
import numpy.testing as nptu
import abipy.data as abidata
import abipy.flowtk as flowtk
import abipy.flowtk.abiphonopy as abiph

from abipy.abio.factories import gs_input
from abipy.core.testing import has_phonopy


def itest_phonopy_flow(fwp, tvars):
    """
    Testing phonopy flow with the scheduler.
    """
    if not has_phonopy():
        raise unittest.SkipTest("This test requires phonopy")

    #print("tvars:\n %s" % str(tvars))
    si_structure = abidata.structure_from_cif("si.cif")

    gsinp = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), kppa=10, ecut=2, spin_mode="unpolarized")
    gsinp["paral_kgb"] = tvars.paral_kgb

    flow = flowtk.Flow(workdir=fwp.workdir, manager=fwp.manager)
    scdims = [2, 2, 2]
    phpy_work = abiph.PhonopyWork.from_gs_input(gsinp, scdims=scdims,
                                                phonopy_kwargs=None, displ_kwargs=None)
    flow.register_work(phpy_work)

    assert hasattr(phpy_work, "phonon")
    assert len(phpy_work.phonopy_tasks) == len(phpy_work)
    assert len(phpy_work.phonopy_tasks) == 1
    assert len(phpy_work.bec_tasks) == 0
    nptu.assert_equal(scdims, phpy_work.scdims)

    # Will remove output files (WFK)
    flow.set_garbage_collector()
    flow.use_smartio()

    flow.build_and_pickle_dump(abivalidate=True)
    scheduler = flow.make_scheduler()
    assert scheduler.start() == 0

    assert not scheduler.exceptions
    flow.show_status()
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()
    assert all(work.finalized for work in flow)

    # The WFK files should have been removed because we called set_garbage_collector
    for task in flow[0]:
        assert not task.outdir.has_abiext("WFK")

    out_filenames = set([os.path.basename(f) for f in phpy_work.outdir.list_filepaths()])
    nmiss = 0
    for f in ["POSCAR", "disp.yaml", "FORCE_SETS", "band.conf", "dos.conf", "band-dos.conf", "README.md"]:
        if f not in out_filenames:
            nmiss += 1
            print("Cannot find %s in work.outdir" % f)
    assert nmiss == 0


def itest_phonopy_gruneisen_flow(fwp, tvars):
    """
    Testing phonopy Gruneisen flow with the scheduler.
    """
    if not has_phonopy():
        raise unittest.SkipTest("This test requires phonopy")

    #print("tvars:\n %s" % str(tvars))
    si_structure = abidata.structure_from_cif("si.cif")

    gsinp = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), kppa=10, ecut=2, spin_mode="unpolarized")
    gsinp["paral_kgb"] = tvars.paral_kgb

    flow = flowtk.Flow(workdir=fwp.workdir, manager=fwp.manager)
    scdims = [2, 2, 2]

    # Grunesein with phonopy
    grun_work = abiph.PhonopyGruneisenWork.from_gs_input(gsinp, voldelta=0.1, scdims=scdims,
                                                         phonopy_kwargs=None, displ_kwargs=None)
    flow.register_work(grun_work)
    assert len(grun_work) == 3
    nptu.assert_equal(scdims, grun_work.scdims)

    # Will remove output files (WFK)
    flow.set_garbage_collector()
    flow.use_smartio()

    flow.build_and_pickle_dump(abivalidate=True)
    scheduler = flow.make_scheduler()
    assert scheduler.start() == 0

    assert not scheduler.exceptions
    flow.show_status()
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()
    # Initialial work + 3 phonopy works.
    assert len(flow) == 4
    assert all(work.finalized for work in flow)

    # The WFK files should have been removed because we called set_garbage_collector
    # FIXME: This does not work because new works that have been created.
    #for task in flow.iflat_tasks():
    #    assert not task.outdir.has_abiext("WFK")

    for work in flow[1:]:
        out_filenames = set([os.path.basename(f) for f in work.outdir.list_filepaths()])
        nmiss = 0
        for f in ["POSCAR", "disp.yaml", "FORCE_SETS", "band.conf", "dos.conf", "band-dos.conf", "README.md"]:
            if f not in out_filenames:
                nmiss += 1
                print("Cannot find %s in work.outdir" % f)
        assert nmiss == 0

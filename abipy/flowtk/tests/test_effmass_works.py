"""Tests for abiphonopy module"""
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
from abipy.abio.factories import gs_input
from abipy.flowtk.effmass_works import EffectiveMassLineWork


class TestEffectiveMassLineWork(AbipyTest):

    def test_work_apu(self):
        """Testing EffectiveMassLineWork."""
        si_structure = abidata.structure_from_cif("si.cif")

        scf_input = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4, spin_mode="unpolarized")
        flow = flowtk.Flow.temporary_flow()

        with self.assertRaises(ValueError):
            work = EffectiveMassLineWork.from_scf_input(scf_input, kpoint=(0, 0, 0), step=0.01, npts=2)

        #red_dirs = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        work = EffectiveMassLineWork.from_scf_input(scf_input, kpoint=(0, 0, 0), step=0.01, npts=9)
                                                    #red_dirs=[1, 0, 0])

        t1 = work[1]
        assert t1.input["kptopt"] == 0
        assert t1.input["iscf"] == -2

        #assert len(work.relax_tasks) == 3
        #assert all(t.input["optcell"] == 3 for t in work)

        flow.register_work(work)
        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

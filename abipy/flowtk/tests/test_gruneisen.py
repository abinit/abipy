"""Tests for abiphonopy module"""
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
from abipy.abio.factories import gs_input
from abipy.flowtk.gruneisen import GruneisenWork


class TestGruneisenWork(AbipyTest):
    """Unit tests for gruneisen module."""

    def test_gruneisen_work(self):
        """Testing GruneisenWork."""
        si_structure = abidata.structure_from_cif("si.cif")

        gs_inp = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4, spin_mode="unpolarized")
        flow = flowtk.Flow.temporary_flow()
        voldelta = gs_inp.structure.volume * 0.02
        work = GruneisenWork.from_gs_input(gs_inp, voldelta, ngqpt=[2, 2, 2], with_becs=False)
        flow.register_work(work)

        assert len(work.relax_tasks) == 3
        assert all(t.input["optcell"] == 3 for t in work)
        assert all(t.input["ionmov"] == 3 for t in work)
        assert all(t.input["ecutsm"] == 0.5 for t in work)
        assert all(t.input["dilatmx"] == 1.05 for t in work)

        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

        with self.assertRaises(ValueError):
            GruneisenWork.from_gs_input(gs_inp, voldelta, ngqpt=[3, 3, 3], with_becs=False)

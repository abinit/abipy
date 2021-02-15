"""Tests for mocks module"""
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
from abipy.abio.factories import gs_input
from abipy.flowtk import mocks


class TestMocks(AbipyTest):
    """Unit tests for mocks module."""
    def test_infinite_flow(self):
        si_structure = abidata.structure_from_cif("si.cif")
        gsinp = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4)

        flow = flowtk.Flow.temporary_flow()
        work = flowtk.Work()
        gstask = work.register_scf_task(gsinp)
        flow.register_work(work)
        flow.allocate()

        mocks.infinite_flow(flow)
        flow.check_status()
        assert (t.status == flow.S_INIT for t in flow)

        mocks.change_task_start(gstask, mocked_status="Error")
        assert gstask.start() == 1 and gstask.status == gstask.S_ERROR

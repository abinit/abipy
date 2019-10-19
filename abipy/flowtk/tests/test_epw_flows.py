"""Tests for eph_flows module"""
from abipy.core.testing import AbipyTest
from abipy.flowtk.eph_flows import GkqPathFlow


class TestEphFlows(AbipyTest):
    """Unit tests for eph_flows module."""

    def test_gkqpath_flows(self):
        """Testing GkqPathFlow."""
        self.skip_if_abinit_not_ge("8.11.0")

        ngkpt = [4, 4, 4]
        gs_inp = self.get_gsinput_alas_ngkpt(ngkpt=ngkpt)

        ngqpt = (2, 2, 2)
        qbounds = [[0.01, 0, 0], [0.02, 0, 0], [0.04, 0, 0], [0.05, 0, 0]]
        qbounds = [[0, 0, 0], [0.5, 0, 0]]

        workdir = self.mkdtemp()
        flow = GkqPathFlow.from_scf_input(workdir, gs_inp, ngqpt, qbounds, ndivsm=2, with_becs=True, ddk_tolerance=None,
                                          test_ft_interpolation=True, prepgkk=0)

        #assert len(work.relax_tasks) == 3
        #assert all(t.input["optcell"] == 3 for t in work)
        #assert all(t.input["ionmov"] == 3 for t in work)
        #assert all(t.input["ecutsm"] == 0.5 for t in work)
        #assert all(t.input["dilatmx"] == 1.05 for t in work)

        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

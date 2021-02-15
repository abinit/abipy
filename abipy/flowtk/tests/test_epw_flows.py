"""Tests for eph_flows module"""

from abipy.core.testing import AbipyTest
from abipy.flowtk.eph_flows import EphPotFlow, GkqPathFlow


class TestEphFlows(AbipyTest):
    """Unit tests for eph_flows module."""

    def test_ephpot_flows(self):
        """Testing EphPotFlow."""
        self.skip_if_abinit_not_ge("8.11.0")

        gs_inp = self.get_gsinput_alas_ngkpt(ngkpt=[4, 4, 4])

        ngqpt = (2, 2, 2)
        qbounds = [[0, 0, 0], [0.5, 0, 0]]

        workdir = self.mkdtemp()
        flow = EphPotFlow.from_scf_input(workdir, gs_inp, ngqpt, qbounds, ndivsm=2, with_becs=True, ddk_tolerance=None,
                                         prepgkk=1)

        #assert len(work.relax_tasks) == 3
        #assert all(t.input["dilatmx"] == 1.05 for t in work)

        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok


    def test_gkqpath_flows(self):
        """Testing GkqPathFlow."""
        self.skip_if_abinit_not_ge("8.11.0")

        gs_inp = self.get_gsinput_alas_ngkpt(ngkpt=[4, 4, 4])
        ngqpt = (2, 2, 2)
        qbounds = [[0, 0, 0], [0.5, 0, 0]]

        workdir = self.mkdtemp()
        flow = GkqPathFlow.from_scf_input(workdir, gs_inp, ngqpt, qbounds, ndivsm=2, with_becs=True, ddk_tolerance=None,
                                          test_ft_interpolation=True, prepgkk=0)

        #assert len(work.relax_tasks) == 3
        #assert all(t.input["dilatmx"] == 1.05 for t in work)

        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

"""Tests for qha module."""

from abipy.core.testing import AbipyTest
from abipy.flowtk.vzsisa import VzsisaFlow


class TestVzsisa(AbipyTest):

    def test_vzsisa_flow(self):
        """Testing VzsisaFlow."""
        workdir = self.mkdtemp()
        ngkpt = [4, 4, 4]
        ngqpt = [2, 2, 2]
        scf_input = self.get_gsinput_alas_ngkpt(ngkpt=ngkpt)

        bo_scales = [0.96, 0.98, 1, 1.02, 1.04]    # EinfVib4(S)
        ph_scales = [1, 1.02, 1.04]                # EinfVib2(D)

        with_becs = False
        with_quad = False
        #with_quad = not structure.has_zero_dynamical_quadrupoles

        flow = VzsisaFlow.from_scf_input(workdir, scf_input, bo_scales, ph_scales, ngqpt,
                                         with_becs, with_quad, edos_ngkpt=None)

        flow.allocate()
        flow.check_status()
        assert len(flow) == 1
        isok, checks = flow.abivalidate_inputs()
        assert isok

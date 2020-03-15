"""Tests for qha module."""

from abipy.core.testing import AbipyTest
from abipy.flowtk.qha import QhaFlow


class TestQha(AbipyTest):

    def test_qhaflow(self):
        """Testing QhaFlow."""
        workdir = self.mkdtemp()
        scf_input = self.get_gsinput_alas_ngkpt(ngkpt=[4, 4, 4])
        v0 = scf_input.structure.volume
        volumes = [0.98 * v0, v0, v0 * 1.02]
        flow = QhaFlow.from_scf_input(workdir, scf_input, volumes, ngqpt=[2, 2, 2],
                                      with_becs=False, edos_ngkpt=[4, 4, 4], metadata={"mpi_id": 123})

        flow.allocate()
        flow.check_status()
        assert len(flow) == 1
        isok, checks = flow.abivalidate_inputs()
        assert isok

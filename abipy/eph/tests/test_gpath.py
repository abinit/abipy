"""Tests for varpeq module."""
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.eph.gpath import GpathFile, GpathRobot


class GpathTest(AbipyTest):

    def test_gpath_file_fixed_k(self):
        """Testing GPATH.nc file with fixed k."""
        with GpathFile(abidata.ref_file("teph4zpr_9o_DS1_GPATH.nc")) as gpath:
            # g(k,q) at fixed k along q-path.
            repr(gpath)
            str(gpath)
            assert gpath.to_string(verbose=2)
            assert not gpath.params
            assert gpath.structure.formula == "Mg1 O1" and len(gpath.structure) == 2
            assert gpath.r.eph_fix_korq == "k"
            assert gpath.r.bstart == 5 and gpath.r.bstop == 8
            assert gpath.r.nk_path == 1
            assert gpath.r.nq_path == 42
            assert len(gpath.ebands_k.kpoints) == 1
            assert len(gpath.ebands_kq.kpoints) == 42
            assert len(gpath.phbands.qpoints) == 42
            self.assert_equal(gpath.r.eph_fix_wavec, (0, 0, 0))

            if self.has_matplotlib():
                assert gpath.plot_g_qpath(show=False)
                with self.assertRaises(ValueError):
                    gpath.plot_g_kpath(show=False)

            # Test jupyter notebook creation
            if self.has_nbformat():
                gpath.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_gpath_file_fixed_q(self):
        """Testing GPATH.nc file with fixed q."""
        with GpathFile(abidata.ref_file("teph4zpr_9o_DS2_GPATH.nc")) as gpath:
            # g(k,q) at fixed q along k-path.
            repr(gpath)
            str(gpath)
            assert gpath.to_string(verbose=2)
            assert not gpath.params
            assert gpath.structure.formula == "Mg1 O1" and len(gpath.structure) == 2
            assert gpath.r.eph_fix_korq == "q"
            assert gpath.r.bstart == 5 and gpath.r.bstop == 8
            assert gpath.r.nk_path == 42
            assert gpath.r.nq_path == 1
            assert len(gpath.ebands_k.kpoints) == 42
            #assert len(gpath.ebands_kq.kpoints) == 1
            assert len(gpath.phbands.qpoints) == 1
            self.assert_equal(gpath.r.eph_fix_wavec, (0.11, 0, 0))

            if self.has_matplotlib():
                assert gpath.plot_g_kpath(show=False)
                with self.assertRaises(ValueError):
                    gpath.plot_g_qpath(show=False)

            # Test jupyter notebook creation
            if self.has_nbformat():
                gpath.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_gpath_robot(self):
        """Testing GpathRobot."""
        files = abidata.ref_files(
            "abinitio_qpath_V1QAVG.nc",
            "interpolated_qpath_V1QAVG.nc",
        )

        with GpathRobot() as robot:
            robot.add_file("one", abidata.ref_file("teph4zpr_9o_DS1_GPATH.nc"))
            robot.add_file("two", abidata.ref_file("teph4zpr_9o_DS1_GPATH.nc"))
            assert len(robot) == 2
            repr(robot); str(robot)
            robot.to_string(verbose=2)

            # Test matplotlib methods
            if self.has_matplotlib():
                assert robot.plot_g_qpath(show=False)

            # Test jupyter notebook creation
            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))

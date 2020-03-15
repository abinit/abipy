"""Tests for v1qavg module."""
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.eph.v1qavg import V1qAvgFile, V1qAvgRobot


class V1qavgTest(AbipyTest):

    def test_v1qavg_file(self):
        """Testing V1qAvgFile."""
        with V1qAvgFile(abidata.ref_file("abinitio_qpath_V1QAVG.nc")) as ncfile:
            repr(ncfile); str(ncfile)
            assert ncfile.to_string(verbose=2)
            assert not ncfile.params
            assert ncfile.structure.formula == "Ga1 P1" and len(ncfile.structure) == 2

            assert ncfile.has_zeff
            assert ncfile.has_dielt
            assert not ncfile.has_quadrupoles
            assert not ncfile.has_efield
            assert ncfile.dvdb_add_lr == 1
            assert ncfile.symv1scf == 0
            assert not ncfile.interpolated
            assert not ncfile.has_maxw
            assert ncfile.qdamp == -1.0

            assert len(ncfile.qpoints) == 12

            # Test matplotlib methods
            if self.has_matplotlib():
                assert ncfile.plot(what_list="all", ispden=0, show=False)
                assert ncfile.plot_maxw(show=False) is None
                assert ncfile.plot_maxw_perts(scale="semilogy", sharey=False, fontsize=8, show=False) is None

            # Test jupyter notebook creation
            if self.has_nbformat():
                ncfile.write_notebook(nbpath=self.get_tmpname(text=True))


class V1qAvgRobotTest(AbipyTest):

    def test_v1qavg_robot(self):
        """Testing V1qAvgRobot."""
        files = abidata.ref_files(
            "abinitio_qpath_V1QAVG.nc",
            "interpolated_qpath_V1QAVG.nc",
        )

        with V1qAvgRobot.from_files(files[0]) as robot:
            robot.add_file("interpolated_v1qavg", files[1])
            assert len(robot) == 2
            repr(robot); str(robot)
            robot.to_string(verbose=2)

            interp_ncfile = robot.abifiles[1]
            assert interp_ncfile.interpolated
            assert interp_ncfile.has_maxw

            # Test matplotlib methods
            if self.has_matplotlib():
                assert robot.plot(ispden=0, vname="v1scf_avg", show=False)
                assert robot.plot_maxw(show=False) is None

                assert interp_ncfile.plot_maxw(show=False)
                assert interp_ncfile.plot_maxw_perts(scale="semilogy", sharey=False, fontsize=8, show=False)

            # Test jupyter notebook creation
            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))

"""Tests for varpeq module."""
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.eph.varpeq import VarpeqFile


#class VarpeqTest(AbipyTest):
#
#    def test_varpeq_file(self):
#        """Testing VarpeqFile."""
#        with VarpeqFile(abidata.ref_file("abinitio_qpath_V1QAVG.nc")) as varpeq:
#            repr(varpeq)
#            str(varpeq)
#            assert varpeq.to_string(verbose=2)
#            assert not varpeq.params
#            assert varpeq.structure.formula == "Ga1 P1" and len(varpeq.structure) == 2
#
#            #print(varpeq.ebands.kpoints.ksampling)
#            for polaron in varpeq.polaron_spin:
#                #print(polaron)
#                #polaron.plot_bz_sampling()
#                #polaron.plot_bqnu_with_phbands(phbands)
#                polaron.plot_ank_with_ebands(ebands)
#
#                if self.has_matplotlib():
#
#            # Test jupyter notebook creation
#            if self.has_nbformat():
#                varpeq.write_notebook(nbpath=self.get_tmpname(text=True))


#class VarpeqRobotTest(AbipyTest):
#
#    def test_varpeq_robot(self):
#        """Testing VarpeqRobot."""
#        files = abidata.ref_files(
#            "abinitio_qpath_V1QAVG.nc",
#            "interpolated_qpath_V1QAVG.nc",
#        )
#
#        with VarpeqRobot.from_files(files[0]) as robot:
#            robot.add_file("interpolated_v1qavg", files[1])
#            assert len(robot) == 2
#            repr(robot); str(robot)
#            robot.to_string(verbose=2)
#
#            interp_ncfile = robot.abifiles[1]
#            assert interp_ncfile.interpolated
#            assert interp_ncfile.has_maxw
#
#            # Test matplotlib methods
#            if self.has_matplotlib():
#                assert robot.plot(ispden=0, vname="v1scf_avg", show=False)
#                assert robot.plot_maxw(show=False) is None
#                assert interp_ncfile.plot_maxw(show=False)
#                assert interp_ncfile.plot_maxw_perts(scale="semilogy", sharey=False, fontsize=8, show=False)
#
#            # Test jupyter notebook creation
#            if self.has_nbformat():
#                robot.write_notebook(nbpath=self.get_tmpname(text=True))

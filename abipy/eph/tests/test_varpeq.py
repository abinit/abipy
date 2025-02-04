"""Tests for varpeq module."""
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.eph.varpeq import VarpeqFile


class VarpeqTest(AbipyTest):

    def test_varpeq_file(self):
        """Testing VarpeqFile."""
        varpeq_filepath = "/Users/giantomassi/git_repos/abinit/_build/tests/POLARON/varpeq6/out_VARPEQ.nc"
        with VarpeqFile(varpeq_filepath) as varpeq:
            repr(varpeq)
            str(varpeq)
            assert varpeq.to_string(verbose=2)
            assert varpeq.structure.formula == "Li1 F1" and len(varpeq.structure) == 2
            params = varpeq.params
            #assert params["nkbz"] ==
            #assert params["ngkpt"] ==

            #print(varpeq.ebands.kpoints.ksampling)
            for polaron in varpeq.polaron_spin:
                print(polaron)
                #assert polaron.spin == 0
                #assert polaron.nstates == 0
                #assert polaron.nb == 0
                #assert polaron.nk == 0
                #assert polaron.nq == 0
                #assert polaron.bstart == 0
                #assert polaron.bstop == 0
                df = polaron.get_final_results_df(with_params=True)
                #polaron.write_a2_bxsf(self, filepath: PathLike, fill_value: float = 0.0) -> None:
                #polaron.write_b2_bxsf(self, filepath: PathLike, fill_value: float = 0.0) -> None:
                if self.has_matplotlib():
                    polaron.plot_scf_cycle(show=False)
                    #polaron.plot_ank_with_ebands(ebands_kpath, ebands_kmesh=None)
                    #polaron.plot_bqnu_with_ddb("in_DDB", with_phdos=True)
                    #polaron.plot_bqnu_with_phbands(phbands_qpath)

            # Test jupyter notebook creation
            #if self.has_nbformat():
            #    varpeq.write_notebook(nbpath=self.get_tmpname(text=True))


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

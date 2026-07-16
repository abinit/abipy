"""Tests for varpeq module."""

import os

import numpy as np
import pytest

from abipy.core.testing import AbipyTest
from abipy.eph.vpq import VpqFile, VpqReader

root = "/Users/giantomassi/git_repos/abinit/_build/tests/tutorespfn_teph4vpq_1-teph4vpq_2-teph4vpq_3-teph4vpq_4-teph4vpq_5-teph4vpq_6-teph4vpq_7-teph4vpq_8-teph4vpq_9-teph4vpq_10"


class ArrayVariable:
    """Minimal array-backed stand-in for a NetCDF variable."""

    def __init__(self, values, dimensions):
        self.values = values
        self.dimensions = dimensions

    def __getitem__(self, index):
        return self.values[index]


class ArrayReader:
    """Minimal reader exposing one array-backed variable."""

    def __init__(self, variable):
        self.variable = variable

    def read_variable(self, varname):
        return self.variable


def test_read_spin_scf_variable_with_optional_hop_dimension():
    """The new hop axis is collapsed while legacy VPQ array shapes are preserved."""
    spin = 1
    variables = {
        "cvflag_spin": ("nsppol", "nstates"),
        "nstep2cv_spin": ("nsppol", "nstates"),
        "scf_hist_spin": ("nsppol", "nstates", "nstep", "six"),
    }

    for dimensions in variables.values():
        shape = tuple(range(2, 2 + len(dimensions)))
        legacy_values = np.arange(np.prod(shape)).reshape(shape)
        legacy_variable = ArrayVariable(legacy_values, dimensions)
        actual = VpqReader.read_spin_scf_variable(ArrayReader(legacy_variable), "unused", spin)
        np.testing.assert_array_equal(actual, legacy_values[spin])

        hop_dimensions = ("nsppol", "hop_nstep", *dimensions[1:])
        hop_shape = (shape[0], 3, *shape[1:])
        hop_values = np.arange(np.prod(hop_shape)).reshape(hop_shape)
        hop_variable = ArrayVariable(hop_values, hop_dimensions)
        actual = VpqReader.read_spin_scf_variable(ArrayReader(hop_variable), "unused", spin)
        np.testing.assert_array_equal(actual, hop_values[spin, -1])


class VarpeqTest(AbipyTest):
    @pytest.mark.xfail(condition=not os.path.exists(root), reason=f"{root=} does not exist")
    def test_varpeq_file(self):
        """Testing VpqFile."""
        filepath = os.path.join(root, "teph4vpq_9o_VPQ.nc")
        with VpqFile(filepath) as vpq:
            repr(vpq)
            str(vpq)
            assert vpq.to_string(verbose=2)
            assert vpq.structure.formula == "Li1 F1" and len(vpq.structure) == 2
            params = vpq.params
            assert params["avg_g"]
            assert params["e_frohl"] == -0.21380923340128977

            # print(vpq.ebands.kpoints.ksampling)
            for polaron in vpq.polaron_spin:
                print(polaron)
                # assert polaron.spin == 0
                # assert polaron.nstates == 0
                # assert polaron.nb == 0
                # assert polaron.nk == 0
                # assert polaron.nq == 0
                # assert polaron.bstart == 0
                # assert polaron.bstop == 0
                df = polaron.get_final_results_df(with_params=True)
                print(df)

                assert df["spgroup"][0] == 225
                self.assert_almost_equal(df["e_frohl"].values, -0.2138092)
                self.assert_equal(df["filter_value"].values, 1)

                # output of bxsf files
                tmp_filepath = self.tmpfileindir("foo.xsf")
                polaron.write_a2_bxsf(tmp_filepath)
                polaron.write_b2_bxsf(tmp_filepath)

                if self.has_matplotlib():
                    polaron.plot_scf_cycle(show=False)
                    # polaron.plot_ank_with_ebands(ebands_kpath, ebands_kmesh=None)
                    # polaron.plot_bqnu_with_ddb("in_DDB", with_phdos=True)
                    # polaron.plot_bqnu_with_phbands(phbands_qpath)

            # Test jupyter notebook creation
            # if self.has_nbformat():
            #    vpq.write_notebook(nbpath=self.get_tmpname(text=True))


# class VarpeqRobotTest(AbipyTest):
#
#    def test_varpeq_robot(self):
#        """Testing VarpeqRobot."""
#        filepaths = abidata.ref_files(
#            "abinitio_qpath_V1QAVG.nc",
#            "interpolated_qpath_V1QAVG.nc",
#        )
#
#        with VarpeqRobot.from_files(filepaths[0]) as robot:
#            robot.add_file("interpolated_v1qavg", filepaths[1])
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

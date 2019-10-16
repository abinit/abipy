# coding: utf-8
"""Tests for optic module."""
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab


class OpticTest(AbipyTest):

    def test_opticfile(self):
        """Testing OpticNcFile"""
        with abilab.abiopen(abidata.ref_file("gaas_121212_OPTIC.nc")) as optic:
            repr(optic); str(optic)
            assert optic.to_string(verbose=2)
            assert optic.structure.formula == "Ga1 As1"
            # TODO: kptopt is not propagated in optic.
            #assert optic.kptopt == 2
            self.assert_almost_equal(optic.broadening, 0.002)
            self.assert_almost_equal(optic.domega, 0.0003)
            self.assert_almost_equal(optic.maxomega, 0.3)
            self.assert_almost_equal(optic.scissor, 0.000)
            self.assert_almost_equal(optic.tolerance, 0.002)
            assert optic.reader.ntemp == 1
            assert optic.params["nspden"] == 1

            assert optic.has_linopt
            assert optic.reader.computed_components["linopt"] == ["xx", "zz"]

            assert optic.has_shg
            assert optic.reader.computed_components["shg"] == ["xyz", "yyy"]

            assert optic.has_leo
            assert optic.reader.computed_components["leo"] == ["xyz"]
            #assert not optic.reader.computed_components["leo2"]

            # Test plot methods
            if self.has_matplotlib():
                assert optic.plot_linear_epsilon(show=False)
                assert optic.plot_linopt(show=False)
                assert optic.plot_shg(show=False)
                assert optic.plot_leo(show=False)

            if self.has_nbformat():
                optic.write_notebook(nbpath=self.get_tmpname(text=True))


class OpticRobotTest(AbipyTest):

    def test_optic_robot(self):
        """Test OpticRobot."""
        files = abidata.ref_files("gaas_444_OPTIC.nc", "gaas_888_OPTIC.nc", "gaas_121212_OPTIC.nc")
        with abilab.OpticRobot.from_files(files) as robot:
            repr(robot); str(robot)
            robot.to_string(verbose=2)
            assert robot.computed_components_intersection["linopt"] == ["xx", "zz"]
            assert robot.computed_components_intersection["shg"] == ["xyz", "yyy"]
            assert robot.computed_components_intersection["leo"] == ["xyz"]
            assert [t[2] for t in robot.sortby("nkpt")] == [10, 60, 182]

            df_params = robot.get_params_dataframe()
            self.assert_equal(df_params["nspden"].values, 1)

            # Test plot methods
            if self.has_matplotlib():
                assert robot.plot_linopt_convergence(show=False)
                assert robot.plot_shg_convergence(show=False)
                assert robot.plot_leo_convergence(show=False)
                assert robot.plot_lattice_convergence(sortby="nkpt", hue="nspden", show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))

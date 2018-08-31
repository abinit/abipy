""""Tests for HIST.nc files."""
from __future__ import division, print_function, unicode_literals

import abipy.data as abidata
from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.dynamics.hist import HistFile, HistRobot
import abipy.core.abinit_units as abu


class HistFileTest(AbipyTest):

    def test_hist_api(self):
        """Testing HistFile API."""
        hist = HistFile(abidata.ref_file("sic_relax_HIST.nc"))
        repr(hist); str(hist)
        hist.to_string(verbose=2)
        assert not hist.params

        an = hist.get_relaxation_analyzer()
        assert hist.num_steps == 7
        assert len(hist.structures) == hist.num_steps
        assert hist.initial_structure is hist.structures[0]
        assert hist.final_structure is hist.structures[-1]
        assert hist.final_structure.composition.reduced_formula == "SiC"
        assert len(hist.final_structure) == hist.reader.natom
        assert len(hist.final_structure) == 2
        assert len(hist.etotals) == hist.num_steps
        self.assert_almost_equal(hist.etotals.to("Ha"), [-10.4914795629442, -10.491527362795, -10.4915307068041,
    -10.4915319277052, -10.4915319344634, -10.491531918083, -10.4915319353756])

        last_stren = [5.01170783044364e-08, 5.01170783044364e-08, 5.01170783046533e-08, 0, 0, 0]
        self.assert_almost_equal(hist.reader.read_value("strten")[-1], last_stren)

        cart_forces_step = hist.reader.read_cart_forces(unit="Ha bohr^-1")
        #self.assert_almost_equal(cart_forces_step[0], [
        #    6.42133418983323e-32, -1.92640025694997e-31, 6.42133418983323e-32,
        #   -6.42133418983323e-32, 1.92640025694997e-31, -6.42133418983323e-32])
        fred = hist.reader.read_reduced_forces()
        assert fred.shape == cart_forces_step.shape

        # Cartesian components of stress tensor (hartree/bohr^3)
        #  sigma(1 1)=  5.01170783E-08  sigma(3 2)=  0.00000000E+00
        #  sigma(2 2)=  5.01170783E-08  sigma(3 1)=  0.00000000E+00
        #  sigma(3 3)=  5.01170783E-08  sigma(2 1)=  0.00000000E+00

        #-Cartesian components of stress tensor (GPa)         [Pressure= -1.4745E-03 GPa]
        #- sigma(1 1)=  1.47449510E-03  sigma(3 2)=  0.00000000E+00
        #- sigma(2 2)=  1.47449510E-03  sigma(3 1)=  0.00000000E+00
        #- sigma(3 3)=  1.47449510E-03  sigma(2 1)=  0.00000000E+00

        cart_stress_tensors, pressures = hist.reader.read_cart_stress_tensors()
        self.assert_almost_equal(pressures[-1], -1.4745E-03)
        self.assert_almost_equal(cart_stress_tensors[-1, 1, 0], 0.0)
        for i in range(3):
            self.assert_almost_equal(cart_stress_tensors[-1, i, i], 5.01170783E-08 * abu.HaBohr3_GPa)

        same_structure = abilab.Structure.from_file(abidata.ref_file("sic_relax_HIST.nc"))
        self.assert_almost_equal(same_structure.frac_coords, hist.final_structure.frac_coords)

        # Test to_xdatcar converter
        xdatcar = hist.to_xdatcar(filepath=None, groupby_type=True)
        assert xdatcar.natoms == [1, 1] and len(xdatcar.structures) == hist.num_steps
        self.assert_almost_equal(xdatcar.structures[0].volume, hist.structures[0].volume)
        self.assert_almost_equal(xdatcar.structures[-1].frac_coords, hist.structures[-1].frac_coords)

        xdatcar_nogroup = hist.to_xdatcar(filepath=None, groupby_type=False)
        assert xdatcar.structures[0] ==  xdatcar_nogroup.structures[0]
        assert xdatcar.structures[-1] ==  xdatcar_nogroup.structures[-1]

        # Test matplotlib plots.
        if self.has_matplotlib():
            assert hist.plot(show=False)
            assert hist.plot_energies(show=False)

        # Test notebook generation.
        if self.has_nbformat():
            hist.write_notebook(nbpath=self.get_tmpname(text=True))

        if self.has_mayavi():
            assert hist.mvplot_trajectories(show=False)
            #assert hist.mvanimate(delay=100)

        hist.close()

    def test_hist_robot(self):
        """Test HistRobot."""
        filepath = abidata.ref_file("sic_relax_HIST.nc")
        with HistRobot.from_files(filepath) as robot:
            robot.add_file("same hist", filepath)
            repr(robot); str(robot)
            assert robot.to_string(verbose=2)

            # From base class
            assert len(robot.get_structure_dataframes()) == 3

            df = robot.get_dataframe()
            assert "alpha" in df

            if self.has_matplotlib():
                what_list = ["energy", "abc", "pressure", "forces"]
                assert robot.gridplot(what=what_list, fontsize=4, show=False)
                assert robot.combiplot(colormap="viridis", show=False)
                assert robot.plot_lattice_convergence(fontsize=10, show=False)
                assert robot.plot_lattice_convergence(what_list=("a", "alpha"), show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))

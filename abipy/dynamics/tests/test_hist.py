""""Tests for HIST.nc files."""
from __future__ import division, print_function, unicode_literals

import abipy.data as abidata
from abipy.core.testing import AbipyTest
from abipy.dynamics.hist import HistFile


class HistFileTest(AbipyTest):

    def test_hist_api(self):
        """Testing HistFile API."""
        hist = HistFile(abidata.ref_file("sic_relax_HIST.nc"))
        repr(hist)
        str(hist)

        assert hist.num_steps == 7
        assert len(hist.structures) == hist.num_steps
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
            self.assert_almost_equal(cart_stress_tensors[-1, i, i], 5.01170783E-08)

        # Test matplotlib plots.
        if self.has_matplotlib():
            hist.plot(show=False)
            hist.plot_energies(show=False)

        # Test notebook generation.
        if self.has_nbformat():
            hist.write_notebook(nbpath=self.get_tmpname(text=True))

        hist.close()

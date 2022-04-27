"""Unit tests for oncvpsp"""

import sys
import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.ppcodes.oncvpsp import OncvOutputParser, psp8_get_densities


def filepath(basename):
    return os.path.join(os.path.dirname(__file__), basename)


class OncvOutputParserTest(AbipyTest):

    def test_nonrelativistic_oxygen_v2(self):
        """
        Parsing the non-relativistic output file produced by ONCVPSPS v2
        """
        # Non-relativistic results
        p = OncvOutputParser(filepath("08_O_nr.out"))
        repr(p); str(p)

        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert p.calc_type == "non-relativistic"
        assert not p.fully_relativistic
        assert p.version == "2.1.1"
        assert p.major_version == 2
        assert p.minor_version == 1
        assert p.patch_level == 1

        assert p.atsym == "O"
        assert p.z == "8.00"
        assert p.iexc == "3"
        assert p.lmax == 1
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 1

        rhov, rhoc, rhom = p.densities["rhoV"], p.densities["rhoC"], p.densities["rhoM"]
        assert rhov.rmesh[0] == 0.0100642
        assert rhov.rmesh[-1] == 3.9647436
        assert rhoc.values[0] == 53.3293576
        assert all(rhom.values == 0.0)

        # Build the plotter
        plotter = p.make_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

        #if self.has_nbformat():

    def test_scalar_relativistic_oxygen_v2(self):
        """
        Parsing the scalar-relativistic output file produced by ONCVPSPS v2
        """
        # Scalar relativistic output
        p = OncvOutputParser(filepath("08_O_sr.out"))
        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert not p.fully_relativistic
        assert p.calc_type == "scalar-relativistic"
        assert p.version == "2.1.1"

        assert p.atsym == "O"
        assert p.z == "8.00"
        assert p.iexc == "3"
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 1

        # Test potentials
        vloc = p.potentials[-1]
        pl0 = {0: -7.4449470, 1: -14.6551019, -1: -9.5661177}

        for l, pot in p.potentials.items():
            assert (pot.rmesh[0], pot.rmesh[-1]) == (0.0099448, 3.9647436)
            str(l)
            assert pot.values[0] == pl0[l]
            assert all(pot.rmesh == vloc.rmesh)

        # Test wavefunctions
        ae_wfs, ps_wfs = p.radial_wfs.ae, p.radial_wfs.ps

        nlk = (1, 0, None)
        ae10, ps10 = ae_wfs[nlk], ps_wfs[nlk]
        assert ae10[0] == (0.009945, -0.092997)
        assert ps10[0] == (0.009945,  0.015273)
        assert ae10[-1] == (3.964744, 0.037697)
        assert ps10[-1] == (3.964744, 0.037694)

        nlk = (2, 1, None)
        ae21, ps21 = ae_wfs[nlk], ps_wfs[nlk]
        assert ae21[0] == (0.009945, 0.001463)
        assert ps21[0] == (0.009945, 0.000396)

        # Test projectors
        prjs = p.projectors
        assert prjs[(1, 0, None)][0] == (0.009945, 0.015274)
        assert prjs[(2, 0, None)][0] == (0.009945, -0.009284)
        assert prjs[(1, 0, None)][-1] == (3.964744, 0.037697)
        assert prjs[(2, 0, None)][-1] == (3.964744, 0.330625)

        assert prjs[(1, 1, None)][0] == (0.009945, 0.000395)
        assert prjs[(2, 1, None)][0] == (0.009945, -0.000282)

        # Test convergence data
        c = p.ene_vs_ecut
        assert c[0].energies[0] == 5.019345
        assert c[0].values[0] == 0.010000
        assert c[0].energies[-1] == 25.317286
        assert c[0].values[-1] == 0.000010
        assert c[1].energies[0] == 19.469226
        assert c[1].values[0] == 0.010000

        # Test log derivatives
        ae0, ps0 = p.atan_logders.ae[0], p.atan_logders.ps[0]
        assert (ae0.energies[0], ae0.values[0]) == (2.000000, 0.706765)
        assert (ps0.energies[0], ps0.values[0]) == (2.000000, 0.703758)
        assert ae0.energies[-1] == -2.000000
        assert ps0.energies[-1] == -2.000000

        ae1, ps1 = p.atan_logders.ae[1], p.atan_logders.ps[1]
        assert (ae1.energies[0], ae1.values[0]) == (2.000000, -2.523018)
        assert ps1.values[0] == -2.521334

        # Build the plotter
        plotter = p.make_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

    def test_fully_relativistic_oxygen_v2(self):
        """
        Parsing the fully-relativistic output file produced by ONCVPSPS v2
        """
        p = OncvOutputParser(filepath("08_O_r.out"))

        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert p.fully_relativistic
        assert p.calc_type == "fully-relativistic"
        assert p.version == "2.1.1"

        assert p.atsym == "O"
        assert p.z == "8.00"
        assert p.iexc == "3"
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 1

        # TODO: Wavefunctions

        # Build the plotter
        plotter = p.make_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

    def test_scalar_relativistic_oxygen_v4(self):
        """
        Parsing the scalar-relativistic output file produced by ONCVPSPS v4
        """
        # Scalar relativistic output
        p = OncvOutputParser(filepath("O_sr_v4.out"))
        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert not p.fully_relativistic
        assert p.calc_type == "scalar-relativistic"
        assert p.version == "4.0.1"

        assert p.atsym == "O"
        assert p.z == "8.00"
        assert p.iexc == "-106131"
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 2

        # Test potentials
        vloc = p.potentials[-1]
        pl0 = {0: -4.6445128, 1: -15.3234007, 2: 20.9698547, -1: -10.0124145}

        for l, pot in p.potentials.items():
            assert (pot.rmesh[0], pot.rmesh[-1]) == (0.0099582, 2.4161798)
            str(l)
            assert pot.values[0] == pl0[l]
            assert all(pot.rmesh == vloc.rmesh)

        # Test wavefunctions
        ae_wfs, ps_wfs = p.radial_wfs.ae, p.radial_wfs.ps

        nlk = (2, 0, None)
        ae20, ps20 = ae_wfs[nlk], ps_wfs[nlk]
        assert ae20[0] == (0.009958, -0.093703)
        assert ps20[0] == (0.009958,  0.013614)
        assert ae20[-1] == (5.998219, 0.002734)
        assert ps20[-1] == (5.998219, 0.002734)

        nlk = (2, 1, None)
        ae21, ps21 = ae_wfs[nlk], ps_wfs[nlk]
        assert ae21[0] == (0.009958, 0.001474)
        assert ps21[0] == (0.009958, 0.000456)

        # Test projectors
        prjs = p.projectors
        assert prjs[(1, 0, None)][0] == (0.009958, 0.090486)
        assert prjs[(2, 0, None)][0] == (0.009958, -0.025921)
        assert prjs[(1, 0, None)][-1] == (1.580056, -0.000000)
        assert prjs[(2, 0, None)][-1] == (1.580056, -0.000000)

        assert prjs[(1, 1, None)][0] == (0.009958, 0.002057)
        assert prjs[(2, 1, None)][0] == (0.009958, -0.000854)

        # Test convergence data
        c = p.ene_vs_ecut
        assert c[0].energies[0] == 12.172858
        assert c[0].values[0] == 0.010000
        assert c[0].energies[-1] == 32.408732
        assert c[0].values[-1] == 0.000010
        assert c[1].energies[0] == 23.772439
        assert c[1].values[0] == 0.010000

        # Test log derivatives
        ae0, ps0 = p.atan_logders.ae[0], p.atan_logders.ps[0]
        assert (ae0.energies[0], ae0.values[0]) == (12.000000, -1.601062)
        assert (ps0.energies[0], ps0.values[0]) == (12.000000, -1.611618)

        assert (ae0.energies[-1], ae0.values[-1]) == (-11.980000, 4.528242)
        assert (ps0.energies[-1], ps0.values[-1]) == (-11.980000, 4.528290)

        ae1, ps1 = p.atan_logders.ae[1], p.atan_logders.ps[1]
        assert (ae1.energies[0], ae1.values[0]) == (12.000000, 1.066272)
        assert ps1.values[0] == 1.286771

        # Build the plotter
        plotter = p.make_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

    def _call_plotter_methods(self, plotter):
        if self.has_matplotlib():
            assert plotter.plot_atan_logders(show=False)
            assert plotter.plot_radial_wfs(show=False)
            assert plotter.plot_projectors(show=False)
            assert plotter.plot_densities(show=False)
            assert plotter.plot_der_densities(order=1, show=False)
            assert plotter.plot_potentials(show=False)
            assert plotter.plot_der_potentials(order=1, show=False)
            assert plotter.plot_ene_vs_ecut(show=False)
            assert plotter.plot_atanlogder_econv(show=False)
            assert plotter.plot_dens_and_pots(show=False)
            assert plotter.plot_waves_and_projs(show=False)
            assert plotter.plot_den_formfact(ecut=40, show=False)

        #if self.has_plotly():

    def test_psp8_get_densities(self):
        """
        Testing get_densities from psp8 format.
        """
        n = psp8_get_densities(filepath("Lu-sp.psp8"),
                               fc_file=sys.stdout,
                               ae_file=sys.stdout, plot=False)
        assert len(n.rmesh) == 600
        self.assert_almost_equal([n.rmesh[0], n.psval[0], n.aeval[0], n.aecore[0]],
            np.fromstring("0.0000000000000E+00  5.9161585718320E-02  3.9966212837901E+03  3.9211427139394E+06", sep=" "))
        self.assert_almost_equal([n.rmesh[-1], n.psval[-1], n.aeval[-1], n.aecore[-1]],
            np.fromstring("5.9900000000000E+00  4.1502186898788E-03  4.1500999839310E-03  4.6618962684673E-06", sep=" "))

        with self.assertRaises(RuntimeError):
            path = os.path.join(abidata.pseudo_dir, "O.psp8")
            psp8_get_densities(path, plot=False)

        # TODO: SOC

"""Unit tests for oncvpsp"""

import sys
import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.core.atom import NlkState
from abipy.ppcodes.oncv_parser import OncvParser
from abipy.ppcodes.oncv_plotter import MultiOncvPlotter, psp8_get_densities


def filepath(basename):
    return os.path.join(abidata.dirpath, "oncv_data", basename)


class OncvOutputParserTest(AbipyTest):

    def test_nonrelativistic_oxygen_v2(self):
        """
        Parsing the non-relativistic output file produced by ONCVPSPS v2
        """
        # Non-relativistic results
        p = OncvParser(filepath("08_O_nr.out"))
        repr(p); str(p)

        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert p.calc_type == "non-relativistic"
        assert not p.relativistic
        assert p.version == "2.1.1"
        assert p.major_version == 2
        assert p.minor_version == 1
        assert p.patch_level == 1

        assert p.atsym == "O"
        assert p.z == 8.00
        assert p.iexc == 3
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
        plotter = p.get_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

        #if self.has_nbformat():

    def test_scalar_relativistic_oxygen_v2(self):
        """
        Parsing the scalar-relativistic output file produced by ONCVPSPS v2
        """
        # Scalar relativistic output
        p = OncvParser(filepath("08_O_sr.out"))
        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert not p.relativistic
        assert p.calc_type == "scalar-relativistic"
        assert p.version == "2.1.1"

        assert p.atsym == "O"
        assert p.z == 8.00
        assert p.iexc == 3
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 1

        #assert p.hints["low"]["ecut"] ==
        #assert p.hints["normal"]["ecut"] ==
        #assert p.hints["high"]["ecut"] ==

        #results = p.get_results()

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
        c = p.kene_vs_ecut
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
        plotter = p.get_plotter()
        repr(plotter); str(plotter)
        assert plotter is not None
        self._call_plotter_methods(plotter)

    def test_relativistic_oxygen_v2(self):
        """
        Parsing the fully-relativistic output file produced by ONCVPSPS v2
        """
        p = OncvParser(filepath("08_O_r.out"))

        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert p.relativistic
        assert p.calc_type == "fully-relativistic"
        assert p.version == "2.1.1"

        assert p.atsym == "O"
        assert p.z == 8.00
        assert p.iexc == 3
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 1
        assert p.rc5 == 1.40
        assert p.rc_l[0] == 1.60
        assert p.rc_l[1] == 1.60

        # TODO: Wavefunctions

        # Build the plotter
        plotter = p.get_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

    def test_scalar_relativistic_oxygen_v4(self):
        """
        Parsing the scalar-relativistic output file produced by ONCVPSPS v4
        """
        # Scalar relativistic output
        p = OncvParser(filepath("O_sr_v4.out"))
        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed

        assert not p.relativistic
        assert p.calc_type == "scalar-relativistic"
        assert p.version == "4.0.1"

        assert p.atsym == "O"
        assert p.z == 8.00
        assert p.iexc == -106131
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 2
        assert p.rc5 == 1.2
        assert p.rc_l[0] == 1.35000
        assert p.rc_l[1] == 1.45000
        assert p.rc_l[2] == 1.25000

        assert p.get_input_str()
        assert p.get_psp8_str()
        assert p.get_upf_str()

        # Calculating optimized projector #   1
        # for l=   0
        nlk = NlkState(n=1, l=0, k=None)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [12.17, 21.67, 27.75, 32.41])

        #Calculating optimized projector #   2
        # for l=   1
        nlk = NlkState(n=2, l=1, k=None)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [22.49, 30.35, 36.39, 41.55])

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
        c = p.kene_vs_ecut
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
        plotter = p.get_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

    def test_relativistic_oxygen_v4(self):
        """
        Parsing the relativistic output file produced by ONCVPSPS v4
        """
        # Fully relativistic output
        p = OncvParser(filepath("O_fr_v4.out"))
        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed
        assert p.is_oncvpsp and not p.is_metapsp

        assert p.relativistic
        assert p.calc_type == "relativistic"
        assert p.version == "4.0.1"

        assert p.atsym == "O"
        assert p.z == 8.00
        assert p.iexc == 3
        assert p.nc == 1
        assert p.nv == 2
        assert p.lmax == 1
        assert p.rc5 == 1.4
        assert p.rc_l[0] == 1.6
        assert p.rc_l[1] == 1.6

        assert p.get_input_str()
        assert p.get_psp8_str()
        assert p.get_upf_str()

        # Calculating optimized projector #   1
        # for l=   0
        nlk = NlkState(n=1, l=0, k=1)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [5.01, 14.64, 21.05, 25.33])

        # Calculating optimized projector #   2
        # for l=   0
        nlk = NlkState(n=2, l=0, k=1)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [4.55, 15.25, 22.10, 26.79])

        #Calculating optimized projector #   1
        # for l=   1
        nlk = NlkState(n=1, l=1, k=1)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [19.49, 24.68, 28.68, 35.11])

        #Calculating optimized projector #   2
        # for l=   1
        nlk = NlkState(n=2, l=1, k=1)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [19.17, 25.12, 29.75, 41.34])

        # Check values associated to Fortran ikap = 2

        #Calculating optimized projector #   1
        # for l=   1
        nlk = NlkState(n=1, l=1, k=2)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [19.50, 24.69, 28.69, 35.19])

        #Calculating optimized projector #   2
        #for l=   1
        nlk = NlkState(n=2, l=1, k=2)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [19.19, 25.13, 29.76, 41.40])

        # Test potentials
        vloc = p.potentials[-1]
        #pl0 = {0: -4.6445128, 1: -15.3234007, 2: 20.9698547, -1: -10.0124145}

        #for l, pot in p.potentials.items():
        #    assert (pot.rmesh[0], pot.rmesh[-1]) == (0.0099582, 2.4161798)
        #    str(l)
        #    assert pot.values[0] == pl0[l]
        #    assert all(pot.rmesh == vloc.rmesh)

        # Test wavefunctions
        ae_wfs, ps_wfs = p.radial_wfs.ae, p.radial_wfs.ps

        # n= 2,  l= 0, kap=-1, all-electron wave function, pseudo w-f
        nlk = NlkState.from_nlkap(n=2, l=0, kap=-1)
        ae20, ps20 = ae_wfs[nlk], ps_wfs[nlk]
        assert ae20[0] == (0.020088, -0.172694)
        assert ps20[0] == (0.020088, 0.030838)
        assert ae20[-1] == (3.976916, 0.037123)
        assert ps20[-1] == (3.976916, 0.037120)

        # scattering, iprj= 2,  l= 0, kap=-1, all-electron wave function, pseudo w-f
        ae_swfs, ps_swfs = p.scattering_wfs.ae, p.scattering_wfs.ps
        nlk = NlkState.from_nlkap(n=2, l=0, kap=-1)
        ae20, ps20 = ae_swfs[nlk], ps_swfs[nlk]
        assert ae20[0] == (0.020088, -0.115970)
        assert ps20[0] == (0.020088, 0.021593)
        assert ae20[-1] == (3.976916, 0.315387)
        assert ps20[-1] == (3.976916, 0.315347)

        #n= 2,  l= 1, kap=-2, all-electron wave function, pseudo w-f
        nlk = NlkState.from_nlkap(n=2, l=1, kap=-2)
        ae21, ps21 = ae_wfs[nlk], ps_wfs[nlk]
        assert ae21[0] == (0.020088, 0.005663)
        assert ps21[0] == (0.020088, 0.001610)
        assert ae21[-1] == (3.976916, 0.099985)
        assert ps21[-1] == (3.976916, 0.099981)

        # Test projectors
        #       l    rmesh        p1           p2
        #!J     0    0.009976     0.069866     0.060936

        prjs = p.projectors

        nlk_1 = NlkState(n=1, l=0, k=1)
        assert prjs[nlk_1][0] == (0.009976, 0.069866)
        assert prjs[nlk_1][-1] == (1.741907, -0.000000)

        nlk_2 = NlkState(n=2, l=0, k=1)
        assert prjs[nlk_2][0] == (0.009976, 0.060936)
        assert prjs[nlk_2][-1] == (1.741907,  0.000000)

        #!J    -1    0.009976     0.001728    -0.000948
        nlk_1 = NlkState(n=1, l=1, k=1)
        nlk_2 = NlkState(n=2, l=1, k=1)
        #print(prjs.keys())
        assert prjs[nlk_1][0] == (0.009976, 0.001728)
        assert prjs[nlk_2][0] == (0.009976, -0.000948)

        #!J     1    0.009976     0.001729    -0.000948
        nlk_1 = NlkState(n=1, l=1, k=2)
        nlk_2 = NlkState(n=2, l=1, k=2)
        assert prjs[nlk_1][0] == (0.009976, 0.001729)
        assert prjs[nlk_2][0] == (0.009976, -0.000948)

        # Test convergence data
        # convergence profiles, (ll=0,lmax), kappa average
        c = p.kene_vs_ecut
        assert c[0].energies[0] == 5.013968
        assert c[0].values[0] == 0.010000
        assert c[0].energies[-1] == 25.333857
        assert c[0].values[-1] == 0.000010
        assert c[1].energies[0] == 19.486256
        assert c[1].values[0] == 0.010000

        # Test log derivatives
        #    log derivativve data for plotting, l= 0
        #    atan(r * ((d psi(r)/dr)/psi(r))), r=  1.62
        #    l, energy, all-electron, pseudopotential

        ae0, ps0 = p.atan_logders.ae[0], p.atan_logders.ps[0]
        assert (ae0.energies[0], ae0.values[0]) == (2.000000, 0.588898)
        assert (ps0.energies[0], ps0.values[0]) == (2.000000, 0.585331)
        assert (ae0.energies[-1], ae0.values[-1]) == (-1.980000, 3.922676)
        assert (ps0.energies[-1], ps0.values[-1]) == (-1.980000, 3.922364)

        # l = -1
        ae1, ps1 = p.atan_logders.ae[-1], p.atan_logders.ps[-1]
        assert (ae1.energies[0], ae1.values[0]) == (2.000000, -2.633391)
        assert (ps1.energies[0], ps1.values[0]) == (2.000000, -2.630644)

        # l = -2
        ae2, ps2 = p.atan_logders.ae[-2], p.atan_logders.ps[-2]
        assert (ae2.energies[0], ae2.values[0]) == (2.000000, -0.713245)
        assert (ps2.energies[0], ps2.values[0]) == (2.000000, -0.537181)

        ## Build the plotter
        plotter = p.get_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

    def test_multi_oncv_plotter(self):
        """Testing MultiOncvPlotter."""
        paths = [filepath("O_fr_v4.out"), filepath("O_sr_v4.out")]
        plotter = MultiOncvPlotter.from_files(paths)
        assert len(plotter) == 2
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
            assert plotter.plot_kene_vs_ecut(show=False)
            assert plotter.plot_atanlogder_econv(show=False)
            assert plotter.plot_den_formfact(ecut=20, show=False)

            if isinstance(plotter, OncvParser) and plotter.parser.is_metapsp:
                assert plotter.plot_vtau(show=False)
                assert plotter.plot_tau(show=False)

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

    def test_scalar_relativistic_cooper_metagga(self):
        """
        Parsing the scalar-relativistic output file produced by metapsp-1.0.1
        """
        # Scalar relativistic output
        p = OncvParser(filepath("29_Cu_m.out"))
        p.scan(verbose=1)
        repr(p); str(p)
        assert p.run_completed
        assert not p.is_oncvpsp and p.is_metapsp

        assert not p.relativistic
        assert p.calc_type == "scalar-relativistic"
        assert p.version == "1.0.1"

        assert p.atsym == "Cu"
        assert p.z == 29.00
        assert p.iexc == 5
        assert p.nc == 3
        assert p.nv == 4
        assert p.lmax == 2
        assert p.rc5 == 1.5
        assert p.rc_l[0] == 1.90000
        assert p.rc_l[1] == 2.10000
        assert p.rc_l[2] == 2.10000

        # In this output, psfile is none
        #assert p.get_input_str()
        #assert p.get_psp8_str()
        #assert p.get_upf_str()

        # Calculating optimized projector #   1
        # for l=   0
        nlk = NlkState(n=1, l=0, k=None)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [12.99, 16.73, 19.91, 22.56])

        #Calculating optimized projector #   2
        # for l=   1
        nlk = NlkState(n=2, l=1, k=None)
        ke = p.kinerr_nlk[nlk]
        self.assert_almost_equal(ke.values_ha, [0.01000, 0.00100, 0.00010, 0.00001])
        self.assert_almost_equal(ke.ecuts, [14.68, 20.8 , 25.39, 29.03])

        # Test potentials
        vloc = p.potentials[-1]
        pl0 = {0: -29.089934, 1: -31.848153, 2: -34.020455, -1:  -27.027998}

        for l, pot in p.potentials.items():
            assert (pot.rmesh[0], pot.rmesh[-1]) == (0.009994, 3.088443)
            str(l)
            assert pot.values[0] == pl0[l]
            assert all(pot.rmesh == vloc.rmesh)

        # Test wavefunctions
        ae_wfs, ps_wfs = p.radial_wfs.ae, p.radial_wfs.ps

        nlk = (3, 2, None)
        ae32, ps32 = ae_wfs[nlk], ps_wfs[nlk]
        assert ae32[0] == (0.009994  , 0.000226)
        assert ps32[0] == (0.009994  , 0.000007)
        assert ae32[-1] == (5.987503 , 0.021621)
        assert ps32[-1] == (5.987503 , 0.021621)

        # Test projectors
        prjs = p.projectors
        assert prjs[(1, 2, None)][0] == (0.009994, 0.000020)
        assert prjs[(2, 2, None)][0] == (0.009994, -0.000010)

        # Test convergence data
        c = p.kene_vs_ecut
        assert c[0].energies[0] == 12.985434
        assert c[0].values[0] == 0.010000
        assert c[0].energies[-1] == 22.558500
        assert c[0].values[-1] == 0.000010
        assert c[1].energies[0] == 19.543935
        assert c[1].values[0] == 0.010000

        # Test log derivatives
        ae0, ps0 = p.atan_logders.ae[0], p.atan_logders.ps[0]
        assert (ae0.energies[0], ae0.values[0]) == (5.000000, 1.148959)
        assert (ps0.energies[0], ps0.values[0]) == (5.000000, 1.089880)

        assert (ae0.energies[-1], ae0.values[-1]) == (-4.980000, 7.624244)
        assert (ps0.energies[-1], ps0.values[-1]) == (-4.980000, 7.624659)

        ae1, ps1 = p.atan_logders.ae[1], p.atan_logders.ps[1]
        assert (ae1.energies[0], ae1.values[0]) == (5.000000, 2.803894)
        assert ps1.values[0] == 3.285840

        # Build the plotter
        plotter = p.get_plotter()
        repr(plotter); str(plotter)
        self._call_plotter_methods(plotter)

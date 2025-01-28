# coding: utf-8
"""Tests for gwr module."""
import numpy as np
#import pymatgen.core.units as pmgu
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.gwr import GwrFile, GwrRobot


class GwrFileTest(AbipyTest):

    def test_gwr_file(self):
        """Testing GwrFile"""
        path = abidata.ref_file("t01o_DS3_GWR.nc")
        with GwrFile(abidata.ref_file(path)) as gwr:
            repr(gwr); str(gwr)
            gwr.to_string(verbose=1)
            assert gwr.structure.formula == "Si2"
            assert len(gwr.sigma_kpoints) == 2
            assert len(gwr.sigma_kpoints) == gwr.nkcalc
            assert gwr.sigma_kpoints[1] == [0.5, 0, 0]
            self.assert_equal(gwr.sigma_kpoints[1].frac_coords, [0.5, 0, 0])
            assert gwr.r.gwr_task == "G0W0"
            assert gwr.r.sig_diago

            self.assert_almost_equal(gwr.ks_dirgaps[0], [2.4542143, 3.4969783])
            self.assert_almost_equal(gwr.qpz0_dirgaps[0], [3.27446, 4.2812], decimal=5)

            # Minimax imaginary tau/omega mesh: !Tabular | # tau, weight(tau), omega, weight(omega)
            ref_data = np.fromstring("""\
2.14325E-01   5.65693E-01   1.17768E-02   2.46619E-02
1.28479E+00   1.69883E+00   4.26019E-02   4.02802E-02
4.04688E+00   4.16331E+00   1.02353E-01   8.69441E-02
1.06328E+01   9.74336E+00   2.42820E-01   2.15922E-01
2.57180E+01   2.20388E+01   6.17146E-01   6.09743E-01
6.02382E+01   5.21455E+01   1.86465E+00   2.44712E+00""", sep=" ")
            ref_data.shape = (6, 4)

            mesh = gwr.minimax_mesh
            assert mesh.ntau == 6
            self.assert_almost_equal(mesh.tau_mesh, ref_data[:, 0], decimal=4)
            self.assert_almost_equal(mesh.tau_wgs, ref_data[:, 1], decimal=4)
            self.assert_almost_equal(mesh.iw_mesh, ref_data[:, 2], decimal=4)
            self.assert_almost_equal(mesh.iw_wgs, ref_data[:, 3], decimal=4)
            self.assert_almost_equal(mesh.min_transition_energy_eV, 3.31673712E-02)
            self.assert_almost_equal(mesh.max_transition_energy_eV, 1.89598634E+00)
            self.assert_almost_equal(mesh.ft_max_err_t2w_cos, 1.69754098E-02)
            self.assert_almost_equal(mesh.ft_max_err_w2t_cos, 5.30899614E-04)
            self.assert_almost_equal(mesh.ft_max_err_t2w_sin, 3.37573093E-01)
            self.assert_almost_equal(mesh.cosft_duality_error, 7.86624507E-04)

            params = gwr.params
            assert params["gwr_ntau"] == 6
            assert params["ecuteps"] == 4
            assert params["ecutsigx"] == 4

            gas_df = gwr.get_dirgaps_dataframe(with_params=True, with_geo=True)
            df = gwr.get_dataframe_sk(spin=0, kpoint=0, with_params=True, with_geo=True)

            if self.has_matplotlib():
                assert mesh.plot_ft_weights(mesh, show=False)
                assert gwr.plot_sigma_imag_axis(kpoint=0, show=False)
                assert gwr.plot_sigma_real_axis(kpoint=0, show=False)
                assert gwr.plot_qps_vs_e0(show=False)
                assert gwr.plot_spectral_functions(show=False)

            if self.has_nbformat():
                assert gwr.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_gwr_robot(self):
        """Testing GwrRobot."""
        path = abidata.ref_file("t01o_DS3_GWR.nc")
        with GwrRobot() as robot:
            robot.add_file("one", path)
            robot.add_file("two", path)
            spin, kpoint, band = 0, 0, 4

            df_sk = robot.get_dataframe_sk(spin, kpoint, with_params=True, ignore_imag=False)
            gaps_df = robot.get_dirgaps_dataframe(sortby="kname", with_params=True)
            df = robot.get_dataframe(sortby="kname", with_params=True, ignore_imag=False)

            if self.has_matplotlib():
                assert robot.plot_selfenergy_conv(spin, kpoint, band, show=False)
                assert robot.plot_qpgaps_convergence(show=False)
                # FIXME
                #assert robot.plot_qpfield_vs_e0("qpe", show=False)
                #assert robot.plot_qpdata_conv_skb(spin, kpoint, band, show=False)

            if self.has_nbformat():
                assert robot.write_notebook(nbpath=self.get_tmpname(text=True))

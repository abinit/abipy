"""Tests for lumi.deltaSCF module"""
import abipy.data as abidata
from abipy.core.testing import AbipyTest
from abipy.lumi.deltaSCF import DeltaSCF

class DeltaSCFTest(AbipyTest):
    
    def test_deltaSCF(self):
        """Testing DeltaSCF"""

        delta_SCF=DeltaSCF.from_four_points_file([abidata.ref_file("site_1_relaxed_gs_out_GSR.nc"),
                                            abidata.ref_file("site_1_unrelaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("site_1_relaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("site_1_unrelaxed_gs_out_GSR.nc")])
        
        self.assert_equal(delta_SCF.natom(),36)
        self.assert_equal(delta_SCF.defect_index('Eu'),0)
        
        self.assert_almost_equal(delta_SCF.delta_r(),0.182233396654827,decimal=5)
        self.assert_almost_equal(delta_SCF.delta_q(),0.8974717558835328,decimal=5)
        self.assert_almost_equal(delta_SCF.effective_mass(),24.254126571027456,decimal=5)

        self.assert_almost_equal(delta_SCF.E_zpl(),1.7086385000075097,decimal=5)
        self.assert_almost_equal(delta_SCF.FWHM_1D(),0.19142848424536402,decimal=5)

        if self.has_matplotlib():
            assert delta_SCF.plot_lineshape_1D_zero_temp(show=False)
            assert delta_SCF.plot_delta_R_distance(defect_symbol='Eu',show=False)
            assert delta_SCF.plot_delta_F_distance(defect_symbol='Eu',show=False)
            assert delta_SCF.displacements_visu(show=False)
            assert delta_SCF.plot_lineshape_1D_zero_temp(show=False)
            assert delta_SCF.draw_displaced_parabolas(show=False)

            nscf_files=[abidata.ref_file("relaxed_gs_nscf_GSR.nc"),
                        abidata.ref_file("unrelaxed_ex_nscf_GSR.nc"),
                        abidata.ref_file("relaxed_ex_nscf_GSR.nc"),
                        abidata.ref_file("unrelaxed_gs_nscf_GSR.nc")]


            assert delta_SCF.plot_four_BandStructures(nscf_files,show=False)




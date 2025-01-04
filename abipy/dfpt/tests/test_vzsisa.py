"""Tests for frozen_phonons"""
import os
import numpy as np
import abipy.data as abidata

from abipy.dfpt.vzsisa import Vzsisa
from abipy.core.testing import AbipyTest


class QhaTest(AbipyTest):

    def test_vzsisa(self):
        """Testing Vzsisa postprocessing tools."""
        # Root points to the directory in the git submodule with the output results.
        root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "Si_v_ZSISA_approximation")

        bo_strains = [96, 98, 100, 102, 104, 106]
        ph_strains = [98, 100, 102, 104, 106] # EinfVib4(D)
        #ph_strains = [96, 98, 100, 102, 104] # EinfVib4(S)
        #ph_strains = [100, 102, 104] # EinfVib2(D)

        gsr_paths = [os.path.join(root, "scale_{:d}_GSR.nc".format(s)) for s in bo_strains]

        ddb_paths = [os.path.join(root, "scale_{:d}_DDB".format(s)) for s in ph_strains]
        anaget_kwargs = {}
        # FIXME
        #FileNotFoundError: [Errno 2] No such file or directory: '/Users/giantomassi/git_repos/abipy/abipy/data/data_v-ZSISA-QHA.git/Si_v_ZSISA_approximation/scale_96_GSR.nc'
        #qha = Vzsisa.from_gsr_ddb_paths(4, gsr_paths, ddb_paths, anaget_kwargs, verbose=1)

        ddb_paths = [os.path.join(root, "scale_{:d}_GSR_DDB".format(s)) for s in bo_strains]
        phdos_paths = [os.path.join(root, "scale_{:d}_PHDOS.nc".format(s)) for s in ph_strains]

        qha = Vzsisa.from_ddb_phdos_files(ddb_paths, phdos_paths)
        tstart, tstop = 0, 800

        # Test basic properties and get methods of qha
        assert qha.bo_nvols == 6
        assert qha.ph_nvols == 5
        assert qha.eos_name == "vinet"
        assert qha.scale_points == "D"
        assert str(qha)
        assert qha.to_string(verbose=1)

        data = qha.fit_energies_vt(qha.bo_volumes, qha.bo_energies[np.newaxis, :].T, tstart=0, tstop=0, num=1)

        self.assert_almost_equal(data.tot_en,
           [[-230.27531394],
           [-230.28382933],
           [-230.28659309],
           [-230.28404059],
           [-230.27657594],
           [-230.26457428]])

        assert data.fits
        self.assert_almost_equal(data.min_en, -230.28659309)
        self.assert_almost_equal(data.min_vol, 40.05665096)
        self.assert_almost_equal(data.temp, 0)

        e2d_v = qha.second_derivative_bo_energy_v(41)
        self.assert_almost_equal(e2d_v, 0.0128727)

        f = qha.get_thermal_expansion_coeff(tstart=0, tstop=1000, num=7, tref=None)
        self.assert_almost_equal(f.values, [-0.0000000e+00,  5.1944952e-07,  7.0431284e-06,  9.4401141e-06,
                                             1.0630764e-05,  1.1414897e-05,  1.2035432e-05])

        vols, fits = qha.vol_E2Vib1(tstop=tstop, tstart=tstart, num=2)
        self.assert_almost_equal(vols, [40.2380733, 40.4418801])

        vols, fits = qha.vol_Einf_Vib1(tstop=tstop, tstart=tstart, num=2)
        self.assert_almost_equal(vols, [40.2445043, 40.4564043])

        vols, fits = qha.vol_Einf_Vib2(tstop=tstop, tstart=tstart, num=2)
        self.assert_almost_equal(vols, [40.2457091, 40.4461412])

        vols, fits = qha.vol_Einf_Vib4(tstop=tstop, tstart=tstart, num=2)
        self.assert_almost_equal(vols, [40.2456922, 40.4467746])
        aas, bbs, ccs = qha.get_abc(volumes=vols)
        self.assert_almost_equal(aas, [3.8466098, 3.8530056])
        self.assert_almost_equal(bbs, [3.8466098, 3.8530056])
        self.assert_almost_equal(ccs, [3.8466098, 3.8530056])

        alphas, betas, gammas = qha.get_angles(vols)
        self.assert_almost_equal(alphas, [60., 60.])
        self.assert_almost_equal(betas, [60., 60.])
        self.assert_almost_equal(gammas, [60., 60.])

        ph_energies = qha.get_free_energy_vt(tstart, tstop, num=2)
        tot_en = qha.phbo_energies[np.newaxis, :].T + ph_energies
        fit = qha.fit_forth(qha.ph_volumes, tot_en, tstart=tstart, tstop=tstop, num=2)
        self.assert_almost_equal(fit.min_vol, [40.24570378, 40.44670049])
        self.assert_almost_equal(fit.min_en, [-230.1654581, -230.5580189])
        self.assert_almost_equal(fit.F2D_V, [0.01425036, 0.01312722])

        vols = qha.vol_E2Vib1_forth(tstart=0, tstop=1000, num=5)
        self.assert_almost_equal(vols, [40.2380458, 40.2411873, 40.3202008, 40.4207685, 40.5276663])

        vols = qha.vol_EinfVib1_forth(tstart=0, tstop=1000, num=5)
        self.assert_almost_equal(vols, [40.2444849, 40.247707 , 40.3291996, 40.4341991, 40.5474182])

        vols = qha.vol_Einf_Vib4_forth(tstart=0, tstop=1000, num=5)
        self.assert_almost_equal(vols, [40.2456751, 40.2432144, 40.3196922, 40.4241022, 40.5411743])

        alphas = qha.get_thermal_expansion_coeff_4th(tstart=0, tstop=1000, num=5, tref=None)
        self.assert_almost_equal(alphas, [-0.0000000e+00,  4.6240312e-06,  9.4381513e-06,  1.1048456e-05, 1.2028183e-05])

        free_energies = qha.get_free_energy_vt(tstart=0, tstop=100, num=2)
        ref_free_energies = [
           [0.12305556, 0.11944085],
           [0.12140287, 0.11786901],
           [0.11976321, 0.11629327],
           [0.11813874, 0.11471864],
           [0.11653109, 0.11314891]
        ]
        self.assert_almost_equal(free_energies, ref_free_energies)

        data = qha.get_thermodynamic_properties(tstart=0, tstop=1000, num=2)
        self.assert_almost_equal(data.tmesh, [   0., 1000.])
        self.assert_almost_equal(data.cv,
            [[0.        , 0.00050313],
             [0.        , 0.0005035 ],
             [0.        , 0.00050386],
             [0.        , 0.0005042 ],
             [0.        , 0.00050453]])

        self.assert_almost_equal(data.free_energy[0], [ 0.12305556, -0.45237198])
        self.assert_almost_equal(data.entropy[0], [0., 0.0009788])
        self.assert_almost_equal(data.zpe, [0.1230556, 0.1214029, 0.1197632, 0.1181387, 0.1165311])

        if self.has_matplotlib():
            assert qha.plot_bo_energies(tstop=tstop, tstart=tstart, num=11, show=False)
            assert qha.plot_vol_vs_t(tstop=tstop, tstart=tstart, num=101, show=False)
            assert qha.plot_abc_vs_t(tstop=tstop, tstart=tstart, num=101, show=False)
            assert qha.plot_abc_vs_t(tstop=tstop, tstart=tstart, num=101, lattice="b", show=False)
            assert qha.plot_thermal_expansion_coeff(tstop=tstop, tstart=tstart ,num=101, show=False)
            assert qha.plot_thermal_expansion_coeff_abc(tstop=tstop, tstart=tstart ,num=101, show=False)
            assert qha.plot_angles_vs_t(tstop=tstop, tstart=tstart, num=101, show=False)
            assert qha.plot_vol_vs_t_4th(tstop=tstop, tstart=tstart, num=101, show=False)
            assert qha.plot_abc_vs_t_4th(tstop=tstop, tstart=tstart, num=101, lattice="a", show=False)
            assert qha.plot_abc_vs_t_4th(tstop=tstop, tstart=tstart, show=False)
            assert qha.plot_thermal_expansion_coeff_4th(tref=293, show=False)
            assert qha.plot_thermal_expansion_coeff_abc_4th(tstop=tstop, tstart=tstart ,num=101, tref=293, show=False)
            assert qha.plot_angles_vs_t_4th(tstop=tstop, tstart=tstart, num=101, angle=3, show=False)

            plotter = qha.get_phdos_plotter()
            plotter.combiplot(show=False)

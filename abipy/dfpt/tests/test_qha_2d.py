"""Tests for QHA_2D"""
import os
import numpy as np
import abipy.data as abidata

from abipy.dfpt.qha_2D import QHA_2D
from abipy.core.testing import AbipyTest


class Qha2dTest(AbipyTest):

    def test_zsisa_approximation(self):

        bo_strains_a = [995, 1000, 1005, 1010, 1015]
        bo_strains_c = [995, 1000, 1005, 1010, 1015]
        phdos_strains_a = [1000, 1005, 1010]
        phdos_strains_c = [1000, 1005, 1010]

        # Root points to the directory in the git submodule with the output results.
        root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_approximation")

        gsr_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_GSR_DDB") for s3 in bo_strains_c] for s1 in bo_strains_a]
        dos_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_PHDOS.nc") for s3 in phdos_strains_c] for s1 in phdos_strains_a]

        bo_strains_ac = [bo_strains_a, bo_strains_c]
        bo_strains_ac = (np.array(bo_strains_ac) - 1000) / 100
        phdos_strains_ac = [phdos_strains_a, phdos_strains_c]
        phdos_strains_ac = (np.array(phdos_strains_ac) - 1000) / 100

        qha = QHA_2D.from_files(gsr_paths, dos_paths, bo_strains_ac, phdos_strains_ac, gsr_file="DDB")

        # Test properties and methods.
        assert qha.pressure == 0

        f_mat = qha.get_vib_free_energies(0, 1000, 2)
        ref_mat = [
            [
                [0.23190045769242368, -1.0257329004129736],
                [0.22879089684005882, -1.0344811855487408],
                [0.0, 0.0],
            ],
            [
                [0.23044275786150215, -1.0297189635203232],
                [0.22735644425413698, -1.038628734783889],
                [0.22427217946933886, -1.0480192052628534],
            ],
            [
                [0.0, 0.0],
                [0.2259076206066778, -1.0429439417530417],
                [0.0, 0.0],
            ],
        ]

        self.assert_almost_equal(f_mat, ref_mat)

        if self.has_matplotlib():
            qha.plot_energies(show=False)
            qha.plot_free_energies(tstop=500, tstart=0, num=6, show=False)
            qha.plot_thermal_expansion(tstop=1000, tstart=0, num=101, show=False)
            qha.plot_lattice(tstop=1000, tstart=0, num=101, show=False)

    def test_qha_2d(self):

        bo_strains_a = [995, 1000, 1005, 1010, 1015]
        bo_strains_c = [995, 1000, 1005, 1010, 1015]
        #bo_strains_a = [1000, 1005, 1010, 1015 , 1020]
        #bo_strains_c = [1000, 1005, 1010, 1015 , 1020]

        # Root points to the directory in the git submodule with the output results.
        root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_QHA")

        #gsr_paths = [[f"scale_{s1}_{s3}/out_GSR.nc" for s3 in bo_strains_c] for s1 in bo_strains_a]
        gsr_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_GSR_DDB") for s3 in bo_strains_c] for s1 in bo_strains_a]
        phdos_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_PHDOS.nc") for s3 in bo_strains_c] for s1 in bo_strains_a]

        bo_strains_ac = [bo_strains_a, bo_strains_c]
        bo_strains_ac = (np.array(bo_strains_ac) - 1000) / 100

        #qha = QHA_2D.from_files(gsr_paths, phdos_paths, bo_strains_ac, phdos_strains_ac, gsr_file="GSR.nc")
        qha = QHA_2D.from_files(gsr_paths, phdos_paths, bo_strains_ac, bo_strains_ac, gsr_file="DDB")

        # Test properties and methods.
        assert qha.pressure == 0

        f_mat = qha.get_vib_free_energies(0, 1000, 2)
        #print(f_mat)
        ref_mat = [
            [
                [ 0.2364734 , -1.0138777 ],
                [ 0.23334375, -1.02191455],
                [ 0.23021055, -1.03048734],
                [ 0.22707712, -1.03958781],
                [ 0.22394442, -1.04917831]
            ],
            [
                [ 0.23500674, -1.0175096 ],
                [ 0.23190046, -1.0257329 ],
                [ 0.2287909 , -1.03448119],
                [ 0.22568001, -1.0437223 ],
                [ 0.22257111, -1.05344296]
            ],
            [
                [ 0.23352563, -1.02131569],
                [ 0.23044276, -1.02971896],
                [ 0.22735644, -1.03862873],
                [ 0.22427218, -1.04801921],
                [ 0.2211813 , -1.05787544]
            ],
            [
                [ 0.23203085, -1.02529729],
                [ 0.22897121, -1.0338748 ],
                [ 0.22590762, -1.04294394],
                [ 0.22284219, -1.05248871],
                [ 0.21977746, -1.06247964]
            ],
            [
                [ 0.2305236 , -1.02944732],
                [ 0.22748688, -1.03820052],
                [ 0.22444611, -1.04742487],
                [ 0.22140302, -1.05711238],
                [ 0.21836046, -1.06723406]
            ]
        ]

        self.assert_almost_equal(f_mat, ref_mat)

        if self.has_matplotlib():
            assert qha.plot_energies(show=False)
            assert qha.plot_free_energies(tstop=500, tstart=0, num=6, show=False)
            assert qha.plot_thermal_expansion(tstop=1000, tstart=0, num=101, show=False)
            assert qha.plot_lattice(tstop=1000, tstart=00, num=101, show=False)

"""Tests for QHA_2D"""
import os
import abipy.data as abidata

from abipy.dfpt.qha_2D import QHA_2D
from abipy.core.testing import AbipyTest


class Qha2dTest(AbipyTest):

    def test_qha_2d(self):

        strains_a = [ 995,1000, 1005, 1010, 1015  ]
        strains_c = [ 995,1000, 1005, 1010, 1015  ]
        strains_a1 = [ 1000, 1005, 1010 ]
        strains_c1 = [ 1000, 1005, 1010 ]

        # Root points to the directory in the git submodule with the output results.
        root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_approximation")

        gsr_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_GSR_DDB") for s3 in strains_c] for s1 in strains_a]
        dos_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_PHDOS.nc") for s3 in strains_c1] for s1 in strains_a1]

        qha = QHA_2D.from_files(gsr_paths, dos_paths, gsr_file="DDB")

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
        #if False:
            qha.plot_energies(show=False)
            qha.plot_free_energies(tstop=500, tstart=0, num=6, show=False)
            qha.plot_thermal_expansion(tstop=1000, tstart=0, num=101, show=False)
            qha.plot_lattice(tstop=1000, tstart=0, num=101, show=False)

"""Tests for qha_general_stress module"""
import os
import abipy.data as abidata

from abipy.dfpt.qha_general_stress import QHA_ZSISA
from abipy.core.testing import AbipyTest


class QhaZSISATest(AbipyTest):

    def test_qha_zsisa(self):

        # Root points to the directory in the git submodule with the output results.
        root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_ThermalStress")

        # Define strain configurations and paths
        strains_a = [1000, 1005, 1010]
        strains_c = [1000, 1005, 1010]

        dos_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_PHDOS.nc") for s3 in strains_c] for s1 in strains_a]
        guess_path = "Relax2o_GSR.nc"
        gsr_BO_paths = os.path.join(root, f"scale_1000_1000/out_GSR_DDB")

        zsisa = QHA_ZSISA.from_files(dos_paths, guess_path, gsr_BO_paths )
        result = zsisa.cal_stress(0.0, pressure=0.0)
        print("Stress calculation result:", result)

        # Make sure zsisa is serializable with pickle as we store an instance in the Abipy Works
        #self.serialize_with_pickle(zsisa)

        #if self.has_matplotlib():
        #if False:

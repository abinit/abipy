"""Tests for qha_general_stress module"""

import os
import abipy.data as abidata
import abipy.core.abinit_units as abu

from abipy.core.testing import AbipyTest
from abipy.electrons.gsr import GsrFile
from abipy.dfpt.qha_general_stress import QHA_ZSISA


class QhaZSISATest(AbipyTest):

    def test_qha_zsisa(self):

        # Root points to the directory in the git submodule with the output results.
        root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_ThermalStress")

        # Define strain configurations and paths
        strains_a = [1000, 1005, 1010]
        strains_c = [1000, 1005, 1010]

        phdos_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_PHDOS.nc") for s3 in strains_c] for s1 in strains_a]
        gsr_guess_path = os.path.join(root, f"find_TEC/Temp_0300_000/Relax2o_GSR.nc")
        gsr_bo_paths = os.path.join(root, f"scale_1000_1000/out_GSR_DDB")
        elastic_BO_paths = os.path.join(root, f"find_TEC/Temp_0300_000/elastic_constant.txt")

        with GsrFile(gsr_guess_path) as gsr:
            structure_guess = gsr.structure
            stress_guess = gsr.cart_stress_tensor * abu.GPa_to_au
            energy_guess = gsr.energy

        zsisa = QHA_ZSISA.from_files(phdos_paths, gsr_bo_paths, verbose=1)
        tdata = zsisa.get_tstress(300.0, 0.0, structure_guess, stress_guess, energy_guess,
                                  mode="TEC", elastic_path=elastic_BO_paths)
        #print("tdata)
        assert tdata.elastic is None

        dtol, gibbs, stress, therm = zsisa.stress_ZSISA_2DOF(300.0, 0.0)

        self.assert_almost_equal(gibbs, -12268.354238255062)
        self.assert_almost_equal(dtol,
                [9.64016894e-09, 9.64014627e-09, 9.27508661e-09, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
        self.assert_almost_equal(stress,
                [4.17248287e-05, 4.17248287e-05, 3.84335822e-05, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
        self.assert_almost_equal(therm,
                [4.8963870998460175e-06, 4.896377032350449e-06, 3.4159028178772595e-06, 0, 0, 0])

        # Root points to the directory in the git submodule with the output results.
        root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_ElasticConstants")

        # Define strain configurations and paths
        strains_a = [1000, 1005, 1010]
        strains_b = [1000, 1005, 1010]
        strains_c = [1000, 1005, 1010]
        strains_d = [1000, 1005, 1010]

        phdos_paths = [[[[os.path.join(root, f"scale_{s1}_{s2}_{s3}_{s4}_1000_1000/out_PHDOS.nc") for s4 in strains_d]
            for s3 in strains_c] for s2 in strains_b] for s1 in strains_a]
        gsr_guess_path = os.path.join(root, f"find_TEC_ECs/Temp_0600_08/Relax2o_GSR.nc")
        gsr_BO_path = os.path.join(root, f"scale_1000_1000_1000_1000_1000_1000/out_GSR.nc")
        elastic_path = os.path.join(root, f"find_TEC_ECs/Temp_0600_08/elastic_constant.txt")

        with GsrFile(gsr_guess_path) as gsr:
            structure_guess = gsr.structure
            stress_guess = gsr.cart_stress_tensor * abu.GPa_to_au

        temp = 600
        pressure = 8.0
        zsisa = QHA_ZSISA.from_files(phdos_paths, gsr_BO_path)
        assert zsisa.dim == (3, 3, 3, 3, 1, 1)

        tdata = zsisa.get_tstress(temp, pressure, structure_guess, stress_guess, energy_guess,
                                  mode='ECs', elastic_path=elastic_path)
        #print(tdata)
        self.assert_almost_equal(tdata.elastic, [
            [201.48244311, 153.66462493, 137.08726476, 0., 0., 0.],
            [153.66462493, 201.48250611, 137.08726776, 0., 0., 0.],
            [137.08725376, 137.08725576, 213.92063478, 0., 0., 0.],
            [0., 0., 0., 28.43553811, 0., 0.],
            [0., 0., 0., 0., 31.676614, 0.],
            [0., 0., 0., 0., 0., 33.597613]]
        )

        dtol, gibbs, stress, therm, elastic = zsisa.stress_ZSISA_3DOF(temp, pressure/abu.HaBohr3_GPa, mode='ECs')

        self.assert_almost_equal(gibbs, -12268.704989002104)
        self.assert_almost_equal(dtol,
            [1.1743820e-09, 1.1743588e-09, 3.9019984e-10, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00])
        self.assert_almost_equal(stress,
            [-0.000218940, -0.000218940, -0.000224898, 0.000000000, 0.000000000, 0.000000000])
        self.assert_almost_equal(therm,
            [4.547207782339523e-06, 4.547201684160713e-06, 1.7056272194618724e-06, 0, 0, 0])

        ecs = [elastic[0,0], elastic[0,1], elastic[0,2], elastic[2,2], elastic[3,3]]
        self.assert_almost_equal(ecs,
            [201.48244311040455, 153.6646249346295, 137.08726476110223, 213.9206347809034, 28.435538105561616])

        # Make sure zsisa is serializable with pickle as we store an instance in the Abipy Works
        #self.serialize_with_pickle(zsisa)

        #if self.has_matplotlib():

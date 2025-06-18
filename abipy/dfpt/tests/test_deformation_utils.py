"""Tests for deformation_utils module"""
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dfpt.deformation_utils import generate_deformations


class DeformationUtilsTest(AbipyTest):

    def test_generate_deformations(self):
        """Testing generate_deformations"""
        eps = 0.005/1.005
        si = self.get_structure("Si")
        structures_dict, inds_6d, spgrp_number = generate_deformations(si, eps)
        assert len(structures_dict) == 3
        assert spgrp_number == 227
        assert inds_6d.shape == (3, 6)
        self.assert_equal(inds_6d[:, 0], [1, 0, 2])
        assert np.all(inds_6d[3:, 1] == 0)

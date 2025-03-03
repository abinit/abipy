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
        structures_dict, strain_inds, spgrp_number = generate_deformations(si, eps)
        assert len(structures_dict) == 3
        assert spgrp_number == 227
        assert strain_inds.shape == (3, 6)
        self.assert_equal(strain_inds[:, 0], [0, 1, 2])
        assert np.all(strain_inds[1:, 1] == 0)

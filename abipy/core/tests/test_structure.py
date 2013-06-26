"""Tests for structure module"""
import numpy as np

from abipy.core.structure import *
from abipy.tests import AbipyTest


#class TestCrystalTools(AbipyTest):
#    "Test crystal tools."
#
#    def test_reciprocal_space(self):
#        """Test crystal tools"""
#        rprimd = np.array([1.,0,0, 0,0.5,0, 0,0.2,1])
#        rprimd.shape = (3,3)
#        rmet, ucvol = metric_vol(rprimd)
#        gprimd = reciprocal_space(rprimd)
#        gmet, bzvol = metric_vol(gprimd)
#        same_rprimd = reciprocal_space(gprimd)
#        self.assert_almost_equal(rprimd, same_rprimd)
#        self.assert_almost_equal(bzvol, (2*np.pi)**3/ucvol)

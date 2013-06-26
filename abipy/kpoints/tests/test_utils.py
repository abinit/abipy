"""Tests kpoints"""

from abipy.kpoints.utils import *
from abipy.tests import *

##########################################################################################

class TestWrapWS(AbipyTest):
    def test_wrap_to_ws(self):
        """Testing wrap_to_ws"""
        self.assertAlmostEqual(wrap_to_ws( 0.5), 0.5)
        self.assertAlmostEqual(wrap_to_ws(-0.5), 0.5)
        self.assertAlmostEqual(wrap_to_ws( 0.2), 0.2)
        self.assertAlmostEqual(wrap_to_ws(-0.3),-0.3)
        self.assertAlmostEqual(wrap_to_ws( 0.7),-0.3)
        self.assertAlmostEqual(wrap_to_ws( 2.3), 0.3)
        self.assertAlmostEqual(wrap_to_ws(-1.2),-0.2)
        #self.assertAlmostEqual( wrap_to_ws([0.5,2.3,-1.2]), np.array([0.5,0.3,-0.2]) )

##########################################################################################

class TestWrapBZ(AbipyTest):
    def test_wrap_to_bz(self):
        """Testing wrap_to_bz"""
        self.assertAlmostEqual(wrap_to_bz( 0.0), 0.0)
        self.assertAlmostEqual(wrap_to_bz( 1.0), 0.0)
        self.assertAlmostEqual(wrap_to_bz( 0.2), 0.2)
        self.assertAlmostEqual(wrap_to_bz(-0.2), 0.8)
        self.assertAlmostEqual(wrap_to_bz( 3.2), 0.2)
        self.assertAlmostEqual(wrap_to_bz(-3.2), 0.8)

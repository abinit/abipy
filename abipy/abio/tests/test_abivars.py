"""Tests for structure module"""
from __future__ import print_function, division, unicode_literals, absolute_import

from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import *
from abipy.core.testing import *
from abipy.abio.abivars import AbinitInputFile

class TestAbinitInputFile(AbipyTest):

    def test_simple_input(self):
        # Reference structure
        s = ("acell 1 1 1 rprim 1 0 0 0 1 0 0 0 1 natom 1"
             "ntypat 1 typat 1 znucl 14 xred 0 0 0")
        inp = AbinitInputFile.from_string(s)

        print(inp)
        si1_structure = inp.structure
        assert inp.ndtset == 1 and len(si1_structure) == 1 and si1_structure.formula == "Si1"

        # Default values for acell (3*1) and rprim (np.eye(3))
        s = "natom 1" "ntypat 1 typat 1 znucl 14 xred 0 0 0"
        inp = AbinitInputFile.from_string(s)
        assert inp.structure == si1_structure

        # Same structure but with * syntax.
        s = "acell 3*1 natom 1 ntypat 1 typat 1 znucl 14 xred 0 0 0"
        inp = AbinitInputFile.from_string(s)
        assert inp.structure == si1_structure

        # xcart instead of xred and more variables using * syntax. 
        s = "acell 1 2*1 natom 1*1 ntypat 1*1 typat 1*1 znucl 1*14 xcart 0 0 0"
        inp = AbinitInputFile.from_string(s)
        assert inp.structure == si1_structure

        # acell and rprimd with unit.
        s = "acell 1 2*1 Bohr rprim 1 0 0 0 1 0 0 0 1 Bohr natom 1*1 ntypat 1*1 typat 1*1 znucl 1*14 xangst 0 0 0"
        inp = AbinitInputFile.from_string(s)
        assert inp.structure == si1_structure

        # TODO Angdeg sqrt()

        #assert 0

    def test_input_with_datasets(self):
        # H2 molecule in a big box
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2 * 1
          acell1 10 10 10  # this is a comment
          acell2 20 20 20
          xcart -0.7 0.0 0.0 0.7 0.0 0.0 
        """
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 2
        s0, s1 = inp.dtsets[0].structure, inp.dtsets[1].structure
        assert s0 != s1  
        assert s1.volume == 8 * s0.volume

        # same input but with global acell
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2 * 1
          acell  10 10 10 ! this is a comment
          acell2 20 20 20 ! bohrs
          xcart -0.7 0.0 0.0 0.7 0.0 0.0 
        """
        inp = AbinitInputFile.from_string(s)
        assert s0 == inp.dtsets[0].structure
        assert s1 == inp.dtsets[1].structure

        # same input in terms of an arithmetic series in acell
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2 * 1
          acell:  10 10 10    
          acell+ 10 10 10 # bohrs
          xcart -0.7 0.0 0.0 0.7 0.0 0.0 
        """
        inp = AbinitInputFile.from_string(s)
        assert s0 == inp.dtsets[0].structure
        assert s1 == inp.dtsets[1].structure

    def test_input_with_serie(self):
        # Test arithmetic and geometric series with ecut.
        s = """
          ndtset 3
          natom 2 ntypat 1 znucl 1 typat 2 * 1
          acell 10 10 10 xcart -0.7 0.0 0.0 0.7 0.0 0.0 
          ecut: 10 ecut+ 5
          pawecutdg: 2 pawecutdg* 3
        """
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 3
        self.assertArrayEqual([dt["ecut"] for dt in inp.dtsets], [10, 15, 20])
        self.assertArrayEqual([dt["pawecutdg"] for dt in inp.dtsets], [2, 6, 18])

        # Test arithmetic series with xcart.
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2 * 1
          acell 1 1 1 
          xcart: -0.7 0.0 0.0 0.7 0.0 0.0 
          xcart+ -0.1 0.0 0.0 0.1 0.0 0.0 
        """
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 2 and inp.structure is None

        s0, s1 = inp.dtsets[0].structure, inp.dtsets[1].structure

        self.assert_almost_equal(s0.cart_coords.ravel() / bohr_to_ang, [-0.7, 0, 0, 0.7, 0, 0]) 
        self.assert_almost_equal(s1.cart_coords.ravel() / bohr_to_ang, [-0.8, 0, 0, 0.8, 0, 0]) 


if __name__ == "__main__":
    import unittest
    unittest.main()

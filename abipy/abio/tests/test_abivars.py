"""Tests for abivars module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import os
import abipy.data as abidata

from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import *
from abipy.core.testing import AbipyTest
from abipy.abio.abivars import AbinitInputFile, AbinitInputParser


class TestAbinitInputParser(AbipyTest):

    def test_static_methods(self):
        """Testing AbinitInputParser static methods."""
        p = AbinitInputParser()
        repr(p); str(p)
        assert p.varname_dtindex("acell1") == ("acell", 1)
        assert p.varname_dtindex("fa1k2") == ("fa1k", 2)

        from math import sqrt
        assert p.eval_abinit_operators(["1/2"]) == [str(1/2)]
        assert p.eval_abinit_operators(["sqrt(3.)"]) == [str(sqrt(3.))]
        assert p.eval_abinit_operators(["+sqrt(3.)"]) == [str(sqrt(3.))]
        assert p.eval_abinit_operators(["-sqrt(3.)"]) == [str(-sqrt(3.))]


class TestAbinitInputFile(AbipyTest):

    def test_simple_input(self):
        # Reference structure
        s = ("acell 1 1 1 rprim 1 0 0 0 1 0 0 0 1 natom 1 "
             "ntypat 1 typat 1 znucl 14 xred 0 0 0 ")
        inp = AbinitInputFile.from_string(s)
        repr(inp); str(inp)
        assert inp.to_string(verbose=1)
        assert inp._repr_html_()
        si1_structure = inp.structure
        assert inp.ndtset == 1 and len(si1_structure) == 1 and si1_structure.formula == "Si1"

        # Default values for acell (3*1) and rprim (np.eye(3))
        s = "natom 1 ntypat 1 typat 1 znucl 14 xred 0 0 0"
        inp = AbinitInputFile.from_string(s)
        assert inp.structure == si1_structure

        # Same structure but with * syntax.
        s = "acell 3*1 natom 1 ntypat 1 typat 1 znucl 14 xred 3*0e0"
        inp = AbinitInputFile.from_string(s)
        assert inp.structure == si1_structure
        assert not inp.has_multi_structures

        # xcart instead of xred and more variables using * syntax.
        s = "acell 1 2*1 natom 1*1 ntypat 1*1 typat 1*1 znucl *14 xcart 0d0 0d0 0d0"
        inp = AbinitInputFile.from_string(s)
        assert inp.structure == si1_structure

        # acell and rprimd with unit.
        s = ("acell 1 2*1 Bohr rprim 1 0 0 0 1 0 0 0 1 Bohr natom 1 "
             "ntypat 1 typat *1 znucl 1*14 xangst 0 0 0")
        with AbinitInputFile.from_string(s) as inp:
            assert inp.structure == si1_structure
            if self.has_nbformat():
                inp.write_notebook(nbpath=self.get_tmpname(text=True))

        # TODO Angdeg sqrt(4) sqrt(4/2)
        #assert 0

    def test_input_with_datasets(self):
        # H2 molecule in a big box
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2*1
          acell1 10 10 10  # this is a comment
          acell2 20 20 20
          xcart -0.7 0.0 0.0 0.7 0.0 0.0
        """
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 2
        dt0 = inp.datasets[0]
        repr(dt0); str(dt0)
        assert dt0.to_string(verbose=2)
        assert dt0._repr_html_()
        s0, s1 = inp.datasets[0].structure, inp.datasets[1].structure
        assert s0 != s1
        assert s1.volume == 8 * s0.volume
        repr(inp); str(inp)

        # same input but with global acell
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2*1
          acell  10 10 10 ! this is a comment
          acell2 20 20 20 ! bohrs
          xcart -0.7 0.0 0.0 0.7 0.0 0.0
        """
        inp = AbinitInputFile.from_string(s)
        assert s0 == inp.datasets[0].structure
        assert s1 == inp.datasets[1].structure
        repr(inp); str(inp)

        d = inp.datasets[1].get_vars()
        assert "natom" not in d and len(d) == 0

        # same input in terms of an arithmetic series in acell
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2*1
          acell:  10 10 10
          acell+ 10 10 10 # bohrs
          xcart -0.7 0.0 0.0 0.7 0.0 0.0
        """
        inp = AbinitInputFile.from_string(s)
        assert s0 == inp.datasets[0].structure
        assert s1 == inp.datasets[1].structure
        str(inp)

    def test_input_with_serie(self):
        # Test arithmetic and geometric series with ecut.
        s = """
          ndtset 3
          natom 2 ntypat 1 znucl 1 typat 2*1
          acell 10 10 10 xcart -0.7 0.0 0.0 0.7 0.0 0.0
          ecut: 10 ecut+ 5
          pawecutdg: 2 pawecutdg* 3
        """
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 3
        self.assertArrayEqual([dt["ecut"] for dt in inp.datasets], [10, 15, 20])
        self.assertArrayEqual([dt["pawecutdg"] for dt in inp.datasets], [2, 6, 18])
        repr(inp); str(inp)

        # Test arithmetic series with xcart.
        s = """
          ndtset 2
          natom 2 ntypat 1 znucl 1 typat 2*1
          acell 1 1 1
          xcart: -0.7 0.0 0.0 0.7 0.0 0.0
          xcart+ -0.1 0.0 0.0 0.1 0.0 0.0
        """
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 2 and inp.structure is None

        s0, s1 = inp.datasets[0].structure, inp.datasets[1].structure
        self.assert_almost_equal(s0.cart_coords.ravel() / bohr_to_ang, [-0.7, 0, 0, 0.7, 0, 0])
        self.assert_almost_equal(s1.cart_coords.ravel() / bohr_to_ang, [-0.8, 0, 0, 0.8, 0, 0])
        str(inp)

    def test_tricky_inputs(self):
        """Testing tricky inputs"""
        s = """\
# define kpt mesh
kptrlatt 2 2 -2
        -2 2 -2
        -2 2  2

#definition of the elementary cell
natom 2
ntypat 2
znucl 31 1*33
typat 1*1 2

acell 3*5.6533 angstrom # expt value

rprim 0   1/2 1/2
      1/2 0   1/2
      1/2 1/2 0

include "foo"
xred 3*0 3*1/4
"""
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 1 and inp.structure is not None and len(inp.structure) == 2
        self.assertArrayEqual(inp.structure[0].frac_coords, [0, 0, 0])
        self.assertArrayEqual(inp.structure[0].specie.symbol, "Ga")
        self.assertArrayEqual(inp.structure[1].frac_coords, [1/4, 1/4, 1/4])
        self.assertArrayEqual(inp.structure[1].specie.symbol, "As")
        mat = 5.6533 * np.array([0, 1/2, 1/2, 1/2, 0, 1/2, 1/2, 1/2, 0])
        mat.shape = (3, 3)
        self.assertArrayEqual(inp.structure[1].lattice.matrix, mat)

        # tutorial/input/tbase2_1
        # 2 datasets with different natom (should use typat[:1] in 2nd dataset)
        s = """\
ndtset 2

#Definition of the unit cell and ecut,
#for which one will have to make a convergence study
acell 10 10 10
ecut 10

   natom1  2             # There are two atoms
   xcart1  -0.7  0.0 0.0 # The starting values of the
            0.7  0.0 0.0 # atomic coordinates

#Second dataset : get the total energy of the isolated atom
   natom2  1             # There is one atom
   xcart2  0.0 0.0 0.0   # The atom is located at the origin

#Definition of the atom types
ntypat 1          # There is only one type of atom
znucl 1           # The keyword "znucl" refers to the atomic number of the

#Definition of the atoms
typat 1 1         # For the first dataset, both numbers will be read,
                  # while for the second dataset, only one number will be read
"""
        inp = AbinitInputFile.from_string(s)
        assert inp.ndtset == 2 and inp.structure is None
        assert len(inp.datasets[0].structure) == 2
        assert len(inp.datasets[1].structure) == 1
        str(inp)

    def test_all_inputs_in_tests(self):
        """
        Try to parse all Abinit input files inside the Abinit `tests` directory.
        Requires $ABINIT_HOME_DIR env variable.
        """
        abi_homedir = os.environ.get("ABINIT_HOME_DIR")
        if abi_homedir is not None:
            #raise self.SkipTest("Environment variable `ABINIT_HOME_DIR` is required for this test.")
            abitests_dir = os.path.join(abi_homedir, "tests")
        else:
            abitests_dir = os.path.join(abidata.dirpath, "refs/si_bse")

        from abipy.abio.abivars import validate_input_parser
        assert os.path.exists(abitests_dir)
        retcode = validate_input_parser(abitests_dir=abitests_dir)
        assert retcode == 0

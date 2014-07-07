"""Tests for structure module"""
from __future__ import print_function, division

import numpy as np
import abipy.data as data

from abipy.core.structure import *
from abipy.core.testing import *
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.finder import SymmetryFinder

class TestStructure(AbipyTest):
    """Unit tests for Structure."""

    def test_structure_from_ncfiles(self):
        """Initialize Structure from Netcdf data files"""

        for filename in data.WFK_NCFILES + data.GSR_NCFILES:
            print("About to read file %s" % filename)
            structure = Structure.from_file(filename)
            print(structure)

            # All nc files produced by ABINIT should have info on the spacegroup.
            self.assertTrue(structure.has_spacegroup)

            # Call pymatgen machinery to get the high-symmetry stars.
            print(structure.hsym_stars)

            if self.which("xcrysden") is not None:
                # Export data in Xcrysden format.
                structure.export(".xsf")

    def test_utils(self):
        """Test utilities for the generation of Abinit inputs."""
        structure = data.structure_from_ucell("MgB2")

        self.serialize_with_pickle(structure)

        pseudos = data.pseudos("12mg.pspnc", "5b.pspnc")
        nval = structure.calc_nvalence(pseudos)
        self.assertEqual(nval, 8)
        shiftk = structure.calc_shiftk()
        self.assert_equal(shiftk, [[0.0, 0.0, 0.5]])

    def test_fphonons(self):
        """ This is not a real test, just to show how to use it ! """
        rprimd = np.array([[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])
        #rprimd = rprimd*6.7468
        rprimd = rprimd*10
        lattice = Lattice(rprimd)
        structure = Structure(lattice, ["Ga", "As"],
                                      [[0, 0, 0], [0.25, 0.25, 0.25]])
        old_structure = structure.copy()


        print(old_structure.lattice._matrix)
        for site in old_structure:
            print(structure.lattice.get_cartesian_coords(site.frac_coords))

        qpoint = [1/2, 1/2,1/2]
        mx_sc = [2, 2, 2]
        scale_matrix = structure.get_smallest_supercell(qpoint, max_supercell=mx_sc)
        print("Scale_matrix = ", scale_matrix)
        #scale_matrix = 2*np.eye(3)
        structure.frozen_phonon(qpoint, 0.1*np.array([[1,1,1], [-1,-1,-1]]), do_real=True, frac_coords=False, max_supercell=mx_sc, scale_matrix = scale_matrix)
        print("Structure = ", structure)

        print(structure.lattice._matrix)
        for site in structure:
            print(structure.lattice.get_cartesian_coords(site.frac_coords))


        # Test symmetries

        real_finder = SymmetryFinder(structure)

        real_symmops = real_finder.get_point_group_operations(cartesian=True)
        real_symrel = real_finder.get_point_group_operations(cartesian=False)

        # # Real space matrix
        # mat = structure.lattice.matrix.T
        # invmat = np.linalg.inv(mat)
        # mat2 = structure.reciprocal_lattice.matrix.T
        # invmat2 = np.linalg.inv(mat2)
        #
        #
        # q1 = np.array([1/2,0,0])
        # qfull1 = np.dot(q1,invmat)
        # qtofind = np.array([1/2,1/2,1/2])
        # qfulltofind = np.dot(qtofind,invmat)
        # print("q1 = ",q1)
        # print("q2 = ",qtofind)
        # # print("nsym = ",len(real_symmops))
        # for symmop in real_symmops:
        #
        #     symcart = symmop.rotation_matrix
        #     print("symcart = ",symcart) # In cartesian coordinates
        #
        #     symrel = np.dot(invmat, np.dot(symcart, mat))
        #
        #     print("symrel = ",symrel) # In reduced coordinates
        #
        #     # Symmetry in reciprocal space:
        #     symrecicart = np.linalg.inv(symcart)
        #
        #     print("symrecicart = ",symrecicart)
        #
        #     symrecirel = np.dot(invmat2, np.dot(symrecicart, mat2))
        #
        #     print("symrecirel = ",symrecirel)
        #
        #     qfull2 = np.dot(symrecirel,q1)
        #     print("Q -> ",qfull2)
        #     # print("symrec = ",symrec)
        #     if(np.allclose(qfull2,-qtofind) or np.allclose(qfull2,qtofind)):
        #       print("I found it !!!!")



        # We should add some checks here


if __name__ == "__main__":
    import unittest
    unittest.main()

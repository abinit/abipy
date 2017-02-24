"""Tests for structure module"""
from __future__ import print_function, division

import numpy as np
import sys
import abipy.data as abidata

from abipy.core.structure import *
from abipy.core.testing import *
from pymatgen.core.lattice import Lattice


class TestStructure(AbipyTest):
    """Unit tests for Structure."""

    def test_structure_from_ncfiles(self):
        """Initialize Structure from Netcdf data files"""

        for filename in abidata.WFK_NCFILES + abidata.GSR_NCFILES:
            print("About to read file %s" % filename)
            structure = Structure.from_file(filename)
            print(structure)
            assert structure.__class__ is Structure

            # All nc files produced by ABINIT should have info on the spacegroup.
            assert structure.has_abi_spacegroup

            # Call pymatgen machinery to get the high-symmetry stars.
            print(structure.hsym_stars)

            geodict = structure.get_dict4frame()
            assert geodict["abispg_num"] is not None

            #if self.which("xcrysden") is not None:
            #    # Export data in Xcrysden format.
            #    structure.export(self.get_tmpname(text=True, suffix=".xsf"))

            if self.has_ase():
                assert structure == Structure.from_ase_atoms(structure.to_ase_atoms())

    def test_utils(self):
        """Test utilities for the generation of Abinit inputs."""
        structure = abidata.structure_from_ucell("MgB2")
        structure.abi_sanitize()
        structure.get_conventional_standard_structure()
        print(structure.abi_string)
        print(structure.spglib_summary(verbose=10))
        #print(structure.__repr_html__())

        self.serialize_with_pickle(structure)

        pseudos = abidata.pseudos("12mg.pspnc", "5b.pspnc")
        assert structure.num_valence_electrons(pseudos) == 8
        self.assert_equal(structure.calc_shiftk() , [[0.0, 0.0, 0.5]])

        # Test notebook generation.
        if self.has_nbformat():
            structure.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_frames_from_structures(self):
        """Testing frames from structures."""
        mgb2 = abidata.structure_from_ucell("MgB2")
        sic = abidata.structure_from_ucell("SiC")
        alas = abidata.structure_from_ucell("AlAs")
        dfs = frames_from_structures([mgb2, sic, alas], index=None, with_spglib=True, cart_coords=False)

        assert dfs.lattice is not None
        assert dfs.coords is not None
        assert dfs.structures is not None
        formulas = [struct.composition.reduced_formula for struct in dfs.structures]
        assert formulas == ["MgB2", "SiC", "AlAs"]

    def test_frozen_phonons(self):
        """ This is not a real test, just to show how to use it ! """
        rprimd = np.array([[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])
        #rprimd = rprimd*6.7468
        rprimd = rprimd *10.60 * 0.529
        lattice = Lattice(rprimd)
        structure = Structure(lattice, ["Ga", "As"],
                                      [[0, 0, 0], [0.25, 0.25, 0.25]])
        old_structure = structure.copy()

        print(old_structure.lattice._matrix)
        for site in old_structure:
            print(structure.lattice.get_cartesian_coords(site.frac_coords))

        #qpoint = [1/2, 1/2,1/2]
        qpoint = [0,0,0]
        mx_sc = [2, 2, 2]
        scale_matrix = structure.get_smallest_supercell(qpoint, max_supercell=mx_sc)
        scale_matrix = 2*np.eye(3)
        #print("Scale_matrix = ", scale_matrix)
        #scale_matrix = 2*np.eye(3)
        natoms = int(np.round(2*np.linalg.det(scale_matrix)))

        structure.write_vib_file(sys.stdout, qpoint, 0.1*np.array([[1,1,1], [1,1,1]]), do_real=True, frac_coords=False, max_supercell=mx_sc, scale_matrix = scale_matrix)
        structure.write_vib_file(sys.stdout, qpoint, 0.1*np.array([[1,1,1], [-1,-1,-1]]), do_real=True, frac_coords=False, max_supercell=mx_sc, scale_matrix = scale_matrix)

        structure.frozen_phonon(qpoint, 0.1*np.array([[1,1,1], [-1,-1,-1]]), do_real=True, frac_coords=False, max_supercell=mx_sc, scale_matrix = scale_matrix)

        print("Structure = ", structure)

        #print(structure.lattice._matrix)
        #for site in structure:
        #    print(structure.lattice.get_cartesian_coords(site.frac_coords))

        # We should add some checks here

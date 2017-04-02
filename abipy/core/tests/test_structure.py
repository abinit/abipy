"""Tests for structure module"""
from __future__ import print_function, division

import numpy as np
import sys
import abipy.data as abidata

from pymatgen.core.lattice import Lattice
from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import *
from abipy.core.testing import AbipyTest



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

            # Export data in Xcrysden format.
            #structure.export(self.get_tmpname(text=True, suffix=".xsf"))
            #visu = structure.visualize("vesta")

            if self.has_ase():
                assert structure == Structure.from_ase_atoms(structure.to_ase_atoms())

    def test_utils(self):
        """Test utilities for the generation of Abinit inputs."""
        # Test as_structure and from/to abivars
        si = Structure.as_structure(abidata.cif_file("si.cif"))
        assert si.formula == "Si2"
        assert si.abi_spacegroup is None and not si.has_abi_spacegroup

        if self.has_matplotlib():
            si.show_bz(show=False)

        assert si is Structure.as_structure(si)
        assert si == Structure.as_structure(si.to_abivars())
        assert si == Structure.from_abivars(si.to_abivars())
        assert len(si.abi_string)
        assert si.reciprocal_lattice == si.lattice.reciprocal_lattice
        #si.si_kptbounds()
        #si.calc_ksampling(nksmall=10)
        #si.calc_ngkpt(nksmall=10)
        #si.calc_shiftk()

        si = Structure.from_material_id("mp-149", api_key="8pkvwRLQSCVbW2Fe")
        assert si.formula == "Si2"

        mgb2 = abidata.structure_from_ucell("MgB2")
        if self.has_ase():
            mgb2.abi_primitive()
        # FIXME: This is buggy
        #print(mgb2.get_sorted_mgb2())
        #assert [site.species_string for site in mgb2.get_sorted_structure()] == ["B", "B", "Mg"]
        mgb2.abi_sanitize()
        mgb2.get_conventional_standard_structure()
        assert len(mgb2.abi_string)
        assert len(mgb2.spglib_summary(verbose=10))
        #print(structure.__repr_html__())

        self.serialize_with_pickle(mgb2)

        pseudos = abidata.pseudos("12mg.pspnc", "5b.pspnc")
        assert mgb2.num_valence_electrons(pseudos) == 8
        assert mgb2.valence_electrons_per_atom(pseudos) == [2, 3, 3]
        self.assert_equal(mgb2.calc_shiftk() , [[0.0, 0.0, 0.5]])

        #bmol = Structure.boxed_molecule(pseudos, cart_coords=[[0, 0, 0], [5, 5, 5]])
        #self.assert_almost_equal(bmol.volume, (10 * bohr_to_ang) ** 3)
        #batom = Structure.boxed_atom(abidata.pseudo("12mg.pspnc"), cart_coords=[1, 2, 3], acell=(10, 20, 30))
        #assert isinstance(batom, Structure)
        #assert len(batom.cart_coords) == 1 and np.all(batom.cart_coords[0] == [1, 2, 3])

        #bcc = Structure.bcc(a=10, species, primitive)
        #fcc = Structure.bcc(a=10, species, primitive)
        #rock = Structure.rocksalt(a=10, species)
        #perov = Structure.ABO3(a=10, species)

        # Test notebook generation.
        if self.has_nbformat():
            mgb2.write_notebook(nbpath=self.get_tmpname(text=True))

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

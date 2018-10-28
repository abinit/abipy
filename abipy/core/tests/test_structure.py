"""Tests for structure module"""
from __future__ import print_function, division, absolute_import, unicode_literals

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
            str(structure)
            structure.to_string(verbose=2)
            assert structure.__class__ is Structure

            # All nc files produced by ABINIT should have info on the spacegroup.
            assert structure.has_abi_spacegroup

            # Call pymatgen machinery to get the high-symmetry stars.
            str(structure.hsym_stars)

            geodict = structure.get_dict4pandas()
            assert geodict["abispg_num"] is not None

            # Export data in Xcrysden format.
            #structure.export(self.get_tmpname(text=True, suffix=".xsf"))
            #visu = structure.visualize(appname="vesta")
            #assert callable(visu)

            if self.has_ase():
                assert structure == Structure.from_ase_atoms(structure.to_ase_atoms())

    def test_utils(self):
        """Test utilities for the generation of Abinit inputs."""
        # Test as_structure and from/to abivars
        si = Structure.as_structure(abidata.cif_file("si.cif"))
        assert si.formula == "Si2"
        assert si.latex_formula == "Si$_{2}$"
        assert si.abi_spacegroup is None and not si.has_abi_spacegroup
        assert "ntypat" in si.to(fmt="abivars")

        spgroup = si.spgset_abi_spacegroup(has_timerev=True)
        assert spgroup is not None
        assert si.has_abi_spacegroup
        assert si.abi_spacegroup.spgid == 227
        kfrac_coords = si.get_kcoords_from_names(["G", "X", "L", "Gamma"])
        self.assert_equal(kfrac_coords,
            ([[0. , 0. , 0. ], [0.5, 0. , 0.5], [0.5, 0.5, 0.5], [0. , 0. , 0. ]]))

        with self.assertRaises(TypeError):
            Structure.as_structure({})
        with self.assertRaises(TypeError):
            Structure.as_structure([])

        si_wfk = Structure.as_structure(abidata.ref_file("si_scf_WFK.nc"))
        assert si_wfk.formula == "Si2"
        si_wfk.print_neighbors(radius=2.5)

        assert si_wfk.has_abi_spacegroup
        # Cannot change spacegroup
        with self.assertRaises(ValueError):
            si_wfk.spgset_abi_spacegroup(has_timerev=True)

        # K and U are equivalent. [5/8, 1/4, 5/8] should return U
        assert si_wfk.findname_in_hsym_stars([3/8, 3/8, 3/4]) == "K"
        assert si_wfk.findname_in_hsym_stars([5/8, 1/4, 5/8]) == "U"

        # TODO: Fix order of atoms in supercells.
        # Test __mul__, __rmul__ (should return Abipy structures)
        assert si_wfk == 1 * si_wfk
        supcell = si_wfk * [2, 2, 2]
        assert len(supcell) == 8 * len(si_wfk) and hasattr(supcell, "abi_string")

        si_abi = Structure.from_file(abidata.ref_file("refs/si_ebands/run.abi"))
        assert si_abi.formula == "Si2"
        self.assert_equal(si_abi.frac_coords, [[0, 0, 0], [0.25, 0.25, 0.25]])

        si_abo = Structure.from_file(abidata.ref_file("refs/si_ebands/run.abo"))
        assert si_abo == si_abi
        assert "ntypat" in si_abi.to(fmt="abivars")

        znse = Structure.from_file(abidata.ref_file("refs/znse_phonons/ZnSe_hex_qpt_DDB"))
        assert len(znse) == 4
        assert znse.formula == "Zn2 Se2"
        self.assert_almost_equal(znse.frac_coords.flat, [
            0.33333333333333,  0.66666666666667, 0.99962203020000,
            0.66666666666667,  0.33333333333333, 0.49962203020000,
            0.33333333333333,  0.66666666666667, 0.62537796980000,
            0.66666666666667,  0.33333333333333, 0.12537796980000])

        from abipy.core.structure import diff_structures
        diff_structures([si_abi, znse], headers=["si_abi", "znse"], fmt="abivars", mode="table")
        diff_structures([si_abi, znse], headers=["si_abi", "znse"], fmt="abivars", mode="diff")

        # From pickle file.
        import pickle
        tmp_path = self.get_tmpname(suffix=".pickle")
        with open(tmp_path, "wb") as fh:
            pickle.dump(znse, fh)
        same_znse = Structure.from_file(tmp_path)
        assert same_znse == znse
        same_znse = Structure.as_structure(tmp_path)
        assert same_znse == znse

        for fmt in ["abivars", "cif", "POSCAR", "json", "xsf", "qe", "siesta", "wannier90"]:
            assert len(znse.convert(fmt=fmt)) > 0

        for fmt in ["abinit", "w90", "siesta"]:
            assert len(znse.get_kpath_input_string(fmt=fmt)) > 0

        oxi_znse = znse.get_oxi_state_decorated()
        assert len(oxi_znse.abi_string)
        from pymatgen.core.periodic_table import Specie
        assert Specie("Zn", 2) in oxi_znse.composition.elements
        assert Specie("Se", -2) in oxi_znse.composition.elements

        system = si.spget_lattice_type()
        assert system == "cubic"

        e = si.spget_equivalent_atoms(printout=True)
        assert len(e.irred_pos) == 1
        self.assert_equal(e.eqmap[0], [0, 1])
        for irr_pos in e.irred_pos:
            assert len(e.eqmap[irr_pos]) > 0
        assert "equivalent_atoms" in e.spgdata

        if self.has_matplotlib():
            assert si.plot_bz(show=False)
            assert si.plot_bz(pmg_path=False, show=False)
            assert si.plot(show=False)
            if sys.version[0:3] > '2.7':
                # pmg broke py compatibility
                assert si.plot_xrd(show=False)

        if self.has_mayavi():
            #assert si.vtkview(show=False)  # Disabled due to (core dumped) on travis
            assert si.mayaview(show=False)

        assert si is Structure.as_structure(si)
        assert si == Structure.as_structure(si.to_abivars())
        assert si == Structure.from_abivars(si.to_abivars())
        assert len(si.abi_string)
        assert si.reciprocal_lattice == si.lattice.reciprocal_lattice

        kptbounds = si.calc_kptbounds()
        ksamp = si.calc_ksampling(nksmall=10)

        shiftk = [[ 0.5,  0.5,  0.5], [ 0.5,  0. ,  0. ], [ 0. ,  0.5,  0. ], [ 0. ,  0. ,  0.5]]
        self.assert_equal(si.calc_ngkpt(nksmall=2), [2, 2, 2])
        self.assert_equal(si.calc_shiftk(), shiftk)
        self.assert_equal(ksamp.ngkpt, [10, 10, 10])
        self.assert_equal(ksamp.shiftk, shiftk)

        lif = Structure.from_abistring("""
acell      7.7030079150    7.7030079150    7.7030079150 Angstrom
rprim      0.0000000000    0.5000000000    0.5000000000
           0.5000000000    0.0000000000    0.5000000000
           0.5000000000    0.5000000000    0.0000000000
natom      2
ntypat     2
typat      1 2
znucl      3 9
xred       0.0000000000    0.0000000000    0.0000000000
           0.5000000000    0.5000000000    0.5000000000
""")
        assert lif.formula == "Li1 F1"
        same = Structure.rocksalt(7.7030079150, ["Li", "F"], units="ang")
        self.assert_almost_equal(lif.lattice.a,  same.lattice.a)

        si = Structure.from_mpid("mp-149")
        assert si.formula == "Si2"

        # Test abiget_spginfo
        d = si.abiget_spginfo(tolsym=None, pre="abi_")
        assert d["abi_spg_symbol"] == "Fd-3m"
        assert d["abi_spg_number"] == 227
        assert d["abi_bravais"] == "Bravais cF (face-center cubic)"

        llzo = Structure.from_file(abidata.cif_file("LLZO_oxi.cif"))
        assert llzo.is_ordered
        d = llzo.abiget_spginfo(tolsym=0.001)
        assert d["spg_number"] == 142

        mgb2_cod = Structure.from_cod_id(1526507, primitive=True)
        assert mgb2_cod.formula == "Mg1 B2"
        assert mgb2_cod.spget_lattice_type() == "hexagonal"

        mgb2 = abidata.structure_from_ucell("MgB2")
        if self.has_ase():
            mgb2.abi_primitive()

        assert [site.species_string for site in mgb2.get_sorted_structure_z()] == ["B", "B", "Mg"]

        s2inds = mgb2.get_symbol2indices()
        self.assert_equal(s2inds["Mg"], [0])
        self.assert_equal(s2inds["B"], [1, 2])

        s2coords = mgb2.get_symbol2coords()
        self.assert_equal(s2coords["Mg"], [[0, 0, 0]])
        self.assert_equal(s2coords["B"],  [[1/3, 2/3, 0.5], [2/3, 1/3, 0.5]])

        # TODO: This part should be tested more carefully
        mgb2.abi_sanitize()
        mgb2.abi_sanitize(primitive_standard=True)
        mgb2.get_conventional_standard_structure()
        assert len(mgb2.abi_string)
        assert len(mgb2.spget_summary(verbose=10))
        #print(structure._repr_html_())

        self.serialize_with_pickle(mgb2)

        pseudos = abidata.pseudos("12mg.pspnc", "5b.pspnc")
        nv = mgb2.num_valence_electrons(pseudos)
        assert nv == 8 and isinstance(nv , int)
        assert mgb2.valence_electrons_per_atom(pseudos) == [2, 3, 3]
        self.assert_equal(mgb2.calc_shiftk() , [[0.0, 0.0, 0.5]])

        bmol = Structure.boxed_molecule(pseudos, cart_coords=[[0, 0, 0], [5, 5, 5]], acell=[10, 10, 10])
        self.assert_almost_equal(bmol.volume, (10 * bohr_to_ang) ** 3)

        # FIXME This is buggy
        #acell = np.array([10, 20, 30])
        #batom = Structure.boxed_atom(abidata.pseudo("12mg.pspnc"), cart_coords=[1, 2, 3], acell=acell)
        #assert isinstance(batom, Structure)
        #assert len(batom.cart_coords) == 1
        #self.assert_equal(batom.cart_coords[0], [1, 2, 3])

        # Function to compute cubic a0 from primitive v0 (depends on struct_type)
        vol2a = {"fcc": lambda vol: (4 * vol) ** (1/3.),
                 "bcc": lambda vol: (2 * vol) ** (1/3.),
                 "zincblende": lambda vol: (4 * vol) ** (1/3.),
                 "rocksalt": lambda vol: (4 * vol) ** (1/3.),
                 "ABO3": lambda vol: vol ** (1/3.),
                 "hH": lambda vol: (4 * vol) ** (1/3.),
                 }

        a = 10
        bcc_prim = Structure.bcc(a, ["Si"], primitive=True)
        assert len(bcc_prim) == 1
        self.assert_almost_equal(a, vol2a["bcc"](bcc_prim.volume))
        bcc_conv = Structure.bcc(a, ["Si"], primitive=False)
        assert len(bcc_conv) == 2
        self.assert_almost_equal(a**3, bcc_conv.volume)
        fcc_prim = Structure.fcc(a, ["Si"], primitive=True)
        assert len(fcc_prim) == 1
        self.assert_almost_equal(a, vol2a["fcc"](fcc_prim.volume))
        fcc_conv = Structure.fcc(a, ["Si"], primitive=False)
        assert len(fcc_conv) == 4
        self.assert_almost_equal(a**3, fcc_conv.volume)
        zns = Structure.zincblende(a / bohr_to_ang, ["Zn", "S"], units="bohr")
        self.assert_almost_equal(a, vol2a["zincblende"](zns.volume))
        rock = Structure.rocksalt(a, ["Na", "Cl"])
        assert len(rock) == 2
        self.assert_almost_equal(a, vol2a["rocksalt"](rock.volume))
        perov = Structure.ABO3(a, ["Ca", "Ti", "O", "O", "O"])
        assert len(perov) == 5
        self.assert_almost_equal(a**3, perov.volume)

        # Test notebook generation.
        if self.has_nbformat():
            assert mgb2.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_dataframes_from_structures(self):
        """Testing dataframes from structures."""
        mgb2 = abidata.structure_from_ucell("MgB2")
        sic = abidata.structure_from_ucell("SiC")
        alas = abidata.structure_from_ucell("AlAs")
        dfs = dataframes_from_structures([mgb2, sic, alas], index=None, with_spglib=True, cart_coords=True)
        dfs = dataframes_from_structures([mgb2, sic, alas], index=None, with_spglib=True, cart_coords=False)

        assert dfs.lattice is not None
        assert dfs.coords is not None
        assert dfs.structures is not None
        formulas = [struct.composition.reduced_formula for struct in dfs.structures]
        assert formulas == ["MgB2", "SiC", "AlAs"]

    def test_frozen_phonon_methods(self):
        """Testing frozen phonon methods (This is not a real test, just to show how to use it!)"""
        rprimd = np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
        #rprimd = rprimd*6.7468
        rprimd = rprimd * 10.60 * 0.529
        lattice = Lattice(rprimd)
        structure = Structure(lattice, ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 0.25]])
        old_structure = structure.copy()

        #print(old_structure.lattice._matrix)
        for site in old_structure:
            print(structure.lattice.get_cartesian_coords(site.frac_coords))

        # TODO: Check all this stuff more carefully
        #qpoint = [0, 0, 0]
        qpoint = [1/2, 1/2, 1/2]
        mx_sc = [2, 2, 2]
        scale_matrix = structure.get_smallest_supercell(qpoint, max_supercell=mx_sc)
        scale_matrix = 2 * np.eye(3)
        #print("Scale_matrix = ", scale_matrix)
        #scale_matrix = 2*np.eye(3)
        natoms = int(np.round(2*np.linalg.det(scale_matrix)))

        structure.write_vib_file(sys.stdout, qpoint, 0.1*np.array([[1, 1, 1], [1, 1, 1]]),
                                 do_real=True, frac_coords=False, max_supercell=mx_sc, scale_matrix=scale_matrix)

        displ = np.array([[1, 1, 1], [-1, -1, -1]])
        structure.write_vib_file(sys.stdout, qpoint, 0.1 * displ,
                                 do_real=True, frac_coords=False, max_supercell=mx_sc, scale_matrix=scale_matrix)

        structure.write_vib_file(sys.stdout, qpoint, 0.1 * displ,
                                 do_real=True, frac_coords=False, max_supercell=mx_sc, scale_matrix=None)

        fp_data = structure.frozen_phonon(qpoint, 0.1 * displ, eta=0.5, frac_coords=False,
                                          max_supercell=mx_sc, scale_matrix=scale_matrix)

        max_displ = np.linalg.norm(displ, axis=1).max()
        self.assertArrayAlmostEqual(fp_data.structure[0].coords,
                                    structure[0].coords + 0.5*displ[0]/max_displ)
        self.assertArrayAlmostEqual(fp_data.structure[8].coords,
                                    structure[1].coords + 0.5*displ[1]/max_displ)

        displ2 = np.array([[1, 0, 0], [0, 1, 1]])

        f2p_data = structure.frozen_2phonon(qpoint, 0.05 * displ, 0.02*displ2, eta=0.5, frac_coords=False,
                                           max_supercell=mx_sc, scale_matrix=scale_matrix)

        d_tot = 0.05*displ+0.02*displ2
        max_displ = np.linalg.norm(d_tot, axis=1).max()
        self.assertArrayAlmostEqual(f2p_data.structure[0].coords,
                                    structure[0].coords + 0.5*d_tot[0]/max_displ)
        self.assertArrayAlmostEqual(f2p_data.structure[8].coords,
                                    structure[1].coords + 0.5*d_tot[1]/max_displ)

        #print("Structure = ", structure)
        #print(structure.lattice._matrix)
        #for site in structure:
        #    print(structure.lattice.get_cartesian_coords(site.frac_coords))

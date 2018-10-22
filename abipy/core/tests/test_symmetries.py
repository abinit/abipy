"""Tests for symmetries module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import abipy.data as abidata

from abipy.core import Structure
from abipy.core.symmetries import LatticeRotation, AbinitSpaceGroup, mati3inv
from abipy.core.testing import AbipyTest
from abipy.abilab import abiopen


class TestSymmetries(AbipyTest):
    """"Test symmetries."""

    def test_mati3inv(self):
        """Testing mati3inv."""
        with self.assertRaises(ValueError):
            mati3inv(np.zeros((3, 3)))
        with self.assertRaises(ValueError):
            mati3inv(np.zeros(3))
        with self.assertRaises(ValueError):
            mati3inv(np.zeros((3, 3), dtype=np.int))

        mat = np.reshape([0, -1, 1, 0, -1, 0, 1, -1, 0], (3, 3))
        mat_it = mati3inv(mat)
        self.assert_equal(np.dot(mat_it.T, mat), np.eye(3, dtype=np.int))
        mat_it = mati3inv(mat, trans=False)
        self.assert_equal(np.dot(mat_it, mat), np.eye(3, dtype=np.int))

    def test_silicon(self):
        """Test silicon space group."""
        structure = Structure.from_file(abidata.ref_file("si_scf_WFK.nc"))

        assert structure.has_abi_spacegroup
        assert structure.abi_spacegroup.is_symmorphic
        spgrp = structure.abi_spacegroup
        repr(spgrp); str(spgrp)
        self.serialize_with_pickle(spgrp, test_eq=True)

        # Classes cover the entire group.
        assert sum(len(cls) for cls in spgrp.groupby_class()) == len(spgrp)

        # Multiplication table.
        mtable = spgrp.mult_table
        for i, oi in enumerate(spgrp):
            for j, oj in enumerate(spgrp):
                ij = mtable[i, j]
                assert ij is not None
                assert oi * oj == spgrp[ij]

        # Operation in the same class have the same trace and determinant.
        for cls in spgrp.groupby_class():
            #print(cls)
            op0 = cls[0]
            ref_trace, ref_det = op0.trace, op0.det
            for op in cls[1:]:
                assert op.trace == ref_trace
                assert op.det == ref_det

        assert spgrp == spgrp
        # FIXME: Temporary disabled spgid is set to 0 in the WFK file.
        #assert spgrp.spgid == 227
        assert spgrp.has_timerev
        assert len(spgrp) == 48 * 2
        assert spgrp.num_spatial_symmetries == 48

        assert spgrp.is_group()
        # TODO
        #si_symrel =
        si_tnons = np.reshape(24 * [0, 0, 0, 0.25, 0.25, 0.25], (48, 3))
        si_symafm = np.ones(48, dtype=np.int)

        self.assert_almost_equal(si_tnons, spgrp.tnons)
        self.assert_almost_equal(si_symafm, spgrp.symafm)
        assert not spgrp.afm_symmops

        for idx, symmop in enumerate(spgrp):
            repr(symmop); str(symmop)
            symmop.to_string(verbose=2)
            assert symmop in spgrp
            assert spgrp.count(symmop) == 1
            assert spgrp.find(symmop) == idx
            if symmop.det == 1: assert symmop.is_proper

        # Test pickle
        self.serialize_with_pickle(spgrp[0], protocols=None, test_eq=True)

        for idx in range(len(spgrp)-1):
            assert spgrp[idx] == spgrp[idx]
            assert spgrp[idx] != spgrp[idx+1]

        for fmop in spgrp.fm_symmops:
            assert fmop.is_fm and not fmop.is_afm

        ucell_coords = np.reshape([site.frac_coords for site in structure], (len(structure), 3))

        err_msg = ""
        for site in structure:
            for symop in spgrp:
                rot_coords = symop.rotate_r(site.frac_coords, in_ucell=True)

                for atom_coords in ucell_coords:
                    #print (atom_coords - rot_coords)
                    if np.allclose(atom_coords,  rot_coords):
                        break
                else:
                    err_msg += "Cannot find symmetrical image of %s\n" % str(rot_coords)

                assert not err_msg

        k1, k2 = [0.5, 0, 0], [0, 0.5, 0]
        ktab = spgrp.symeq(k1, k2, atol=None)
        assert ktab.isym != -1
        self.assert_equal(spgrp[ktab.isym].rotate_k(k1) - np.array(k2), ktab.g0)

        k1, k2 = [0.0, 0, 0], [0, 0.5, 0]
        ktab = spgrp.symeq(k1, k2, atol=None)
        assert ktab.isym == -1

        # Test little group with Gamma point.
        lg_gamma = spgrp.find_little_group(kpoint=[0, 0, 0])
        assert len(lg_gamma) == len(spgrp)
        repr(lg_gamma); str(lg_gamma)
        assert lg_gamma.is_symmorphic and not lg_gamma.on_bz_border
        for o1, (o2, g0) in zip(spgrp, lg_gamma.iter_symmop_g0()):
            assert o1 == o2
            assert np.all(g0 == 0)

        # Little group with X point.
        #    'G' : (0.000, 0.000, 0.000),
        #    'X' : (0.500, 0.000, 0.500),
        #    'W' : (0.500, 0.250, 0.750),
        #    'L' : (0.500, 0.500, 0.500),
        #    'K' : (0.375, 0.375, 0.750),
        #    'U' : (0.625, 0.250, 0.625)

        lg_x = spgrp.find_little_group(kpoint=[0.5, 0, 0.5])
        assert len(lg_x) == 32
        repr(lg_x); str(lg_x)
        assert lg_x.is_symmorphic and lg_x.on_bz_border

        # This is just to test from_structure but one should always try to init from file.
        other_spgroup = AbinitSpaceGroup.from_structure(structure, has_timerev=True)
        assert other_spgroup.has_timerev
        assert other_spgroup.spgid == 227
        assert len(other_spgroup) == 48 * 2


class LatticeRotationTest(AbipyTest):

    def test_base(self):
        """Testing LatticeRotation."""
        E = LatticeRotation([1, 0, 0, 0, 1, 0, 0, 0, 1])
        I = LatticeRotation([-1,  0,  0, 0, -1,  0, 0,  0, -1])

        assert E.isE and E.is_proper and E.inverse() == E
        assert E.name == "1+"
        assert I.isI and not I.is_proper and I.inverse() == I
        assert I.name == "-2-"

        # Test Basic operations
        assert E != I
        assert +E == E
        assert -I == E
        assert E * I == I
        assert E ** 2 == E
        assert E ** -1 == E
        assert I ** 0 == E
        assert I ** -1 == I
        assert I ** 3 == I
        assert I ** -3 == I

        # Test pickle.
        self.serialize_with_pickle([E, I])


class BilbaoPointGroupTest(AbipyTest):

    def test_database(self):
        """Testing BilbaoPointGroup database."""
        from abipy.core.symmetries import bilbao_ptgroup, sch_symbols
        for sch_symbol in sch_symbols:
            print(sch_symbol)
            ptg = bilbao_ptgroup(sch_symbol)
            repr(ptg); str(ptg)
            assert len(ptg.to_string())
            #for irrep_name in ptg.irrep_names: ptg.show_irrep(irrep_name)
            assert ptg.auto_test() == 0


class LittleGroupTest(AbipyTest):

    def test_silicon_little_group(self):
        """Testing little group in Silicon."""
        with abiopen(abidata.ref_file("si_scf_WFK.nc")) as wfk_file:

            spgrp = wfk_file.structure.abi_spacegroup
            assert spgrp is not None
            repr(spgrp); str(spgrp)

            kpoints = [[0,0,0],
                       [0.5, 0, 0],
                       [1/3, 1/3, 1/3],
                       [1/4,1/4,0],
                      ]

            for ik, kpoint in enumerate(kpoints):
                ltk = spgrp.find_little_group(kpoint)
                repr(ltk); str(ltk)
                #if ik == 1:
                #    assert ltk.onborder_and_nonsymmorphic
                #else:
                #    assert not ltk.onborder_and_nonsymmorphic

                #wfk_file.classify_ebands(0, kpoint, bands_range=range(0,5))


# reduced_symmetry_matrices =
#  1, 0, 0,
#  0, 1, 0,
#  0, 0, 1,
#  -1, 0, 0,
#  0, -1, 0,
#  0, 0, -1,
#  0, -1, 1,
#  0, -1, 0,
#  1, -1, 0,
#  0, 1, -1,
#  0, 1, 0,
#  -1, 1, 0,
#  -1, 0, 0,
#  -1, 0, 1,
#  -1, 1, 0,
#  1, 0, 0,
#  1, 0, -1,
#  1, -1, 0,
#  0, 1, -1,
#  1, 0, -1,
#  0, 0, -1,
#  0, -1, 1,
#  -1, 0, 1,
#  0, 0, 1,
#  -1, 0, 0,
#  -1, 1, 0,
#  -1, 0, 1,
#  1, 0, 0,
#  1, -1, 0,
#  1, 0, -1,
#  0, -1, 1,
#  1, -1, 0,
#  0, -1, 0,
#  0, 1, -1,
#  -1, 1, 0,
#  0, 1, 0,
#  1, 0, 0,
#  0, 0, 1,
#  0, 1, 0,
#  -1, 0, 0,
#  0, 0, -1,
#  0, -1, 0,
#  0, 1, -1,
#  0, 0, -1,
#  1, 0, -1,
#  0, -1, 1,
#  0, 0, 1,
#  -1, 0, 1,
#  -1, 0, 1,
#  -1, 1, 0,
#  -1, 0, 0,
#  1, 0, -1,
#  1, -1, 0,
#  1, 0, 0,
#  0, -1, 0,
#  1, -1, 0,
#  0, -1, 1,
#  0, 1, 0,
#  -1, 1, 0,
#  0, 1, -1,
#  1, 0, -1,
#  0, 0, -1,
#  0, 1, -1,
#  -1, 0, 1,
#  0, 0, 1,
#  0, -1, 1,
#  0, 1, 0,
#  0, 0, 1,
#  1, 0, 0,
#  0, -1, 0,
#  0, 0, -1,
#  -1, 0, 0,
#  1, 0, -1,
#  0, 1, -1,
#  0, 0, -1,
#  -1, 0, 1,
#  0, -1, 1,
#  0, 0, 1,
#  0, -1, 0,
#  0, -1, 1,
#  1, -1, 0,
#  0, 1, 0,
#  0, 1, -1,
#  -1, 1, 0,
#  -1, 0, 1,
#  -1, 0, 0,
#  -1, 1, 0,
#  1, 0, -1,
#  1, 0, 0,
#  1, -1, 0,
#  0, 1, 0,
#  1, 0, 0,
#  0, 0, 1,
#  0, -1, 0,
#  -1, 0, 0,
#  0, 0, -1,
#  0, 0, -1,
#  0, 1, -1,
#  1, 0, -1,
#  0, 0, 1,
#  0, -1, 1,
#  -1, 0, 1,
#  1, -1, 0,
#  0, -1, 1,
#  0, -1, 0,
#  -1, 1, 0,
#  0, 1, -1,
#  0, 1, 0,
#  0, 0, 1,
#  1, 0, 0,
#  0, 1, 0,
#  0, 0, -1,
#  -1, 0, 0,
#  0, -1, 0,
#  -1, 1, 0,
#  -1, 0, 0,
#  -1, 0, 1,
#  1, -1, 0,
#  1, 0, 0,
#  1, 0, -1,
#  0, 0, 1,
#  0, 1, 0,
#  1, 0, 0,
#  0, 0, -1,
#  0, -1, 0,
#  -1, 0, 0,
#  1, -1, 0,
#  0, -1, 0,
#  0, -1, 1,
#  -1, 1, 0,
#  0, 1, 0,
#  0, 1, -1,
#  0, 0, -1,
#  1, 0, -1,
#  0, 1, -1,
#  0, 0, 1,
#  -1, 0, 1,
#  0, -1, 1,
#  -1, 1, 0,
#  -1, 0, 1,
#  -1, 0, 0,
#  1, -1, 0,
#  1, 0, -1,
#  1, 0, 0 ;

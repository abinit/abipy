"""Tests for symmetries module"""
import numpy as np
import abipy.data as data

from abipy.core import Structure
from abipy.core.symmetries import *
from abipy.core.testing import *

class TestSymmetries(AbipyTest):
    """"Test symmetries."""

    def test_silicon(self):
        """Test silicon space group."""
        structure = Structure.from_file(data.ref_file("si_nscf_WFK-etsf.nc"))

        spgrp = structure.spacegroup
        print(spgrp)

        self.assertTrue(spgrp.spgid == 227)
        self.assertTrue(spgrp.has_timerev)
        self.assertTrue(len(spgrp) == 48 * 2)
        self.assertTrue(spgrp.num_spatial_symmetries == 48)

        self.assertTrue(spgrp.is_group())
        # TODO
        #si_symrel = 
        si_tnons = np.reshape(24 * [0, 0, 0, 0.25, 0.25, 0.25], (48, 3))
        si_symafm = np.ones(48, dtype=np.int)

        self.assert_almost_equal(si_tnons, spgrp.tnons)
        self.assert_almost_equal(si_symafm, spgrp.symafm)

        for (idx, symmop) in enumerate(spgrp):
            self.assertTrue(symmop in spgrp)
            self.assertTrue(spgrp.count(symmop) == 1)
            self.assertTrue(spgrp.find(symmop) == idx)
            self.assertTrue(abs(symmop.det) == 1)

        for idx in range(len(spgrp)-1):
            self.assertTrue(spgrp[idx] == spgrp[idx])
            self.assertTrue(spgrp[idx] != spgrp[idx+1])

        for fmop in spgrp.fm_symmops:
            self.assertTrue(fmop.is_fm)

        ucell_coords = np.reshape([site.frac_coords for site in structure], (len(structure), 3))

        for site in structure:
            for symop in spgrp:
                rot_coords = symop.rotate_r(site.frac_coords, in_ucell=True)

                found = False
                for atom_coords in ucell_coords:
                    #print (atom_coords - rot_coords)
                    if np.allclose(atom_coords,  rot_coords):
                        found = True
                        break

                self.assertTrue(found)

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

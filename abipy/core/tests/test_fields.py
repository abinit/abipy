#!/usr/bin/env python
"""Tests for core.field module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import abipy.data as abidata

from abipy.core.fields import *
from abipy.core.testing import AbipyTest

from abipy.core import Density
from abipy.iotools import *


class TestScalarField(AbipyTest):
    """Unit tests for ScalarField."""

    def test_base_class(self):
        """Testing ScalarField base class."""
        structure = abidata.structure_from_ucell("Si")
        nspinor, nsppol, nspden = 1, 1, 1
        xyz_shape = (2, 3, 4)
        datar = np.zeros((nsppol,) + xyz_shape)

        field = ScalarField(nspinor, nsppol, nspden, datar, structure, iorder="c")

        repr(field); str(field)
        assert field.datar.ndim ==  4
        assert len(field) == nsppol
        assert field.shape == datar.shape
        assert datar.shape == (field.nspden, field.nx, field.ny, field.nz)
        assert field.nx == 2 and field.ny == 3 and field.nz == 4
        assert field.nz == field.mesh.nz
        assert field.is_collinear

        self.assert_equal(field.mean(space="r"), field.datar.mean(axis=0))
        self.assert_equal(field.mean(space="g"), field.datag.mean(axis=0))
        self.assert_equal(field.std(space="r"), field.datar.std(axis=0))
        self.assert_equal(field.std(space="g"), field.datag.std(axis=0))

        field.export(self.get_tmpname(text=True, suffix=".xsf"))
        visu = field.visualize("vesta")
        assert callable(visu)

        #assert field.datar_xyz.ndim == 4
        #assert field.datar_xyz.shape[-3:], xyz_shape)
        other = field + field
        assert other.nspden == field.nspden  and np.all(other.datar == 2 * field.datar)
        other = other - field
        assert other.nspden == field.nspden
        self.assert_almost_equal(other.datar, field.datar)
        self.assert_almost_equal(other.datag, field.datag)
        other = -field
        assert other.structure == field.structure
        self.assert_almost_equal(other.datar, -field.datar)

        with self.assertRaises(TypeError):
            field + "hello"

        nspinor_ncoll, nsppol_ncoll, nspden_ncoll = 1, 1, 4
        datar_ncoll = np.zeros((nspden_ncoll,) + xyz_shape)
        field_ncoll = ScalarField(nspinor_ncoll, nsppol_ncoll, nspden_ncoll, datar_ncoll, structure, iorder="c")
        with self.assertRaises(ValueError):
            field + field_ncoll

        with self.assertRaises(ValueError):
            field._check_space("foo")

    def test_silicon_density(self):
        """Testing density object (spin unpolarized)."""
        si_den = Density.from_file(abidata.ref_file("si_DEN.nc"))
        repr(si_den); str(si_den)
        assert si_den.nspinor == 1 and si_den.nsppol == 1 and si_den.nspden == 1
        assert si_den.is_collinear
        assert si_den.structure.formula == "Si2"
        assert si_den.mesh.shape == (18, 18, 18)

        # Read data directly from file.
        with ETSF_Reader(abidata.ref_file("si_DEN.nc")) as r:
            nelect_file = r.read_value("number_of_electrons")

        ne = 8
        assert ne == nelect_file
        self.assert_almost_equal(si_den.get_nelect(), ne)
        self.assert_almost_equal(si_den.total_rhor.sum() * si_den.structure.volume / si_den.mesh.size, ne)
        self.assert_almost_equal(si_den.total_rhog[0, 0, 0] * si_den.structure.volume, ne)

        totden = si_den.total_rhor_as_density()
        self.assert_equal(totden.datar.flatten(), si_den.total_rhor.flatten())

        other = si_den - si_den
        assert other.nspden == si_den.nspden
        self.assert_equal(other.datar, 0)

        magfield = si_den.magnetization_field
        assert magfield.shape == si_den.mesh.shape
        self.assert_equal(magfield, 0)
        assert si_den.magnetization == 0

        nup, ndown = si_den.nelect_updown
        assert nup == ndown
        self.assert_almost_equal(nup, ne / 2)
        self.assert_almost_equal(si_den.zeta, 0)

        # Export to chgcar and re-read it.
        chgcar_path = self.get_tmpname(text=True)
        chgcar = si_den.to_chgcar(filename=chgcar_path)
        assert hasattr(chgcar, "structure")
        assert not chgcar.is_spin_polarized

        poscar_path = self.get_tmpname(text=True)
        si_den.structure.to(fmt="poscar", filename=poscar_path)

        same_den = Density.from_chgcar_poscar(chgcar_path, poscar_path)
        self.assert_almost_equal(same_den.datar, si_den.datar)

        # Export data in xsf format.
        visu = si_den.export(".xsf")
        assert callable(visu)

        # Test CUBE files.
        tmp_cubefile = self.get_tmpname(text=True)
        si_den.export_to_cube(tmp_cubefile, spin="total")
        total_den = Density.from_cube(tmp_cubefile)

        assert total_den.structure == si_den.structure
        assert abs(total_den.get_nelect().sum() - ne) < 1e-3

        # TODO@DW
        #den.ae_core_density_on_mesh(cls, valence_density, structure, rhoc_files, maxr=2.0, nelec=None,
        #                        method='mesh3d_dist_gridpoints', small_dist_mesh=(8, 8, 8), small_dist_factor=1.5):

    def test_ni_density(self):
        """Testing density object (spin polarized, collinear)."""
        ni_den = Density.from_file(abidata.ref_file("ni_666k_DEN.nc"))

        repr(ni_den); str(ni_den)
        assert ni_den.nspinor == 1 and ni_den.nsppol == 2 and ni_den.nspden == 2
        assert ni_den.is_collinear
        assert ni_den.structure.formula == "Ni1"
        assert ni_den.mesh.shape == (27, 27, 27)
        ne = 18
        self.assert_almost_equal(ni_den.get_nelect(), ne)
        self.assert_almost_equal(ni_den.total_rhor.sum() * ni_den.structure.volume / ni_den.mesh.size, ne)
        self.assert_almost_equal(ni_den.total_rhog[0, 0, 0] * ni_den.structure.volume, ne)

        totden = ni_den.total_rhor_as_density()
        self.assert_equal(totden.datar.flatten(), ni_den.total_rhor.flatten())

        other = ni_den - ni_den
        assert other.nspden == ni_den.nspden
        self.assert_equal(other.datar, 0)

        magfield = ni_den.magnetization_field
        assert magfield.shape == ni_den.mesh.shape
        self.assert_equal(magfield, ni_den.datar[0] - ni_den.datar[1])
        self.assert_almost_equal(ni_den.magnetization, 0.650144, decimal=5)

        nup, ndown = ni_den.nelect_updown
        self.assert_almost_equal(nup, 9.32507195)
        self.assert_almost_equal(ndown, 8.674928)
        self.assert_almost_equal(ni_den.zeta[3, 3, 3], 0.31311881970587324)

        # Export data in xsf format.
        visu = ni_den.export(".xsf")
        assert callable(visu)

        # Test CUBE files.
        tmp_cubefile = self.get_tmpname(text=True)
        ni_den.export_to_cube(tmp_cubefile, spin="total")
        total_den = Density.from_cube(tmp_cubefile)

        assert total_den.structure == ni_den.structure
        assert abs(total_den.get_nelect().sum() - ne) < 1e-3

        # Export to chgcar and re-read it.
        chgcar_path = self.get_tmpname(text=True)
        chgcar = ni_den.to_chgcar(filename=chgcar_path)
        assert hasattr(chgcar, "structure")
        assert chgcar.is_spin_polarized

        poscar_path = self.get_tmpname(text=True)
        ni_den.structure.to(fmt="poscar", filename=poscar_path)

        same_den = Density.from_chgcar_poscar(chgcar_path, poscar_path)
        assert same_den.nspinor == 1 and same_den.nsppol == 2 and same_den.nspden == 2
        self.assert_almost_equal(same_den.datar, ni_den.datar, decimal=5)

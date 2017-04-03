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

        print(field)
        assert field.datar.ndim ==  4
        assert len(field) == nsppol
        assert field.shape == datar.shape
        assert datar.shape == (field.nspden, field.nx, field.ny, field.nz)
        assert field.nx == 2 and field.ny == 3 and field.nz == 4
        assert field.nz == field.mesh.nz
        assert field.is_collinear

        assert np.all(field.mean(space="r") == field.datar.mean(axis=0))
        assert np.all(field.mean(space="g") == field.datag.mean(axis=0))

        assert np.all(field.std(space="r") == field.datar.std(axis=0))
        assert np.all(field.std(space="g") == field.datag.std(axis=0))

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

    def test_silicon_density(self):
        """Testing density object."""
        density = Density.from_file(abidata.ref_file("si_DEN.nc"))
        repr(density)
        str(density)
        assert density.nspinor == 1 and density.nsppol == 1 and density.nspden == 1
        assert density.is_collinear
        assert density.structure.formula == "Si2"
        self.assert_almost_equal(density.get_nelect(), 8)
        self.assert_almost_equal(density.total_rhor.sum() * density.structure.volume / density.mesh.size, 8)
        self.assert_almost_equal(density.total_rhog[0, 0, 0] * density.structure.volume, 8)

        totden = density.total_rhor_as_density()
        self.assert_equal(totden.datar.flatten(), density.total_rhor.flatten())

        other = density - density
        assert other.nspden == density.nspden
        self.assert_equal(other.datar, 0)

        magfield = density.magnetization_field
        assert magfield.shape == density.mesh.shape
        self.assert_equal(magfield, 0)
        assert density.magnetization == 0

        nup, ndown = density.nelect_updown
        assert nup == ndown
        self.assert_almost_equal(nup, 8 / 2)
        self.assert_almost_equal(density.zeta, 0)

        # Export to chgcar
        chgcar_path = self.get_tmpname(text=True)
        chgcar = density.to_chgcar(filename=chgcar_path)
        assert hasattr(chgcar, "structure")
        assert not chgcar.is_spin_polarized

        #same_density = Density.from_chgcar_file(chgcar_path)
        #self.assert_almost_equal(same_density.datar, density.datar)

        # TODO@DW
        #den.ae_core_density_on_mesh(cls, valence_density, structure, rhoc_files, maxr=2.0, nelec=None,
        #                        method='mesh3d_dist_gridpoints', small_dist_mesh=(8, 8, 8), small_dist_factor=1.5):

    def test_ncread_density(self):
        """Read density from NC example data files"""
        assert abidata.DEN_NCFILES

        # This section is mainly used to spot possible problems with nspden, nsppol, nspinor ....
        for path in abidata.DEN_NCFILES:
            print("Reading DEN file %s " % path)

            # Read data directly from file.
            with ETSF_Reader(path) as r:
                nelect_file = r.read_value("number_of_electrons")

            # Compute nelect from data.
            den = Density.from_file(path)
            print(den)
            structure = den.structure
            rhor_tot = den.total_rhor
            rhog_tot = den.total_rhog
            nelect_calc = den.get_nelect().sum()

            # Diff between nelect computed and the one written on file.
            self.assert_almost_equal(nelect_calc, nelect_file)
            self.assert_almost_equal(rhog_tot[0,0,0] * structure.volume, nelect_file)

            if self.which("xcrysden") is not None:
                # Export data in xsf format.
                visu = den.export(".xsf")
                assert callable(visu)

            # Test CUBE files.
            tmp_cubefile = self.get_tmpname(text=True)
            den.export_to_cube(tmp_cubefile)
            total_den = Density.from_cube(tmp_cubefile)

            assert total_den.structure == den.structure
            assert abs(total_den.get_nelect().sum() - nelect_file) < 1e-3

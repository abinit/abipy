#!/usr/bin/env python
"""Tests for core.field module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import os
import pymatgen.core.units as pmgu
import abipy.data as abidata

from pymatgen.core.units import bohr_to_angstrom
from abipy.core.fields import _Field, FieldReader, Density, VxcPotential, VhartreePotential, VhxcPotential
from abipy.core.fields import core_density_from_file
from abipy.core.testing import AbipyTest
from abipy.iotools import *


class TestScalarField(AbipyTest):
    """Unit tests for _Field."""

    def test_base_class(self):
        """Testing _Field base class."""
        structure = abidata.structure_from_ucell("Si")
        nspinor, nsppol, nspden = 1, 1, 1
        xyz_shape = (2, 3, 4)
        datar = np.random.rand(*((nsppol,) + xyz_shape))

        field = _Field(nspinor, nsppol, nspden, datar, structure, iorder="c")

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
        #visu = field.visualize(appname="vesta")
        #assert callable(visu)

        # Field "algebra"
        #assert field.datar_xyz.ndim == 4
        #assert field.datar_xyz.shape[-3:], xyz_shape)
        other = field + field
        assert other.nspden == field.nspden  and np.all(other.datar == 2 * field.datar)
        other = other - field
        assert other.nspden == field.nspden
        self.assert_almost_equal(other.datar, field.datar)
        self.assert_almost_equal(other.datag, field.datag)
        other = -field
        self.assert_equal(other.datar, -field.datar)
        self.assert_equal(other.datag, -field.datag)
        assert other.structure == field.structure

        other = field * np.pi
        self.assert_equal(other.datar, field.datar * np.pi)
        self.assert_almost_equal(other.datag, field.datag * np.pi)
        self.assert_equal((np.pi * field).datar, other.datar)

        other = - field / np.pi
        self.assert_equal(other.datar, - field.datar / np.pi)
        other = abs(field)
        self.assert_equal(other.datar, abs(field.datar))
        #other = np.pi / field
        #self.assert_equal(other.datar, np.pi / field.datar)

        other = field + 3
        other = (field / field)
        self.assert_equal(other.datar, 1.0)

        with self.assertRaises(TypeError):
            field + "hello"

        nspinor_ncoll, nsppol_ncoll, nspden_ncoll = 1, 1, 4
        datar_ncoll = np.zeros((nspden_ncoll,) + xyz_shape)
        field_ncoll = _Field(nspinor_ncoll, nsppol_ncoll, nspden_ncoll, datar_ncoll, structure, iorder="c")
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
        assert si_den.is_density_like
        assert not si_den.is_potential_like

        self.assert_almost_equal(si_den.mesh.fft_g2r(si_den.datag), si_den.datar)

        # Read data directly from file.
        with FieldReader(abidata.ref_file("si_DEN.nc")) as r:
            nelect_file = r.read_value("number_of_electrons")
            same_si_den = r.read_field()
            assert np.all(same_si_den.datar == si_den.datar)

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

        df = si_den.integrate_in_spheres(rcut_symbol=None, out=True)
        assert "frac_coords" in df
        self.assert_almost_equal(df["ntot"].values, 2 * [2.010537])
        self.assert_almost_equal(df["rsph_ang"].values, 2 * [1.11])
        df = si_den.integrate_in_spheres(rcut_symbol=2, out=False)

        if self.has_matplotlib():
            assert si_den.plot_line(0, 1, num=1000, show=False)
            assert si_den.plot_line([0, 0, 0], [1, 0, 0], num=1000, cartesian=True, show=False)

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

        # Test creation of AE core density. Use low parameters to reduce time
        rhoc = {"Si": core_density_from_file(os.path.join(abidata.pseudo_dir, "Si.fc"))}
        core_den_1 = Density.ae_core_density_on_mesh(si_den, si_den.structure, rhoc, maxr=1.5,
                                                     method='get_sites_in_sphere', small_dist_mesh=(6, 6, 6))
        core_den_2 = Density.ae_core_density_on_mesh(si_den, si_den.structure, rhoc, maxr=1.5,
                                                     method='mesh3d_dist_gridpoints', small_dist_mesh=(6, 6, 6))
        self.assertAlmostEqual(np.sum(core_den_1.datar) * si_den.mesh.dv, 20, delta=0.5)
        self.assertArrayAlmostEqual(core_den_1.datar, core_den_2.datar)
        with self.assertRaises(ValueError):
            Density.ae_core_density_on_mesh(si_den, si_den.structure, rhoc, maxr=1, nelec=20, tol=0.001,
                                            method='get_sites_in_sphere', small_dist_mesh=(2,2,2))


    def test_ni_density(self):
        """Testing density object (spin polarized, collinear)."""
        ni_den = Density.from_file(abidata.ref_file("ni_666k_DEN.nc"))
        repr(ni_den); str(ni_den)
        ni_den.to_string(verbose=1)

        self.assert_almost_equal(ni_den.mesh.vectors, pmgu.bohr_to_angstrom *
            np.reshape([0.0000000, 3.3259180, 3.3259180,
                        3.3259180, 0.0000000, 3.3259180,
                        3.3259180, 3.3259180, 0.0000000], (3, 3)))

        assert ni_den.nspinor == 1 and ni_den.nsppol == 2 and ni_den.nspden == 2
        assert ni_den.is_collinear
        assert ni_den.structure.formula == "Ni1"
        assert ni_den.mesh.shape == (27, 27, 27)
        assert ni_den.is_density_like
        assert not ni_den.is_potential_like
        ne = 18
        self.assert_almost_equal(ni_den.get_nelect(), ne)
        self.assert_almost_equal(ni_den.total_rhor.sum() * ni_den.structure.volume / ni_den.mesh.size, ne)
        self.assert_almost_equal(ni_den.total_rhog[0, 0, 0] * ni_den.structure.volume, ne)

        totden = ni_den.total_rhor_as_density()
        self.assert_equal(totden.datar.flatten(), ni_den.total_rhor.flatten())
        self.assert_almost_equal(ni_den.mesh.fft_g2r(ni_den.datag), ni_den.datar)

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

        df = ni_den.integrate_in_spheres(rcut_symbol=None, out=True)
        assert "frac_coords" in df
        self.assert_almost_equal(df["nup"].values, 8.9778673)
        self.assert_almost_equal(df["ndown"].values, 8.2989892)
        self.assert_almost_equal(df["mz"].values, 0.6788782)
        self.assert_almost_equal(df["rsph_ang"].values, 1.24)
        self.assert_almost_equal(df["nup"].values + df["ndown"].values, df["ntot"].values)

        if self.has_matplotlib():
            assert ni_den.plot_line([0, 0, 0],  [1, 1, 1], num=1000, show=False)
            assert ni_den.plot_line_neighbors(site_index=0, radius=1, num=50, max_nn=10) is None
            assert ni_den.plot_line_neighbors(site_index=0, radius=3, num=50, max_nn=10, show=False)

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

    def test_potentials(self):
        """Testing Potentials."""
        vxc = VxcPotential.from_file(abidata.ref_file("ni_666k_VXC.nc"))
        assert vxc.nsppol == 2 and vxc.nspden == 2
        assert vxc.is_collinear
        assert not vxc.is_density_like
        assert vxc.is_potential_like
        assert vxc.datar.dtype == np.float
        fact = pmgu.Ha_to_eV / pmgu.bohr_to_angstrom ** 3
        self.assert_almost_equal(vxc.datar[0, 0, 0, 0], -2.40411892342838 * fact)
        self.assert_almost_equal(vxc.datar[0, 0, 0, 1], -2.31753083824603 * fact)

        if self.has_matplotlib():
            assert vxc.plot_line([0, 0, 0],  [1, 1, 1], num=1000, show=False)

        vh = VhartreePotential.from_file(abidata.ref_file("ni_666k_VHA.nc"))
        assert vh.nsppol == 2 and vh.nspden == 2
        assert vh.is_collinear
        assert not vh.is_density_like
        assert vh.is_potential_like

        vhxc = VhxcPotential.from_file(abidata.ref_file("ni_666k_VHXC.nc"))
        assert vhxc.nsppol == 2 and vhxc.nspden == 2
        assert vhxc.is_collinear
        assert not vhxc.is_density_like
        assert vhxc.is_potential_like

        # VH + VXC = VHXC
        self.assert_almost_equal((vh + 1.0 * vxc).datar, vhxc.datar)

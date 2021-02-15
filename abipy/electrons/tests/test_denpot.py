# coding: utf-8
"""Tests for DEN/POT files."""
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.electrons.denpot import DensityNcFile


class DensityNcFileTest(AbipyTest):

    def test_silicon_density_ncfile(self):
        """Testing si_DEN.nc"""
        with DensityNcFile(abidata.ref_file("si_DEN.nc")) as denc:
            repr(denc); str(denc)
            assert denc.structure.formula == "Si2"
            assert denc.ebands.structure == denc.structure
            assert str(denc.xc) == "LDA_XC_TETER93"
            assert denc.params["nsppol"] == 1

            assert np.all(denc.reader.read_ngfft3() == np.array([18, 18, 18]))

            # Kpoint sampling
            # kptopt1 1 ngkpt1 8 8 8 nshiftk1 1 shiftk1   0 0 0
            assert denc.ebands.kpoints.is_mpmesh
            assert not denc.ebands.kpoints.is_path

            ksamp = denc.ebands.kpoints.ksampling
            str(ksamp); repr(ksamp)

            self.assert_equal(ksamp.mpdivs, [8, 8, 8])
            self.assert_equal(ksamp.kptrlatt_orig, 8 * np.eye(3, 3))
            self.assert_equal(ksamp.shifts_orig, 0)

            self.assert_equal(ksamp.kptrlatt, 8 * np.eye(3, 3))
            self.assert_almost_equal(ksamp.shifts.flatten(), [0.0, 0.0, 0.0])
            assert ksamp.kptopt == 1
            assert ksamp.is_mesh

            density = denc.density
            assert density.nspinor == 1 and density.nsppol == 1 and density.nspden == 1
            self.assert_almost_equal(density.get_nelect(), 8)

            # Test converters.
            denc.write_chgcar(filename=self.get_tmpname(text=True))
            denc.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            denc.write_cube(filename=self.get_tmpname(text=True), spin="total")

            if self.has_matplotlib():
                assert denc.ebands.plot(show=False)

            # Test ipython notebooks.
            if self.has_nbformat():
                denc.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_nickel_density_ncfile(self):
        """Testing ni_666k_DEN.nc"""
        with abilab.abiopen(abidata.ref_file("ni_666k_DEN.nc")) as denc:
            repr(denc); str(denc)
            assert denc.to_string(verbose=2)
            assert denc.structure.formula == "Ni1"
            assert denc.ebands.structure == denc.structure
            assert str(denc.xc) == "PBE"

            # Kpoint sampling
            # ngkpt 6 6 6 nshiftk  4
            # shiftk   1/2 1/2 1/2 1/2 0.0 0.0 0.0 1/2 0.0 0.0 0.0 1/2
            assert not denc.ebands.kpoints.is_mpmesh
            assert not denc.ebands.kpoints.is_path
            ksamp = denc.ebands.kpoints.ksampling
            repr(ksamp); str(ksamp)
            assert ksamp.mpdivs is None

            self.assert_equal(ksamp.kptrlatt_orig, 6 * np.eye(3, 3))
            self.assert_equal(ksamp.shifts_orig.flatten(), [
                1/2, 1/2, 1/2, 1/2, 0.0, 0.0, 0.0, 1/2, 0.0, 0.0, 0.0, 1/2])

            self.assert_equal(ksamp.kptrlatt.flatten(), [6, -6, 6, -6, 6, 6, -6, -6, 6])
            self.assert_almost_equal(ksamp.shifts.flatten(), [0.5, 0.5, 0.5])

            assert ksamp.kptopt == 1
            assert ksamp.is_mesh

            density = denc.density
            assert density.nspinor == 1 and density.nsppol == 2 and density.nspden == 2
            assert density.is_density_like
            assert not density.is_potential_like
            self.assert_almost_equal(density.get_nelect(), 18)

            # Test converters.
            denc.write_chgcar(filename=self.get_tmpname(text=True))
            denc.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            denc.write_cube(filename=self.get_tmpname(text=True), spin="total")

            if self.has_matplotlib():
                assert denc.ebands.plot(show=False)

            # Test ipython notebooks.
            if self.has_nbformat():
                denc.write_notebook(nbpath=self.get_tmpname(text=True))


class VxcNcFileTest(AbipyTest):

    def test_vhartree_ncfile(self):
        """Testing VXC netcdf file."""
        with abilab.abiopen(abidata.ref_file("ni_666k_VHA.nc")) as ncfile:
            repr(ncfile); str(ncfile)
            assert ncfile.ebands.nsppol == 2 and ncfile.ebands.nspden == 2
            assert ncfile.structure.formula == "Ni1"
            assert ncfile.ebands.structure == ncfile.structure
            assert str(ncfile.xc) == "PBE"

            # Kpoint sampling
            #ngkpt   6 6 6  nshiftk  4
            #shiftk  1/2 1/2 1/2 1/2 0.0 0.0 0.0 1/2 0.0 0.0 0.0 1/2
            assert not ncfile.ebands.kpoints.is_mpmesh
            assert not ncfile.ebands.kpoints.is_path

            vh = ncfile.vh
            assert vh.netcdf_name == "vhartree"
            assert vh.nspinor == 1 and vh.nsppol == 2 and vh.nspden == 2
            assert vh.is_potential_like

            ## Test converters.
            #ncfile.write_chgcar(filename=self.get_tmpname(text=True))
            #ncfile.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            #ncfile.write_cube(filename=self.get_tmpname(text=True), spin="total")

            if self.has_matplotlib():
                assert ncfile.ebands.plot(show=False)

            # Test ipython notebooks.
            if self.has_nbformat():
                ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_vxc_ncfile(self):
        """Testing VXC netcdf file."""
        with abilab.abiopen(abidata.ref_file("ni_666k_VXC.nc")) as ncfile:
            repr(ncfile); str(ncfile)
            assert ncfile.ebands.nsppol == 2 and ncfile.ebands.nspden == 2
            assert ncfile.structure.formula == "Ni1"
            assert ncfile.ebands.structure == ncfile.structure
            assert str(ncfile.xc) == "PBE"

            # Kpoint sampling
            #ngkpt   6 6 6  nshiftk  4
            #shiftk  1/2 1/2 1/2 1/2 0.0 0.0 0.0 1/2 0.0 0.0 0.0 1/2
            assert not ncfile.ebands.kpoints.is_mpmesh
            assert not ncfile.ebands.kpoints.is_path

            vxc = ncfile.vxc
            assert vxc.netcdf_name == "exchange_correlation_potential"
            assert vxc.nspinor == 1 and vxc.nsppol == 2 and vxc.nspden == 2
            assert vxc.is_potential_like
            assert not vxc.is_density_like

            ## Test converters.
            #ncfile.write_chgcar(filename=self.get_tmpname(text=True))
            #ncfile.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            #ncfile.write_cube(filename=self.get_tmpname(text=True), spin="total")

            if self.has_matplotlib():
                assert ncfile.ebands.plot(show=False)

            # Test ipython notebooks.
            if self.has_nbformat():
                ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_vhxc_ncfile(self):
        """Testing VHXC netcdf file."""
        with abilab.abiopen(abidata.ref_file("ni_666k_VHXC.nc")) as ncfile:
            repr(ncfile); str(ncfile)
            assert ncfile.ebands.nsppol == 2 and ncfile.ebands.nspden == 2
            assert ncfile.structure.formula == "Ni1"
            assert ncfile.ebands.structure == ncfile.structure
            assert str(ncfile.xc) == "PBE"

            # Kpoint sampling
            #ngkpt   6 6 6  nshiftk  4
            #shiftk  1/2 1/2 1/2 1/2 0.0 0.0 0.0 1/2 0.0 0.0 0.0 1/2
            assert not ncfile.ebands.kpoints.is_mpmesh
            assert not ncfile.ebands.kpoints.is_path

            vhxc = ncfile.vhxc
            assert vhxc.netcdf_name == "vhxc"
            assert vhxc.netcdf_name == "vhxc"
            assert vhxc.nspinor == 1 and vhxc.nsppol == 2 and vhxc.nspden == 2
            assert vhxc.is_potential_like
            assert not vhxc.is_density_like

            ## Test converters.
            #ncfile.write_chgcar(filename=self.get_tmpname(text=True))
            #ncfile.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            #ncfile.write_cube(filename=self.get_tmpname(text=True), spin="total")

            if self.has_matplotlib():
                assert ncfile.ebands.plot(show=False)

            # Test ipython notebooks.
            if self.has_nbformat():
                ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_vpot_ncfile(self):
        """Testing POT netcdf file."""
        with abilab.abiopen(abidata.ref_file("ni_666k_POT.nc")) as ncfile:
            repr(ncfile); str(ncfile)
            assert ncfile.ebands.nsppol == 2 and ncfile.ebands.nspden == 2
            assert ncfile.structure.formula == "Ni1"
            assert ncfile.ebands.structure == ncfile.structure
            assert str(ncfile.xc) == "PBE"

            vks = ncfile.vks
            assert vks.netcdf_name == "vtrial"
            assert vks.nspinor == 1 and vks.nsppol == 2 and vks.nspden == 2
            assert vks.is_potential_like
            assert not vks.is_density_like

            ## Test converters.
            #ncfile.write_chgcar(filename=self.get_tmpname(text=True))
            #ncfile.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            #ncfile.write_cube(filename=self.get_tmpname(text=True), spin="total")

            if self.has_matplotlib():
                assert ncfile.ebands.plot(show=False)

            # Test ipython notebooks.
            if self.has_nbformat():
                ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

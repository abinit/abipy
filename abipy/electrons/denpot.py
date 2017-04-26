# coding: utf-8
"""Density/potential files in netcdf/fortran format."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from monty.string import marquee
from monty.termcolor import cprint
from monty.functools import lazy_property
from abipy.core.mixins import (AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter,
    AbinitFortranFile, CubeFile)
from abipy.core.fields import FieldReader
from abipy.electrons.ebands import ElectronsReader
#from abipy.tools import duck
#from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "DensityNcFile",
    "VhartreeNcFile",
    "VxcNcFile",
    "VhxcNcFile",
    "PotNcFile",
]


class _DenPotNcReader(ElectronsReader, FieldReader):
    """Object used to read data from density/potential files in netcdf format."""


class _NcFileWithField(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    Base class providing commong methods for netcdf files with density/potential
    """
    _field_name = None

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(_NcFileWithField, self).__init__(filepath)
        self.reader = _DenPotNcReader(filepath)

    @lazy_property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    @lazy_property
    def xc(self):
        """:class:`XcFunc` object with info on the exchange-correlation functional."""
        return self.reader.read_abinit_xcfunc()

    @lazy_property
    def _field(self):
        """
        The field object provided by the subclass.
        Methods of this base class should use self._field to implement
        logic common to the sub-classes.
        """
        return getattr(self, self.__class__._field_name)

    def close(self):
        self.reader.close()

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(str(self.structure))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))
        app("XC functional: %s" % str(self.xc))

        # Add info on the field
        app(marquee(self._field.__class__.__name__, mark="="))
        app(str(self._field))

        return "\n".join(lines)

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("fig = ncfile.ebands.kpoints.plot()"),
            nbv.new_code_cell("fig = ncfile.ebands.plot()"),
            nbv.new_code_cell("fig = ncfile.ebands.get_edos().plot()"),
            nbv.new_code_cell("#cube = ncfile.write_cube(filename=None)"),
            nbv.new_code_cell("#xsf_path = ncfile.write_xsf(filename=None"),
            nbv.new_code_cell("#chgcar = ncfile.write_chgcar(filename=None"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class DensityNcFile(_NcFileWithField):
    """
    Netcdf File containing the electronic density.

    Usage example:

    .. code-block:: python

        with DensityNcFile("foo_DEN.nc") as ncfile:
            ncfile.density
            ncfile.ebands.plot()
    """
    _field_name = "density"

    @lazy_property
    def density(self):
        return self.reader.read_density()

    def to_string(self, verbose=0):
        s = super(DensityNcFile, self).to_string(verbose=verbose)

        # Add density related stuff.
        lines = [" "]
        app = lines.append
        app("Integrated electronic and magnetization densities in atomic spheres:")
        df = self.density.integrate_in_spheres(rcut_symbol=None, out=False)
        app(str(df))
        app("Total magnetization from unit cell integration: %s" % str(self.density.magnetization))

        return s + "\n".join(lines)

    def write_chgcar(self, filename=None):
        """
        Write density in CHGCAR format. Return :class:`ChgCar` instance.
        """
        if filename is None:
            filename = self.basename.replace(".nc", "_CHGCAR")
            cprint("Writing density in CHGCAR format to file: %s" % filename, "yellow")

        return self.density.to_chgcar(filename=filename)

    def write_xsf(self, filename=None):
        """
        Write density in XSF format.
        """
        if filename is None:
            filename = self.basename.replace(".nc", ".xsf")
            cprint("Writing density in XSF format to file: %s" % filename, "yellow")

        return self.density.export(filename)

    def write_cube(self, filename=None, spin="total"):
        if filename is None:
            filename = self.basename.replace(".nc", ".cube")
            cprint("Writing density in CUBE format to file: %s" % filename, "yellow")

        return self.density.export_to_cube(filename, spin=spin)


class VhartreeNcFile(_NcFileWithField):
    _field_name = "vh"

    @lazy_property
    def vh(self):
        """Hartree potential."""
        return self.reader.read_vh()


class VxcNcFile(_NcFileWithField):
    _field_name = "vxc"

    @lazy_property
    def vxc(self):
        """XC potential."""
        return self.reader.read_vxc()


class VhxcNcFile(_NcFileWithField):
    _field_name = "vhxc"

    @lazy_property
    def vhxc(self):
        """Hartree + XC potential."""
        return self.reader.read_vhxc()


class PotNcFile(_NcFileWithField):
    _field_name = "vks"

    @lazy_property
    def vks(self):
        """Hartree + XC potential + sum of local pseudo-potential terms."""
        return self.reader.read_vks()


class DensityFortranFile(AbinitFortranFile):
    """
    Class representing the _DEN fortran file containing density.
    Provides methods to run Cut3D and convert data into different formats.
    """

    def _convert(self, cut3d_input, workdir=None):
        """
        Internal function to run a conversion using cut3d.
        """
        workdir = tempfile.mkdtemp() if workdir is None else workdir

        # local import to avoid circular references
        from abipy.flowtk import Cut3D
        outfile, converted_file = Cut3D().cut3d(cut3d_input, workdir)

        return converted_file

    def get_cube(self, out_filepath, workdir=None):
        """
        Runs cut3d to convert the density to the cube format

        Args:
            out_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            workdir: directory where cut3d is executed.

        Returns:
            (CubeFile) the converted file as a CubeFile object.
        """
        from abipy.abio.inputs import Cut3DInput
        return CubeFile(self._convert(cut3d_input=Cut3DInput.den_to_cube(self.filepath, out_filepath),
                                      workdir=workdir))

    def get_xsf(self, out_filepath, shift=None, workdir=None):
        """
        Runs cut3d to convert the density to the xsf format

        Args:
            out_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            shift: a list of three integers defining the shift along the x, y, z axis. None if no shift is required.
            workdir: directory where cut3d is executed.

        Returns:
            (string) path to the converted file.
        """
        from abipy.abio.inputs import Cut3DInput
        return self._convert(cut3d_input=Cut3DInput.den_to_xsf(self.filepath, output_filepath=out_filepath, shift=shift),
                             workdir=workdir)

    def get_tecplot(self, out_filepath, workdir=None):
        """
        Runs cut3d to convert the density to the tecplot format

        Args:
            out_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            workdir: directory where cut3d is executed.

        Returns:
            (string) path to the converted file.
        """
        from abipy.abio.inputs import Cut3DInput
        return self._convert(cut3d_input=Cut3DInput.den_to_tecplot(self.filepath, out_filepath), workdir=workdir)

    def get_molekel(self, out_filepath, workdir=None):
        """
        Runs cut3d to convert the density to the molekel format

        Args:
            out_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            workdir: directory where cut3d is executed.

        Returns:
            (string) path to the converted file.
        """
        from abipy.abio.inputs import Cut3DInput
        return self._convert(cut3d_input=Cut3DInput.den_to_molekel(self.filepath, out_filepath), workdir=workdir)

    def get_3d_indexed(self, out_filepath, workdir=None):
        """
        Runs cut3d to convert the density to the 3D indexed format

        Args:
            out_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            workdir: directory where cut3d is executed.

        Returns:
            (string) path to the converted file.
        """
        from abipy.abio.inputs import Cut3DInput
        return self._convert(cut3d_input=Cut3DInput.den_to_3d_indexed(self.filepath, out_filepath), workdir=workdir)

    def get_3d_formatted(self, out_filepath, workdir=None):
        """
        Runs cut3d to convert the density to the 3D formatted format

        Args:
            out_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            workdir: directory where cut3d is executed.

        Returns:
            (string) path to the converted file.
        """
        from abipy.abio.inputs import Cut3DInput
        return self._convert(cut3d_input=Cut3DInput.den_to_3d_indexed(self.filepath, out_filepath), workdir=workdir)

    def get_hirshfeld(self, structure, all_el_dens_paths=None, fhi_all_el_path=None, workdir=None):
        """
        Runs cut3d to get the Hirshfeld charges. Requires all-electron density files, so at least one source
        should be specified (all_el_dens_paths or fhi_all_el_path)

        Args:
            structure: a Structure object representing the structure used to generate the density
            all_el_dens_paths: a list of paths to the all-electron density files correspinding to the elements defined
                in the abinit input.
            fhi_all_el_path: path to the folder containing the fhi all-electron density files that will be used
                to automatically determine the path of the required densities.
            workdir: directory where cut3d is executed.

        Returns:
            (HirshfeldCharges) the calculated Hirshfeld charges.
        """
        if all_el_dens_paths is None and fhi_all_el_path is None:
            raise ValueError("At least one source of all electron densities should be provided.")
        if all_el_dens_paths is not None and fhi_all_el_path is not None:
            raise ValueError("all_el_dens_paths and fhi_all_el_path are mutually exclusive.")

        # local import to avoid circular references
        from abipy.flowtk import Cut3D
        from abipy.abio.inputs import Cut3DInput

        if all_el_dens_paths is not None:
            cut3d_input = Cut3DInput.hirshfeld(self.filepath, all_el_dens_paths)
        else:
            cut3d_input = Cut3DInput.hirshfeld_from_fhi_path(self.filepath, structure, fhi_all_el_path)

        workdir = tempfile.mkdtemp() if workdir is None else workdir

        cut3d = Cut3D()
        outfile, converted_file = cut3d.cut3d(cut3d_input, workdir)

        from abipy.electrons.charges import HirshfeldCharges
        return HirshfeldCharges.from_cut3d_outfile(structure=structure, filepath=cut3d.stdout_fname)

# coding: utf-8
"""Density/potential files in netcdf/fortran format."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import tempfile
import numpy as np

from monty.string import marquee
from monty.termcolor import cprint
from monty.functools import lazy_property
from abipy.core.mixins import (AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter,
    AbinitFortranFile, CubeFile)
from abipy.flowtk import Cut3D
from abipy.core.fields import FieldReader
from abipy.abio.inputs import Cut3DInput
from abipy.electrons.ebands import ElectronsReader


__all__ = [
    "DensityNcFile",
    "VhartreeNcFile",
    "VxcNcFile",
    "VhxcNcFile",
    "PotNcFile",
]


class Cut3dDenPotNcFile(AbinitNcFile, Has_Structure):
    """
    Netcdf file with structure and density/potential produced by CUT3d
    Unlike _NcFileWithField subclasses, this object does not contain an electronic band-structure
    and it's mainly used to convert from Fortran DEN/POT to netcdf.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Cut3dDenPotNcFile
    """
    def __init__(self, filepath):
        super(Cut3dDenPotNcFile, self).__init__(filepath)
        self.reader = FieldReader(filepath)
        self.field = self.reader.read_field()

    @property
    def structure(self):
        """|Structure| object."""
        return self.field.structure

    def close(self):
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}


class _DenPotNcReader(ElectronsReader, FieldReader):
    """Object used to read data from density/potential files in netcdf format."""


class _NcFileWithField(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    Base class providing commong methods for netcdf files with density/potential
    """
    field_name = None

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(_NcFileWithField, self).__init__(filepath)
        self.reader = _DenPotNcReader(filepath)

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def xc(self):
        """:class:`XcFunc` object with info on the exchange-correlation functional."""
        return self.reader.read_abinit_xcfunc()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @lazy_property
    def field(self):
        """
        The field object provided by the subclass.
        Methods of this base class should use self.field to implement
        logic common to the sub-classes.
        """
        return getattr(self, self.__class__.field_name)

    def close(self):
        """Close the file."""
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
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))
        app("XC functional: %s" % str(self.xc))

        # Add info on the field
        app(marquee(self.field.__class__.__name__, mark="="))
        app(str(self.field))

        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.structure.plot(show=False)
        yield self.ebands.plot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.ebands.kpoints.plot();"),
            nbv.new_code_cell("ncfile.ebands.plot();"),
            nbv.new_code_cell("ncfile.ebands.get_edos().plot();"),
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

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: DensityNcFile
    """
    field_name = "density"

    @lazy_property
    def density(self):
        """Density object."""
        return self.reader.read_density()

    def to_string(self, verbose=0):
        """String representation."""
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
        Write density in XSF format (xcrysden_)
        """
        if filename is None:
            filename = self.basename.replace(".nc", ".xsf")
            cprint("Writing density in XSF format to file: %s" % filename, "yellow")

        return self.density.export(filename)

    def write_cube(self, filename=None, spin="total"):
        """Write density in CUBE format."""
        if filename is None:
            filename = self.basename.replace(".nc", ".cube")
            cprint("Writing density in CUBE format to file: %s" % filename, "yellow")

        return self.density.export_to_cube(filename, spin=spin)


class VhartreeNcFile(_NcFileWithField):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VhartreeNcFile
    """
    field_name = "vh"

    @lazy_property
    def vh(self):
        """Hartree potential."""
        return self.reader.read_vh()


class VxcNcFile(_NcFileWithField):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VxcNcFile
    """
    field_name = "vxc"

    @lazy_property
    def vxc(self):
        """XC potential."""
        return self.reader.read_vxc()


class VhxcNcFile(_NcFileWithField):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VhxcNcFile
    """
    field_name = "vhxc"

    @lazy_property
    def vhxc(self):
        """Hartree + XC potential."""
        return self.reader.read_vhxc()


class PotNcFile(_NcFileWithField):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PotNcFile
    """
    field_name = "vks"

    @lazy_property
    def vks(self):
        """Hartree + XC potential + sum of local pseudo-potential terms."""
        return self.reader.read_vks()


class DensityFortranFile(AbinitFortranFile):
    """
    Class representing the _DEN fortran file containing density.
    Provides methods to run Cut3D and convert data into different formats.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: DensityFortranFile
    """

    def _convert(self, cut3d_input, workdir=None):
        """
        Internal function to run a conversion using cut3d.
        """
        workdir = os.path.abspath(tempfile.mkdtemp() if workdir is None else workdir)
        outfile, converted_file = Cut3D().cut3d(cut3d_input, workdir)

        return converted_file

    def get_cube(self, out_filepath, workdir=None):
        """
        Runs cut3d to convert the density to the cube format

        Args:
            out_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on the compiler.
            workdir: directory where cut3d is executed.

        Returns:
            (CubeFile) the converted file as a CubeFile object.
        """
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
        return self._convert(cut3d_input=Cut3DInput.den_to_xsf(self.filepath,
                             output_filepath=out_filepath, shift=shift), workdir=workdir)

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
        return self._convert(cut3d_input=Cut3DInput.den_to_tecplot(self.filepath, out_filepath),
                             workdir=workdir)

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
        return self._convert(cut3d_input=Cut3DInput.den_to_molekel(self.filepath, out_filepath),
                             workdir=workdir)

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
        return self._convert(cut3d_input=Cut3DInput.den_to_3d_indexed(self.filepath, out_filepath),
                             workdir=workdir)

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
        return self._convert(cut3d_input=Cut3DInput.den_to_3d_indexed(self.filepath, out_filepath),
                             workdir=workdir)

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

        if all_el_dens_paths is not None:
            cut3d_input = Cut3DInput.hirshfeld(self.filepath, all_el_dens_paths)
        else:
            cut3d_input = Cut3DInput.hirshfeld_from_fhi_path(self.filepath, structure, fhi_all_el_path)

        workdir = os.path.abspath(tempfile.mkdtemp() if workdir is None else workdir)

        cut3d = Cut3D()
        outfile, converted_file = cut3d.cut3d(cut3d_input, workdir)

        from abipy.electrons.charges import HirshfeldCharges
        return HirshfeldCharges.from_cut3d_outfile(structure=structure, filepath=cut3d.stdout_fname)

    def get_density(self, workdir=None):
        """
        Invoke cut3d to produce a netcdf file with the density, read the file and return Density object.

        Args:
            workdir: directory in which cut3d is executed.
        """
        workdir = os.path.abspath(tempfile.mkdtemp() if workdir is None else workdir)
        output_filepath = os.path.join(workdir, "field_CUT3DDENPOT.nc")
        # FIXME Converters with nspden > 1 won't work since cut3d asks for the ispden index.
        cut3d_input = Cut3DInput(infile_path=self.filepath, output_filepath=output_filepath,
                                 options=[15, output_filepath, 0, 0])

        outfile, _ = Cut3D().cut3d(cut3d_input, workdir)
        with Cut3dDenPotNcFile(output_filepath) as nc:
            assert nc.field.is_density_like and nc.field.netcdf_name == "density"
            return nc.field

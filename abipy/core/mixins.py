# coding: utf-8
"""This module ..."""
from __future__ import print_function, division, unicode_literals

import abc
import os
import six
import collections

from time import ctime
from monty.os.path import which
from monty.string import is_string
from monty.functools import lazy_property
from pymatgen.io.abinit.events import EventsParser
from pymatgen.io.abinit.abiinspect import GroundStateScfCycle, D2DEScfCycle, Relaxation
from pymatgen.io.abinit.abitimer import AbinitTimerParser


__all__ = [
    "AbinitNcFile",
    "Has_Structure",
    "Has_ElectronBands",
    "Has_PhononBands",
]

@six.add_metaclass(abc.ABCMeta)
class _File(object):
    """
    Abstract base class defining the methods that must be implemented 
    by the concrete classes representing the different files produced by ABINIT.
    """
    def __init__(self, filepath):
        self._filepath = os.path.abspath(filepath)

    def __repr__(self):
        return "<%s, %s>" % (self.__class__.__name__, self.relpath)

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a string."""
        if isinstance(filepath, cls):
            return filepath

        try:
            return cls(filepath)
        except:
            import traceback
            msg = traceback.format_exc()
            msg += "\n Perhaps the subclass %s must redefine the classmethod from_file\n" % cls
            raise ValueError(msg)

    @property
    def filepath(self):
        """Absolute path of the file."""
        return self._filepath

    @property
    def relpath(self):
        """Relative path."""
        try:
            return os.path.relpath(self.filepath)
        except OSError:
            # current working directory may not be defined!
            return self.filepath

    @property
    def basename(self):
        """Basename of the file."""
        return os.path.basename(self.filepath)

    @property
    def filetype(self):
        """String defining the filetype."""
        return self.__class__.__name__

    def filestat(self):
        """Dictionary with file metadata"""
        return get_filestat(self.filepath)

    @abc.abstractmethod
    def close(self):
        """Close file."""
    
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement. It automatically closes the file."""
        self.close()


class TextFile(_File):
    def __enter__(self):
        # Open the file
        self._file
        return self

    def __iter__(self):
        return iter(self._file)

    @lazy_property
    def _file(self):
        """File object open in read-only mode."""
        return open(self.filepath, mode="rt")

    def close(self):
        """Close the file."""
        try:
            self._file.close()
        except:
            pass

    def seek(self, offset, whence=0):
        """Set the file's current position, like stdio's fseek()."""
        self._file.seek(offset, whence)


class AbinitTextFile(TextFile):
    """Class for the ABINIT main output file and the log file."""

    @lazy_property
    def events(self):
        """List of ABINIT events reported in the file."""
        return EventsParser().parse(self.filepath)

    @lazy_property
    def timer_data(self):
        """Timer data."""
        return AbinitTimerParser().parse(self.filepath)


class AbinitInputFile(TextFile):
    """Class representing the input file."""

class AbinitLogFile(AbinitTextFile):
    """Class representing the log file."""


class AbinitOutputFile(AbinitTextFile):
    """Class representing the main output file."""

    def next_gs_scf_cycle(self):
        """Return the next :class:`GroundStateScfCycle` in the file. None if not found."""
        return GroundStateScfCycle.from_stream(self)

    def next_d2de_scf_cycle(self):
        """:class:`GroundStateScfCycle` with information on the GS iterations. None if not found."""
        return D2DEScfCycle.from_stream(self)

    def compare_gs_scf_cycles(self, others, show=True):
        """
        Produce and returns a list of `matplotlib` figure comparing the GS self-consistent 
        cycle in self with the ones in others.

        Args:
            others: list of `AbinitOutputFile` objects or strings with paths to output files.
            show: True to diplay plots.
        """
        for i, other in enumerate(others):
            if is_string(other): others[i] = self.__class__.from_file(other)
        
        fig, figures = None, []
        while True:
            cycle = self.next_gs_scf_cycle()
            if cycle is None: break 

            fig = cycle.plot(show=False)
            for i, other in enumerate(others):
                other_cycle = other.next_gs_scf_cycle()
                if other_cycle is None: break
                last = (i == len(others) - 1)
                fig = other_cycle.plot(fig=fig, show=show and last)
                if last: figures.append(fig)

        self.seek(0)
        for other in others: other.seek(0)
        return figures

    def compare_d2de_scf_cycles(self, others, show=True):
        """
        Produce and returns a `matplotlib` figure comparing the DFPT self-consistent 
        cycle in self with the ones in others.
                                                                                                   
        Args:
            others: list of `AbinitOutputFile` objects or strings with paths to output files.
            show: True to diplay plots.
        """
        for i, other in enumerate(others):
            if is_string(other): others[i] = self.__class__.from_file(other)

        fig, figures = None, []
        while True:
            cycle = self.next_d2de_scf_cycle()
            if cycle is None: break 

            fig = cycle.plot(show=False)
            for i, other in enumerate(others):
                other_cycle = other.next_d2de_scf_cycle()
                if other_cycle is None: break
                last = (i == len(others) - 1)
                fig = other_cycle.plot(fig=fig, show=show and last)
                if last: figures.append(fig)

        self.seek(0)
        for other in others: other.seek(0)
        return figures
                                                                                                   

@six.add_metaclass(abc.ABCMeta)
class AbinitNcFile(_File):
    """
    Abstract class representing a Netcdf file with data saved
    according to the ETSF-IO specifications (when available).
    """
    def ncdump(self, *nc_args, **nc_kwargs):
        """Returns a string with the output of ncdump."""
        return NcDumper(*nc_args, **nc_kwargs).dump(self.filepath)


@six.add_metaclass(abc.ABCMeta)
class Has_Structure(object):
    """Mixin class for :class:`AbinitNcFile` containing crystallographic data."""

    @abc.abstractproperty
    def structure(self):
        """Returns the :class:`Structure` object."""

    def show_bz(self):
        """
        Gives the plot (as a matplotlib object) of the symmetry line path in the Brillouin Zone.
        """
        return self.structure.hsym_kpath.get_kpath_plot()

    def export_structure(self, filepath):
        """
        Export the structure on file.

        returns:
            Instance of :class:`Visualizer`
        """
        return self.structure.export(filepath)

    def visualize_structure_with(self, visu_name):
        """
        Visualize the crystalline structure with the specified visualizer.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        from abipy.iotools.visualizer import Visualizer
        visu = Visualizer.from_name(visu_name)

        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                return self.export_structure(ext)
            except visu.Error:
                pass
        else:
            raise visu.Error("Don't know how to export data for visu_name %s" % visu_name)


@six.add_metaclass(abc.ABCMeta)
class Has_ElectronBands(object):
    """Mixin class for :class:`AbinitNcFile` containing electron data."""

    @abc.abstractproperty
    def ebands(self):
        """Returns the :class:`ElectronBands` object."""

    @property
    def nsppol(self):
        """Number of spin polarizations"""
        return self.ebands.nsppol

    #@property
    #def nspinor(self):
    #    """Number of spinors"""
    #    return self.ebands.nspinor

    #@property
    #def nspden(self):
    #    """Number of indepedendent spin-density components."""
    #    return self.ebands.nspden

    @property
    def mband(self):
        """Maximum number of bands."""
        return self.ebands.mband

    @property
    def nband(self):
        """Maximum number of bands."""
        return self.ebands.nband

    @property
    def nelect(self):
        """Number of elecrons per unit cell"""
        return self.ebands.nelect

    @property
    def kpoints(self):
        """Iterable with the Kpoints."""
        return self.ebands.kpoints

    @property
    def nkpts(self):
        """Number of k-points."""
        return len(self.kpoints)

    def plot_ebands(self, **kwargs):
        """Plot the electron energy bands. See the :func:`ElectronBands.plot` for the signature."""
        return self.ebands.plot(**kwargs)

    def plot_ebands_with_edos(self, dos, **kwargs):
        return self.ebands.plot_with_edos(dos, **kwargs)


@six.add_metaclass(abc.ABCMeta)
class Has_PhononBands(object):
    """Mixin class for :class:`AbinitNcFile` containing phonon data."""

    @abc.abstractproperty
    def phbands(self):
        """Returns the :class:`PhononBands` object."""

    def plot_phbands(self, **kwargs):
        """
        Plot the electron energy bands. See the :func:`PhononBands.plot` for the signature.""
        """
        return self.phbands.plot(**kwargs)

    #def plot_phbands_with_phdos(self, phdos, **kwargs):
    #    return self.phbands.plot_with_phdos(phdos, **kwargs)


class NcDumper(object):
    """Wrapper object for the ncdump tool."""

    def __init__(self, *nc_args, **nc_kwargs):
        """
        Args:
            nc_args: Arguments passed to ncdump.
            nc_kwargs: Keyword arguments passed to ncdump
        """
        self.nc_args = nc_args
        self.nc_kwargs = nc_kwargs

        self.ncdump = which("ncdump")

    def dump(self, filepath):
        """Returns a string with the output of ncdump."""
        if self.ncdump is None:
            return "Cannot find ncdump tool in PATH"
        else:
            from subprocess import check_output
            return check_output(["ncdump", filepath])


_ABBREVS = [
    (1 << 50, 'Pb'),
    (1 << 40, 'Tb'),
    (1 << 30, 'Gb'),
    (1 << 20, 'Mb'),
    (1 << 10, 'kb'),
    (1, 'b'),
]


def size2str(size):
    """Convert size to string with units."""
    for factor, suffix in _ABBREVS:
        if size > factor:
            break
    return "%.2f " % (size / factor) + suffix


def get_filestat(filepath):
    stat = os.stat(filepath)
    return collections.OrderedDict([
        ("Name", os.path.basename(filepath)),
        ("Directory", os.path.dirname(filepath)),
        ("Size", size2str(stat.st_size)),
        ("Access Time", ctime(stat.st_atime)),
        ("Modification Time", ctime(stat.st_mtime)),
        ("Change Time", ctime(stat.st_ctime)),
    ])

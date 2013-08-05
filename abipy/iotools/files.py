"""This module ..."""
from __future__ import division, print_function

import abc
import os
import collections

from time import ctime
from pymatgen.io.abinitio import EventParser
from abipy.tools import which

__all__ = [
    "AbinitNcFile",
]


class AbinitFile(object):
    """
    Abstract base class defining the methods that must be  implemented 
    by the concrete classes representing the different files produced by ABINIT.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, filepath):
        self._filepath = os.path.abspath(filepath)

    def __repr__(self):
        return "<%s at %s, filepath = %s>" % (self.__class__.__name__, id(self), self.filepath)

    def __str__(self):
        return repr(self)

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a string."""
        if isinstance(filepath, cls):
            return filepath

        try:
            return cls(filepath)
        except Exception as exc:
            msg = "%s\n Perhaps the subclass %s must redefine the classmethod from_file\n" % (str(exc), cls)
            raise ValueError(msg)

    @property
    def filepath(self):
        """Absolute path of the file."""
        return self._filepath

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

    #@abc.abstractmethod
    #def close(self):


class AbinitTextFile(AbinitFile):
    """Class for the ABINIT main output file and the log file."""

    @property
    def events(self):
        """List of ABINIT events reported in the file."""
        try:
            return self._events
        except AttributeError:
            self._events = EventParser().parse(self.filepath)
            return self._events

    @property
    def timer_data(self):
        """Timer data."""
        try:
            return self._timer_data

        except AttributeError:
            from abipy.htc.abitimer import AbinitTimerParser
            parser = AbinitTimerParser()
            parser.read(self.filepath)
            self._timer_data = parser
            return self._timer_data


class AbinitOutputFile(AbinitTextFile):
    """Class representing the main output file."""


class AbinitLogFile(AbinitTextFile):
    """Class representing the log file."""


class AbinitNcFile(AbinitFile):
    """
    Abstract class representing a Netcdf file with data saved
    according to the ETSF-IO specifications (when available).
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_structure(self):
        """Returns the `Structure` object."""

    def ncdump(self, *nc_args, **nc_kwargs):
        """Returns a string with the output of ncdump."""
        return NcDumper(*nc_args, **nc_kwargs).dump(self.filepath)


class NcDumper(object):
    """Wrapper object for the ncdump tool."""

    def __init__(self, *nc_args, **nc_kwargs):
        """
        Args:
            nc_args:
                Arguments passed to ncdump.
            nc_kwargs:
                Keyword arguments passed to ncdump
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

            cmd = ["ncdump", filepath]
            return check_output(cmd)


_ABBREVS = [
    (1 << 50L, 'Pb'),
    (1 << 40L, 'Tb'),
    (1 << 30L, 'Gb'),
    (1 << 20L, 'Mb'),
    (1 << 10L, 'kb'),
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

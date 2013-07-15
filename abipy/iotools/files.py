"""This module ..."""
from __future__ import division, print_function

import abc
import os

from abipy.tools import which

__all__ = [
    "AbinitNcFile",
]

class NcDumper(object):
    """This object wraps the ncdump tool"""

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


class AbinitNcFile(object):
    """
    Abstract base class defining the methods that must be 
    implemented by the concrete classes representing the Netcdf file 
    produces by ABINIT.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, filepath):
        self._filepath = os.path.abspath(filepath)

    def __repr__(self):
        return "<%s at %s, filetype = %s>" % (self.__class__.__name__, id(self), self.filetype)

    def __str__(self):
        return self.summary

    @property
    def filetype(self):
        """String defining the filetype."""
        return self.__class__.__name__

    #@abstractclassmethod
    def from_ncfile(cls, filepath):
        """Initialize the object from a Netcdf file"""
        try:
            return cls(filepath)
        except:
            raise ValueError("Subclass must define the class method from_ncfile")

    @property
    def filepath(self):
        """Absolute path of the file."""
        return self._filepath

    @property
    def basename(self):
        """Basename of the file"""
        return os.path.basename(self.filepath)

    @abc.abstractmethod
    def get_structure(self):
        """Returns the `Structure` object."""

    #@abc.abstractproperty
    #def summary(self):
    #    """String summarizing the most important properties."""

    #@abc.abstractproperty
    #def fileinfo(self):
    #    """File metadata"""

    def ncdump(self, *nc_args, **nc_kwargs):
        """Returns a string with the output of ncdump."""
        return NcDumper(*nc_args, **nc_kwargs).dump(self.filepath)

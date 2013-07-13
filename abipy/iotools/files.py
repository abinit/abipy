"""This module ..."""
from __future__ import division, print_function

import abc
import os

from abipy.tools import which

__all__ = [
    "AbinitNcFile",
]

class AbinitNcFile(object):
    """
    Abstract base class defining the methods that must be 
    implemented by the concrete classes representing the Netcdf file 
    produces by ABINIT.
    """
    __metaclass__ = abc.ABCMeta

    def __repr__(self):
        return "<%s at %s, filetype = %s>" % (self.__class__.__name__, id(self), self.filetype)

    def __str__(self):
        return self.summary

    @property
    def filetype(self):
        """String defining the filetype."""
        return self.__class__.__name__

    #@abstractclassmethod
    def from_ncfile(cls, path):
        """Initialize the object from a Netcdf file"""
        try:
            return cls(path)
        except:
            raise ValueError("Subclass must define the class method from_ncfile")

    @abc.abstractproperty
    def filepath(self):
        """Absolute path of the file."""

    @property
    def basename(self):
        """Basename of the file"""
        return os.path.basename(self.filepath)

    #@abc.abstractproperty
    #def summary(self):
    #    """String summarizing the most important properties."""

    #@abc.abstractproperty
    #def fileinfo(self):
    #    """File metadata"""

    def ncdump(self, *nc_args, **nc_kwargs):
        """Returns a string with the output of ncdump."""
        ncdump = which("ncdump")
        if ncdump is None:
            return "Cannot find ncdump tool in PATH"
        else:
            from subprocess import check_output
            cmd = ["ncdump", self.filepath]
            return check_output(cmd)

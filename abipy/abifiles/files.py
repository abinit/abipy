"""This module ..."""
from __future__ import division, print_function

from abipy.waves.wfkfile import WFK_File

import abc

__all__ = [
    "abiopen",
]

class AbstractNcFile(object):
    """
    Abstract base class defining the methods that must be 
    implemented by the concrete classes.
    """
    __metaclass__ = abc.ABCMeta

    def __repr__(self):
        return "<%s at %s, filetype = %s>" % (
            self.__class__.__name__, id(self), self.filetype)

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
    def summary(self):
        """String summarizing the most important properties."""

    @abc.abstractproperty
    def filepath(self):
        """Absolute path of the file."""


def abiopen(filename):
    """
    Factory function that returns the appropriate object
    from the extension of filename.

    Args:
        filename:
            string with the filename.
    """
    ext2class = {
        #"GSR.nc": GSR_File
        #"PHBST.nc": PHBST_File,
        "WFK.nc": WFK_File,
    }

    ext = os.path.split("_")[-1]
    try:
        klass = ext2class[ext]
    except KeyError:
        raise KeyError("Unsupported extension %s" % ext)

    return klass.from_ncfile(path)

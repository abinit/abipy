"""abipy is a set of tools for the analysis of ABINIT results."""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Setup everything
#-----------------------------------------------------------------------------

# Don't forget to also update setup.py when this changes!
import sys
if sys.version[0:3] < '2.7':
    raise ImportError('Python Version 2.7 or above is required for abipy.')
del sys

#-----------------------------------------------------------------------------
# Setup the top level names
#-----------------------------------------------------------------------------
import abipy.core.constants 

from abipy.core import release
from abipy.profile import abipy_env
from abipy.core import *
from abipy.kpoints import *
from abipy.electrons import *
from abipy.waves import *
from abipy.htc import Launcher, MassLauncher, AbinitInput

from abipy.tests import get_reference_file, get_ncfile

# Release data
__author__ = ''
for author, email in release.authors.itervalues():
  __author__ += author + ' <' + email + '>\n'
del author, email

__license__  = release.license
__version__  = release.version


from abipy.waves import WFK_File
from abipy.electrons import SIGRES_File
from abipy.phonons import PHBST_File
from abipy.iotools.files import AbinitNcFile

def ncfile_subclass_from_filename(filename):
    ext2ncfile = {
        "SIGRES.nc": SIGRES_File,
        "WFK-etsf.nc": WFK_File,
        "MDF.nc" : MDF_File,
        #"GSR.nc": GSR_File
        #"PHDOS.nc": PHDOS_File,
        #"PHBST.nc": PHBST_File,
    }

    ext = filename.split("_")[-1]
    try:
        return ext2ncfile[ext]
    except KeyError:
        raise KeyError("No ncfile subclass has been registered for extension %s" % ext)

def abiopen(filepath):
    """
    Factory function that returns the appropriate `AbinitNcFile` object

    Args:
        filepath:
            string with the filename or `AbinitNcFile` instance
    """
    if isinstance(filepath, AbinitNcFile):
        return filepath

    # Assume string.
    ncfile = ncfile_subclass_from_filename(filepath)
    return ncfile.from_ncfile(filepath)

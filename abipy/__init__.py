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
#import abipy.core.constants 

from abipy.core import release
from abipy.profile import abipy_env
from abipy.core import *
from abipy.electrons import *
from abipy.phonons import *
from abipy.waves import *
from abipy.htc import Launcher, MassLauncher, AbinitInput

# Release data
__author__ = ''
for author, email in release.authors.itervalues():
  __author__ += author + ' <' + email + '>\n'
del author, email

__license__  = release.license
__version__  = release.version


def abifile_subclass_from_filename(filename):
    from abipy.iotools.files import AbinitFile, AbinitLogFile, AbinitOutputFile
    from abipy.electrons import SIGRES_File, GSR_File
    from abipy.waves import WFK_File
    #from abipy.phonons import PHDOS_File, PHBST_File

    ext2ncfile = {
        "SIGRES.nc": SIGRES_File,
        "WFK-etsf.nc": WFK_File,
        "MDF.nc" : MDF_File,
        "GSR.nc": GSR_File
        #"PHDOS.nc": PHDOS_File,
        #"PHBST.nc": PHBST_File,
    }

    #if filename.endswith(".abi"):
    #    return AbinitInputFile
                                                                                        
    if filename.endswith(".abo"):
        return AbinitOutputFile
    
    if filename.endswith(".abl"):
        return AbinitLogFile

    # CIF files.
    if filename.endswith(".cif"):
        from abipy.core.structure import Structure
        return Structure.from_file(filename)

    ext = filename.split("_")[-1]
    try:
        return ext2ncfile[ext]
    except KeyError:
        raise KeyError("No class has been registered for extension %s" % ext)


def abiopen(filepath):
    """
    Factory function that opens any file supported by abipy.

    Args:
        filepath:
            string with the filename. 
    """
    cls = abifile_subclass_from_filename(filepath)
    return cls.from_file(filepath)

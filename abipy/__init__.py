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

# Support for tests
from abipy.tests import * 

# Release data
__author__ = ''
for author, email in release.authors.itervalues():
  __author__ += author + ' <' + email + '>\n'
del author, email

__license__  = release.license
__version__  = release.version


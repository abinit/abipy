"""abipy is a set of tools for the analysis of ABINIT results."""

#-----------------------------------------------------------------------------
# Setup the top level names
#-----------------------------------------------------------------------------

from abipy.core import release

# Release data
__author__ = ''
for author, email in release.authors.values():
    __author__ += author + ' <' + email + '>\n'
del author, email

__license__ = release.license
__version__ = release.version

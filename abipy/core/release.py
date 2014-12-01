"""Release data for the abipy project."""
#from __future__ import print_function, division, unicode_literals

# Name of the package for release purposes.  This is the name which labels
# the tarballs and RPMs made by distutils, so it's best to lowercase it.
name = 'abipy'

# version information.  An empty _version_extra corresponds to a full
# release.  'dev' as a _version_extra string means this is a development version
_version_major = 0
_version_minor = 1
_version_micro = '0'  # use '' for first of series, number for 1 and above
#_version_extra = 'dev'
_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro: _ver.append(_version_micro)
if _version_extra: _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

version = __version__  # backwards compatibility name

description = "Python modules and scripts to analyze the results of ab-initio calculation performed with ABINIT"

long_description = \
    """
    The goal of abipy is to create a comprehensive environment for
    interactive and exploratory computing.  To support this goal,
    abipy provides:

    * A set of pythons objects to store the results of the calculations.

    * A set of scripts for performing common tasks such as plotting the results
      of the calculation.

    The latest development version is always available from site <https://github.com/gmatteo/abipy>
    """

license = 'GPL'

authors = {'Matteo': ('Matteo Giantomassi', 'gmatteo at gmail.com'),
           'Antonius': ('Gabriel Antonius', 'gabriel.antonius at gmail.com'),
}

author = 'The ABINIT group'

author_email = 'gmatteo at gmail com'

url = 'https://github.com/gmatteo/abipy'

#download_url = 'https://github.com/gmatteo/abipy'

platforms = ['Linux', 'darwin']

keywords = ["ABINIT", "ab initio", "first principles"]

classifiers=[
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    #"Programming Language :: Python :: 3.2",
    #"Programming Language :: Python :: 3.3",
    #"Programming Language :: Python :: 3.4",
    #"Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    #"License :: OSI Approved :: MIT License",
    #"Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Information Analysis",
    "Topic :: Scientific/Engineering :: Physics",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Software Development :: Libraries :: Python Modules"
],

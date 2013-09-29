#!/usr/bin/env python
"""Setup script for abipy."""
import sys

# This check is also made in abipy/__init__, don't forget to update both when
# changing Python version requirements.
if sys.version[0:3] < '2.7':
    sys.stderr.write("abipy requires Python version 2.7 or above. Exiting.")
    sys.exit(1)

#-------------------------------------------------------------------------------
import os
import sys

from glob import glob
from setuptools import find_packages

#-------------------------------------------------------------------------------
# Useful globals and utility functions
#-------------------------------------------------------------------------------

# A few handy globals
isfile = os.path.isfile
pjoin = os.path.join

# A little utility we'll need below, since glob() does NOT allow you to do exclusion on multiple endings!
def file_doesnt_end_with(test, endings):
    """
    Returns true if test is a file and its name does NOT end with any
    of the strings listed in endings.
    """
    if not isfile(test):
        return False
    for e in endings:
        if test.endswith(e):
            return False
    return True

#---------------------------------------------------------------------------
# Basic project information
#---------------------------------------------------------------------------

# release.py contains version, authors, license, url, keywords, etc.
release_file = pjoin('abipy','core','release.py')

with open(release_file) as f:
    code = compile(f.read(), release_file, 'exec')
    exec(code)

#---------------------------------------------------------------------------
# Find packages
#---------------------------------------------------------------------------
#
#def find_packages():
#    """
#    Find all of abipy's packages.
#    """
#    return find_packages(exclude=())


#---------------------------------------------------------------------------
# Find package data
#---------------------------------------------------------------------------

def find_package_data():
    "Find abipy's package_data."
    # This is not enough for these things to appear in an sdist.
    # We need to muck with the MANIFEST to get this to work
    package_data = {
        'abipy.data' : ['*','pseudos/*','runs/*','cifs/*','benchmarks/*'],
        'abipy.data.runs' : ['data_*/outdata/*'],
        'abipy.gui.awx' : ['images/*'],
    }
    return package_data

def find_exclude_package_data():
    package_data = {
        'abipy.data' : ['pseudos','runs','cifs','benchmarks','runs/data_*'],
    }
    return package_data



#---------------------------------------------------------------------------
# Find scripts
#---------------------------------------------------------------------------

def find_scripts():
    """Find abipy scripts."""
    scripts = []
    #
    # All python files in abipy/scripts
    pyfiles = glob(pjoin('abipy','scripts',"*.py"))
    scripts.extend(pyfiles)
    return scripts


def get_long_desc():
    with open("README.rst") as f:
        long_desc = f.read()
        #ind = long_desc.find("\n")
        #long_desc = long_desc[ind + 1:]
        return long_desc


#-----------------------------------------------------------------------------
# Function definitions
#-----------------------------------------------------------------------------

def cleanup():
    """Clean up the junk left around by the build process."""

    if "develop" not in sys.argv:
        import shutil
        try:
            shutil.rmtree('abipy.egg-info')
        except:
            try:
                os.unlink('abipy.egg-info')
            except:
                pass

# List of external packages we rely on.
# Note setup install will download them from Pypi if they are not avaiable.
install_requires = [
    "pymatgen>=2.7.3",
    "matplotlib>=1.1",
    "scipy>=0.10",
    "numpy>=1.5",
    "netCDF4",
    "PyDispatcher",
    "pyyaml>=3.1.0",
    #"wxPython",
    #"wxmplot",
    #"termcolor",
]

#---------------------------------------------------------------------------
# Find all the packages, package data, and data_files
#---------------------------------------------------------------------------

# Get the set of packages to be included.
#all_packages = find_packages(exclude=())
my_packages = find_packages(exclude=())

my_scripts = find_scripts()

my_package_data = find_package_data()
my_excl_package_data = find_exclude_package_data()
#data_files = find_data_files()

# Create a dict with the basic information
# This dict is eventually passed to setup after additional keys are added.
setup_args = dict(
      name             = name,
      version          = version,
      description      = description,
      long_description = long_description,
      author           = author,
      author_email     = author_email,
      url              = url,
      license          = license,
      platforms        = platforms,
      keywords         = keywords,
      install_requires = install_requires,
      packages         = my_packages,
      package_data     = my_package_data,
      exclude_package_data = my_excl_package_data,
      scripts          = my_scripts,
      #download_url     = download_url,
      )

if __name__ == "__main__":
    from setuptools import setup
    setup(**setup_args)
    # Create the abipyrc file
    #execfile('abipy/profile.py')
    cleanup()

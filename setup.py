#!/usr/bin/env python
"""Setup script for abipy."""
from __future__ import print_function

import sys
import os
import shutil
import numpy as np

from glob import glob
from setuptools import find_packages, setup, Extension
from setuptools.command.install import install

# This check is also made in abipy/__init__, don't forget to update both when
# changing Python version requirements.
if sys.version[0:3] < '2.7':
    sys.stderr.write("abipy requires Python version 2.7 or above. Exiting.")
    sys.exit(1)

class MyInstall(install):
    """
    Extends setuptools install so that we can create the abipy configuration files in ~/abinit/.abipy
    """
    def run(self):
        install.run(self)
        self.__write_user_cfgfiles()

    def __write_user_cfgfiles(self, overwrite=False): 
        """
        This function creates the abipy configurations files.
        
            overwrite:
                True if pre-existent files can be overwritten.
        """
        conf_dir = os.path.join(os.getenv("HOME"), ".abinit", "abipy")
        print("Creating configuration files in %s ... " % conf_dir, end="")

        if not os.path.exists(conf_dir): os.makedirs(conf_dir)

        # TODO: add taskmanager and scheduler
        path_data =[
            ("nodecounter", "-1"),
            #("taskmanager.yml", lambda: fh: fh.write(get_taskmanager_template()))
            #("scheduler.yml", lambda: fh: fh.write(get_scheduler_template()))))
        ]
        path_data = [(os.path.join(conf_dir, f), _) for (f, _) in path_data]

        count = 0
        for p, data in path_data:
            exist = os.path.exists(p) 
            if exist and not overwrite: continue
            # keep a backup copy if we are going to overwrite the file.
            if exist: shutil.move(p, p + ".bkp") 
            with open(p, "w") as fh:
                fh.write(data)
                count += 1

        print("(Wrote %d configuration files)" % count)


# Install ipython with notebook support.
with_ipython = True
if '--with-ipython' in sys.argv:
    with_ipython = True
    sys.argv.remove('--with-ipython')

with_cython = True
try:
    from Cython.Distutils import build_ext
except ImportError:
    with_cython = False

ext_modules = []

# Disable cython for the time being.
with_cython = False
if with_cython:
    #define_macros = [("CYTHON_TRACE", "1")]
    ext_modules += [
        Extension("abipy.extensions.klib", ["abipy/extensions/klib.pyx"], include_dirs=[np.get_include()])
    ]
    cmdclass.update({'build_ext': build_ext})

else:
    ext_modules += [
        Extension("abipy.extensions.klib", ["abipy/extensions/klib.c"], include_dirs=[np.get_include()])
    ]

ext_modules = []

#-------------------------------------------------------------------------------
# Useful globals and utility functions
#-------------------------------------------------------------------------------

# A little utility we'll need below, since glob() does NOT allow you to do exclusion on multiple endings!
def file_doesnt_end_with(test, endings):
    """
    Returns true if test is a file and its name does NOT end with any
    of the strings listed in endings.
    """
    if not os.path.isfile(test):
        return False
    for e in endings:
        if test.endswith(e):
            return False
    return True

#---------------------------------------------------------------------------
# Basic project information
#---------------------------------------------------------------------------

# release.py contains version, authors, license, url, keywords, etc.
release_file = os.path.join('abipy', 'core', 'release.py')

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
    """Find abipy's package_data."""
    # This is not enough for these things to appear in an sdist.
    # We need to muck with the MANIFEST to get this to work
    package_data = {
        'abipy.data': ['*', 'pseudos/*', 'runs/*', 'cifs/*', 'benchmarks/*'],
        'abipy.data.runs': ['data_*/outdata/*'],
        'abipy.htc': ['anaddb_vars.json'],
        'abipy.gui.awx': ['images/*'],
    }
    return package_data


def find_exclude_package_data():
    package_data = {
        'abipy.data': ['pseudos', 'runs', 'cifs', 'benchmarks', 'runs/data_*'],
    }
    return package_data


#---------------------------------------------------------------------------
# Find scripts
#---------------------------------------------------------------------------

def find_scripts():
    """Find abipy scripts."""
    scripts = []
    # All python files in abipy/scripts
    pyfiles = glob(os.path.join('abipy', 'scripts', "*.py"))
    scripts.extend(pyfiles)
    return scripts


def get_long_desc():
    with open("README.rst") as fh:
        return fh.read()


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
# Note setup install will download them from Pypi if they are not available.
install_requires = [
    "termcolor>=1.1.0",
    "apscheduler>=2.1.1",
    "PyDispatcher",
    "numpy",
    #"numpy>=1.8",  # We need this one for the ufuncs
    "scipy>=0.10",
    #"matplotlib>=1.1",
    "pyyaml>=3.1.0",
    "netCDF4",
    "pymatgen>=2.9.0",
    #"fabric",
    #"paramiko",
    "wxmplot",
    #"asciitable",
    #"psutil",
]

if with_ipython:
    install_requires += [
        "ipython>=1.1.0",
        "pyzmq",     # for the notebook
        "jinja2",    
    ]

#if with_cython:
#    install_requires += [
#        "cython",
#    ]

#print("install_requires\n", install_requires)


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
      cmdclass={'install': MyInstall},
      ext_modules=ext_modules,
      )

if __name__ == "__main__":
    setup(**setup_args)
    # Create the abipyrc file
    #execfile('abipy/profile.py')
    cleanup()

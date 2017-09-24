#!/usr/bin/env python
"""Setup script for abipy."""
from __future__ import print_function

import sys
import os
import shutil
import numpy as np

from glob import glob
from setuptools import find_packages, setup, Extension

# This check is also made in abipy/__init__, don't forget to update both when
# changing Python version requirements.
if sys.version[0:3] < '2.7':
    sys.stderr.write("abipy requires Python version 2.7 or above. Exiting.")
    sys.exit(1)

# Install ipython with notebook support.
with_ipython = False
#with_ipython = True
if '--with-ipython' in sys.argv:
    with_ipython = True
    sys.argv.remove('--with-ipython')

with_cython = True
try:
    from Cython.Distutils import build_ext
except ImportError:
    with_cython = False

cmdclass = {}
ext_modules = []

# Disable cython for the time being.
with_cexts = False
if with_cexts:
    with_cython = False
    if with_cython:
        import numpy as np
        #define_macros = [("CYTHON_TRACE", "1")]
        ext_modules += [
            Extension("abipy.extensions.klib", ["abipy/extensions/klib.pyx"], include_dirs=[np.get_include()])
        ]
        cmdclass.update({'build_ext': build_ext})

    else:
        ext_modules += [
            Extension("abipy.extensions.klib", ["abipy/extensions/klib.c"], include_dirs=[np.get_include()])
        ]

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
release_file = os.path.join('abipy','core','release.py')

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
    #top = os.path.join("abipy", "data", "refs")
    #ref_files = {}
    #for root, dirs, files in os.walk(top):
    #    root = root.replace("/", ".")
    #    ref_files[root] = [os.path.join(root, f) for f in files]
    #print(ref_files)

    # This is not enough for these things to appear in an sdist.
    # We need to muck with the MANIFEST to get this to work
    package_data = {
        'abipy.data': ["cifs/*.cif", "pseudos/*", "hgh_pseudos/*", "runs/*", "managers/*", "refs/*.nc", "variables/*"],
        'abipy.data.refs' : [
            "al_g0w0_spfunc/*",
            "alas_nl_dfpt/*",
            "alas_phonons/*",
            "mgb2_fatbands/*",
            "ni_ebands/*",
            "si_bse/*",
            "si_bse_kpoints/*",
            "si_ebands/*",
            "si_g0w0/*",
            "sio2_screening/*",
            "znse_phonons/*",
        ],
        'abipy.htc': ["*.json"],
        'abipy.gui.awx' : ['images/*'],
        'abipy.lessons': ["*.man"],
    }

    #package_data.update(ref_files)
    return package_data


def find_exclude_package_data():
    package_data = {
        'abipy.data': ["managers", 'benchmarks', 'runs/flow_*', 'runs/gspert'],
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
        try:
            shutil.rmtree('abipy.egg-info')
        except (IOError, OSError):
            try:
                os.unlink('abipy.egg-info')
            except:
                pass

# List of external packages we rely on.
# Note setup install will download them from Pypi if they are not available.
install_requires = [
    "six",
    "prettytable",
    "tabulate",
    "apscheduler==2.1.0",
    "pydispatcher>=2.0.5",
    "tqdm",
    "html2text",
    "pyyaml>=3.11",
    "pandas",
    "numpy>=1.9",
    "scipy>=0.14",
    "spglib",
    "pymatgen>=4.7.2",
    "netCDF4",
    "matplotlib>=1.5",
    "seaborn",
]

if with_ipython:
    install_requires += [
        "ipython",
        "jupyter",
        "nbformat",
    ]

#if with_cython:
#    install_requires += [
#        "cython",
#    ]

with_wxpython = False
if with_wxpython:
    install_requires += [
        "wxmplot",
        "wxpython",
    ]

#print("install_requires\n", install_requires)


#---------------------------------------------------------------------------
# Find all the packages, package data, and data_files
#---------------------------------------------------------------------------

# Create a dict with the basic information
# This dict is eventually passed to setup after additional keys are added.
setup_args = dict(
      name=name,
      version=version,
      description=description,
      long_description=long_description,
      author=author,
      author_email=author_email,
      maintainer=maintainer,
      maintainer_email=maintainer_email,
      url=url,
      license=license,
      platforms=platforms,
      keywords=keywords,
      classifiers=classifiers,
      install_requires=install_requires,
      packages=find_packages(exclude=()),
      package_data=find_package_data(),
      exclude_package_data=find_exclude_package_data(),
      scripts=find_scripts(),
      download_url=download_url,
      ext_modules=ext_modules,
      )

if __name__ == "__main__":
    setup(**setup_args)

    print("""
Please read the following if you are about to use abipy for the first time:

[1]
    abipy needs to know about the cluster/computer you are running on. This information
    is provided via the manager.yml and scheduler.yml files. These files must be located
    in ~/.abinit/abipy or in the working directory in which you execute the flow.
    Examples are provided in abipy/data/managers

[2]
    If you are completely new to abipy you may want to start from the abipy lessons.
    The simplest way is to move to an empty directory, start an ipython session and type:

    In [1]: from abipy.lessons.lesson_kpoint_convergence import Lesson()

followed by:

    In [2]: Lesson()

This will print the lessons documentation with further instructions.
Have fun!
""")

    cleanup()

#!/usr/bin/env python
# flake8: noqa
"""Setup script for AbiPy."""

import sys
import os
import shutil

from glob import glob
from setuptools import find_packages, setup

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
        'abipy.panels': [
            "assets/img/*",
        ],
        'abipy.htc': [
            "protocols/*.yml",
        ],
        'abipy.data': [
            "cifs/*.cif",
            "pseudos/*",
            "hgh_pseudos/*",
            "runs/*",
            "managers/*",
            "refs/*.nc",
            "refs/*.log",
            "refs/*.abo",
        ],
        'abipy.data.refs': [
            "al_eph/*",
            "al_g0w0_spfunc/*",
            "alas_nl_dfpt/*",
            "alas_phonons/*",
            #"diamond_sigeph/*",
            "gaas_optic/*",
            "mgb2_fatbands/*",
            "ni_ebands/*",
            "si_bse/*",
            "si_bse_kpoints/*",
            "si_ebands/*",
            "si_g0w0/*",
            "sio2_screening/*",
            "znse_phonons/*",
        ],
        #'abipy.gui.awx': ['images/*'],
    }

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
        return f.read()


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
            except Exception:
                pass


# List of external packages we rely on.
# Note setup install will download them from Pypi if they are not available.
#with open("requirements.txt", "rt") as fh:
#    install_requires = [s.strip() for s in fh]

install_requires = [
    "monty",
    "tabulate",
    #"apscheduler==2.1.0",
    "apscheduler",
    "pydispatcher>=2.0.5",
    "tqdm",
    "pyyaml>=3.11",
    "pandas",
    "numpy",
    "scipy",
    "spglib",
    "pymatgen>=2022.0.14",
    "netCDF4",
    "matplotlib",
    "seaborn",
    "plotly",
    "ipython",
    "chart-studio",
    #pydantic,
    #pymongo,
    #panel,
]

with_wxpython = False
if with_wxpython:
    install_requires += [
        "wxmplot",
        "wxpython",
    ]


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
      long_description_content_type="text/x-rst",
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
Please read the following if you are about to use AbiPy for the first time:

AbiPy needs to know about the cluster/computer you are running on.
This information is provided via two Yaml configuration files: manager.yml and scheduler.yml.
These files must be located either in ~/.abinit/abipy or in the working directory in which you execute the flow.
Examples are provided in abipy/data/managers.
See also the HTML page:

    http://abinit.github.io/abipy/workflows/manager_examples.html

TIPS:

    1) Issue `rehash` in the shell if the AbiPy scripts cannot be found after the installation
    2) Use `abicheck.py --with-flow` to validate the final configuration before running large calculations.

Have fun!
""")
    cleanup()

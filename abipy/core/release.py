# coding: utf-8
"""Release data for the AbiPy project."""
from __future__ import print_function, division, unicode_literals

from collections import OrderedDict

# Name of the package for release purposes.  This is the name which labels
# the tarballs and RPMs made by distutils, so it's best to lowercase it.
name = 'abipy'

# version information.  An empty _version_extra corresponds to a full
# release.  'dev' as a _version_extra string means this is a development version
_version_major = 0
_version_minor = 7
_version_micro = '0'  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
#_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro: _ver.append(_version_micro)
if _version_extra: _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

version = __version__  # backwards compatibility name

# The minimum Abinit version compatible with AbiPy
min_abinit_version = "8.0.8"

description = "Python package to automate ABINIT calculations and analyze the results."

long_description = """\
AbiPy is a Python library to analyze the results produced by `ABINIT <https://www.abinit.org>`_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
AbiPy also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with `Pymatgen <http://www.pymatgen.org>`_ allowing users to
benefit from the different tools and python objects available in the pymatgen ecosystem.

AbiPy can be used in conjunction with  `matplotlib <http://matplotlib.org>`_, `pandas <http://pandas.pydata.org>`_,
`ipython <https://ipython.org/index.html>`_ and `jupyter <http://jupyter.org/>`_
thus providing a powerful and user-friendly environment for data analysis and visualization.
Check out the list of plotting scripts available in our :doc:`gallery </examples/index>`.
To learn more about the integration between jupyter and AbiPy, visit our collection of `notebooks
<http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_
and the `AbiPy lessons <http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/lessons/index.ipynb>`_.

The latest development version is always available from <https://github.com/abinit/abipy>
"""

license = 'GPL'

author = 'M. Giantomassi and the AbiPy group'
author_email = 'matteo.giantomassi@uclouvain.be'
maintainer = "Matteo Giantomassi"
maintainer_email = author_email
authors = OrderedDict([
    ('Matteo', ('M. Giantomassi', 'nobody@nowhere')),
    ('Michiel', ('M. J. van Setten', 'nobody@nowhere')),
    ('Guido', ('G. Petretto', 'nobody@nowhere')),
    ('Henrique', ('H. Miranda', 'nobody@nowhere')),
])

url = "https://github.com/abinit/abipy"
download_url = "https://github.com/abinit/abipy"
platforms = ['Linux', 'darwin']
keywords = ["ABINIT", "ab-initio", "density-function-theory", "first-principles", "electronic-structure", "pymatgen"]
classifiers=[
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Information Analysis",
    "Topic :: Scientific/Engineering :: Physics",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Software Development :: Libraries :: Python Modules",
]

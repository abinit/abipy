.. :Repository: https://github.com/gmatteo/abipy
.. :Author: Matteo Giantomassi (http://github.com/gmatteo)

.. image:: https://pypip.in/v/abipy/badge.png
        :target: https://pypi.python.org/pypi/abipy

.. image:: https://img.shields.io/travis/gmatteo/abipy/master.svg    
        :target: https://travis-ci.org/gmatteo/abipy

.. image:: https://img.shields.io/badge/license-GPL-blue.svg


About
=====

Abipy is a Python library for the analysis of the results produced by ABINIT (http://www.abinit.org), 
an open-source program for the ab-initio calculation of the physical properties of materials 
within Density Functional Theory and Maby-Body perturbation theory.
It also provides tools to generate Abinit input files and workflows to automate 
ab-initio calculations and typical convergence studies. 
Abipy is interfaced with Pymatgen (http://www.pymatgen.org) and this allows user to 
benefit from the different tools and python objects available in the pymatgen ecosystem

Note that the majority of the post-processing tools require output files in netcdf format so
we strongly suggest to compile Abinit with netcdf support 
(Use `--with_trio_flavor="netcdf-fallback"` at configure time to activate the internal netcdf library,
to link Abinit against an external netcdf library see the configuration examples 
provided by abiconfig (https://github.com/abinit/abiconfig).

Abipy supports both Python 2.7 as well as Python >= 3.4.
Note however that Python 2.7 is more intensively tested than py3k especially at the level of workflows 
so we recommend py2.7 if you plan to run automatic calculations with abipy.

Documentation available at `abipy docpage`_.
Check out the the list of plotting scripts available in the 
`examples/plot <http://pythonhosted.org/abipy/examples/plot/index.html>`_ directory.
To learn more about abipy, you can visit our example collection of `ipython notebooks
<http://nbviewer.ipython.org/github/gmatteo/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_
and the
`abipy lessons <http://nbviewer.ipython.org/github/gmatteo/abipy/blob/master/abipy/examples/notebooks/lessons/index.ipynb>`_.

Abipy is free to use. However, we also welcome your help to improve this library by making your own contributions.  
Please report any bugs and issues at AbiPy's `Github page <https://github.com/gmatteo/abipy>`_.

Getting AbiPy
=============

Stable version
--------------

The version at the Python Package Index (PyPI) is always the latest stable release
that can be installed with:: 

    pip install abipy

Note, however, that you may need to install pymatgen and other critical dependencies manually.
For this reason, we strongly suggest to install the required python packages through one 
of the following python distributions::

    - `Anaconda <https://continuum.io/downloads>`_
    - `Canopy <https://www.enthought.com/products/canopy>`_

We routinely use anaconda to test new developments with multiple versions of Python and multiple virtual environments.
The anaconda distribution already provides the most critical dependencies (`matplotlib`, `scipy`, `numpy`) 
in the form of precompiled packages and `netcdf4` can be easily installed with::

    conda install netcdf4

Additional information on the steps required to install abipy with anaconda are available
in the `anaconda howto <http://pythonhosted.org/abipy/users/howto_anaconda.html>`_.

We are also working with the [Spack](https://github.com/LLNL/spack) community
to provide packages for abipy and Abinit in order to faciliate the installation on large supercomputing centers.

Advanced users who need to compile a local version of the python interpreter and install the abipy dependencies
manually can consult this `howto <http://pythonhosted.org/abipy/users/howto_compile_python_and_bootstrap_pip.html>`_.

Developmental version
---------------------

The developmental version is at the abipy's `Github repo <https://github.com/gmatteo/abipy>`_. 
The Github version include test files for complete unit testing. 
After cloning the source, type::

    python setup.py install

or::

    python setup.py develop

to install the package in developmental mode.
This is the recommended approach, especially if you are planning to implement new features.

To run the suite of unittests, make sure you have py.test (nose) installed and then just type::

    py.test 

or::

    nosetests

in the abipy root directory. 
Unit tests require two additional packages that can be installed with

   `pip install nose-exclude scripttest``

Requirements
============

All required dependencies should be automatically taken care of if you install abipy using conda (pip). 
Otherwise, these packages should be available on `PyPI <http://pypi.python.org>`_.

Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

1. wxPython - For the GUI 
2. wxmplot

The directory `abipy.gui.demos` contains demos that can be used to test the installation 
(run the script `runall.py` to have an overview of the different graphical interfaces).

Using abipy
===========

Basic usage
-----------

There are a variety of ways to use abipy, and most of them are illustrated in the `abipy/examples` directory.
Below is a brief description of the different directories found there:

  * plot - scripts showing how to produce plots with matplotlib

  * notebooks - juptyer notebooks 
    (use `jupyter notebook FILE` to open the notebook in your browser, usa `conda install jupyter` to install
    the package)

The directory `abipy/data/runs` contains python scripts that can be used to automate typical ab-initio calculations.

Examples of the basic capabilities can be found in the 
`example page <http://pythonhosted.org/abipy/examples/index.html>`_ of the  official documentation.

If the examples stops with the error message::
    
    "ValueError: unknown locale: UTF-8"

add the following line to your `.bashrc` file inside your home (`.profile` if MacOSx)::

    export LC_ALL=C

reload the environment with `source ~/.bashrc` and rerun.

Advanced Usage
--------------

Users are strongly encouraged to explore the detailed `api docs <http://pythonhosted.org/abipy/api/index.html>`_.

License
=======

Abipy is released under the GNU GPL License. The terms of the license are as follows::

    abipy is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2.1 of the License, or
    (at your option) any later version.

    abipy is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with abipy.  
    If not, see <http://www.gnu.org/licenses/>.


.. _`abipy docpage` : http://pythonhosted.org/abipy

:Repository: https://github.com/gmatteo/abipy
:Author: Matteo Giantomassi (http://github.com/gmatteo)

.. image:: https://pypip.in/v/abipy/badge.png
        :target: https://pypi.python.org/pypi/abipy

.. image:: https://coveralls.io/repos/gmatteo/abipy/badge.svg
  :target: https://coveralls.io/r/gmatteo/abipy


About
=====
Abipy is an open-source library for the analysis of the results produced by ABINIT (http://www.abinit.org), 
an open-source program for the ab-initio calculation of the physical properties of materials 
within Density Functional Theory (DFT).
Abipy is written in Python and is designed with the philosophy that you should be able to create 
simple plots with just a few commands.

Documentation available at http://pythonhosted.org/abipy/
Check out the the list of plotting scripts available in the 
`examples/plot <http://pythonhosted.org/abipy/examples/plot/index.html>`_ directory.
To learn more about abipy, you can visit our example collection of ipython `notebooks 
<http://nbviewer.ipython.org/github/gmatteo/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_
and the abipy 
`lessons <http://nbviewer.ipython.org/github/gmatteo/abipy/blob/master/abipy/examples/notebooks/lessons/index.ipynb>`_.

Abipy is free to use. However, we also welcome your help to improve this library by making your own contributions.  
These contributions can be in the form of additional tools or modules you develop, or even simple things 
such as bug reports. 
Please report any bugs and issues at abipy's `Github page <https://github.com/gmatteo/abipy>`_. 

Getting abipy
=============

Stable version
--------------

The version at the Python Package Index (PyPI) is always the latest stable
release that will be hopefully, be relatively bug-free. 
The easiest way to install abipy is to use easy_install or pip, as follows::

    easy_install abipy

or::

    pip install abipy


**Note**: You may need to install pymatgen before installing abipy as abipy depends on pymatgen 
Besides, abipy required additional dependencies such as netcdf4 and wxpython for the graphical interface.
Users who want to use the graphical interface are suggested to install the required python packages (wxpython)
through one of the following python distributions::

    #. `Canopy <https://www.enthought.com/products/canopy>`_

    #. `Anaconda <http://continuum.io/downloads`_


Developmental version
---------------------

The developmental version is at the abipy's `Github repo <https://github.com/gmatteo/abipy>`_. 
The Github version include test files for complete unit testing. 
After cloning the source, you can type::

    python setup.py install

or to install the package in developmental mode::

    python setup.py develop

To run the suite of unittests, make sure you have nose installed and then just type::

    nosetests

in the abipy root directory.


Requirements
============

All required dependencies should be automatically taken care of if you install abipy using easy_install or pip. 
Otherwise, these packages should be available on `PyPI <http://pypi.python.org>`_.

  #. Python 2.7 required (Python 3.0+ not supported) 
  #. pymatgen
  #. numpy 
  #. matplotlib 
  #. scipy 
  #. netCDF4
  #. pyYaml

The following packages are much easier to install with easy_install:

  #. pyYaml
  #. pyCifRW
  #. pyhull
  #. PyDispatcher

for netcdf4 and hdf see http://www.unidata.ucar.edu/software/netcdf/docs/build_default.html


pyhull:
 export CC=gcc 
 and optionally:
 export CFLAGS='-O3 -g' 


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

  * notebooks - ipython notebooks 
    (use `ipython notebook FILE` to open the notebook in your browser)

The directory `abipy/data/runs` contains scripts that can be used to automate typical ab-initio calculations.

Examples of the basic capabilities can be found in the 
`example page http://pythonhosted.org/abipy/examples/index.html`_ of the  official documentation.

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

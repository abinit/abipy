
Abipy is an open-source library for analyzing the results produced by ABINIT (http://www.abinit.org), 
an open-source program for the ab-initio calculation of the physical properties of materials within Density Functional Theory (DFT).
Abipy is written in Python and is designed with the philosophy that you should be able to create simple plots with just a few commands.

Abipy is free to use. However, we also welcome your help to improve this library by making your own contributions.  
These contributions can be in the form of additional tools or modules you develop, or even simple things such as bug reports. 
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

Developmental version
---------------------

The developmental version is at the abipy's `Github repo <https://github.com/gmatteo/abipy>`_. 
The developmental version is likely to be more unstable, but may contain new features. 
The Github version include test files as well for complete unit testing. 
After cloning the source, you can type::

    python setup.py install

or to install the package in developmental mode::

    python setup.py develop

Requirements
============

All required dependencies should be automatically taken care of if you
install abipy using easy_install or pip. 
Otherwise, these packages should be available on `PyPI <http://pypi.python.org>`_.

1. Python 2.7+ required. 

2. pymatgen 2.7.3+

3. numpy - For array, matrix and other numerical manipulations. 

4. matplotlib 1.1+

5. scipy 0.10+

6. netCDF4

Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

1. nose - For complete unittesting.

Using abipy
===========

Basic usage
-----------

Here are some quick examples of the core capabilities and objects:

.. code-block:: pycon

    >>> import abipy as ab
    >>>

The above example illustrates only the most basic capabilities of abipy.

.. note:: Examples

    A good way to explore the functionality of abipy is to look at examples.
    We have created a `Github wiki page <https://github.com/gmatteo/abipy/wiki>`_ 
    to allow users to share their Github gists (essentially mini git repos of scripts)
    performing various kinds of functions with abipy. 
    Please feel free to check them out and we welcome your contributions as well!

Advanced Usage
--------------

Users are strongly encouraged to explore the detailed `api docs <http://pythonhosted.org/abipy/api/index..html>`_.

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

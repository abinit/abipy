.. image:: _static/abipy_logo.png
   :width: 300 px
   :alt: abipy
   :align: center


===================
Abipy Documentation
===================

.. image:: https://pypip.in/v/abipy/badge.png
        :target: https://pypi.python.org/pypi/abipy

.. image:: https://img.shields.io/travis/gmatteo/abipy/master.svg
        :target: https://travis-ci.org/gmatteo/abipy

.. image:: https://img.shields.io/badge/license-GPL-blue.svg


.. htmlonly::

    :Release: |version|
    :Date: |today|

AbiPy is a Python library to analyze the results produced by `ABINIT <http://www.abinit.org>`_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
It also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with `Pymatgen <http://www.pymatgen.org>`_ and this allows users to
benefit from the different tools and python objects available in the pymatgen ecosystem

AbiPy supports both Python 2.7 as well as Python >= 3.4.
Note however that Python 2.7 is more intensively tested than py3k especially at the level of workflows
so we still recommend py2.7 if you plan to run automatic calculations with AbiPy.

AbiPy can be used in conjunction with
`ipython <https://ipython.org/index.html>`_ and `jupyter <http://jupyter.org/>`_
thus providing a powerful still user-friendly environment for data analysis and visualization.


Getting Started
===============

.. toctree::
   :maxdepth: 1

   users/index.rst

Basic Usage
===========

Matplotlib Examples
===================

.. toctree::
   :maxdepth: 1

   examples/index.rst

.. toctree::
   :maxdepth: 1

   users/index.rst

API documentation
=================

.. toctree::
   :maxdepth: 1

   api/index.rst
   devel/index.rst

License
=======

AbiPy is released under the GPL License. The terms of the license are as follows:

.. literalinclude:: ../LICENSE

Indices and tables
==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`

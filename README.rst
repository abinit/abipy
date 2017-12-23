.. :Repository: https://github.com/abinit/abipy
.. :Author: Matteo Giantomassi (http://github.com/abinit)

.. image:: https://badge.fury.io/py/abipy.svg
        :target: https://badge.fury.io/py/abipy

.. image:: https://travis-ci.org/abinit/abipy.svg?branch=develop
        :target: https://travis-ci.org/abinit/abipy

.. image:: https://coveralls.io/repos/github/abinit/abipy/badge.svg?branch=develop
        :target: https://coveralls.io/github/abinit/abipy?branch=develop

.. image:: https://img.shields.io/badge/license-GPL-blue.svg


About
=====

AbiPy is a Python library to analyze the results produced by `ABINIT <https://www.abinit.org>`_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
It also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with `Pymatgen <http://www.pymatgen.org>`_ and this allows users to
benefit from the different tools and python objects available in the pymatgen ecosystem

Official documentation of the stable version available at the `abipy docpage`_.

AbiPy can be used in conjunction with `matplotlib <http://matplotlib.org>`_, `pandas <http://pandas.pydata.org>`_,
`ipython <https://ipython.org/index.html>`_ and `jupyter <http://jupyter.org/>`_
thus providing a powerful and user-friendly environment for data analysis and visualization.
Check out the list of plotting scripts available in our
`examples/plot <http://pythonhosted.org/abipy/examples/plot/index.html>`_ gallery.
To learn more about the integration between jupyter and AbiPy, visit our collection of `notebooks
<http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_
and the
`AbiPy lessons <https://nbviewer.jupyter.org/github/abinit/abitutorials/tree/master/abitutorials/index.ipynb>`_.

AbiPy supports both Python 2.7 as well as Python >= 3.4.
Python 2.7 is more intensively tested than py3k especially at the level of workflows
so we still recommend py2.7 if you plan to run automatic calculations with AbiPy.

Note that the majority of the post-processing tools available in AbiPy require output files in
``netcdf`` format so we strongly suggest to compile Abinit with netcdf support
(use ``--with_trio_flavor="netcdf-fallback"`` at configure time to activate the internal netcdf library,
to link Abinit against an external netcdf library please consult the configuration examples
provided by `abiconfig <https://github.com/abinit/abiconfig>`_.

AbiPy is free to use. However, we also welcome your help to improve this library by making your own contributions.
Please report any bugs and issues at AbiPy's `Github page <https://github.com/abinit/abipy>`_.

Getting AbiPy
=============

Stable version
--------------

The version at the Python Package Index (PyPI) is always the latest stable release
that can be installed with::

    pip install abipy

Note that you may need to install pymatgen and other critical dependencies manually.
In this case, please consult the detailed installation instructions provided by the
`pymatgen howto <http://pymatgen.org/index.html#standard-install>`_ to install pymatgen 
and then follow the instructions in `our howto <http://pythonhosted.org/abipy/installation.html>`_.

The installation process is greatly simplified if you install the required 
python packages through one of the following python distributions:

  * `Anaconda <https://continuum.io/downloads>`_

  * `Canopy <https://www.enthought.com/products/canopy>`_

We routinely use ``conda`` to test new developments with multiple versions of Python and multiple virtual environments.
The anaconda distribution already provides the most critical dependencies (``matplotlib``, ``scipy``, ``numpy``)
in the form of pre-compiled packages and ``netcdf4`` can be easily installed with::

    conda install netcdf4

Additional information on the steps required to install AbiPy with anaconda are available
in the `anaconda howto <http://pythonhosted.org/abipy/installation.html>`_.
We are also working with the `Spack <https://github.com/LLNL/spack>`_ community
to provide packages for AbiPy and Abinit in order to facilitate the installation on large supercomputing centers.

---------------------
Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

``ipython``

    Required to interact with the AbiPy/Pymatgen objects in the ipython shell
    (strongly recommended, already provided by ``conda``).

``jupyter`` and ``nbformat``

    Required to generate jupyter notebooks.
    Install these two packages with ``conda install jupyter nbformat`` or use ``pip``.
    Recommended but you will also need a web browser to open the notebook.

``wxPython`` and ``wxmplot`` for the GUI

    Use ``conda install wxpython``
    The directory ``abipy.gui.demos`` contains demos that can be used to test the installation.
    of the GUI (run the script ``runall.py`` to have an overview of the different graphical interfaces).

Developmental version
---------------------

Getting the developmental version of AbiPy is easy. You can clone it from the 
`github repository <https://github.com/abinit/abipy>`_ using this command:

.. code-block:: console

   $ git clone https://github.com/abinit/abipy

After cloning the repository, type::

    $ python setup.py install

or alternately::

    $ python setup.py develop

to install the package in developmental mode 
(this is the recommended approach, especially if you are planning to implement new features).

The documentation of the **developmental** version is hosted on `github pages <http://abinit.github.io/abipy>`_.

The Github version include test files for complete unit testing.
To run the suite of unit tests, make sure you have ``pytest`` (recommended) 
or ``nose`` installed and then just type::

    $ pytest

or::

    $ nosetests

in the AbiPy root directory.
Unit tests require two additional packages that can be installed with::

   $ pip install nose-exclude scripttest

Note that several unit tests check the integration between AbiPy and Abinit. 
In order to run the tests, you need a working set of Abinit executables and  
a ``manager.yml`` configuration file.
A pre-compiled sequential version of Abinit for Linux and OSx can be installed directly from the anaconda cloud with::

    $ conda install abinit -c abinit

For further information on the syntax of the configuration file, please consult the 
`workflows <http://pythonhosted.org/abipy/workflows.html>`_ section.

Contributing to AbiPy is relatively easy.
Just send us a `pull request <https://help.github.com/articles/using-pull-requests/>`_.
When you send your request, make ``develop`` the destination branch on the repository
AbiPy uses the `Git Flow <http://nvie.com/posts/a-successful-git-branching-model/>`_ branching model.
The ``develop`` branch contains the latest contributions, and ``master`` is always tagged and points
to the latest stable release.

Using AbiPy
===========

Basic usage
-----------

There are a variety of ways to use AbiPy, and most of them are illustrated in the ``abipy/examples`` directory.
Below is a brief description of the different directories found there:

  * `plot <http://pythonhosted.org/abipy/examples/plot/index.html>`_

    scripts showing how to produce plots with ``matplotlib``

  * `notebooks <http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_

    jupyter notebooks
    (use ``jupyter notebook FILE`` to open the notebook in your browser,
    use ``conda install jupyter`` to install the package)

The directory ``abipy/data/runs`` contains python scripts that can be used to automate typical ab-initio calculations.

Command line tools
------------------

The following scripts can be invoked directly from the terminal:

  * ``abiopen.py``
  * ``abistruct.py``
  * ``abicomp.py``
  * ``abiview.py``
  * ``abicheck.py``

For further information, please consult the `official documentation <http://pythonhosted.org/abipy/scripts.html>`_.

Advanced Usage
--------------

Users are strongly encouraged to explore the detailed `api docs <http://pythonhosted.org/abipy/api/index.html>`_.

License
=======

AbiPy is released under the GNU GPL license. For more details see the LICENSE file.

.. _`abipy docpage` : http://pythonhosted.org/abipy

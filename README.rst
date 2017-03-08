.. :Repository: https://github.com/abinit/abipy
.. :Author: Matteo Giantomassi (http://github.com/abinit)

.. image:: https://pypip.in/v/abipy/badge.png
        :target: https://pypi.python.org/pypi/abipy

.. image:: https://img.shields.io/travis/gmatteo/abipy/master.svg
        :target: https://travis-ci.org/gmatteo/abipy

.. image:: https://img.shields.io/badge/license-GPL-blue.svg


About
=====

AbiPy is a Python library to analyze the results produced by `ABINIT <http://www.abinit.org>`_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
It also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with `Pymatgen <http://www.pymatgen.org>`_ and this allows users to
benefit from the different tools and python objects available in the pymatgen ecosystem





Official documentation of the stable version available at `abipy docpage`_.

AbiPy can be used in conjunction with `matplotlib <http://matplotlib.org>`_,
`ipython <https://ipython.org/index.html>`_ and `jupyter <http://jupyter.org/>`_
thus providing a powerful still user-friendly environment for data analysis and visualization.
Check out the list of plotting scripts available in the
`examples/plot <http://pythonhosted.org/abipy/examples/plot/index.html>`_ directory.
To learn more about the integration between jupyter and AbiPy, visit our collection of `notebooks
<http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_
and the
`AbiPy lessons <http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/lessons/index.ipynb>`_.

AbiPy supports both Python 2.7 as well as Python >= 3.4.
Note however that Python 2.7 is more intensively tested than py3k especially at the level of workflows
so we still recommend py2.7 if you plan to run automatic calculations with AbiPy.

Note, however, that the majority of the post-processing tools available in AbiPy require output files in
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
For this reason, we strongly suggest to install the required python packages through one
of the following python distributions:

  * `Anaconda <https://continuum.io/downloads>`_

  * `Canopy <https://www.enthought.com/products/canopy>`_

We routinely use ``conda`` to test new developments with multiple versions of Python and multiple virtual environments.
The anaconda distribution already provides the most critical dependencies (``matplotlib``, ``scipy``, ``numpy``)
in the form of pre-compiled packages and ``netcdf4`` can be easily installed with::

    conda install netcdf4

Additional information on the steps required to install AbiPy with anaconda are available
in the `anaconda howto <http://pythonhosted.org/abipy/users/howto_anaconda.html>`_.

We are also working with the `Spack <https://github.com/LLNL/spack>`_ community
to provide packages for AbiPy and Abinit in order to facilitate the installation on large supercomputing centers.

Advanced users who need to compile a local version of the python interpreter and install the AbiPy dependencies
manually can consult this `howto <http://pythonhosted.org/abipy/users/howto_compile_python_and_bootstrap_pip.html>`_.

Developmental version
---------------------

The developmental version is at the AbiPy's `Github repo <https://github.com/abinit/abipy>`_.
The Github version include test files for complete unit testing.
After cloning the source, type::

    $ python setup.py install

or::

    $ python setup.py develop

to install the package in developmental mode (this is the recommended approach, especially if you are
planning to implement new features).

To run the suite of unit tests, make sure you have py.test (nose) installed and then just type::

    $ py.test

or::

    $ nosetests

in the AbiPy root directory.
Unit tests require two additional packages that can be installed with::

   $ pip install nose-exclude scripttest

Contributing to AbiPy is relatively easy.
Just send us a `pull request <https://help.github.com/articles/using-pull-requests/>`_.
When you send your request, make ``develop`` the destination branch on the repository
AbiPy uses the `Git Flow <http://nvie.com/posts/a-successful-git-branching-model/>`_ branching model.
The ``develop`` branch contains the latest contributions, and ``master`` is always tagged and points
to the latest stable release.

Requirements
============

All required dependencies should be automatically taken care of if you are using ``conda``

Otherwise, these packages should be available on `PyPI <http://pypi.python.org>`_.

Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

  * ``matplotlib``

  * ``ipython``

    Required to interact with the AbiPy/Pymatgen objects in the ipython shell.

  * ``jupyter`` and ``nbformat``

    Required to generate jupyter notebooks to analyze data.
    Install these two packages with `conda install jupyter nbformat`

  * ``wxPython`` and ``wxmplot`` for the GUI

    Use ``conda install wxpython``
    The directory ``abipy.gui.demos`` contains demos that can be used to test the installation.
    of the GUI (run the script ``runall.py`` to have an overview of the different graphical interfaces).

Using AbiPy
===========

Basic usage
-----------

There are a variety of ways to use AbiPy, and most of them are illustrated in the ``abipy/examples`` directory.
Below is a brief description of the different directories found there:

  * ``plot``

    scripts showing how to produce plots with ``matplotlib``

  * ``notebooks``

    jupyter notebooks
    (use ``jupyter notebook FILE`` to open the notebook in your browser,
    use ``conda install jupyter`` to install the package)


The directory ``abipy/data/runs`` contains python scripts that can be used to automate typical ab-initio calculations.

The following scripts can be invoked directly from the terminal:

  * `abiopen.py`

    Script to open outputs file produced by Abinit (usually in netcdf format but
    other files are supported as well). By default the script starts an interactive ipython
    session so that one can interact with the file and call its methods.
    Alternatively, it is possible to generate automatically a jupyter notebook to execute code.

  * `abistruct.py`

    Script to analyze/export/visualize the crystal structure saved in the netcdf files produced by ABINIT.

  * `abicomp.py`

    Script to analyze/compare results stored in multiple netcdf files.
    By default the script displays the results/plots in the shell.
    Use `--ipython` to start an ipython terminal or `-nb` to generate a jupyter notebook.

  * `abicheck.py`

    This script checks that the environment on the local machine is properly configured.


Examples of the basic capabilities can be found in the
`example page <http://pythonhosted.org/abipy/examples/index.html>`_ of the  official documentation.

If the examples stops with the error message::

    "ValueError: unknown locale: UTF-8"

add the following line to your ``.bashrc`` file inside your ``$HOME`` (``.profile`` if MacOSx)::

    export LC_ALL=C

reload the environment with ``source ~/.bashrc`` and rerun.

Advanced Usage
--------------

Users are strongly encouraged to explore the detailed `api docs <http://pythonhosted.org/abipy/api/index.html>`_.

License
=======

AbiPy is released under an GNU GPL license. For more details see the LICENSE file.

.. _`abipy docpage` : http://pythonhosted.org/abipy

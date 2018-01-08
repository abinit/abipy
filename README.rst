.. :Repository: https://github.com/abinit/abipy
.. :Author: Matteo Giantomassi (http://github.com/abinit)

.. list-table::
    :stub-columns: 1
    :widths: 10 90

    * - Package
      - |pypi-version| |download-with-anaconda| |supported-versions|
    * - Continuous Integration
      - |travis-status| |coverage-status| 
    * - Documentation
      - |launch-nbviewer| |launch-binder| 


About
=====

AbiPy is a Python_ library to analyze the results produced by Abinit_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
It also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with pymatgen_ and this allows users to
benefit from the different tools and python objects available in the pymatgen ecosystem.

The official documentation of the stable version available at the `abipy docpage`_.

AbiPy can be used in conjunction with matplotlib_, pandas_, scipy_, seaborn_, ipython_ and jupyter_ notebooks
thus providing a powerful and user-friendly environment for data analysis and visualization.
Check out the list of plotting scripts available in our
`examples/plot <http://pythonhosted.org/abipy/examples/plot/index.html>`_ gallery.
To learn more about the integration between jupyter_ and AbiPy, visit `our collection of notebooks
<https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb>`_

.. AbiPy supports both Python 2.7 as well as Python >= 3.4.
.. Python 2.7 is more intensively tested than py3k especially at the level of workflows
.. so we still recommend py2.7 if you plan to run automatic calculations with AbiPy.

Note that the majority of the post-processing tools available in AbiPy require output files in
netcdf_ format so we strongly suggest to compile Abinit with netcdf support
(use ``--with_trio_flavor="netcdf-fallback"`` at configure time to activate the internal netcdf library,
to link Abinit against an external netcdf library please consult the configuration examples provided by abiconfig_).

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
python packages through `Anaconda <https://continuum.io/downloads>`_.
We routinely use conda_ to test new developments with multiple versions of Python and multiple virtual environments.
The anaconda distribution already provides the most critical dependencies (matplotlib_, scipy_, numpy_)
in the form of pre-compiled packages and netcdf4-python_ can be easily installed with::

    conda install netcdf4

Additional information on the steps required to install AbiPy with anaconda are available
in the `anaconda howto <http://pythonhosted.org/abipy/installation.html>`_.
We are also working with the spack_ community
to provide packages for AbiPy and Abinit in order to facilitate the installation on large supercomputing centers.

---------------------
Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

``ipython``

    Required to interact with the AbiPy/Pymatgen objects in the ipython_ shell
    (strongly recommended, already provided by conda_).

``jupyter`` and ``nbformat``

    Required to generate jupyter_ notebooks.
    Install these two packages with ``conda install jupyter nbformat`` or use ``pip``.
    Recommended but you will also need a web browser to open the notebook.

.. ``wxPython`` and ``wxmplot`` for the GUI
..    Use ``conda install wxpython``
..    The directory ``abipy.gui.demos`` contains demos that can be used to test the installation.
..    of the GUI (run the script ``runall.py`` to have an overview of the different graphical interfaces).

Developmental version
---------------------

Getting the developmental version of AbiPy is easy. 
You can clone it from the `github repository <https://github.com/abinit/abipy>`_ using this command::

    git clone https://github.com/abinit/abipy

After cloning the repository, type::

    python setup.py install

or alternately::

    python setup.py develop

to install the package in developmental mode 
(this is the recommended approach, especially if you are planning to implement new features).

The documentation of the **developmental** version is hosted on `github pages <http://abinit.github.io/abipy>`_.

The Github version include test files for complete unit testing.
To run the suite of unit tests, make sure you have pytest_ installed and then type::

    pytest

in the AbiPy root directory.
Unit tests require two additional packages that can be installed with::

   $ pip install nose-exclude scripttest

Note that several unit tests check the integration between AbiPy and Abinit. 
In order to run the tests, you need a working set of Abinit executables and  a ``manager.yml`` configuration file.
A pre-compiled sequential version of Abinit for Linux and OSx can be installed directly from the 
`abinit channel <https://anaconda.org/abinit>`_ with::

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

There are a variety of ways to use AbiPy, and most of them are illustrated in the :file:`abipy/examples` directory.
Below is a brief description of the different directories found there:

  * `plot <http://pythonhosted.org/abipy/examples/plot/index.html>`_

    scripts showing how to produce plots with matplotlib_

  * `notebooks <http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_

    jupyter notebooks
    (use ``jupyter notebook FILE`` to open the notebook in your browser,
    use ``conda install jupyter`` to install the package)

The directory :file:`abipy/examples/flows` contains python scripts that can be used 
to automate typical ab-initio calculations.

Command line tools
------------------

The following scripts can be invoked directly from the terminal:

* ``abicheck.py``
* ``abidoc.py``
* ``abiopen.py``
* ``abistruct.py``
* ``abicomp.py``
* ``abinp.py``
* ``abirun.py``
* ``abiview.py``

For further information, please consult the `official documentation <http://pythonhosted.org/abipy/scripts.html>`_.

Advanced Usage
--------------

Users are strongly encouraged to explore the detailed `API docs <http://pythonhosted.org/abipy/api/index.html>`_.

License
=======

AbiPy is released under the GNU GPL license. For more details see the LICENSE file.

.. _Python: http://www.python.org/
.. _Abinit: https://www.abinit.org
.. _pymatgen: http://www.pymatgen.org
.. _`abipy docpage` : http://pythonhosted.org/abipy
.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _scipy: https://www.scipy.org/
.. _seaborn: https://seaborn.pydata.org/
.. _ipython: https://ipython.org/index.html
.. _jupyter: http://jupyter.org/
.. _netcdf: https://www.unidata.ucar.edu/software/netcdf/docs/faq.html#whatisit
.. _abiconfig: https://github.com/abinit/abiconfig
.. _conda: https://conda.io/docs/
.. _netcdf4-python: http://unidata.github.io/netcdf4-python/
.. _spack: https://github.com/LLNL/spack
.. _pytest: https://docs.pytest.org/en/latest/contents.html
.. _numpy: http://www.numpy.org/


.. |pypi-version| image:: https://badge.fury.io/py/abipy.svg
    :alt: PyPi version
    :target: https://badge.fury.io/py/abipy

.. |travis-status| image:: https://travis-ci.org/abinit/abipy.svg?branch=develop
    :alt: Travis status
    :target: https://travis-ci.org/abinit/abipy

.. |coverage-status| image:: https://coveralls.io/repos/github/abinit/abipy/badge.svg?branch=develop
    :alt: Coverage status
    :target: https://coveralls.io/github/abinit/abipy?branch=develop

.. |download-with-anaconda| image:: https://anaconda.org/abinit/abipy/badges/installer/conda.svg   
    :alt: Download with Anaconda
    :target: https://conda.anaconda.org/abinit

.. |launch-binder| image:: https://mybinder.org/badge.svg 
    :alt: Launch binder
    :target: https://mybinder.org/v2/gh/abinit/abipy/develop

.. |launch-nbviewer| image:: https://img.shields.io/badge/render-nbviewer-orange.svg
    :alt: Launch nbviewer
    :target: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/abipy.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/abipy

.. |requires| image:: https://requires.io/github/abinit/abipy/requirements.svg?branch=develop
     :target: https://requires.io/github/abinit/abipy/requirements/?branch=develop
     :alt: Requirements Status

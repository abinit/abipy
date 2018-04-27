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
      - |docs-stable| |docs-devel| |launch-nbviewer| |launch-binder| 

About
=====

AbiPy is a Python_ library to analyze the results produced by Abinit_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
It also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with pymatgen_ and this allows users to
benefit from the different tools and python objects available in the pymatgen ecosystem.

The official documentation of the stable version is available at the `abipy docpage`_.
while the documentation of the **developmental** version is hosted on `github pages <http://abinit.github.io/abipy>`_.

AbiPy can be used in conjunction with matplotlib_, pandas_, scipy_, seaborn_, ipython_ and jupyter_ notebooks
thus providing a powerful and user-friendly environment for data analysis and visualization.
Check out our `gallery of plotting scripts <http://abinit.github.io/abipy/gallery/index.html>`_
and the `gallery of AbiPy workflows <http://abinit.github.io/abipy/flow_gallery/index.html>`_.

To learn more about the integration between jupyter_ and AbiPy, visit `our collection of notebooks
<https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb>`_
or click the **Launch Binder** badge to start a Docker image with Abinit, AbiPy and all the other python dependencies
required to run the code inside the jupyter notebooks.
The notebook will be opened in your browser after building.

AbiPy is free to use. However, we also welcome your help to improve this library by making your own contributions.
Please report any bugs and issues at AbiPy's `Github page <https://github.com/abinit/abipy>`_.

.. important::

    Note that the majority of the post-processing tools available in AbiPy require output files in
    netcdf_ format so we **strongly** suggest to compile Abinit with netcdf support
    (use ``--with_trio_flavor="netcdf-fallback"`` at configure time to activate the internal netcdf library,
    to link Abinit against an external netcdf library please consult the configuration examples provided by abiconfig_).

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
and then follow the instructions in `our howto <http://abinit.github.io/abipy/installation>`_.

The installation process is greatly simplified if you install the required 
python packages through `Anaconda <https://continuum.io/downloads>`_.
We routinely use conda_ to test new developments with multiple Python versions and multiple virtual environments.
The anaconda distribution already provides the most critical dependencies (matplotlib_, scipy_, numpy_, netcdf4-python_)
in the form of pre-compiled packages that can be easily installed with e.g.::

    conda install numpy scipy netcdf4

To install AbiPy with conda, download the `minconda installer <https://conda.io/miniconda.html>`_
(select python3.6 and the version corresponding to your operating system).
Create a new conda_ environment (let's call it ``abipy3.6``) based on python3.6 with::

    conda create --name abipy3.6 python=3.6

and activate it with::

    source activate abipy3.6

You should see the name of the conda environment in the shell prompt.

Now add ``conda-forge``, ``matsci`` and ``abinit`` to your conda channels with::

    conda config --add channels conda-forge
    conda config --add channels matsci
    conda config --add channels abinit

These are the channels from which we will download pymatgen, abipy and abinit.

Finally, install AbiPy from the abinit-channel_ with::

    conda install abipy -c abinit

One of the big advantages of conda over pip is that conda can also install 
libraries and executables written in Fortran.
For example, one can easily install abinit inside the conda environment with::

    conda install abinit -c abinit
    abinit --version

Additional information on the steps required to install AbiPy with anaconda are available
in the `anaconda howto <http://abinit.github.io/abipy/installation#anaconda-howto>`_.

We are also collaborating with the spack_ community
to provide packages for AbiPy and Abinit in order to facilitate the installation on large supercomputing centers.

Developmental version
---------------------

Getting the developmental version of AbiPy is easy. 
Clone the `github repository <https://github.com/abinit/abipy>`_ with::

    git clone https://github.com/abinit/abipy

For pip, use::

    pip install -r requirements.txt
    pip install -r requirements-optional.txt

If you are using conda_,  create a new environment (``abipy3.6``) based on python3.6 with::

    conda create -n abipy3.6 python=3.6
    source activate abipy3.6

Add ``conda-forge``, ``matsci`` and ``abinit`` to your channels with::

    conda config --add channels conda-forge
    conda config --add channels matsci
    conda config --add channels abinit

and install the AbiPy dependencies with::

    conda install --file ./requirements.txt
    conda install --file ./requirements-optional.txt

Once the requirements have been installed (either with pip or conda), execute::

    python setup.py install

or alternately::

    python setup.py develop

to install the package in developmental mode 
(this is the recommended approach, especially if you are planning to implement new features).

The Github version include test files for complete unit testing.
To run the suite of unit tests, make sure you have pytest_ installed and then type::

    pytest

in the AbiPy root directory.
Unit tests require ``scripttest`` that can be installed with::

    pip install scripttest

Note that several unit tests check the integration between AbiPy and Abinit. 
In order to run the tests, you need a working set of Abinit executables and  a ``manager.yml`` configuration file.
A pre-compiled sequential version of Abinit for Linux and OSx can be installed directly from the 
abinit-channel_ with::

    conda install abinit -c abinit

For further information on the syntax of the configuration file, please consult the 
`workflows <http://abinit.github.io/abipy/workflows/taskmanager.html>`_ section.

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

  * `examples/plot <http://abinit.github.io/abipy/gallery/index.html>`_

    Scripts showing how to read data from netcdf files and produce plots with matplotlib_

  * `examples/flows <http://abinit.github.io/abipy/flow_gallery/index.html>`_.

    Scripts showing how to generate an AbiPy flow, run the calculation and use ipython to analyze the data.

Additional jupyter notebooks with the Abinit tutorials written with AbiPy are available in the
`abitutorial repository <https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb>`_

Users are strongly encouraged to explore the detailed `API docs <http://abinit.github.io/abipy/api/index.html>`_.

Command line tools
------------------

The following scripts can be invoked directly from the terminal:

* ``abiopen.py``    Open file inside ipython.
* ``abistruct.py``  Swiss knife to operate on structures.
* ``abiview.py``    Visualize results from file.
* ``abicomp.py``    Compare results extracted from multiple files.
* ``abicheck.py``   Validate integration between AbiPy and Abinit
* ``abirun.py``     Execute AbiPy flow from terminal.
* ``abidoc.py``     Document Abinit input variables and Abipy configuration files.
* ``abinp.py``      Build input files (simplified interface for the AbiPy factory functions).

Use ``SCRIPT --help`` to get the list of supported commands and 
``SCRIPT COMMAND --help`` to get the documentation for ``COMMAND``.
For further information, please consult the `official documentation <http://abinit.github.io/abipy/scripts/index.html>`_.

License
=======

AbiPy is released under the GNU GPL license. For more details see the LICENSE file.

.. _Python: http://www.python.org/
.. _Abinit: https://www.abinit.org
.. _abinit-channel: https://anaconda.org/abinit
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

.. |docs-stable| image:: https://img.shields.io/badge/docs-stable_version-blue.svg
     :alt: Documentation stable version
     :target: http://pythonhosted.org/abipy/

.. |docs-devel| image:: https://img.shields.io/badge/docs-devel_version-ff69b4.svg
     :alt: Documentation development version
     :target: http://abinit.github.io/abipy

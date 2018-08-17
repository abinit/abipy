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
      - |docs-github| |launch-nbviewer| |launch-binder| 

About
=====

AbiPy is a Python_ library to analyze the results produced by Abinit_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
It also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with pymatgen_ and this allows users to
benefit from the different tools and python objects available in the pymatgen ecosystem.

The official documentation of the stable version is available at the `abipy docpage`_,
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
python packages through `Anaconda <https://continuum.io/downloads>`_ (or conda). 
See `Installing conda`_ to install conda itself.
We routinely use conda_ to test new developments with multiple Python versions and multiple virtual environments.
The anaconda distribution already provides the most critical dependencies (matplotlib_, scipy_, numpy_, netcdf4-python_)
in the form of pre-compiled packages that can be easily installed with e.g.::

    conda install numpy scipy netcdf4

Create a new conda_ environment (let's call it ``abienv``) based on python3.6 with::

    conda create --name abienv python=3.6

and activate it with::

    source activate abienv

You should see the name of the conda environment in the shell prompt.

Now add ``conda-forge``, ``matsci`` and ``abinit`` to your conda channels with::

    conda config --add channels conda-forge
    conda config --add channels matsci
    conda config --add channels abinit

These are the channels from which we will download pymatgen, abipy and abinit.

Finally, install AbiPy from the abinit-channel_ with::

    conda install abipy -c abinit

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

If you are using conda_ (see `Installing conda`_ to install conda itself),  create a new environment (``abienv``) 
based on python3.6 with::

    conda create -n abienv python=3.6
    source activate abienv

Add ``conda-forge``, ``matsci`` and ``abinit`` to your channels with::

    conda config --add channels conda-forge
    conda config --add channels matsci
    conda config --add channels abinit

and install the AbiPy dependencies with::

    conda install --file ./requirements.txt
    conda install --file ./requirements-optional.txt

The second command is needed for Jupyter only.
Once the requirements have been installed (either with pip or conda), execute::

    python setup.py install

or alternately::

    python setup.py develop

to install the package in developmental mode. 
This is the recommended approach, especially if you are planning to implement new features.

Note, however, that the developmental version of AbiPy is kept in sync with the
developmental version of pymatgen thus ```python setup.py develop``` may 
try to download new versions from the PyPi portal and then fail with e.g. the error message::

    ...
    processing dependencies for abipy==0.6.0.dev0
    error: scipy 1.0.0 is installed but scipy>=1.0.1 is required by {'pymatgen'}

due to inconsistent dependencies.
To solve the problem, use conda to update scipy to a version >= 1.0.1 with::

    conda install "scipy>=1.0.1"

then issue again python setup.py develop. If this fails, supposing you were upgrading abipy inside an already existing conda environment, try to restart by creating from scratch a fresh conda environment, see above.

Use::

    conda info pymatgen

to display information about the installed version of pymatgen.

Also note that the BLAS/Lapack libraries provided by conda have multithreading support activated by default.
Each process will try to use all of the cores on your machine, which quickly overloads things 
if there are multiple processes running. 
(Also, this is a shared machine, so it is just rude behavior in general).
To disable multithreading, add these lines to your ~/.bash_profile::

    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1

and then activate these settings with::

    source ~/.bash_profile

The Github version include test files for complete unit testing.
To run the suite of unit tests, make sure you have pytest_ installed and then type::

    pytest

in the AbiPy root directory. A quicker check might be obtained with:: 

    pytest abipy/core/tests -v

Unit tests require ``scripttest`` that can be installed with::

    pip install scripttest

Two tests rely on the availability of a 
`pymatgen PMG_MAPI_KEY <http://pymatgen.org/usage.html#setting-the-pmg-mapi-key-in-the-config-file>` in ~/.pmgrc.yaml.

Note that several unit tests check the integration between AbiPy and Abinit. 
In order to run the tests, you will need a working set of Abinit executables and  a ``manager.yml`` configuration file.

Contributing to AbiPy is relatively easy.
Just send us a `pull request <https://help.github.com/articles/using-pull-requests/>`_.
When you send your request, make ``develop`` the destination branch on the repository
AbiPy uses the `Git Flow <http://nvie.com/posts/a-successful-git-branching-model/>`_ branching model.
The ``develop`` branch contains the latest contributions, and ``master`` is always tagged and points
to the latest stable release.


Installing Abinit
=================

One of the big advantages of conda over pip is that conda can also install
libraries and executables written in Fortran.
A pre-compiled sequential version of Abinit for Linux and OSx can be installed directly from the
abinit-channel_ with::

    conda install abinit -c abinit

Otherwise, follow the usual abinit installation instructions, and make sure abinit can be run with the command::

    abinit --version


Configuration files for Abipy
=============================

In order to run the Abipy tests, you will need a ``manager.yml`` configuration file.
For a detailed description of the syntax used in this configuration file
please consult the `TaskManager documentation <http://abinit.github.io/abipy/workflows/taskmanager.html>`_.

At this stage, for the purpose of checking the installation, you might 
take the ``shell_nompi_manager.yml`` file from the ``abipy/data/managers`` directory
of this repository, and copy it with new name ``manager.yml`` to your `$HOME/.abinit/abipy` directory.
Open this file and make sure that the ``pre_run`` section contains the shell commands
needed to setup the environment before launching Abinit (e.g. Abinit is in $PATH), unless it is available from the environment (e.g. conda).

To complete the configuration files for Abipy, you might also copy the ``simple_scheduler.yml`` file from the same directory, 
and copy it with name ``scheduler.yml``. Modifications are needed if you are developer.

Checking the installation
=========================

Now open the python interpreter and import the following three modules
to check that the python installation is OK::

    import spglib
    import pymatgen
    from abipy import abilab

then quit the interpreter.

The Abinit executables are placed inside the anaconda directory associated to the ``abienv`` environment::

    which abinit
    /Users/gmatteo/anaconda3/envs/abienv/bin/abinit

To perform a basic validation of the build, execute::

    abinit -b

Abinit should echo miscellaneous information, starting with::

    DATA TYPE INFORMATION: 
    REAL:      Data type name: REAL(DP) 
               Kind value:      8
               Precision:      15

and ending with:: 

    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Default optimizations:
      --- None ---


    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If successful, one can start to use the AbiPy scripts from the command line to analyze the output results.
Execute::

    abicheck.py

You should see (with minor changes)::

    $ abicheck.py
    AbiPy Manager:
    [Qadapter 0]
    ShellAdapter:localhost
    Hardware:
       num_nodes: 2, sockets_per_node: 1, cores_per_socket: 2, mem_per_node 4096,
    Qadapter selected: 0

    Abinitbuild:
    Abinit Build Information:
        Abinit version: 8.8.2
        MPI: True, MPI-IO: True, OpenMP: False
        Netcdf: True

    Abipy Scheduler:
    PyFlowScheduler, Pid: 19379
    Scheduler options: {'weeks': 0, 'days': 0, 'hours': 0, 'minutes': 0, 'seconds': 5}

    Installed packages:
    Package         Version
    --------------  ---------
    system          Darwin
    python_version  3.6.5
    numpy           1.14.3
    scipy           1.1.0
    netCDF4         1.4.0
    apscheduler     2.1.0
    pydispatch      2.0.5
    yaml            3.12
    pymatgen        2018.6.11


    Abipy requirements are properly configured

If the script fails with the error message::

    Abinit executable does not support netcdf
    Abipy requires Abinit version >= 8.0.8 but got 0.0.0

it means that your environment is not property configured or that there's a problem
with the binary executable.
In this case, look at the files produced in the temporary directory of the flow.
The script reports the name of the directory, something like::

    CRITICAL:pymatgen.io.abinit.tasks:Error while executing /var/folders/89/47k8wfdj11x035svqf8qnl4m0000gn/T/tmp28xi4dy1/job.sh

Check the `job.sh` script for possible typos, then search for possible error messages in `run.err`.

The last test consists in executing a small calculation with AbiPy and Abinit.
Inside the shell, execute::

    abicheck.py --with-flow

to run a GS + NSCF band structure calculation for Si.
If the software stack is properly configured, the output should end with::

    Work #0: <BandStructureWork, node_id=313436, workdir=../../../../var/folders/89/47k8wfdj11x035svqf8qnl4m0000gn/T/tmpygixwf9a/w0>, Finalized=True
      Finalized works are not shown. Use verbose > 0 to force output.

    all_ok reached

    Submitted on: Sat Jul 28 09:14:28 2018
    Completed on: Sat Jul 28 09:14:38 2018
    Elapsed time: 0:00:10.030767
    Flow completed successfully

    Calling flow.finalize()...

    Work #0: <BandStructureWork, node_id=313436, workdir=../../../../var/folders/89/47k8wfdj11x035svqf8qnl4m0000gn/T/tmpygixwf9a/w0>, Finalized=True
      Finalized works are not shown. Use verbose > 0 to force output.

    all_ok reached


    Test flow completed successfully

Great, if you've reached this part it means that you've installed AbiPy and Abinit on your machine!
We can finally start to run the scripts in this repo or use one of the AbiPy script to analyze  the results.


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
`abitutorial repository <https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb>`_.

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

For further information, please consult the `scripts docs <http://abinit.github.io/abipy/scripts/index.html>`_ section.


Installing conda
================

A brief install guide, in case you have not yet used conda ... For a more extensive description, see our
`Anaconda Howto <http://abinit.github.io/abipy/installation#anaconda-howto>`_.

Download the `miniconda installer <https://conda.io/miniconda.html>`_.
Select python3.6 and the version corresponding to your operating system.

As an example, if you are a Linux user, download and install `miniconda` on your local machine with::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

while for MacOSx use::

    curl -o https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

Answer ``yes`` to the question::

    Do you wish the installer to prepend the Miniconda3 install location
    to PATH in your /home/gmatteo/.bashrc ? [yes|no]
    [no] >>> yes

Source your ``.bashrc`` file to activate the changes done by ``miniconda`` to your ``$PATH``::

    source ~/.bashrc

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

.. |docs-github| image:: https://img.shields.io/badge/docs-ff69b4.svg
     :alt: AbiPy Documentation
     :target: http://abinit.github.io/abipy

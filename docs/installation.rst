=============
Getting AbiPy
=============

.. contents::
   :backlinks: top

--------------
Stable version
--------------

The version on the `Python Package Index <https://pypi.org/project/abipy/>`_ (PyPI) is always
the latest **stable** release. You can install it in user mode with::

    pip install abipy --user

You may need to install some optional dependencies manually.
If so, follow the detailed installation instructions in the
`pymatgen howto <https://pymatgen.org/installation.html>`_.

Installation is much simpler if you get the required
packages through the Anaconda_ distribution.
We routinely use conda_ to test new developments with multiple Python versions and virtual environments.
Anaconda already provides the most critical dependencies (numpy_, scipy_, matplotlib_, netcdf4-python_)
as pre-compiled packages that can be installed with, for example::

    conda install numpy scipy netcdf4

For more details on installing AbiPy with Anaconda, see the :ref:`anaconda_howto` section
and the `conda-based-install <http://pymatgen.org/installation.html#conda-based-install>`_
section of the pymatgen_ documentation.

We are also working with the spack_ community
to provide AbiPy and Abinit packages that make installation easier on large supercomputing centers.

Advanced users who need to compile their own Python interpreter and install the AbiPy dependencies
manually can consult the :ref:`howto_compile_python_and_bootstrap_pip` section.

---------------------
Optional dependencies
---------------------

The following libraries are needed only for certain features:

ipython_

    Required to interact with AbiPy/pymatgen objects in the ipython shell
    (strongly recommended, already provided by conda_).

jupyter_ and nbformat_

    Required to generate jupyter notebooks (recommended).
    Install both packages with ``conda install jupyter nbformat`` or with pip_.
    You will also need a web browser to open the notebooks.

.. _anaconda_howto:

--------------
Anaconda Howto
--------------

Download the Anaconda installer for your OS from the `official website <https://www.continuum.io/downloads>`_.
If you are installing Anaconda on a cluster, you may find it convenient to download the installer
directly from the terminal with ``wget``.

Run the bash script in the terminal and follow the instructions on screen.
By default, the installer creates an ``anaconda`` directory in your home
and adds a line to your ``.bashrc`` to make the Anaconda executables available.
Once the installation is complete, activate the ``base`` environment with::

    source ~/anaconda/bin/activate base

The output of ``which python`` should now show that you are using the Python interpreter provided by Anaconda.

Use the conda_ command-line interface to install packages that are not included in the official distribution.
For example, you can install ``pyyaml`` and ``netcdf4`` with::

    conda install pyyaml netcdf4

If a package is not available in the official conda repository, you can
download it from one of the conda channels, or fall back to ``pip install`` if no conda package exists.

Fortunately, some conda channels provide all the dependencies needed by AbiPy.
Add ``conda-forge`` to your channels with::

    conda config --add channels conda-forge

This is the channel from which pymatgen, AbiPy and Abinit will be downloaded.

Finally, install AbiPy with::

    conda install abipy

To check the installation, open the ipython_ shell and type::

    # make sure spglib library works
    import spglib

    # make sure pymatgen is installed
    import pymatgen

    from abipy import abilab

conda_ can also create separate environments with different
versions of the Python interpreter or of other libraries.
This is very useful for keeping different versions and branches apart.
More information is available on the `official conda website <http://conda.pydata.org/docs/test-drive.html>`_.

.. _developmental_version:

---------------------
Developmental version
---------------------

Getting the developmental version of AbiPy is easy.
Clone it from our `GitHub repository <https://github.com/abinit/abipy>`_ with::

    git clone https://github.com/abinit/abipy

then, inside the repository, type::

    python -m pip install .

or, to install the package in development (editable) mode::

    python -m pip install --editable .

Development mode is the recommended approach if you plan to implement new features.
In this case, you may prefer to fork AbiPy on GitHub first and then clone your fork,
so that you can push changes to your fork and later get them merged into the main branch.

The documentation of the **developmental** version is hosted on `GitHub Pages <http://abinit.github.io/abipy>`_.

The GitHub version includes the test files needed for complete unit testing.
To run the test suite, make sure pytest_ is installed and run::

    pytest

in the AbiPy root directory.

Several unit tests check the integration between AbiPy and Abinit.
To run them, you need a working set of Abinit executables and
a ``manager.yml`` configuration file.
For the syntax of the configuration file, see the :ref:`taskmanager` section.

A pre-compiled sequential version of Abinit for Linux and macOS can be installed from the abinit-channel_ with::

    conda install abinit -c conda-forge

Examples of configuration files for compiling Abinit on clusters are available
in the abiconfig_ package.

Contributing to AbiPy is easy: just send us a `pull request <https://help.github.com/articles/using-pull-requests/>`_.
When you open the request, choose ``develop`` as the destination branch.
AbiPy uses the `Git Flow <http://nvie.com/posts/a-successful-git-branching-model/>`_ branching model:
the ``develop`` branch contains the latest contributions, while ``master`` is always tagged and points
to the latest stable release.

If you share your developments, please take the time to write unit tests that cover at least
the basic functionality of your code.

.. _installing_without_internet_access:

----------------------------------
Installing without internet access
----------------------------------

This section explains how to set up a virtual environment with AbiPy on a cluster that cannot access the internet.
The idea is to create a virtual environment with AbiPy on a machine with internet access, copy the required files
to the offline cluster, and then perform an offline installation there.
We use Conda for the Python installation, since it reduces the risk of incompatibilities,
and pip for the packages, since it offers a convenient syntax for offline installation.

First, you need Conda on the machine with internet access.
If it is not available by default, follow the :ref:`instructions for installing Conda <anaconda_howto>`.
Then create a conda virtual environment with a given Python version, for example 3.12::

    conda create --name abienv python=3.12
    conda activate abienv

Install AbiPy in this environment, write the list of installed packages to ``requirements.txt``,
and download all the wheels (``.whl`` files) into a ``packages/`` folder::

    pip install abipy
    pip list --format=freeze > requirements.txt
    pip download -r requirements.txt -d packages/

Next, copy ``requirements.txt``, the ``packages/`` folder and the Miniconda installer to the offline cluster.
You may need to use ``scp`` from a computer that can reach both machines.
The Miniconda installer is needed only if Conda is not already available on the offline cluster.
From a computer that can access both locations, run::

    scp -r connected_cluster:/file/and/folder/location/* .
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    scp -r requirements.txt packages/ Miniconda3-latest-Linux-x86_64.sh disconnected_cluster:/desired/location/

If Conda is not available on the offline cluster,
follow the :ref:`Conda installation instructions <anaconda_howto>` there as well.
Then create an **offline** virtual environment on the offline cluster::

    conda create --name abienv --offline python=3.12
    conda activate abienv

and install the packages from the local folder with::

    pip install --no-index --find-links=packages/ -r requirements.txt

At this step, the installation may fail because of missing or incompatible packages.
Some of these issues can be solved by repeating the steps above (except for the environment creation)
for the packages reported as missing or incompatible, updating ``requirements.txt`` and ``packages/``,
and trying again.
Once you see::

        Successfully installed abipy-x.y.z

you can quickly test the installation by running ``python`` followed by ``import abipy``.


.. _howto_compile_python_and_bootstrap_pip:

---------------
Troubleshooting
---------------

^^^^^^^^^^^^^^^^^^^^^
unknown locale: UTF-8
^^^^^^^^^^^^^^^^^^^^^

If Python stops with the error message::

    "ValueError: unknown locale: UTF-8"

add the following line to the ``.bashrc`` file in your ``$HOME`` (``.profile`` on macOS)::

    export LC_ALL=C

then reload the environment with ``source ~/.bashrc`` and run the code again.

^^^^^^^^^^^^^^^^^^^^
netcdf does not work
^^^^^^^^^^^^^^^^^^^^

The version of hdf5 installed by conda may not be compatible with python-netcdf.
Try the hdf5/netcdf4 libraries provided by conda-forge::

    conda uninstall hdf4 hdf5
    conda config --add channels conda-forge
    conda install netcdf4

These packages are known to work on macOS::

    conda list hdf4
    hdf4                      4.2.12                        0    conda-forge
    conda list hdf5
    hdf5                      1.8.17                        9    conda-forge
    conda list netcdf4
    netcdf4                   1.2.7               np112py36_0    conda-forge

^^^^^^^^^^^^^^^^^^^
UnicodeDecodeError
^^^^^^^^^^^^^^^^^^^

If Python 2.7 raises `UnicodeDecodeError: 'ascii' codec can't decode byte ...`
when opening files with abiopen, add

.. code-block:: python

    import sys
    reload(sys)
    sys.setdefaultencoding("utf8")

at the beginning of your script.

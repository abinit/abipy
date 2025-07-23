=============
Getting AbiPy
=============

.. contents::
   :backlinks: top

--------------
Stable version
--------------

The version at the `Python Package Index <https://pypi.org/project/abipy/>`_ (PyPI) is always
the latest **stable** release that can be installed in user mode with::

    pip install abipy --user

Note that you may need to install some optional dependencies manually.
In this case, please consult the detailed installation instructions provided in the
`pymatgen howto <https://pymatgen.org/installation.html>`_ to install these optional packages.

The installation process is greatly simplified if you install the required
packages through the Anaconda_ distribution.
We routinely use conda_ to test new developments with multiple Python versions and multiple virtual environments.
The anaconda distribution already provides the most critical dependencies (numpy_, scipy_, matplotlib_, netcdf4-python_)
in the form of pre-compiled packages that can be easily installed with e.g.::

    conda install numpy scipy netcdf4

Additional information on the steps required to install AbiPy with anaconda
are available in the :ref:`anaconda_howto` howto as well as in the
`conda-based-install <http://pymatgen.org/installation.html#conda-based-install>`_
section of the pymatgen_ documentation.

We are also working with the spack_ community
to provide packages for AbiPy and Abinit in order to facilitate the installation on large supercomputing centers.

Advanced users who need to compile a local version of the python interpreter and install the AbiPy dependencies
manually can consult the :ref:`howto_compile_python_and_bootstrap_pip` document.

---------------------
Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

ipython_

    Required to interact with the AbiPy/Pymatgen objects in the ipython shell
    (strongly recommended, already provided by conda_).

jupyter_ and nbformat_

    Required to generate jupyter notebooks.
    Install these two packages with ``conda install jupyter nbformat`` or use pip_.
    To use ``jupyter`` you will also need a web browser to open the notebook.
    (recommended)

.. _anaconda_howto:

--------------
Anaconda Howto
--------------

Download the anaconda installer from the `official web-site <https://www.continuum.io/downloads>`_.
by choosing the version that matches your OS.
You may want to use the ``wget`` utility to download the anaconda script directly from the terminal
(useful if you are installing anaconda on a cluster).

Run the bash script in the terminal and follow the instructions on screen.
By default, the installer creates the ``anaconda`` directory in your home.
Anaconda will add one line to your ``.bashrc`` to enable access to the anaconda executables.
Once the installation is completed, execute::

    source ~/anaconda/bin/activate base

to activate the ``base`` environment.
The output of ``which python`` should show that you are using the python interpreter provided by anaconda.

Use the conda_ command-line interface to install the packages not included in the official distribution.
For example, you can install ``pyyaml`` and ``netcdf4`` with::

    conda install pyyaml netcdf4

Remember that if a package is not available in the official conda repository, you can always
download the package from one of the conda channels or use ``pip install`` if no conda package is available.

Fortunately there are conda channels providing all dependencies needed by AbiPy.
Now add ``conda-forge`` to your conda channels with::

    conda config --add channels conda-forge

This is the channel from which we will download pymatgen, abipy and abinit.

Finally, install AbiPy with::

    conda install abipy

Once you have completed the installation of AbiPy and pymatgen, open the ipython_ shell and type::

    # make sure spglib library works
    import spglib

    # make sure pymatgen is installed
    import pymatgen

    from abipy import abilab

to check the installation.

Note that one can use conda_ to create different environments with different
versions of the python interpreter or different libraries.
Further information are available on the `conda official website <http://conda.pydata.org/docs/test-drive.html>`_.
Using different environments is very useful to keep different versions and branches separate.

.. _developmental_version:

---------------------
Developmental version
---------------------

Getting the developmental version of AbiPy is easy.
You can clone it from our  `github repository <https://github.com/abinit/abipy>`_ using::

    git clone https://github.com/abinit/abipy

After cloning the repository, type::

    python setup.py install

or alternately::

    python setup.py develop

to install the package in developmental mode
(Develop mode is the recommended approach if you are planning to implement new features.
In this case you may also opt to first fork AbiPy on Git and then clone your own fork.
This will allow you to push any changes to you own fork and also get them merged in the main branch).

The documentation of the **developmental** version is hosted on `github pages <http://abinit.github.io/abipy>`_.

The Github version include test files for complete unit testing.
To run the suite of unit tests, make sure you have pytest_ installed and issue::

    pytest

in the AbiPy root directory.

Note that several unit tests check the integration between AbiPy and Abinit.
In order to run the tests, you need a working set of Abinit executables and
a ``manager.yml`` configuration file.
For further information on the syntax of the configuration file, please consult the :ref:`taskmanager` section.

A pre-compiled sequential version of Abinit for Linux and OSx can be installed directly from the abinit-channel_ with::

    conda install abinit -c conda-forge

Examples of configuration files to configure and compile Abinit on clusters can be found
in the abiconfig_ package.

Contributing to AbiPy is relatively easy.
Just send us a `pull request <https://help.github.com/articles/using-pull-requests/>`_.
When you send your request, make ``develop`` the destination branch on the repository
AbiPy uses the `Git Flow <http://nvie.com/posts/a-successful-git-branching-model/>`_ branching model.
The ``develop`` branch contains the latest contributions, and ``master`` is always tagged and points
to the latest stable release.

If you choose to share your developments please take some time to develop some unit tests of at least the
basic functionalities of your code

.. _installing_without_internet_access:

----------------------------------
Installing without internet access
----------------------------------

Here, it is described how to set up a virtual environment with AbiPy on a cluster that cannot reach out to the internet.
One first creates a virtual environment with AbiPy on a cluster/computer with access, then ports the required files
to the cluster without access, and performs an offline installation.
We use Conda for the Python installation and pip for the packages, as the former reduces the odds that incompatibilities arise,
while the latter provides convenient syntax for offline package installation.

One first needs Conda on the cluster with access.
If not available by default, follow the :ref:`instructions for installing Conda <anaconda_howto>`.
Next, set up a conda virtual environment with a designated Python version, for example 3.12::

    conda create --name abienv python=3.12
    conda activate abienv

We then install AbiPy in this virtual environment, followed by creating requirements.txt, and creating
a folder packages/ containing all the wheels (.whl format)::

    pip install abipy
    pip list --format=freeze > requirements.txt
    pip download -r requirements.txt -d packages/

Next, the .txt file, the folder, and the miniconda installer must be forwarded to the cluster without internet access.
You may have to use a computer that has access to both locations with the scp command.
If the offline cluster does not have Conda preinstalled, the Miniconda executable must be ported so that
an offline Conda installation can be performed.
Thus, from a computer that can access both locations, execute::

    scp -r connected_cluster:/file/and/folder/location/* .
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    scp -r requirements.txt packages/ Miniconda3-latest-Linux-x86_64.sh disconnected_cluster:/desired/location/

If conda is not available on the cluster that cannot access the internet, 
follow the :ref:`Conda installation instructions <anaconda_howto>` once more.
Next, one can set up an **offline** virtual environment on the cluster without internet access::

    conda create --name abienv --offline python=3.12
    conda activate abienv

At this step, AbiPy might fail to install due to missing/incompatible packages.
Some of these issues may be solved by repeating the above steps (excluding the environment creation) 
for packages that are listed as missing/incompatible during the installation procedure, by updating the requirements.txt and packages/ and trying to install again.
Upon reading::

        Successfully installed abipy-x.y.z

You can quickly test your installation by running ``python`` followed by ``import abipy``.


.. _howto_compile_python_and_bootstrap_pip:

---------------
Troubleshooting
---------------

^^^^^^^^^^^^^^^^^^^^^
unknown locale: UTF-8
^^^^^^^^^^^^^^^^^^^^^

If python stops with the error message::

    "ValueError: unknown locale: UTF-8"

add the following line to your ``.bashrc`` file inside your ``$HOME`` (``.profile`` if MacOSx)::

    export LC_ALL=C

reload the environment with ``source ~/.bashrc`` and rerun the code.

^^^^^^^^^^^^^^^^^^^^
netcdf does not work
^^^^^^^^^^^^^^^^^^^^

The version of hdf5 installed by conda may not be compatible with python netcdf.
Try the hdf5/netcdf4 libraries provided by conda forge::

    conda uninstall hdf4 hdf5
    conda config --add channels conda-forge
    conda install netcdf4

These packages are known to work on MacOsX::

    conda list hdf4
    hdf4                      4.2.12                        0    conda-forge
    conda list hdf5
    hdf5                      1.8.17                        9    conda-forge
    conda list netcdf4
    netcdf4                   1.2.7               np112py36_0    conda-forge

^^^^^^^^^^^^^^^^^^^
UnicodeDecodeError
^^^^^^^^^^^^^^^^^^^

Python2.7 raises an `UnicodeDecodeError: 'ascii' codec can't decode byte ...`
when trying to open files with abiopen. Add

.. code-block:: python

    import sys
    reload(sys)
    sys.setdefaultencoding("utf8")

at the beginning of your script.

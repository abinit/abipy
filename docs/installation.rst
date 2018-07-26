=============
Getting AbiPy
=============

.. contents::
   :backlinks: top

--------------
Stable version
--------------

The version at the `Python Package Index <https://pypi.python.org/pypi/abipy>`_ (PyPI) is always 
the latest **stable** release that can be installed with::

    pip install abipy

Note that you may need to install pymatgen_ and other critical dependencies manually.
In this case, please consult the detailed installation instructions provided in the
`pymatgen howto <http://pymatgen.org/index.html#standard-install>`_ to install pymatgen 
and then follow the instructions in the :ref:`netcdf4_installation` section.

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
Choose the version that matches your OS and select python3.6.
You may want to use the ``wget`` utility to download the anaconda script directly from the terminal
(useful if you are installing anaconda on a cluster).

Run the bash script in the terminal and follow the instructions.
By default, the installer creates the ``anaconda`` directory in your home.
Anaconda will add one line to your ``.bashrc`` to enable access to the anaconda executables.
Once the installation is completed, execute::

    source ~/anaconda/bin/activate root

to activate the ``root`` environment.
The output of ``which python`` should show that you are using the python interpreter provided by anaconda.

Use the conda_ command-line interface to install the packages not included in the official distribution.
For example, you can install ``pyyaml`` and ``netcdf4`` with::

    conda install pyyaml netcdf4

Remember that if a package is not available in the official conda repository, you can always
download the package from one of the conda channels or use ``pip install`` if no conda package is available.

Fortunately there are conda channels providing all dependencies needed by AbiPy.
To install the pymatgen_ package from the matsci_ channel, use::

    conda install pymatgen --channel matsci

then install Abipy from the abinit-channel_ with::

    conda install abipy --channel abinit

Visit `materials.sh <http://materials.sh>`_ for instructions on how to use the
matsci channel to install pymatgen and other packages.

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

    conda install abinit -c abinit

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

.. _howto_compile_python_and_bootstrap_pip:

-------------------------------------
How to compile the Python interpreter
-------------------------------------

This section discusses how to install a local version of the python interpreter as well
as the most important dependencies needed by AbiPy.
This approach may be needed if you want to use AbiPy on a machine (e.g. a cluster)
in which you don't have root privileges and the version of the python interpreter is too old 
or if for some reasons you prefer not to use ``anaconda``.
In this case you cannot use a `virtual environment <https://virtualenv.pypa.io/en/latest/>`_ 
on top of the preexisting python library.

First of all, you have to create a new directory containing your python interpreter
as well as as the libraries and the other executables needed by AbiPy.
Let's assume we decided to create this directory inside ``$HOME`` and let's call it ``local``::

    mkdir $HOME/local

Now change your ``~/.bashrc`` file and add the following three lines::

    export PATH=$HOME/local/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
    export C_INCLUDE_PATH=$HOME/include/:$C_INCLUDE_PATH

so that other scripts and tools will know where to find the new binaries/libraries/include files they need.

Get the python tarball from the `python official site <https://www.python.org>`_ and unpack it.
Configure the package with the ``--prefix`` option and compile the code
(use the ``-j`` option to speedup the compilation with threads)::

    ./configure --prefix=$HOME/local
    make -j4

If you plan to use graphical tools you need to make sure that the ``Tkinter`` graphical backends 
is installed and functional at the time of compilation of python, see below.

At the end, you should get the list of modules that could not be built because
your system does not provide the required libraries.
The installation should be OK for AbiPy if you get::

    Python build finished, but the necessary bits to build these modules were not found:
    _sqlite3           bsddb185           dl              
    imageop            sunaudiodev                        
    To find the necessary bits, look in setup.py in detect_modules() for the module's name.

If, on the other hand, python has been built without ``bz2`` or ``_tkinter`` you are in trouble 
because these packages are required.

``bz2`` is more fundamental than ``_tkinter`` because it is used to compress/uncompress files.
AbiPy/Pymatgen won't work without ``bz2`` and you have to install the ``bzip`` library with the C headers.
The source code is available from `bzip.org <www.bzip.org>`_
See also `this post <http://stackoverflow.com/questions/12806122/missing-python-bz2-module>`_ on stackoverflow.

``Tkinter`` is less important than ``bz2`` but without it you won't be able to use the ``matplotlib`` graphical back-end.
If you want ``matplotlib`` with the Tk back-end, you have to install Tk/Tcl. 
Get the tarball from the `tcl.tk site <www.tcl.tk/software/tcltk/downloads.html>`_, configure 
with ``--prefix`` and ``make && make install`` as usual.
Then reconfigure python. 

Once you have solved the problem with the missing modules, you can run the tests with::

    make test 

and install the python interpreter with::

    make install

Now we have our python interpreter installed in ``$HOME/local``::

    which python 
    $HOME/local/bin/python

but we still need to install ``easy_install`` and ``pip`` so that we can automatically 
download and install other python packages.

To install ``easy_install``::

    wget https://bootstrap.pypa.io/ez_setup.py -O - | python

    which easy_install
    $HOME/local/bin/easy_install

For more info, consult the `setuptools page <https://pypi.python.org/pypi/setuptools>`_

Now use ``easy_install`` to install ``pip``::

    easy_install pip

    # Upgrade setuptools with
    pip install setuptools --upgrade

Henceforth we can start to use ``pip`` to install the python modules.
Start with ``cython`` and ``numpy``::

    pip install cython 
    pip install numpy

The installation of ``scipy`` is more complicated due to the need for the BLAS and LAPACK libraries.
Try first::

    pip install scipy

If the installer does not find ``BLAS/LAPACK`` in your system, consult the
`scipy documentation <http://www.scipy.org/scipylib/building/linux.html#id1>`_.

.. _netcdf4_installation:

---------------------------------------------------
How to install HDF5/Netcdf4 and the python bindings
---------------------------------------------------

Obtain the latest ``HDF5`` software from the `official hd5 web-site <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_.
Configure the package with ``--enable-hl --enable-shared`` and the ``--prefix`` option as usual.
Build and install with::

    make
    make install

Finally define the environment variable ``$HDF5_DIR`` with::

    export HDF5_DIR=$HOME/local

Get the latest stable netCDF-C release from `this page <http://www.unidata.ucar.edu/downloads/netcdf/index.jsp>`_.
Configure with::

    configure --prefix=$HOME/local --enable-netcdf-4 --enable-shared \
      CPPFLAGS="-I$HDF5_DIR/include" LDFLAGS="-L$HDF5_DIR/lib"

Build and install with ``make && make install``
Define the environment variable ``$NETCDF4_DIR``::

    export NETCDF4_DIR=$HOME/local

Now we can download and install the python interface with::

    pip install netcdf4

You may want to consult the official `netcdf4-python documentation <http://unidata.github.io/netcdf4-python>`_.

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

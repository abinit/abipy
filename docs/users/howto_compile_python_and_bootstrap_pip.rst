.. _howto_compile_python_and_bootstrap_pip:

*********************************************************
How to compile the Python interpreter and bootstrap `pip` 
*********************************************************

This document discusses how to install a local version of the python interpreter as well
as the most important dependencies needed by AbiPy.

This approach may be needed if you want to use AbiPy on a machine (e.g. a cluster)
in which you don't have root privileges and the version of the python interpreter is too old.
In this case you cannot use a `virtual environment <https://virtualenv.pypa.io/en/latest/>`_ 
on top of the preexisting python library.

How to compile the Python interpreter
=====================================

First of all, you have to create a new directory to contain your python interpreter
as well as as the libraries and the other executables needed by AbiPy.
Let's assume we decided to create this directory inside `$HOME` and let's call it `local`::

    mkdir $HOME/local

Now change you `~/.bashrc` file by adding the following lines::

    export PATH=$HOME/local/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
    export C_INCLUDE_PATH=$HOME/include/:$C_INCLUDE_PATH

so that other scripts and tools will know where to find the new binaries/libraries/include files they need.

Get the python tarball from the `official site <https://www.python.org>`_ and unpack it.
Configure the package with the `prefix` option and compile the code
(use the `-j` option to speedup the compilation with threads)::

    $ ./configure --prefix=$HOME/local
    $ make -j4

At the end, you should get the list of modules that could not be built because 
your system does not provide the required libraries.
The installation should be OK for AbiPy if you get::

    Python build finished, but the necessary bits to build these modules were not found:
    _sqlite3           bsddb185           dl              
    imageop            sunaudiodev                        
    To find the necessary bits, look in setup.py in detect_modules() for the module's name.

If, on the other hand, python has been built without `bz2` or `_tkinter` you are in trouble 
because these packages are needed.

`bz2` is more fundamental than `_tkinter` because it is used to compress/uncompress files.
AbiPy won't work without `bz2` and you have to install the `bzip` library with the C headers.
The source code is available `from bzip.org <www.bzip.org>`_
See also `this post <http://stackoverflow.com/questions/12806122/missing-python-bz2-module>`_ on stackoverflow.

`Tkinter` is less important than `bz2` but without it you won't be able to use the `matplotlib` graphical back-end.
If you want `matplotlib` with the Tk back-end, you have to install Tk/Tcl. 
Get the tarball from the `tcl.tk site <www.tcl.tk/software/tcltk/downloads.html>`_, configure with `--prefix` and 
`make && make install` as usual.
Then reconfigure python. 

Once you have solved the problem with the missing modules, you can run the tests with::

    $ make test 

and install the `python` interpreter with::

    $ make install

Now we have our python interpreter installed in `$HOME/local`::

    $ which python 
    $HOME/local/bin/python

but we still need to install `easy_install` and `pip` so that we can automatically 
download and install other python packages.

To install `easy_install`::

    $ wget https://bootstrap.pypa.io/ez_setup.py -O - | python

    $ which easy_install
    $HOME/local/bin/easy_install

See also https://pypi.python.org/pypi/setuptools

Now use `easy_install` to install `pip`::

    $ easy_install pip

    # Upgrade setuptools with
    $ pip install setuptools --upgrade

Henceforth we can start to use `pip` to install the python modules.
Start with `cython` and `numpy`::

    $ pip install cython 
    $ pip install numpy

The installation of `scipy` is more complicated due to the need for the BLAS and LAPACK libraries.
Try first::

    $ pip install scipy

If the installer does not find BLAS/LAPACK in your system, consult the
`scipy documentation <http://www.scipy.org/scipylib/building/linux.html#id1>`_.


How to install HDF5/Netcdf4 and the python bindings
===================================================

Obtain the latest HDF5 Software from the `official web-site <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_.
Configure the package with `--enable-hl --enable-shared` and the `--prefix` option as usual.
Build and install with:: 

    make 
    make install

Finally define the environment variable `$HDF5_DIR` with::

    export HDF5_DIR=$HOME/local

Get the latest stable netCDF-C release from `here <http://www.unidata.ucar.edu/downloads/netcdf/index.jsp>`_.
Configure with::

    configure --prefix=$HOME/local --enable-netcdf-4 --enable-shared \
      CPPFLAGS="-I$HDF5_DIR/include" LDFLAGS="-L$HDF5_DIR/lib"

Build and install with `make && make install`
Define the environment variable `$NETCDF4_DIR`::

    export HDF5_DIR=$HOME/local

Now we can download and install the python interface with::

    pip install netcdf4

You may want to consult the official `netcdf4-python documentation <http://unidata.github.io/netcdf4-python>`_.

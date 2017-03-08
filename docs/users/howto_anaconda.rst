.. _howto_anaconda:

*************************************************************************
How to use Anaconda to install Python and the most important dependencies
*************************************************************************

Download the anaconda installer from the `official web-site <https://www.continuum.io/downloads>`_.
Choose the version that matches your OS and select python2.7.
You may want to use the ``wget`` utility to download the anaconda script directly from the terminal
(useful if you are installing anaconda on a cluster).

Run the bash script in the terminal and follow the instructions.
By default, the installer creates the ``anaconda`` directory in your home.
Once the installation is completed, execute::

    $ source ~/anaconda/bin/activate root

to activate the ``root`` environment.
The output of ``which python`` should show that you are using the python interpreter provided by anaconda.

Use the ``conda`` command-line interface to install the packages not included in the official distribution.
For example, you can install ``pyyaml`` and ``netcdf4`` by simply typing::

    $ conda install pyyaml, netcdf4

Remember that if a package is not available in the official conda repository, you can always
use ``pip install`` or download the package from one of the conda channels.
For example, if you encounter problems while installing the spacegroup library
with ``pip install pyspglib``, you can install the pre-compiled library from the ``jochym`` channel with::

    $ conda install -c jochym pyspglib

Now you can install ``pymatgen`` and ``abipy``.
Use::

    $ pip install pymatgen
    $ pip install abipy

for the stable version.
If you want to use the developmental version, clone (or fork) the repositories on github
and install the packages with::

    $ python setup.py install

Once you have completed the installation of abipy and pymatgen, open the ``ipython`` shell and type::

    from abipy import abilab

to check the installation.

Optionally, you may want to execute::

    $ conda install wxpython

to install the ``wxpython`` graphical toolkit required for the GUIs.

Note that one can use ``conda`` to create different enviroments with different
versions of the python interpreter or different libraries.
Further information are available on the
`official website <http://conda.pydata.org/docs/test-drive.html>`_.

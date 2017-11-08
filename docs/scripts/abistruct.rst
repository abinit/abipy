.. _abistruct:

^^^^^^^^^^^^^^^^
``abistruct.py``
^^^^^^^^^^^^^^^^

This script is useful if you want to analyze/export/visualize the crystalline structure 
stored in one the ``netcdf`` files produced by ABINIT (other formats are supported e.g. 
``cif`` files or ``POSCAR``, see also docs below).
Also in this case, it is possible to analyze the structure object either inside ``ipython`` or
``jupyter`` notebooks.

The syntax of the command is::

    $ abistruct.py command FILE [options]

The documentation for a given ``command`` is accessible with::

    $ abistruct.py command --help 

Use e.g.:: 

    $ abistruct.py spglib --help

to get the options supported by the ``spglib`` command.

The ``convert`` command is quite useful if you need to convert the crystalline structure
from one format to another one.
For example, one can read a ``cif`` file and print the corresponding Abinit variables with::

    $ abistruct.py convert CIF

To get the list of options, use

.. command-output:: abistruct.py --help

.. argparse::
   :ref: abipy.scripts.abistruct.get_parser
   :prog: abistruct.py

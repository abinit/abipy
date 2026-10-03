.. _abistruct-script:

^^^^^^^^^^^^^^^^
``abistruct.py``
^^^^^^^^^^^^^^^^

This script reads a |Structure| object from a file and performs predefined operations
depending on the ``COMMAND`` and ``options`` given on the command line.
The syntax is::

    abistruct.py COMMAND FILE [options]

where ``FILE`` is any file from which AbiPy can extract a Structure object. This includes
most of the netcdf output files, Abinit input and output files in text format,
and the other formats supported by pymatgen_, e.g. CIF_ files, POSCAR_, etc.

The documentation of a given ``COMMAND`` is available with::

    abistruct.py COMMAND --help

For example::

    $ abistruct.py spglib --help

lists the options supported by the spglib_ command.

The ``convert`` command is useful for converting the crystalline structure
from one format to another.
For example, to read a CIF_ file and print the corresponding Abinit variables, use::

    $ abistruct.py convert si.cif

.. note::

    The script can fetch data from the Materials Project database and
    the COD_ database.
    To access the Materials Project database, register on
    https://www.materialsproject.org to get your personal access token.
    Then create a `.pmgrc.yaml` configuration file in your $HOME and add your token with the line::

        PMG_MAPI_KEY: your_token_goes_here

You can analyze the structure object either in a jupyter_ notebook, e.g.::

    abistruct.py notebook si.cif

or directly in the ipython_ shell with::

    abistruct.py ipython si.cif

Several other commands are available. To get the full list, use:

.. command-output:: abistruct.py --help

Complete command line reference:

.. argparse::
   :ref: abipy.scripts.abistruct.get_parser
   :prog: abistruct.py

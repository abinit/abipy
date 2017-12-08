.. _abistruct:

^^^^^^^^^^^^^^^^
``abistruct.py``
^^^^^^^^^^^^^^^^

This script reads a ``Structure`` object from file and performs predefined operations
depending on the ``command`` and the ``options`` specified on the command line.
The syntax of the command is::

    $ abistruct.py command FILE [options]

where ``FILE`` is any file from which AbiPy can extract a Structure object (this includes
the majority of the nectdf output files, Abinit input and output files in text format
as well as other formats supported by pymatgen e.g. ``cif`` files, ``POSCAR`` etc.

The documentation for a given ``command`` is accessible with::

    $ abistruct.py command --help 

Use e.g.:: 

    $ abistruct.py spglib --help

to get the list of options supported by the ``spglib`` command.

The ``convert`` command is quite useful if you need to convert the crystalline structure
from one format to another one.
For example, one can read a cif file and print the corresponding Abinit variables with::

    $ abistruct.py convert si.cif


.. NOTE::

    The script can fetch data from the materials project database and 
    the `COD database http://www.crystallography.net/cod`_
    To access the materials project database, please register on 
    <https://www.materialsproject.org> to get your personal access token.
    Then create a `.pymrc.yaml` configuration file inside your $HOME and add your token with the line::

        PMG_MAPI_KEY: your_token_goes_here

It is possible to analyze the structure object either a ``jupter`` notebook with e.g.::

    $ abistruct.py jupter si.cif

or directly inside the ``ipython`` shell with::

    $ abistruct.py notebook si.cif

Several other options are available. To get the full list, use

.. command-output:: abistruct.py --help

Complete command line reference

.. argparse::
   :ref: abipy.scripts.abistruct.get_parser
   :prog: abistruct.py

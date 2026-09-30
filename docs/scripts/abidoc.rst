.. _abidoc.py:

^^^^^^^^^^^^^
``abidoc.py``
^^^^^^^^^^^^^

This script provides a command line interface to the Abinit documentation.
For example, to get the documentation of the ``ecut`` input variable, use::

    abidoc.py man ecut

For the full list of commands, use:

.. command-output:: abidoc.py --help

Complete command line reference:

.. argparse::
   :ref: abipy.scripts.abidoc.get_parser
   :prog: abidoc.py

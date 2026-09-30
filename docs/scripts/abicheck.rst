.. _abicheck.py:

^^^^^^^^^^^^^^^
``abicheck.py``
^^^^^^^^^^^^^^^

This script checks that the options in ``manager.yml`` and ``scheduler.yml``,
as well as the environment on the local machine, are properly configured.
See the :ref:`taskmanager` documentation for a detailed description of these YAML_ files.

.. command-output:: abicheck.py --no-colors

Use ``abicheck.py --with-flow`` to run a small AbiPy flow and
check the interface with the Abinit executables.

Complete command line reference:

.. argparse::
   :ref: abipy.scripts.abicheck.get_parser
   :prog: abicheck.py

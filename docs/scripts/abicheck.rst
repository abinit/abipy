.. _abicheck.py:

^^^^^^^^^^^^^^^
``abicheck.py``
^^^^^^^^^^^^^^^

This script checks that the options specified in ``manager.yml``, ``scheduler.yml``,
and the environment on the local machine are properly configured.
Please consult the documentation on the :ref:`taskmanager` for a more detailed description of these YAML_ files.

.. command-output:: abicheck.py --no-colors

The command ``abicheck.py --with-flow`` can be used to run a small AbiPy flow in order to
check the interface with the Abinit executables.

Complete command line reference

.. argparse::
   :ref: abipy.scripts.abicheck.get_parser
   :prog: abicheck.py

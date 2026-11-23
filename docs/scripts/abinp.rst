.. _abinp-script:

^^^^^^^^^^^^
``abinp.py``
^^^^^^^^^^^^

This script provides a simplified interface to the AbiPy API for building input files.
It is especially useful for newcomers who are not yet familiar with the programmatic interface
for building workflows:
``abinp.py`` can automatically generate input files from any
file that provides the crystalline structure of the system, and the generated output can then be customized.

Other commands operate directly on existing input files
and can be used to get data directly from Abinit.

.. argparse::
   :ref: abipy.scripts.abinp.get_parser
   :prog: abinp.py

.. _abinp.py:

^^^^^^^^^^^^
``abinp.py``
^^^^^^^^^^^^

This script provides a simplified interface to the AbiPy API for constructing input files.
It is a useful tool especially for newcomers who are not familiar with the programmatic interface
for building workflows.
In this case, indeed, one can use ``abinp.py`` to generate automatically input files from 
file providing the crystalline structure of the system and then customize the generated output.

There are also commands operating on input files directly. 
These options could be useful to get data directly from Abinit.

.. argparse::
   :ref: abipy.scripts.abinp.get_parser
   :prog: abinp.py

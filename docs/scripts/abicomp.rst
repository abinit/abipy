.. _abicomp.py:

^^^^^^^^^^^^^^
``abicomp.py``
^^^^^^^^^^^^^^

This script compares results stored in multiple netcdf_ files.
For example, one can compare the crystalline structure used in different calculations
or compare the electronic bands stored in two or more netcdf_ files (e.g. GSR.nc_ or ``WFK.nc``).
By default the script displays the results/plots directly within the shell.
Use the command::

    abicomp.py struct out1_GSR.nc out2_GSR.nc

to compare the structures reported in two ``GSR.nc`` files and print the result to screen.
Use ``--ipython`` to start an ipython_ terminal or ``-nb`` to generate a jupyter_ notebook, e.g.::

    abicomp.py ebands out1_GSR.nc out2_GSR.nc -nb

.. command-output:: abicomp.py --help

Complete command line reference

.. argparse::
   :ref: abipy.scripts.abicomp.get_parser
   :prog: abicomp.py

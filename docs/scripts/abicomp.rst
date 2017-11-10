.. _abicomp:

^^^^^^^^^^^^^^
``abicomp.py``
^^^^^^^^^^^^^^

This script is used to analyze/compare results stored in multiple ``netcdf`` files.
For example, one can compare the crystalline structure used in different calculations
or compare the electronic bands stored in two or more ``netcdf`` files (e.g. ``GSR.nc`` or ``WFK.nc``).
By default the script displays the results/plots directly within the shell.
Use the command::

    $ abicomp.py struct out1_GSR.nc out2_GSR.nc

to compare the structures reported in two ``GSR`` files and print the result to screen.
Use ``--ipython`` to start an ``ipython`` terminal or ``-nb`` to generate a ``jupyter`` notebook, e.g.::

    $ abicomp.py ebands out1_GSR.nc out2_GSR.nc -nb

.. command-output:: abicomp.py --help

.. argparse::
   :ref: abipy.scripts.abicomp.get_parser
   :prog: abicomp.py

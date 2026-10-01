.. _abicomp-script:

^^^^^^^^^^^^^^
``abicomp.py``
^^^^^^^^^^^^^^

This script compares the results stored in multiple netcdf_ files.
For example, you can compare the crystalline structures used in different calculations
or the electronic bands stored in two or more netcdf_ files (e.g. GSR.nc_ or ``WFK.nc``).

Depending on COMMAND, ``abicomp`` either starts an ``ipython`` session in which you can interact
with the ``robot``, or prints the results to the screen.

For instance, the command::

    abicomp.py structure out1_GSR.nc out2_GSR.nc

compares the crystalline structures stored in two ``GSR.nc`` files and prints the result to the screen, while::

    abicomp.py gsr out*_GSR.nc

starts an ipython session.

Use the ``-p`` option if you just want to print information about the files without opening an ipython session, e.g.::

    abicomp.py gsr out1_GSR.nc out2_GSR.nc -p

The ``-nb`` option automatically generates a jupyter_ notebook, e.g.::

    abicomp.py gsr out1_GSR.nc out2_GSR.nc -nb

Finally, use ``-e`` (``--expose``) to generate matplotlib plots automatically::

    abicomp.py gsr out1_GSR.nc out2_GSR.nc -e -sns=poster

The seaborn_ plot style and settings can be changed from the command line with the `-sns` option.

.. command-output:: abicomp.py --help

Complete command line reference:

.. argparse::
   :ref: abipy.scripts.abicomp.get_parser
   :prog: abicomp.py

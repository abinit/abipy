.. _abicomp.py:

^^^^^^^^^^^^^^
``abicomp.py``
^^^^^^^^^^^^^^

This script compares output results stored in multiple netcdf_ files.
For example, one can compare the crystalline structure used in different calculations
or compare the electronic bands stored in two or more netcdf_ files (e.g. GSR.nc_ or ``WFK.nc``).

Depending on COMMAND, ``abicomp`` either starts an ``ipython`` session so that the user can interact 
with the ``robot`` or print the results to screen.

For instance, the command::

    abicomp.py structure out1_GSR.nc out2_GSR.nc

compares the crystalline structures reported in two ``GSR.nc`` files and print the result to screen while::

    abicomp.py gsr out*_GSR.nc

starts an ipython session.

Use the ``-p`` option if you just want to get information on the file without opening it, e.g.::

    abicomp.py gsr out1_GSR.nc out2_GSR.nc -p

It is possible to generate automatically a jupyter_ notebook with the ``-nb`` option e.g.::

    abicomp.py gsr out1_GSR.nc out2_GSR.nc -nb

Finally, use ``-e`` (``--expose``) to generate matplotlib plots automatically::

    abicomp.py gsr out1_GSR.nc out2_GSR.nc -e -sns=poster

seaborn_ plot style and settings can be changed from the command line interface with the `-sns` option 

.. command-output:: abicomp.py --help

Complete command line reference

.. argparse::
   :ref: abipy.scripts.abicomp.get_parser
   :prog: abicomp.py

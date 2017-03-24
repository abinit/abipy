==================
Command line tools
==================

.. _abiopen:

^^^^^^^^^^^^^^
``abiopen.py``
^^^^^^^^^^^^^^

AbiPy provides python objects associated to several Abinit output files.
These objects implement methods to analyze and plot the results.
The examples in our :doc:`gallery </examples/index>` use this API to plot data with ``matplotlib``.

The ``abiopen.py`` script provides a handy interface to the AbiPy objects.
It can be used to open Abinit files directly in the ``ipython``` shell or in a ``jupyter`` notebook and interact with
the associated object (called ``abifile`` in the ``ipython`` terminal).
The syntax of the script is::

    $ abiopen.py FILE [options]

where ``FILE`` is one of the files supported by AbiPy (usually in ``netcdf`` format but other 
files are supported as well).
By default ``abiopen`` starts the ``ipython`` terminal.
Alternatively, it is possible to generate automatically a ``jupyter`` notebook with the ``-nb`` option e.g.::

    $ abiopen.py out_FATBANDS.nc -nb

which will produce a notebook to visualize the electronic fatbands produced with ``prtdos 3`` inside a web browser.

Use the ``-p`` option if you just want to get information on the file without opening it, e.g.::

    $ abiopen.py out_GSR.nc -p

``abiopen.py`` uses the file extension to decide what to do with the file and the type
of python object that should be instantiated.
The list of supported file extensions is obtained with:

.. command-output:: abiopen.py --help

.. WARNING::

    AbiPy uses the ``.abi`` extension for Abinit input files, ``.abo`` for output files and ``.log`` for log files.
    Try to follow this convention in your calculations to facilitate the integration with AbiPy.

.. _abistruct:

^^^^^^^^^^^^^^^^
``abistruct.py``
^^^^^^^^^^^^^^^^

This script is useful if you want to analyze/export/visualize the crystalline structure 
stored in one the ``netcdf`` files produced by ABINIT (other formats are supported e.g. 
``cif`` files or ``POSCAR``, see also docs below).
Also in this case, it is possible to analyze the structure object either inside ``ipython`` or
``jupyter`` notebooks.

The syntax of the command is::

    $ abistruct.py command FILE [options]

The documentation for a given ``command`` is accessible with::

    $ abistruct.py command --help 

Use e.g.:: 

    $ abistruct.py spglib --help

to get the options supported by the ``spglib`` command.

The ``convert`` command is quite useful if you need to convert the crystalline structure
from one format to another one.
For example, one can read a ``cif`` file and print the corresponding Abinit variables with::

    $ abistruct.py convert CIF abivars

To get the list of options, use

.. command-output:: abistruct.py --help


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

.. _abidoc:

^^^^^^^^^^^^^
``abidoc.py``
^^^^^^^^^^^^^

This script provides a command line interface to the Abinit documentation.
For example, the documentation for the ``ecut`` input variable can be obtained with::

    $ abidoc.py man ecut

For the full list of commands supported use:

.. command-output:: abidoc.py --help

.. _abicheck:

^^^^^^^^^^^^^^^
``abicheck.py``
^^^^^^^^^^^^^^^

This script checks that the options specified in ``manager.yml``, ``scheduler.yml``,
and the environment on the local machine are properly configured.
Please consult the documentation on the :ref:`taskmanager` for a more detailed description of these YAML files.

.. command-output:: abicheck.py --no-colors

.. _abirun:

^^^^^^^^^^^^^
``abirun.py``
^^^^^^^^^^^^^

This script allows the user to submit the calculations contained in the AbiPy Flow 
(for further detail, consult the :ref:`taskmanager` documentation).
It provides a command line interface as well as a graphical interface based on ``wxpython``.

.. command-output:: abirun.py --help

.. command-output:: abirun.py doc_scheduler

.. command-output:: abirun.py . doc_manager

At the time of writing (|today|), AbiPy supports the following resource managers:

    * ``bluegene``
    * ``moab``
    * ``pbspro``
    * ``sge``
    * ``shell``
    * ``slurm``
    * ``torque``

To obtain the list of options supported by a particular resource manager e.g. ``slurm``::

    $ abirun.py . doc_manager slurm


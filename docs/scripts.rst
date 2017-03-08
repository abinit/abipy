==================
Command line tools
==================

.. _abiopen:

^^^^^^^^^^^^^^
``abiopen.py``
^^^^^^^^^^^^^^

This script opens Abinit outputs file (usually in ``netcdf`` format but other files are supported as well). 
The syntax of the command is::

    $ abiopen.py command FILE [options]

The documentation for a given `command` is accessible with::

    $ abiopen.py command --help 

e.g. ``abiopen.py spglib --help``.

By default the script starts an interactive ``ipython`` session so that the use can interact with the file 
(called ``abifile`` in the ``ipython`` terminal) and have access to its methods.
Alternatively, it is possible to generate automatically a ``jupyter`` notebook with the ``-nb`` option
to execute code and visualize the results inside a web browser.
``abiopen.py`` employs the file extension to decide what to do with the file and the type
of python object that should be instantiated.

The list of supported file extensions is obtained with:

.. command-output:: abiopen.py --help

.. WARNING::

    AbiPy uses the ``.abi`` extension for Abinit input files, ``.abo`` for output files 
    and ``.log`` for log files.
    Try to follow this convention in your calculations to facilitate the integration with AbiPy.


.. _abistruct:

^^^^^^^^^^^^^^^^
``abistruct.py``
^^^^^^^^^^^^^^^^

This script is useful if you want to analyze/export/visualize the crystalline structure 
reported in one the ``netcdf`` files produced by ABINIT (other formats are supported, see docs belows).
Also in this case, it is possible to analyze the structure object either inside ``ipython`` or
``jupyter`` notebooks.

To get the list of options, use

.. command-output:: abistruct.py --help


.. _abicomp:

^^^^^^^^^^^^^^
``abicomp.py``
^^^^^^^^^^^^^^

This script is used to analyze/compare results stored in multiple ``netcdf`` files.
For example, one can compare the crystalline structure used in different calculations
or compare the electronic bands stored in two or more ``netcdf`` files. 
By default the script displays the results/plots directly within the shell.
Use ``--ipython`` to start an ``ipython`` terminal or ``-nb`` to generate a ``jupyter`` notebook.

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


^^^^^^^^^^^^^^^
``abicheck.py``
^^^^^^^^^^^^^^^

This script checks that the options in ``manager.yml``, ``scheduler.yml``,
and the environment on the local machine are properly configured.

.. command-output:: abicheck.py --no-colors

^^^^^^^^^^^^^
``abirun.py``
^^^^^^^^^^^^^

This script allows the user to submit the calculations contained in the `Flow`.
It provides a command line interface as well as a graphical interface based on wxpython.

.. command-output:: abirun.py --help

.. command-output:: abirun.py doc_scheduler

.. command-output:: abirun.py . doc_manager

To obtain the list of options supported by a particular resource manager use e.g. slurm, use::

    $ abirun.py . doc_manager slurm

At the time of writing (|date|), AbiPy supports the following resource managers:

    . bluegene
    . moab
    . pbspro
    . sge
    . shell
    . slurm
    . torque

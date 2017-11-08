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


.. argparse::
   :ref: abipy.scripts.abirun.get_parser
   :prog: abirun.py

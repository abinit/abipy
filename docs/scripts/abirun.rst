.. _abirun-script:

^^^^^^^^^^^^^
``abirun.py``
^^^^^^^^^^^^^

This script submits the calculations contained in an AbiPy Flow
(for more details, see the :ref:`taskmanager` documentation).

.. command-output:: abirun.py --help

.. command-output:: abirun.py doc_scheduler

.. command-output:: abirun.py . doc_manager

At the time of writing (|today|), AbiPy supports the following resource managers:

* ``shell``
* pbspro_
* slurm_
* IBM loadleveler_
* moab_
* sge_
* torque_

To get the list of options supported by a particular resource manager, e.g. ``slurm``, use::

    abirun.py . doc_manager slurm

Complete command line reference:

.. argparse::
   :ref: abipy.scripts.abirun.get_parser
   :prog: abirun.py

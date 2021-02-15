.. _abirun.py:

^^^^^^^^^^^^^
``abirun.py``
^^^^^^^^^^^^^

This script allows the user to submit the calculations contained in the AbiPy Flow 
(for further detail, consult the :ref:`taskmanager` documentation).

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

To obtain the list of options supported by a particular resource manager e.g. ``slurm``::

    abirun.py . doc_manager slurm

Complete command line reference

.. argparse::
   :ref: abipy.scripts.abirun.get_parser
   :prog: abirun.py

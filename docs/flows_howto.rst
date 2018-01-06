.. _flows-howto:

************
Flows How-To
************

This is a list of Frequently Asked Questions about the AbiPy flows and the :ref:`abirun.py` script. 
Feel free to suggest new entries!

.. important::

    The execution of the flows require the proper configuration of ``manager.yml`` and,
    optionally, ``scheduler.yml``.
    Please consult the documentation at ....

Suggestions:

* Start with the examples available in examples/flows before embarking on large scale calculations.
* Make sure the Abinit executable compiled on the machine can be executed both on the frontend 
  and the compute node (ask your sysadmin)
* If you are running on clusters in which the architecture of the compute node is completely different
  from the one available on the frontend, use ``shell_runner``
* Use the ``debug`` command

Do not:

* Change manually the input files and the submission scripts
* Submit jobs manually when the scheduler is running 
* Use a too small delay for the scheduler 


.. contents::
   :backlinks: top

Get all the TaskManager options
-------------------------------

The :ref:`abidoc.py` script provides three commands to get the documentation
for the options supported in ``manager.yml`` and ``scheduler.yml``.

Use::

    abidoc.py manager

to document all the options supported by |TaskManager| and::

    abidoc.py scheduler

for the scheduler options.

If your environment is properly configured, you should be able to get
information about the Abinit version used by the AbiPy with::

    abidoc.py abibuild

.. code-block:: bash

    Abinit Build Information:
    Abinit version: 8.7.2
    MPI: True, MPI-IO: True, OpenMP: False
    Netcdf: True

    Use --verbose for additional info

.. important::

    Netcdf support must be activated in Abinit as AbiPy might use
    these files to extract data and/or fix runtime errors.

At this point, you can try to run a small flow for testing purpose with::

    abicheck.py --with-flow

Reduce the number of files produced by the Flow 
-----------------------------------------------

When running many calculations, 
Use ``prtwf -1`` to tell Abinit to produce the wavefunction file, only
if SCF cycle didn't converged so that AbiPy can reuse the file to restart the calculation.

Note that it's possibile to use::

    flow.use_smartio()

to activate this mode for all tasks that are not supposed to produce WFK files
for their children.

Extend tasks/works with specialized code
----------------------------------------

Remember that pickle_ does not support classes defined inside scripts. 
If you need to subclass one of the AbiPy Tasks/Works/Flows, define the subclass 
in a separated python module and import the module inside your script.

Kill a scheduler running in background
--------------------------------------

Use the official API::

    abirun.py FLOWDIR cancel

to cancel all jobs of the flow that are in queue and kill the scheduler.

Compare multiple output files
-----------------------------

The :ref:`abicomp.py` script

Try to understand why a task failed
------------------------------------

There are several reasons why a task could fail.
Some of these reasons could be related to hardaware failure, disk quota, 
OS or resource manager errors.
others are related to Abinit-specific errors.

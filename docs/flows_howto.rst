.. _flows-howto:

************
Flows How-To
************

This is a list of FAQs about the AbiPy flows and the :ref:`abirun.py` script. 
Feel free to suggest new entries!

.. important::

    The execution of the flows require the proper configuration of ``manager.yml`` and,
    optionally, ``scheduler.yml``.
    Please consult the documentation available via the abidoc.py script. See FAQ below.

Suggestions:

* Start with the examples available in examples/flows before embarking on large scale calculations.
* Make sure the Abinit executable compiled on the machine can be executed both on the front end 
  and the compute node (ask your sysadmin)
* If you are running on clusters in which the architecture of the compute node is completely different
  from the one available on the front end, use ``shell_runner``
* Use the ``debug`` command

Do not:

* Change manually the input files and the submission scripts
* Submit jobs manually when the scheduler is running 
* Use a too small delay for the scheduler 


.. contents::
   :backlinks: top

How to get all the TaskManager options
--------------------------------------

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

How to limit the number of cores used by the scheduler
------------------------------------------------------

Add the following options to scheduler.yml

.. code-block:: yaml

    # Limit on the number of jobs that can be present in the queue. (DEFAULT: 200)
    max_njobs_inqueue: 2

    # Maximum number of cores that can be used by the scheduler.
    max_ncores_used: 4

How to reduce the number of files produced by the Flow 
------------------------------------------------------

When running many calculations, 
Use ``prtwf -1`` to tell Abinit to produce the wavefunction file only
if SCF cycle didn't converged so that AbiPy can reuse the file to restart the calculation.

Note that it's possible to use::

    flow.use_smartio()

to activate this mode for all tasks that are not supposed to produce WFK files
for their children.

How to extend tasks/works with specialized code
-----------------------------------------------

Remember that pickle_ does not support classes defined inside scripts (`__main__`).
This means that `abirun.py` will likely raise an exception when trying to 
reconstruct the object from the pickle file:

.. code-block:: python

    AttributeError: Cannot get attribute 'MyWork' on <module '__main__' 

If you need to subclass one of the AbiPy Tasks/Works/Flows, define the subclass 
in a separated python module and import the module inside your script.
We suggest to create a python module in the AbiPy package e.g. `abipy/flowtk/my_works.py`
in order to have an absolute import that allows one to use

.. code-block:: python

    from abipy.flowtk.my_works import MyWork

in the script without worrying about relative paths and relative imports.


Kill a scheduler running in background
--------------------------------------

Use the official API::

    abirun.py FLOWDIR cancel

to cancel all jobs of the flow that are in queue and kill the scheduler.

Compare multiple output files
-----------------------------

The :ref:`abicomp.py` script

Try to understand why a task failed
-----------------------------------

There are several reasons why a task could fail.
Some of these reasons could be related to hardware failure, disk quota, 
OS errors or resource manager errors.
Others are related to Abinit-specific errors.

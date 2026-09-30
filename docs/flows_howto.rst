.. _flows-howto:

************
Flows How-To
************

This page collects frequently asked questions about AbiPy flows and the abirun.py_ script.
Feel free to suggest new entries!

.. important::

    Running a flow requires properly configured ``manager.yml`` and ``scheduler.yml`` files.
    The available options are documented by the abidoc.py script (see the FAQ below).

Suggestions:

* Start from the examples in `~abipy/examples/flows` before embarking on large-scale calculations.
* Make sure the Abinit executable compiled on your machine runs both on the front end
  and on the compute nodes (ask your sysadmin).
* If the compute nodes have a different architecture from
  the front end, use ``shell_runner``.
* Use the ``abirun.py FLOWDIR debug`` command to investigate problems.

Please DO NOT:

* Manually change input files or submission scripts while the scheduler is running.
* Manually submit jobs while the scheduler is running.
* Use a very short delay for the scheduler, as this may overload the front end.


.. contents::
   :backlinks: top

How to get all the TaskManager options
--------------------------------------

The abidoc.py_ script provides three commands to document
the options supported in ``manager.yml`` and ``scheduler.yml``.

Use::

    abidoc.py manager

to list all the options supported by the |TaskManager|, and::

    abidoc.py scheduler

for the scheduler options.

If your environment is properly configured, you can get
information about the Abinit version used by AbiPy with::

    abidoc.py abibuild

.. code-block:: bash

    Abinit Build Information:
    Abinit version: 8.7.2
    MPI: True, MPI-IO: True, OpenMP: False
    Netcdf: True

    Use --verbose for additional info

.. important::

    Netcdf support must be activated in Abinit, since AbiPy uses
    netcdf files to extract data and to fix runtime errors.

You can then run a small test flow with::

    abicheck.py --with-flow

How to limit the number of cores used by the scheduler
------------------------------------------------------

Add the following options to `~/.abinit/abipy/scheduler.yml`:

.. code-block:: yaml

    # Limit on the number of jobs that can be present in the queue. (DEFAULT: 200)
    max_njobs_inqueue: 2

    # Maximum number of cores that can be used by the scheduler.
    max_ncores_used: 4

How to reduce the number of files produced by the Flow
------------------------------------------------------

When running many calculations, use ``prtwf -1`` so that Abinit writes the wavefunction file only
if the SCF cycle did not converge. AbiPy can then use this file to restart the calculation.

You can also call::

    flow.use_smartio()

to activate this mode for all the tasks that are not expected to produce WFK files for their children.

How to extend tasks/works with specialized code
-----------------------------------------------

Remember that pickle_ does not support classes defined inside scripts (`__main__`).
As a consequence, `abirun.py` will raise the following exception when trying to
reconstruct the object from the pickle file:

.. code-block:: python

    AttributeError: Cannot get attribute 'MyWork' on <module '__main__'

If you need to subclass one of the AbiPy Tasks/Works/Flows, define the subclass
in a separate Python module and import it in your script.
We suggest creating the module inside the AbiPy package, e.g. `abipy/flowtk/my_works.py`,
so that it has an absolute import path and you can write

.. code-block:: python

    from abipy.flowtk.my_works import MyWork

in your script without worrying about relative paths and relative imports.


Kill a scheduler running in background
--------------------------------------

Use::

    abirun.py FLOWDIR cancel

to cancel all the jobs of the flow that are in the queue and kill the scheduler.

Compare multiple output files
-----------------------------

Use the abicomp.py_ script.

Try to understand why a task failed
-----------------------------------

A task can fail for several reasons.
Some are related to hardware failures, disk quotas, OS errors or resource manager errors,
while others are caused by Abinit-specific errors.

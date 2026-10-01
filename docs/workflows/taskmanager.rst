.. contents::
   :backlinks: top

.. _taskmanager:

^^^^^^^^^^^
TaskManager
^^^^^^^^^^^

Besides post-processing tools and a programmatic interface to generate input files,
AbiPy provides a Pythonic API to run small Abinit tasks directly or to submit calculations on supercomputing clusters.
This section explains how to create the configuration files needed to interface AbiPy with Abinit.

We assume that Abinit is already available on your machine and that you know how to configure
your environment so that the operating system can find and execute Abinit.
In other words, we assume you know how to set the ``$PATH`` and ``$LD_LIBRARY_PATH`` (``$DYLD_LIBRARY_PATH`` on macOS)
environment variables, load modules with ``module load``, run MPI applications with ``mpirun``, etc.

.. IMPORTANT::

    Before proceeding with the rest of the tutorial, make sure that you can run Abinit interactively
    with simple input files and that it works as expected.
    It is also a very good idea to run the Abinit test suite with the `runtests.py script <https://asciinema.org/a/40324>`_
    before starting production calculations.

.. TIP::

    A pre-compiled sequential version of Abinit for Linux and macOS can be installed directly
    from the conda-forge channel with::

        conda install abinit -c conda-forge

--------------------------------
How to configure the TaskManager
--------------------------------

The ``TaskManager`` takes care of task submission.
This includes creating the submission script,
initializing the environment and optimizing the parallel execution
(number of MPI processes, number of OpenMP threads, automatic parallelization with the Abinit ``autoparal`` feature).

AbiPy reads the information needed to create the correct ``TaskManager`` for a specific cluster (or personal computer)
from the ``manager.yml`` configuration file.
This file is written in YAML_, a human-readable data serialization language commonly used for configuration files.
A good introduction to the YAML syntax is available `here <http://yaml.org/spec/1.1/#id857168>`_;
see also this `reference card <http://www.yaml.org/refcard.html>`_,
and experiment with the syntax using a `YAML validator <https://yamline.com/validator/>`_.

By default, AbiPy looks for a ``manager.yml`` file first in the current working directory, i.e.
the directory in which you run your script, and then in ``$HOME/.abinit/abipy``.
If no file is found, the code aborts immediately.

The ``TaskManager`` also needs to know the type of queueing system available on the cluster,
the list of queues and their specifications.
In AbiPy, queueing systems (resource managers) are supported via ``qadapters``.
At the time of writing (|today|), AbiPy provides ``qadapters`` for the following resource managers:

* ``shell``
* pbspro_
* slurm_
* IBM loadleveler_
* moab_
* sge_
* torque_

Manager configuration files for typical cases are available in ``~abipy/data/managers``.

We first discuss how to configure AbiPy on a personal computer and then look at the more
complicated case in which calculations must be submitted to a queue.

-----------------------------------
TaskManager for a personal computer
-----------------------------------

Let's start from the simplest case, i.e. a personal computer on which we can execute
applications directly from the shell (``qtype: shell``).
Here the configuration file is relatively simple, because we can run Abinit
directly without generating a script and submitting it to a resource manager.
In its simplest form, the ``manager.yml`` file consists of a list of ``qadapters``:

.. code-block:: yaml

    qadapters:
        -  # qadapter_0
        -  # qadapter_1

Each item in the ``qadapters`` list is a YAML dictionary with the following sub-dictionaries:

``queue``
    Dictionary with the name of the queue and optional parameters
    used to build and customize the header of the submission script.

``job``
    Dictionary with the options used to set up the environment before submitting the job.

``limits``
    Dictionary with the constraints that must be satisfied to run with this ``qadapter``.

``hardware``
    Dictionary with information on the hardware available on this queue.
    It is used by the Abinit ``autoparal`` feature to optimize parallel execution.

The ``qadapter`` is therefore responsible for all the interactions with a specific
queue management system (shell, Slurm, PBS, etc.), including the format of the
submission script, job submission and job management.

.. NOTE::

    Multiple ``qadapters`` are useful on clusters with several queues,
    but we postpone the discussion of this rather technical point.
    For now, we use a ``manager.yml`` with a single adapter.

A typical configuration file used on a laptop to run jobs via the shell is:

.. code-block:: yaml

    qadapters: # List of `qadapters` objects  (just one in this simplified example)

    -  priority: 1
       queue:
            qtype: shell        # "Submit" jobs via the shell.
            qname: localhost    # "Submit" to the localhost queue
                                # (it's a fake queue in this case)

        job:
            pre_run: "export PATH=$HOME/git_repos/abinit/build_gcc/src/98_main:$PATH"
            mpi_runner: "mpirun"

        limits:
            timelimit: 1:00:00   #  Time-limit for each task.
            max_cores: 2         #  Max number of cores that can be used by a single task.

        hardware:
            num_nodes: 1
            sockets_per_node: 1
            cores_per_socket: 2
            mem_per_node: 4 Gb

The ``job`` section is the most critical one, in particular the ``pre_run`` option,
which is executed by the shell script before invoking Abinit.
In this example, the Abinit executable is not in the default search path,
so the directory containing the Abinit executables must be prepended to ``$PATH``.
Adapt ``pre_run`` to your Abinit installation and make sure that ``mpirun`` is also in ``$PATH``.
If you do not use a parallel version of Abinit, simply set ``mpi_runner: null``
(``null`` is the YAML_ equivalent of Python's ``None``).
This approach also lets you safely switch between multiple Abinit versions.

Copy this example and adapt the ``hardware`` and ``limits`` sections to
your machine. In particular, make sure that ``max_cores`` does not exceed the number of physical cores
available on your computer.
Save the file in the current working directory and run the abicheck.py_ script provided by AbiPy.
If everything is configured properly, you should see something like this in the terminal:

.. command-output:: abicheck.py --no-colors

This message tells us that everything is in place and we can finally run our first calculation.

.. note::

    This laptop has 1 socket with 2 CPUs and 4 Gb of memory in total, so we don't want to run
    Abinit tasks with more than 2 CPUs. This is why ``max_cores`` is set to 2.
    The ``timelimit`` option is ignored with ``qtype: shell``, but it becomes
    important when submitting jobs on a cluster: the value is used to generate the submission script,
    and Abinit uses it to exit from iterative algorithms (e.g. the SCF cycle) before the time limit is reached
    and write files from which the calculation can be restarted.

The directory ``~abipy/data/runs`` contains Python scripts that generate workflows for typical ab-initio calculations.
Here we focus on configuring the manager and executing the flow, so we don't discuss how to
generate input files and create Flow objects in Python.
This topic is covered in detail in our collection of `jupyter notebooks
<http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_.

Let's start from the simplest example, the ``run_si_ebands.py`` script, which generates
a flow to compute the band structure of silicon at the Kohn-Sham level:
a GS calculation to get the density, followed by an NSCF run along a k-path in the first Brillouin zone.

Go to ``~abipy/data/runs`` and execute ``run_si_ebands.py`` to generate the flow::

    cd ~abipy/data/runs
    ./run_si_ebands.py

You should now have a ``flow_si_ebands`` directory with the following structure:

.. code-block:: console

    tree flow_si_ebands/

    flow_si_ebands/
    ├── __AbinitFlow__.pickle
    ├── indata
    ├── outdata
    ├── tmpdata
    └── w0
    ├── indata
    ├── outdata
    ├── t0
    │   ├── indata
    │   ├── job.sh
    │   ├── outdata
    │   ├── run.abi
    │   ├── run.files
    │   └── tmpdata
    ├── t1
    │   ├── indata
    │   ├── job.sh
    │   ├── outdata
    │   ├── run.abi
    │   ├── run.files
    │   └── tmpdata
    └── tmpdata

    15 directories, 7 files

``w0/`` contains the input files of the first workflow (the only one in this example).
``w0/t0/`` and ``w0/t1/`` contain the input files needed for the SCF and the NSCF run, respectively.

Note that all the task directories (``w0/t0``, ``w0/t1``) have the same structure:

   * ``run.abi``: Abinit input file.
   * ``run.files``: Abinit files file.
   * ``job.sh``: Submission/shell script.
   * ``outdata``: Directory with output data files.
   * ``indata``: Directory with input data files.
   * ``tmpdata``: Directory with temporary files.

.. DANGER::

   ``__AbinitFlow__.pickle`` is the pickle file used to save the status of the `Flow`. Don't touch it!

The ``job.sh`` script has been generated by the ``TaskManager`` from the information in ``manager.yml``.
Since we are using ``qtype: shell``, it is a simple shell script that executes the code directly.
The script becomes more complex when jobs are submitted to a resource manager on a cluster.

We usually interact with an AbiPy flow via the :ref:`abirun.py <abirun-script>` script, whose syntax is::

     abirun.py FLOWDIR command [options]

where ``FLOWDIR`` is the directory containing the flow and ``command`` is the action to perform.
Use ``abirun.py --help`` to get the list of available commands and ``abirun.py COMMAND --help`` to see the
options supported by ``COMMAND``.

``abirun.py`` reconstructs the Python Flow from the ``__AbinitFlow__.pickle`` file in ``FLOWDIR``
and calls the methods of the object according to the command-line options.

Use::

    abirun.py flow_si_ebands status

to get a summary of the status of the tasks, and::

    abirun.py flow_si_ebands deps

to print the dependencies of the tasks in textual format.

.. code-block:: console

    <ScfTask, node_id=75244, workdir=flow_si_ebands/w0/t0>

    <NscfTask, node_id=75245, workdir=flow_si_ebands/w0/t1>
      +--<ScfTask, node_id=75244, workdir=flow_si_ebands/w0/t0>

.. TIP::

    Alternatively, use ``abirun.py flow_si_ebands networkx``
    to visualize the dependencies with the networkx_ package.

Here we have a flow with one work (``w0``) containing two tasks.
The second task (``w0/t1``) depends on the first one, a ``ScfTask``;
more specifically, ``w0/t1`` needs the density file produced by ``w0/t0``.
This means that ``w0/t1`` cannot be executed until the first task has completed.
AbiPy is aware of this dependency and uses it to manage the submission and execution
of the flow.

Tasks can be launched with two commands: ``single`` and ``rapid``.
The ``single`` command executes the first task of the flow in the ``READY`` state, i.e. the first task
whose dependencies have been fulfilled.
``rapid``, on the other hand, submits **all the tasks** of the flow in the ``READY`` state.
Let's run the flow with the ``rapid`` command:

.. code-block:: console

    abirun.py flow_si_ebands rapid

    Running on gmac2 -- system Darwin -- Python 2.7.12 -- abirun-0.1.0
    Number of tasks launched: 1

    Work #0: <BandStructureWork, node_id=75239, workdir=flow_si_ebands/w0>, Finalized=False
    +--------+-------------+-----------------+--------------+------------+----------+-----------------+----------+-----------+
    | Task   | Status      | Queue           | MPI|Omp|Gb   | Warn|Com   | Class    | Sub|Rest|Corr   | Time     |   Node_ID |
    +========+=============+=================+==============+============+==========+=================+==========+===========+
    | w0_t0  | Submitted   | 71573@localhost | 2|  1|2.0    | 1|  0      | ScfTask  | (1, 0, 0)       | 0:00:00Q |     75240 |
    +--------+-------------+-----------------+--------------+------------+----------+-----------------+----------+-----------+
    | w0_t1  | Initialized | None            | 1|  1|2.0    | NA|NA      | NscfTask | (0, 0, 0)       | None     |     75241 |
    +--------+-------------+-----------------+--------------+------------+----------+-----------------+----------+-----------+


What happened here?
The ``rapid`` command tried to execute all the ``READY`` tasks but, since the second task depends
on the first one, only the first task was submitted.
Note that the SCF task (``w0_t0``) has been submitted with 2 MPI processes.
Before submitting a task, AbiPy
invokes Abinit to get all the parallel configurations compatible with the limits
specified by the user (e.g. ``max_cores``), selects an "optimal" configuration according
to a given policy, and submits the task with the optimized parameters.
Since no other task can be executed at this point, the script exits,
and we have to wait for the SCF task to complete before running the second part of the flow.

At each iteration, :ref:`abirun.py <abirun-script>` prints a table with the status of the tasks.
The columns have the following meaning:

``Queue``
    String of the form ``JobID @ QueueName``, where JobID is the process identifier when running in the shell,
    or the job ID assigned by the resource manager (e.g. slurm) when submitting to a queue.
``MPI``
    Number of MPI processes. This value is obtained automatically by calling Abinit in ``autoparal`` mode
    and cannot exceed ``max_ncpus``.
``OMP``
    Number of OpenMP threads.
``Gb``
    Memory requested in Gb. Meaningless when ``qtype: shell``.
``Warn``
    Number of warning messages found in the log file.
``Com``
    Number of comments found in the log file.
``Sub``
    Number of submissions. It can be > 1 if AbiPy encounters a problem and resubmits the task
    with different parameters, without performing any operation that changes the physics of the system.
``Rest``
    Number of restarts. AbiPy restarts the job if convergence has not been reached.
``Corr``
    Number of corrections performed by AbiPy to fix runtime errors.
    These operations can change the physics of the system.
``Time``
    Time spent in the queue (if the string ends with Q) or running time (if the string ends with R).
``Node_ID``
    Identifier used by AbiPy for each node of the flow.

.. NOTE::
     When jobs are submitted through the shell, there is almost no difference between
     job submission and job execution. The situation is completely different with
     a resource manager, where the task gets a priority value and enters the queue.

If you run ``status`` again, you should see that the first task has completed.
We can now run ``rapid`` again to launch the |NscfTask|.
The second task won't take long; if you run ``status`` once more, you should see that the entire flow
has completed successfully.

To understand in more detail what happened, use the ``history`` command to get
the list of operations performed by AbiPy on each task:

.. code-block:: console

    abirun.py flow_si_ebands history

    ==============================================================================================================================
    =================================== <ScfTask, node_id=75244, workdir=flow_si_ebands/w0/t0> ===================================
    ==============================================================================================================================
    [Mon Mar  6 21:46:00 2017] Status changed to Ready. msg: Status set to Ready
    [Mon Mar  6 21:46:00 2017] Setting input variables: {'max_ncpus': 2, 'autoparal': 1}
    [Mon Mar  6 21:46:00 2017] Old values: {'max_ncpus': None, 'autoparal': None}
    [Mon Mar  6 21:46:00 2017] Setting input variables: {'npband': 1, 'bandpp': 1, 'npimage': 1, 'npspinor': 1, 'npfft': 1, 'npkpt': 2}
    [Mon Mar  6 21:46:00 2017] Old values: {'npband': None, 'npfft': None, 'npkpt': None, 'npimage': None, 'npspinor': None, 'bandpp': None}
    [Mon Mar  6 21:46:00 2017] Status changed to Initialized. msg: finished autoparallel run
    [Mon Mar  6 21:46:00 2017] Submitted with MPI=2, Omp=1, Memproc=2.0 [Gb] submitted to queue
    [Mon Mar  6 21:46:15 2017] Task completed status set to ok based on abiout
    [Mon Mar  6 21:46:15 2017] Finalized set to True

    =============================================================================================================================
    ================================== <NscfTask, node_id=75245, workdir=flow_si_ebands/w0/t1> ==================================
    =============================================================================================================================
    [Mon Mar  6 21:46:15 2017] Status changed to Ready. msg: Status set to Ready
    [Mon Mar  6 21:46:15 2017] Adding connecting vars {u'irdden': 1}
    [Mon Mar  6 21:46:15 2017] Setting input variables: {u'irdden': 1}
    [Mon Mar  6 21:46:15 2017] Old values: {u'irdden': None}
    [Mon Mar  6 21:46:15 2017] Setting input variables: {'max_ncpus': 2, 'autoparal': 1}
    [Mon Mar  6 21:46:15 2017] Old values: {'max_ncpus': None, 'autoparal': None}
    [Mon Mar  6 21:46:15 2017] Setting input variables: {'npband': 1, 'bandpp': 1, 'npimage': 1, 'npspinor': 1, 'npfft': 1, 'npkpt': 2}
    [Mon Mar  6 21:46:15 2017] Old values: {'npband': None, 'npfft': None, 'npkpt': None, 'npimage': None, 'npspinor': None, 'bandpp': None}
    [Mon Mar  6 21:46:15 2017] Status changed to Initialized. msg: finished autoparallel run
    [Mon Mar  6 21:46:15 2017] Submitted with MPI=2, Omp=1, Memproc=2.0 [Gb] submitted to queue
    [Mon Mar  6 21:49:48 2017] Task completed status set to ok based on abiout
    [Mon Mar  6 21:49:48 2017] Finalized set to True


A closer look at the logs reveals that, before submitting the first task, AbiPy executed
Abinit in ``autoparal`` mode to get the list of possible parallel configurations, and only then submitted the calculation.
AbiPy then monitors the output files produced by the task to understand what is happening.
When the first task completes, the status of the second task is automatically changed to ``READY``,
the ``irdden`` input variable is added to the input file of the second task, and a symbolic link to
the ``DEN`` file produced by ``w0/t0`` is created in the ``indata`` directory of ``w0/t1``.
Another autoparal run is executed for the NSCF calculation, and the second task is finally submitted.

The command line interface is very flexible and sometimes it is the only tool available.
In some cases, however, we would like a global view of what is happening.
The command::

    $ abirun.py flow_si_ebands notebook

generates a jupyter_ notebook with predefined Python code that
displays a graphical representation of the status of the flow in a web browser
(requires jupyter_, nbformat_ and, of course, a web browser).

Expert users may prefer::

    $ abirun.py flow_si_ebands ipython

to open the flow in the ipython_ shell and access its API directly.

Once ``manager.yml`` is properly configured, you can
use AbiPy objects to invoke Abinit and perform useful operations.
For example, the |AbinitInput| object can give you the list of k-points in the IBZ,
the list of independent DFPT perturbations, the parallel configurations reported by ``autoparal``, etc.

This programmatic interface can be used in scripts to simplify the creation of input files and workflows.
For example, you can call Abinit to get the list of perturbations for each q-point in the IBZ and then
automatically generate all the input files for the DFPT calculations (this is actually how
the AbiPy factory functions generate DFPT workflows).

``manager.yml`` is also used to invoke other executables (``anaddb``, ``optic``, ``mrgddb``, etc.),
thus providing an interface between Python and the Fortran executables.
Thanks to this interface, relatively simple ab-initio calculations can be performed directly in AbiPy.
For instance, you can open a ``DDB`` file in a jupyter notebook, call ``anaddb`` to compute
the phonon frequencies, and plot the DOS and the phonon band structure with matplotlib_.

.. TIP::

    The command::

        abirun.py . doc_manager

    gives the full documentation of the entries of ``manager.yml``.

.. command-output:: abirun.py . doc_manager

.. _scheduler:

------------------------------
How to configure the scheduler
------------------------------

In the previous example, we ran a simple band structure calculation for silicon in a few seconds
on a laptop, but more complex flows may take hours or even days to complete.
In such cases, the ``single`` and ``rapid`` commands are not practical, because we would have
to monitor the flow and re-run ``abirun.py`` every time a new task becomes ``READY``.
It is much easier to delegate all this repetitive work to a ``python scheduler``:
a process that runs in the background, submits tasks automatically and performs the actions
needed to complete the flow.

The scheduler parameters are declared in the YAML_ file ``scheduler.yml``.
As before, AbiPy looks first in the working directory and then in ``$HOME/.abinit/abipy``.
Create a ``scheduler.yml`` file in the working directory by copying the example below:

.. code-block:: yaml

    seconds: 5   # number of seconds to wait.
    #minutes: 0  # number of minutes to wait.
    #hours: 0    # number of hours to wait.

This file tells the scheduler to wake up every 5 seconds, inspect the status of the tasks
in the flow, and perform the actions needed to reach completion.

.. IMPORTANT::

    Remember to set the time interval to a reasonable value.
    A short interval increases the submission rate, but it also increases the CPU load
    and the pressure on the hardware and the resource manager.
    A very long interval, on the other hand, can reduce throughput, especially
    when submitting many small jobs.

We are now ready to run our first calculation with the scheduler.
To make things more interesting, we execute a slightly more complex flow that computes
the G0W0 corrections to the direct band gap of silicon at the Gamma point.
The flow consists of the following six tasks:

- 0: Ground-state calculation to get the density.
- 1: NSCF calculation with several empty states.
- 2: Calculation of the screening using the WFK produced by task 1.
- 3-4-5: Evaluation of the self-energy matrix elements with different values of nband,
  using the WFK produced by task 1 and the SCR file produced by task 2.

Generate the flow with::

    ./run_si_g0w0.py

and let the scheduler manage the submission with::

     abirun.py flow_si_g0w0 scheduler

You should see the following output in the terminal:

.. code-block:: console

    abirun.py flow_si_g0w0 scheduler

    Abipy Scheduler:
    PyFlowScheduler, Pid: 72038
    Scheduler options: {'seconds': 10, 'hours': 0, 'weeks': 0, 'minutes': 0, 'days': 0}

``Pid`` is the process identifier of the scheduler (also saved in the ``_PyFlowScheduler.pid`` file).

.. IMPORTANT::

    A ``_PyFlowScheduler.pid`` file in ``FLOWDIR`` means that a scheduler is running the flow.
    There must be only one scheduler associated with a given flow.

The scheduler makes AbiPy flows much more powerful, since
complex ab-initio workflows can be automated with little effort: write
a Python script that implements the flow, run it with
``abirun.py FLOWDIR scheduler``, and finally analyze the results with the AbiPy/pymatgen tools.
Even complex convergence studies for G0W0 calculations can be implemented along these lines,
as shown in this `video <https://youtu.be/M9C6iqJsvJI>`_.
Sooner or later, however, flows become too large or too computationally expensive
to run on a personal computer, and we have to move to a supercomputing center.
The next section explains how to configure AbiPy to run on a cluster with a queue management system.

.. TIP::

    Use ``abirun.py . doc_scheduler`` to get the full list of options supported by the scheduler.

.. command-output:: abirun.py doc_scheduler

.. _abipy-on-cluster:

------------------------------
Configuring AbiPy on a cluster
------------------------------

In this section, we discuss how to configure the manager to run flows on a cluster.
The configuration depends on the specific queue management system (Slurm, PBS, etc.), so
we assume that you are already familiar with job submission and know which options
must be specified in the submission script for your job to be accepted
and executed by the resource manager (username, queue name, memory, ...).

Let's assume that our computing center uses slurm_ and that jobs must be submitted to the ``default_queue`` partition.
In the best case, the system administrator already provides
an ``Abinit module`` that can be loaded with ``module load`` before invoking the code.
To make things a bit more challenging, however, we assume that we compiled our own version of Abinit
in the build directory ``${HOME}/git_repos/abinit/build_impi`` using the following two modules
provided by the system administrator::

    compiler/intel/composerxe/2013_sp1.1.106
    intelmpi

In this case, we have to configure the environment carefully, because the Slurm submission
script must load the modules and modify ``$PATH`` so that our version of Abinit can be found.
A ``manager.yml`` with a single ``qadapter`` looks like this:

.. code-block:: yaml

    qadapters:
      - priority: 1

        queue:
           qtype: slurm
           qname: default_queue
           qparams: # Slurm options added to job.sh
              mail_type: FAIL
              mail_user: john@doe

        job:
            modules:
                - compiler/intel/composerxe/2013_sp1.1.106
                - intelmpi
            shell_env:
                 PATH: ${HOME}/git_repos/abinit/build_impi/src/98_main:$PATH
            pre_run:
               - ulimit -s unlimited
            mpi_runner: mpirun

        limits:
           timelimit: 0:20:0
           max_cores: 16
           min_mem_per_proc: 1Gb

        hardware:
            num_nodes: 120
            sockets_per_node: 2
            cores_per_socket: 8
            mem_per_node: 64Gb

.. TIP::

    The command::

        abirun.py FLOWDIR doc_manager script

    prints the submission script that AbiPy will generate at runtime.

Let's discuss the options in more detail, starting from the ``queue`` section:

``qtype``
    String specifying the resource manager. This option tells AbiPy which ``qadapter`` to use to generate
    and submit the submission scripts and to kill jobs in the queue, as well as how to interpret the other options passed by the user.

``qname``
    Name of the submission queue (string, MANDATORY).

``qparams``
    Dictionary with the parameters passed to the resource manager.
    We use the *normalized* version of the options, i.e. dashes in the official parameter name
    are replaced by underscores, e.g. ``--mail-type`` becomes ``mail_type``.
    For the list of supported options, use the ``doc_manager`` command.
    Use ``qverbatim`` to pass additional options that are not included in the template.

Note that we do not specify the number of cores in ``qparams``, since AbiPy finds
an appropriate value at runtime.

The ``job`` section is the most critical one, because it defines how to set up the environment
and how to run the code.
The ``modules`` entry lists the modules to load, while ``shell_env`` modifies the
``$PATH`` environment variable so that the OS can find our Abinit executable.

.. IMPORTANT::

    Some resource managers execute your ``.bashrc`` before loading the new modules.

We also increase the stack size with ``ulimit`` before running the code, and we run Abinit
with the ``mpirun`` provided by the modules.

The ``limits`` section defines the constraints that must be satisfied to run on this queue,
while ``hardware`` describes the hardware available on this queue.
Every job has a ``timelimit`` of 20 minutes, cannot use more than ``max_cores`` cores,
and the first submission requests 1 Gb of memory per process.
The actual number of cores is determined at runtime by calling Abinit in ``autoparal`` mode
to get all the parallel configurations up to ``max_cores``.
If a job is killed because of insufficient memory, AbiPy resubmits the task with more memory,
up to the maximum given by ``mem_per_node``.

``limits`` supports more advanced options as well, and new options
will be added over time.

To get the complete list of options supported by the Slurm ``qadapter``, use:

.. command-output:: abirun.py . doc_manager slurm

.. IMPORTANT::

    If you need to cancel all tasks that have been submitted to the resource manager, use::

        abirun.py FLOWDIR cancel

    The script asks for confirmation before killing all the jobs belonging to the flow.

Once ``manager.yml`` is properly configured for your cluster, you can
use the scheduler to automate job submission.
Your flows will most likely take hours or even days to complete and, in principle,
you would need to keep an active connection to the machine to keep the scheduler alive
(if your session expires, all the subprocesses launched from your terminal,
including the Python scheduler, are killed).
Fortunately, the standard Unix tool ``nohup`` comes to the rescue.

For long-running jobs, we strongly recommend starting the scheduler with::

     nohup abirun.py FLOWDIR scheduler > sched.stdout 2> sched.stderr &

This command runs the scheduler in the background and redirects ``stdout`` and ``stderr``
to ``sched.stdout`` and ``sched.stderr``, respectively.
The process identifier of the scheduler is saved in the ``_PyFlowScheduler.pid`` file in ``FLOWDIR``,
which is removed automatically when the scheduler completes.
Thanks to ``nohup``, we can close the session, let the scheduler work overnight,
and reconnect the next day to collect the results.

.. IMPORTANT::

    Use ``abirun.py FLOWDIR cancel`` to cancel the jobs of a flow managed by
    a scheduler. AbiPy detects that a scheduler is attached to the flow,
    cancels the jobs of the flow and kills the scheduler as well.


.. _inspecting-the-flow:

-------------------
Inspecting the Flow
-------------------

:ref:`abirun.py <abirun-script>` also provides tools to analyze the results of the flow at runtime.
The simplest command is::

    abirun.py FLOWDIR tail

which is similar to Unix ``tail``, but a bit smarter:
it prints only the final part of the output files
of the tasks that are ``RUNNING``.

If matplotlib_ is installed, you may want to use::

    $ abirun.py FLOWDIR inspect

Several AbiPy tasks provide an ``inspect`` method that produces matplotlib figures
with data extracted from the output files.
For example, a ``GsTask`` plots the evolution of the ground-state SCF cycle.
The ``inspect`` command of :ref:`abirun.py <abirun-script>` simply loops over the tasks of the flow and
calls the ``inspect`` method of each one.

The command::

    abirun.py FLOWDIR inputs

prints the input files of the tasks (use ``--nids`` to select a subset of
tasks or, alternatively, replace ``FLOWDIR`` with the ``FLOWDIR/w0/t0`` syntax).

The command::

    abirun.py FLOWDIR listext EXTENSION

prints a table with the nodes of the flow that have produced an Abinit output file with the given
extension. For example::

    abirun.py FLOWDIR listext GSR.nc

shows the nodes of the flow that have produced a GSR.nc_ file.

The command::

    abirun.py FLOWDIR notebook

generates a jupyter_ notebook with predefined Python code that
displays a graphical representation of the status of the flow in a web browser
(requires jupyter_, nbformat_ and, of course, a web browser).

Expert users may prefer::

    abirun.py FLOWDIR ipython

to open the flow in the ipython_ shell and access its API directly.


.. _event-handlers:

--------------
Event handlers
--------------

An event handler is an action executed in response to a particular event.
AbiPy tasks come with built-in event handlers that are executed
to fix typical Abinit runtime errors.

To list the event handlers installed in a given flow, use::

    abirun.py FLOWDIR handlers

The ``--verbose`` option gives a more detailed description of the actions performed
by the event handlers:

.. code-block:: console

    abirun.py FLOWDIR handlers --verbose

    List of event handlers installed:
    event name = !DilatmxError
    event documentation:

    This Error occurs in variable cell calculations when the increase in the
    unit cell volume is too large.

    handler documentation:

    Handle DilatmxError. Abinit produces a netcdf file with the last structure before aborting
    The handler changes the structure in the input with the last configuration and modify the value of dilatmx.

    event name = !TolSymError
    event documentation:

    Class of errors raised by Abinit when it cannot detect the symmetries of the system.
    The handler assumes the structure makes sense and the error is just due to numerical inaccuracies.
    We increase the value of tolsym in the input file (default 1-8) so that Abinit can find the space group
    and re-symmetrize the input structure.

    handler documentation:

    Increase the value of tolsym in the input file.

    event name = !MemanaError
    event documentation:

    Class of errors raised by the memory analyzer.
    (the section that estimates the memory requirements from the input parameters).

    handler documentation:

    Set mem_test to 0 to bypass the memory check.

    event name = !MemoryError
    event documentation:

    This error occurs when a checked allocation fails in Abinit
    The only way to go is to increase memory

    handler documentation:

    Handle MemoryError. Increase the resources requirements

.. NOTE::

    New error handlers will be added in future versions of AbiPy/Abinit.
    Please let us know if you need handlers for errors that commonly occur in your calculations.

.. _flow-troubeshooting:

---------------
Troubleshooting
---------------

Two :ref:`abirun.py <abirun-script>` commands are especially useful when something goes wrong: ``events`` and ``debug``.

To print the Abinit events (warnings, errors, comments) found in the log files of the tasks, use::

    abirun.py FLOWDIR events

To analyze error files and log files for possible error messages, use::

    abirun.py FLOWDIR debug

By default, these commands analyze the entire flow, so the terminal output can be very verbose.
To focus on a particular task, e.g. ``w0/t1``, use::

    abirun.py FLOWDIR/w0/t1 events

To select all the tasks in a work directory, e.g. ``w0``, use::

    abirun.py FLOWDIR/w0 events

To select an arbitrary subset of nodes of the flow, use::

    abirun.py FLOWDIR events -nids=12,13,16

where ``nids`` is a list of AbiPy node identifiers.

.. TIP::

    ``abirun.py events --help`` is your best friend.

.. command-output:: abirun.py events --help

To get information about the Abinit executable called by AbiPy, use::

    abirun.py abibuild

or the verbose variant::

    abirun.py abibuild --verbose

TODO: How to reset tasks

.. _task_policy:

----------
TaskPolicy
----------

At this point, you may wonder why all these parameters are needed in the configuration file.
The reason is that, before submitting a job to a resource manager, AbiPy uses the autoparal
feature of ABINIT to get all the possible parallel configurations with ``ncpus <= max_cores``.
Based on these results, AbiPy selects the "optimal" one and updates the ABINIT input file
and the submission script accordingly.
This is a very useful feature, especially for calculations with ``paral_kgb=1``, which require
setting ``npkpt``, ``npfft``, ``npband``, etc.
If more than one ``QueueAdapter`` is specified, AbiPy first computes all the possible
configurations and then selects the "optimal" ``QueueAdapter`` according to a given policy.

In some cases, you may want to impose constraints on the "optimal" configuration.
For example, you may want to select only the configurations whose parallel efficiency is greater than 0.7
and whose number of MPI processes is divisible by 4.
Such constraints can be enforced via the ``condition`` dictionary, whose syntax is similar to
the one used in mongodb_.

.. code-block:: yaml

    policy:
        autoparal: 1
        max_ncpus: 10
        condition: {$and: [ {"efficiency": {$gt: 0.7}}, {"tot_ncpus": {$divisible: 4}} ]}

The parallel efficiency is defined as :math:`\epsilon = \dfrac{T_1}{N T_N}`, where :math:`N` is the number
of MPI processes and :math:`T_j` is the wall time needed to complete the calculation with :math:`j` MPI processes.
For perfect scaling, :math:`\epsilon` equals one.
The parallel speedup with :math:`N` processes is given by :math:`S = T_1 / T_N`.
Note that ``autoparal = 1`` automatically modifies both the ``job.sh`` script and the input file
so that the job runs in parallel with the optimal configuration.
For example, you can use ``paral_kgb = 1`` in GS calculations, and AbiPy will automatically set
``npband``, ``npfft``, ``npkpt``, ... for you!
If no configuration satisfies the given condition, AbiPy uses the configuration
with the highest parallel speedup (not necessarily the most efficient one).

``policy``
    This section controls the automatic parallelization of the run. Here, AbiPy uses
    the ``autoparal`` capabilities of Abinit to determine an optimal configuration with
    at **most** ``max_ncpus`` MPI processes. Setting ``autoparal`` to 0 disables automatic parallelization.
    Other values of ``autoparal`` are not supported.

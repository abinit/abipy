.. contents::
   :backlinks: top

.. _taskmanager:

^^^^^^^^^^^
TaskManager
^^^^^^^^^^^

Besides post-processing tools and a programmatic interface to generate input files,
AbiPy also provides a pythonic API to execute small Abinit tasks directly or submit calculations on supercomputing clusters.
This section discusses how to create the configuration files required to interface AbiPy with Abinit.

We assume that Abinit is already available on your machine and that you know how to configure
your environment so that the operating system can load and execute Abinit.
In other words, we assume that you know how to set the ``$PATH`` and ``$LD_LIBRARY_PATH`` (``$DYLD_LIBRARY_PATH`` on Mac) 
environment variables, load modules with ``module load``, run MPI applications with ``mpirun``, etc.

.. IMPORTANT:: 

    Please make sure that you can execute Abinit interactively with simple input files and 
    that it works as expected before proceeding with the rest of the tutorial.
    It's also a very good idea to run the Abinit test suite with the `runtest.py script <https://asciinema.org/a/40324>`_ 
    before running production calculations.

.. TIP::

    A pre-compiled sequential version of Abinit for Linux and OSx can be installed directly 
    from the abinit-channel_ on the anaconda cloud with::

        conda install abinit --channel abinit

--------------------------------
How to configure the TaskManager
--------------------------------

The ``TaskManager`` takes care of task submission. 
This includes the creation of the submission script,
the initialization of the environment as well as the optimization of the parallel algorithms
(number of MPI processes, number of OpenMP threads, automatic parallelization with Abinit ``autoparal`` feature). 

AbiPy obtains the information needed to create the correct ``TaskManager`` for a specific cluster (personal computer)
from the ``manager.yml`` configuration file.
The file is written in YAML_ a human-readable data serialization language commonly used for configuration files
(a good introduction to the YAML syntax can be found `here <http://yaml.org/spec/1.1/#id857168>`_.
See also this `reference card <http://www.yaml.org/refcard.html>`_)

By default, AbiPy looks for a ``manager.yml`` file in the current working directory i.e.
the directory in which you execute your script in first and then inside ``$HOME/.abinit/abipy``.
If no file is found, the code aborts immediately.

An important piece of information for the ``TaskManager`` is the type of queueing system available on the cluster,
the list of queues and their specifications. 
In AbiPy queueing systems or resource managers are supported via ``quadapters``.
At the time of writing (|today|), AbiPy provides ``qadapters`` for the following resource managers:

* ``shell``
* pbspro_
* slurm_
* IBM loadleveler_
* moab_
* sge_
* torque_

Manager configuration files for typical cases are available inside ``~abipy/data/managers``.

We first discuss how to configure AbiPy on a personal computer and then we look at the more
complicated case in which the calculation must be submitted to a queue.

-----------------------------------
TaskManager for a personal computer
-----------------------------------

Let's start from the simplest case i.e. a personal computer in which we can execute 
applications directly from the shell (``qtype: shell``).
In this case, the configuration file is relatively easy because we can run Abinit
directly without having to generate and submit a script to the resource manager.
In its simplest form, the ``manager.yml`` file consists of a list of ``qadapters``:

.. code-block:: yaml

    qadapters:
        -  # qadapter_0
        -  # qadapter_1

Each item in the ``qadapters`` list is essentially a YAML dictionary with the following sub-dictionaries:

``queue``
    Dictionary with the name of the queue and optional parameters 
    used to build and customize the header of the submission script.

``job``
    Dictionary with the options used to prepare the environment before submitting the job.

``limits``
    Dictionary with the constraints that must be fulfilled in order to run with this ``qadapter``.

``hardware``
    Dictionary with information on the hardware available on this particular queue.
    Used by Abinit ``autoparal`` to optimize parallel execution.

The ``qadapter`` is therefore responsible for all interactions with a specific 
queue management system (shell, Slurm, PBS, etc), including handling all details 
of queue script format as well as queue submission and management.

.. NOTE::

    Multiple ``qadapters`` are useful if you are running on a cluster with different queues 
    but we post-pone the discussion of this rather technical point.
    For the time being, we use a ``manager.yml`` with a single adapter. 

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

The ``job`` section is the most critical one, in particular the ``pre_run`` option
that will be executed by the shell script before invoking Abinit. 
In this case Abinit is not installed by default (the executable is not already in the path).
The directory where the Abinit executables are located hence have to be prepended to the original ``$PATH`` variable.
Change ``pre_run`` according to your Abinit installation and make sure that ``mpirun`` is also in ``$PATH``.
If you don't use a parallel version of Abinit, just set ``mpi_runner: null`` 
(``null`` is the YAML_ version of the Python ``None``). 
Note this approach also allows you to safely use multiple versions.

Copy this example and change the entries in the ``hardware`` and the ``limits`` section according to
your machine, in particular make sure that ``max_cores`` is not greater than the number of physical cores
available on your personal computer.
Save the file in the current working directory and run the :ref:`abicheck.py` script provided by AbiPy.
If everything is configured properly, you should see something like this in the terminal.

.. command-output:: abicheck.py --no-colors

This message tells us that everything is in place and we can finally run our first calculation.

.. note:

    This laptop has 1 socket with 2 CPUs and 4 Gb of memory in total, hence I don't want to run
    Abinit tasks with more than 2 CPUs. This is the reason why ``max_cores`` is set to 2.
    The ``timelimit`` option is not used when you are using ``qname: shell``, but it becomes 
    important when you submit jobs on a cluster because this value is used to generate the submission script
    and Abinit will use this value to exit from iterative algorithms e.g. the SCF cycle before the timeline 
    and produce files from which it can then restart.

The directory ``~abipy/data/runs`` contains python scripts to generate workflows for typical ab-initio calculations.
Here we focus on the configuration of the manager and the execution of the flow so we don't discuss how to 
generate input files and create Flow objects in python.
This topic is covered in more detail in our collection of `jupyter notebooks
<http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_

Let's start from the simplest example i.e. the ``run_si_ebands.py`` script that generates 
a flow to compute the band structure of silicon at the Kohn-Sham level 
(GS calculation to get the density followed by a NSCF run along a k-path in the first Brillouin zone).

Cd to ``~abipy/data/runs`` and execute ``run_si_ebands.py`` to generate the flow::

    cd ~abipy/data/runs
    ./run_si_ebands.py

At this point, you should have a directory ``flow_si_ebands`` with the following structure:

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

``w0/`` is the directory containing the input files of the first workflow (well, we have only one workflow in our example).
``w0/t0/`` and ``w0/t1/`` contain the input files need to run the SCF and the NSC run, respectively.

You might have noticed that each task directory (``w0/t0``, ``w0/t1``) presents the same structure:
    
   * ``run.abi``: Abinit input file.
   * ``run.files``: Abinit files file.
   * ``job.sh``: Submission/shell script.
   * ``outdata``: Directory with output data files.
   * ``indata``: Directory with input data files.
   * ``tmpdata``: Directory with temporary files.

.. DANGER::

   ``__AbinitFlow__.pickle`` is the pickle file used to save the status of the `Flow`. Don't touch it! 

The ``job.sh`` script has been generated by the ``TaskManager`` using the information provided by ``manager.yml``.
In this case it is a simple shell script that executes the code directly as we are using ``qtype: shell``. 
The script will get more complicated when we start to submit jobs on a cluster with a resource manager.

We usually interact with the AbiPy flow via the :ref:`abirun.py` script whose syntax is::

     abirun.py FLOWDIR command [options]

where ``FLOWDIR`` is the directory containing the flow and ``command`` defines the action to perform 
(use ``abirun.py --help`` to get the list of possible commands).

``abirun.py`` reconstructs the python Flow from the pickle file ``__AbinitFlow__.pickle`` located in ``FLOWDIR``
and invokes the methods of the object depending on the options passed via the command line.

Use::

    abirun.py flow_si_ebands status

to get a summary with the status of the different tasks and::

    abirun.py flow_si_ebands deps

to print the dependencies of the tasks in textual format.

.. code-block:: console

    <ScfTask, node_id=75244, workdir=flow_si_ebands/w0/t0>

    <NscfTask, node_id=75245, workdir=flow_si_ebands/w0/t1>
      +--<ScfTask, node_id=75244, workdir=flow_si_ebands/w0/t0>

.. TIP:: 

    Alternatively one can use ``abirun.py flow_si_ebands networkx``
    to visualize the connections with the networkx_ package.

In this case, we have a flow with one work (``w0``) that contains two tasks. 
The second task (``w0/t1``)  depends on first one that is a ``ScfTask``, 
more specifically ``w0/t1`` depends on the density file produced by ``w0/t0``.
This means that ``w0/t1`` cannot be executed/submitted until we have completed the first task. 
AbiPy is aware of this dependency and will use this information to manage the submission/execution
of our flow.

There are two commands that can be used to launch tasks: ``single`` and ``rapid``.
The ``single`` command executes the first task in the flow that is in the ``READY`` state that is a task
whose dependencies have been fulfilled. 
``rapid``, on the other hand, submits **all tasks** of the flow that are in the ``READY`` state.
Let's try to run the flow with the ``rapid`` command...

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


What's happening here?
The ``rapid`` command tried to execute all tasks that are ``READY`` but since the second task depends 
on the first one only the first task gets submitted.
Note that the SCF task (``w0_t0``) has been submitted with 2 MPI processes. 
Before submitting the task, indeed, AbiPy
invokes Abinit to get all the possible parallel configurations compatible within the limits 
specified by the user (e.g. ``max_cores``), select an "optimal" configuration according 
to some policy and then submit the task with the optimized parameters.
At this point, there's no other task that can be executed, the script exits
and we have to wait for the SCF task before running the second part of the flow.

At each iteration, :ref:`abirun.py` prints a table with the status of the different tasks.
The meaning of the columns is as follows:

``Queue`` 
    String in the form ``JobID @ QueueName`` where JobID is the process identifier if we are in the shell
    or the job ID assigned by the resource manager (e.g. slurm) if we are submitting to a queue.
``MPI`` 
    Number of MPI processes used. This value is obtained automatically by calling Abinit in ``autoparal mode``, 
    cannot exceed ``max_ncpus``.
``OMP`` 
    Number of OpenMP threads.
``Gb`` 
    Memory requested in Gb. Meaningless when ``qtype: shell``.
``Warn`` 
    Number of warning messages found in the log file.
``Com`` 
    Number of comments found in the log file.
``Sub``  
    Number of submissions. It can be > 1 if AbiPy encounters a problem and resubmit the task 
    with different parameters without performing any operation that can change the physics of the system).
``Rest``
    Number of restarts. AbiPy can restart the job if convergence has not been reached.
``Corr``
    Number of corrections performed by AbiPy to fix runtime errors. 
    These operations can change the physics of the system.
``Time``
    Time spent in the queue (if string ends with Q) or running time (if string ends with R).
``Node_ID``
    Node identifier used by AbiPy to identify each node of the flow.

.. NOTE:: 
     When the submission is done through the shell there's almost no difference between 
     job submission and job execution. The scenario is completely different if you are submitting 
     jobs to a resource manager because the task will get a priority value and will enter the queue.

If you execute ``status`` again, you should see that the first task is completed.
We can thus run ``rapid`` again to launch the |NscfTask|.
The second task won't take long and if you issue ``status`` again, you should see that the entire flow
completed successfully.

To understand what happened in more detail, use the ``history`` command to get 
the list of operations performed by AbiPy on each task.

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


A closer inspection of the logs reveal that before submitting the first task, python has executed
Abinit in ``autoparal`` mode to get the list of possible parallel configuration and the calculation is then submitted.
At this point, AbiPy starts to look at the output files produced by the task to understand  what's happening.
When the first task completes, the status of the second task is automatically changed to ``READY``, 
the ``irdden`` input variable is added to the input file of the second task and a symbolic link to
the ``DEN`` file produced by ``w0/t0`` is created in the ``indata`` directory of ``w0/t1``.
Another auto-parallel run is executed for the NSCF calculation and the second task is finally submitted.

The command line interface is very flexible and sometimes it's the only tool available.
However, there are cases in which we would like to have a global view of what's happening.
The command::

    $ abirun.py flow_si_ebands notebook

generates a jupyter_ notebook with pre-defined python code that can be executed 
to get a graphical representation of the status of our flow inside a web browser
(requires jupyter_, nbformat_ and, obviously, a web browser).

Expert users may want to use::

    $ abirun.py flow_si_ebands ipython

to open the flow in the ipython_ shell to have direct access to the API provided by the flow.

Once ``manager.yml`` is properly configured, it is possible 
to use the AbiPy objects to invoke Abinit and perform useful operations.
For example, one can use the |AbinitInput| object to get the list of k-points in the IBZ,
the list of independent DFPT perturbations, the possible parallel configurations reported by ``autoparal`` etc.

This programmatic interface can be used in scripts to facilitate the creation of input files and workflows.
For example, one can call Abinit to get the list of perturbations for each q-point in the IBZ and then
generate automatically all the input files for DFPT calculations (actually this is the approach used to
generated DFPT workflows in the AbiPy factory functions).

Note that ``manager.yml`` is also used to invoke other executables (``anaddb``, ``optic``, ``mrgddb``, etcetera)
thus creating some sort of interface between the python language and the Fortran executables.
Thanks to this interface, one can perform relatively simple ab-initio calculations directly in AbiPy.
For instance one can open a ``DDB`` file in a jupyter notebook, call ``anaddb`` to compute 
the phonon frequencies and plot the DOS and the phonon band structure with matplotlib_.

.. TIP::

        abirun.py . doc_manager

    gives the full documentation for the different entries of ``manager.yml``.

.. command-output:: abirun.py . doc_manager

.. _scheduler:

------------------------------
How to configure the scheduler
------------------------------

In the previous example, we ran a simple band structure calculation for silicon in a few seconds 
on a laptop but one might have more complicated flows requiring hours or even days to complete.
For such cases, the ``single`` and ``rapid`` commands are not handy because we are supposed 
to monitor the evolution of the flow and re-run ``abirun.py`` when a new task is ``READY``.
In these cases, it is much easier to delegate all the repetitive work to a ``python scheduler``,
a process that runs in the background, submits tasks automatically and performs the actions
required to complete the flow.

The parameters for the scheduler are declared in the YAML_ file ``scheduler.yml``.
Also in this case, AbiPy will look first in the working directory and then inside ``$HOME/.abinit/abipy``.
Create a ``scheduler.yml`` in the working directory by copying the example below:

.. code-block:: yaml

    seconds: 5   # number of seconds to wait.
    #minutes: 0  # number of minutes to wait.
    #hours: 0    # number of hours to wait.

This file tells the scheduler to wake up every 5 seconds, inspect the status of the tasks
in the flow and perform the actions required for reach completion

.. IMPORTANT::

    Remember to set the time interval to a reasonable value.
    A small value leads to an increase of the submission rate but it also increases the CPU load 
    and the pressure on the hardware and on the resource manager.
    A too large time interval can have a detrimental effect on the throughput, especially 
    if you are submitting many small jobs.

At this point, we are ready to run our first calculation with the scheduler.
To make things more interesting, we execute a slightly more complicated flow that computes
the G0W0 corrections to the direct band gap of silicon at the Gamma point.
The flow consists of the following six tasks:

- 0: Ground state calculation to get the density.
- 1: NSCF calculation with several empty states. 
- 2: Calculation of the screening using the WFK produced by task 2.
- 3-4-5: Evaluation of the Self-Energy matrix elements with different values of nband 
  using the WFK produced by task 2 and the SCR file produced by task 3

Generate the flow with::

    ./run_si_g0w0.py

and let the scheduler manage the submission with::

     abirun.py flow_si_g0w0 scheduler

You should see the following output on the terminal

.. code-block:: console

    abirun.py flow_si_ebands scheduler

    Abipy Scheduler:
    PyFlowScheduler, Pid: 72038
    Scheduler options: {'seconds': 10, 'hours': 0, 'weeks': 0, 'minutes': 0, 'days': 0}

``Pid`` is the process identifier associated the scheduler (also saved in in the ``_PyFlowScheduler.pid`` file).

.. IMPORTANT:: 

    A ``_PyFlowScheduler.pid`` file in ``FLOWDIR`` means that there's a scheduler running the flow.
    Note that there must be only one scheduler associated to a given flow.

As you can easily understand the scheduler brings additional power to the AbiPy flow because
it is possible to automate complicated ab-initio workflows with little effort: write
a script that implements the flow in python and save it to disk, run it with 
``abirun.py FLOWDIR scheduler`` and finally use the AbiPy/Pymatgen tools to analyze the final results.
Even complicated convergence studies for G0W0 calculations can be implemented along these lines
as shown by this `video <https://youtu.be/M9C6iqJsvJI>`_.
The only problem is that at a certain point our flow will become too big or too computational expensive
that cannot be executed on a personal computer anymore and we have to move to a supercomputing center.
The next section discusses how to configure AbiPy to run on a cluster with a queue management system.

.. TIP:: 
    
    Use ``abirun.py . doc_scheduler`` to get the full list of options supported by the scheduler.

.. command-output:: abirun.py doc_scheduler

.. _abipy-on-cluster:

------------------------------
Configuring AbiPy on a cluster
------------------------------

In this section we discuss how to configure the manager to run flows on a cluster.
The configuration depends on specific queue management system (Slurm, PBS, etc) hence
we assume that you are already familiar with job submissions and you know the options 
that mush be specified in the submission script in order to have your job accepted 
and executed by the management system (username, name of the queue, memory ...)

Let's assume that our computing center uses slurm_ and our jobs must be submitted to the ``default_queue`` partition.
In the best case, the system administrator of our cluster (or you create one yourself) already provides 
an ``Abinit module`` that can be loaded directly with ``module load`` before invoking the code.
To make things a little bit more difficult, however, we assume the we had to compile our own version of Abinit
inside the build directory ``${HOME}/git_repos/abinit/build_impi`` using the following two modules
already installed by the system administrator::

    compiler/intel/composerxe/2013_sp1.1.106
    intelmpi

In this case, we have to be careful with the configuration of our environment because the Slurm submission
script should load the modules and modify our ``$PATH`` so that our version of Abinit can be found.
A ``manager.yml`` with a single ``qadapter`` looks like:

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

    abirun.py FLOWDIR doc_manager script

    prints to screen the submission script that will be generated by AbiPy at runtime.

Let's discuss the different options in more detail. Let's start from the ``queue`` section:

``qtype`` 
    String specifying the resource manager. This option tells AbiPy which ``qadapter`` to use to generate the submission
    script, submit them, kill jobs in the queue and how to interpret the other options passed by the user. 

``qname``
    Name of the submission queue (string, MANDATORY)

``qparams`` 
    Dictionary with the parameters passed to the resource manager. 
    We use the *normalized* version of the options i.e. dashes in the official name of the parameter
    are replaced by underscores e.g. ``--mail-type`` becomes ``mail_type``.
    For the list of supported options use the ``doc_manager`` command.
    Use ``qverbatim`` to pass additional options that are not included in the template.

Note that we are not specifying the number of cores in ``qparams`` because AbiPy will find an appropriate value
at run-time.

The ``job`` section is the most critical one because it defines how to configure the environment
before executing the application and how to run the code.
The ``modules`` entry specifies the list of modules to load, ``shell_env`` allows us to modify the 
``$PATH`` environment variables so that the OS can find our Abinit executable.

.. IMPORTANT::

    Various resource managers will first execute your ``.bashrc`` before starting to load the new modules.

We also increase the size of the stack with ``ulimit`` before running the code and we run Abinit 
with the ``mpirun`` provided by the modules.

The ``limits`` section defines the constraints that must be fulfilled in order to run on this queue
while ``hardware`` is a dictionary with info on the hardware available on this queue.
Every job will have a ``timelimit`` of 20 minutes, cannot use more that ``max_cores`` cores,
and the first job submission will request 1 Gb of memory.
Note that the actual number of cores will be determined at runtime by calling Abinit in ``autoparal`` mode
to get all parallel configurations up to ``max_cores``.
If the job is killed due to insufficient memory, AbiPy will resubmit the task with increased resources
and it will stop when it reaches the maximum amount given by ``mem_per_node``.

Note that there are more advances options supported by ``limits`` and other options
will be added as time goes by.

The get the complete list of options supported by the Slurm ``qadapter`` use:

.. command-output:: abirun.py . doc_manager slurm

.. IMPORTANT::

    If you need to cancel all tasks that have been submitted to the resource manager, use::

        abirun.py FLOWDIR cancel

    Note that the script will ask for confirmation before killing all the jobs belonging to the flow.

Once you have a ``manager.yml`` properly configured for your cluster, you can start
to use the scheduler to automate job submission.
Very likely your flows will require hours or even days to complete and, in principle, 
you should maintain an active connection to the machine in order to keep your scheduler alive
(if your session expires, all subprocesses launched within your terminal, 
including the python scheduler, will be automatically killed).
Fortunately there is a standard Unix tool called ``nohup`` that comes to our rescue.

For long-running jobs, we strongly suggest to start the scheduler with::

     nohup abirun.py FLOWDIR scheduler > sched.stdout 2> sched.stderr &

This command executes the scheduler in background and redirects the ``stdout`` and ``stderr``
to ``sched.log`` and ``sched.err``, respectively.
The process identifier of the scheduler is saved in the ``_PyFlowScheduler.pid`` file inside ``FLOWDIR``
and this file is removed automatically when the scheduler completes its execution.
Thanks to the ``nohup`` command, we can close our session, let the scheduler work overnight
and reconnect the day after to collect our data.

.. IMPORTANT:: 

    Use ``abirun.py FLOWDIR cancel`` to cancel the jobs of a flow that is being executed by
    a scheduler. AbiPy will detect that there is a scheduler already attached to the flow 
    and will cancel the jobs of the flow and kill the scheduler as well.


.. _inspecting-the-flow:

-------------------
Inspecting the Flow
-------------------

:ref:`abirun.py` also provides tools to analyze the results of the flow at runtime.
The simplest command is::

    abirun.py FLOWDIR tail

that is the analogous of Unix tail but a little bit more smarter in the 
sense that ``abirun.py`` will only print to screen the final part of the output files
of the tasks that are ``RUNNING``.

If you have matplotlib_ installed, you may want to use::

    $ abirun.py FLOWDIR inspect

Several AbiPy tasks, indeed, provide an ``inspect`` method producing matplotlib figures
with data extracted from the output files. 
For example, a ``GsTask`` prints the evolution of the ground-state SCF cycle.
The inspect command of :ref:`abirun.py` just loops over the tasks of the flow and 
calls the ``inspect`` method on each of them.

The command::

    abirun.py FLOWDIR inputs

prints the input files of the different tasks (can use ``--nids`` to select a subset of
tasks or, alternatively, replace ``FLOWDIR`` with the ``FLOWDIR/w0/t0`` syntax)

The command::

    abirun.py FLOWDIR listext EXTENSION

prints a table with the nodes of the flow who have produced an Abinit output file with the given 
extension. Use e.g.::

    abirun.py FLOWDIR listext GSR.nc

to show the nodes of the flow who have produced a GSR.nc_ file.

The command::

    abirun.py FLOWDIR notebook

generates a jupyter_ notebook with pre-defined python code that can be executed 
to get a graphical representation of the status of the flow inside a web browser
(requires jupyter_, nbformat_ and, obviously, a web browser).

Expert users may want to use::

    abirun.py FLOWDIR ipython

to open the flow in the ipython_ shell to have direct access to the API provided by the flow.


.. _event-handlers:

--------------
Event handlers
--------------

An event handler is an action that is executed in response of a particular event.
The AbiPy tasks are equipped with built-in events handlers that are be executed 
to fix typical Abinit runtime errors.

To list the event handlers installed in a given flow use::

    abirun.py FLOWDIR handlers

The ``--verbose`` option produces a more detailed description of the action performed
by the event handlers.

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

    New error handlers will be added in the new versions of Abipy/Abinit.
    Please, let us know if you need handlers for errors commonly occuring in your calculations. 

.. _flow-troubeshooting:

---------------
Troubleshooting
---------------

There are two :ref:`abirun.py` commands that are very useful especially if something goes wrong: ``events`` and ``debug``.

To print the Abinit events (Warnings, Errors, Comments) found in the log files of the different tasks use::

    abirun.py FLOWDIR events

To analyze error files and log files for possible error messages, use::

    abirun.py FLOWDIR debug

By default, these commands will analyze the entire flow so the output on the terminal can be very verbose.
If you are interested in a particular task e.g. ``w0/t1`` use the syntax::

    abirun.py FLOWDIR/w0/t1 events

to select all the tasks in a work directory e.g. ``w0`` use::

    abirun.py FLOWDIR/w0 events

to select an arbitrary subset of nodes of the flow use the syntax::

    abirun.py FLOWDIR events -nids=12,13,16

where ``nids`` is a list of AbiPy node identifiers.

.. TIP:: 

    ``abirun.py events --help`` is your best friend

.. command-output:: abirun.py events --help 

To get information on the Abinit executable called by AbiPy, use::

    abirun.py abibuild

or the verbose variant::

    abirun.py abibuild --verbose 

TODO: How to reset tasks 

.. _task_policy:

----------
TaskPolicy
----------

At this point, you may wonder why we need to specify all these parameters in the configuration file.
The reason is that, before submitting a job to a resource manager, AbiPy will use the autoparal 
feature of ABINIT to get all the possible parallel configurations with ``ncpus <= max_cores``. 
On the basis of these results, AbiPy selects the "optimal" one, and changes the ABINIT input file 
and the submission script accordingly .
(this is a very useful feature, especially for calculations done with ``paral_kgb=1`` that require 
the specification of ``npkpt``, ``npfft``, ``npband``, etc).
If more than one ``QueueAdapter`` is specified, AbiPy will first compute all the possible 
configuration and then select the "optimal" ``QueueAdapter`` according to some kind of policy

In some cases, you may want to enforce some constraint on the "optimal" configuration. 
For example, you may want to select only those configurations whose parallel efficiency is greater than 0.7 
and whose number of MPI nodes is divisible by 4. 
One can easily enforce this constraint via the ``condition`` dictionary whose syntax is similar to 
the one used in mongodb_.

.. code-block:: yaml

    policy:
        autoparal: 1
        max_ncpus: 10
        condition: {$and: [ {"efficiency": {$gt: 0.7}}, {"tot_ncpus": {$divisible: 4}} ]}

The parallel efficiency is defined as $\epsilon = \dfrac{T_1}{T_N * N}$ where $N$ is the number 
of MPI processes and $T_j$ is the wall time needed to complete the calculation with $j$ MPI processes. 
For a perfect scaling implementation $\epsilon$ is equal to one.
The parallel speedup with N processors is given by $S = T_N / T_1$.
Note that ``autoparal = 1`` will automatically change your ``job.sh`` script as well as the input file 
so that we can run the job in parallel with the optimal configuration required by the user. 
For example, you can use ``paral_kgb = 1`` in GS calculations and AbiPy will automatically set the values 
of ``npband``, ``npfft``, ``npkpt`` ... for you! 
Note that if no configuration fulfills the given condition, AbiPy will use the optimal configuration 
that leads to the highest parallel speedup (not necessarily the most efficient one).

``policy`` 
    This section governs the automatic parallelization of the run: in this case AbiPy will use 
    the ``autoparal`` capabilities of Abinit to determine an optimal configuration with 
    **maximum** ``max_ncpus`` MPI nodes. Setting ``autoparal`` to 0 disable the automatic parallelization. 
    Other values of autoparal are not supported.

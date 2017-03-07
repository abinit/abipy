.. _running_calculations:

====================
Running calculations
====================

Besides post-processing tools and a programmatic interface to generate input files,
AbiPy also provides a pythonic interface to execute small Abinit tasks or submit calculations on supercomputing centers.
This section discusses how to create the configuration files required to interface AbiPy with Abinit.

We assume that Abinit is already available on your machine and that you know how to configure
your environment: set the `$PATH` and `$LD_LIBRARY_PATH`(`DYLD_LIBRARY_PATH` on Mac) 
environment variables, load modules with `module load`, etc.)
so that the operating system can load and execute code. 

.. IMPORTANT:: Please make sure that you can execute Abinit interactively with simple input files and 
               that the code works as expected before proceeding with the rest of the tutorial.


--------------------------------
How to configure the TaskManager
--------------------------------

The `TaskManager` is responsible for task submission 
(creation of the submission script, initialization of the environment) as well as for the 
optimization of the parameters used in parallel runs 
(number of MPI processes, number of OpenMP threads, automatic parallelization with ABINIT `autoparal` feature). 

AbiPy knows how to run/submit the code with the correct environment and the appropriate syntax
thanks to the options given in the `manager.yml` configuration file.
The configuration file is written in `YAML <https://en.wikipedia.org/wiki/YAML>`_,
a human-readable data serialization language commonly used for configuration files
(a good introduction to the YAML syntax can be found `here <http://yaml.org/spec/1.1/#id857168>`_.
See also this `reference card <http://www.yaml.org/refcard.html>`_)

By default, AbiPy looks for a `manager.yml` file in the current working directory i.e.
the directory in which you execute your script and then inside `$HOME/.abinit/abipy`.
If no `manager.yml` is found, the code aborts immediately.
Configuration files for typical cases are available in `~abipy/data/managers`.
We first discuss how to configure Abipy on a personal computer and then we look at the more
complicated cases in which the calculation must be submitted to a queue.

-----------------------------------
TaskManager for a personal computer
-----------------------------------

Let's start from the simplest case i.e. a personal computer in which we can execute Abinit directly from the shell.
In this case, the configuration is relatively easy because we can execute the code 
directly without having to submit a script to the resource manager to allocate resources (memory, MPI ranks, threads ...)
`manager.yml` contains a list of `QueueAdapters` objects. 
Each `QueueAdapter` is responsible for all interactions with a specific queue management system (Slurm, PBS, bash, etc).
This includes handling all details of queue script format as well as queue submission and management.

The configuration file I use on my laptop is

.. code-block:: yaml

    qadapters: # List of `qadapters` objects  (just one in this simplified example)

	- queue:
	    qtype: shell        # "Submit" jobs via the shell.
	    qname: localhost    # "Submit" to the `localhost` queue. (fake queue in this case)

	  priority: 1

	  job:
	    pre_run: "export PATH=$HOME/git_repos/abinit_eph/build_gcc/src/98_main:$PATH"
	    mpi_runner: mpirun

	  limits: 
	    timelimit: 1:00:00   #  Time-limit for each task.
	    max_cores: 2         #  Max number of cores that can be used by a single task.

	  # Hardware specification, used by autoparal to optimize parallel execution.
	  hardware:  
	     num_nodes: 1
	     sockets_per_node: 1
	     cores_per_socket: 2
	     mem_per_node: 4 Gb

For the sake of brevity, we give a brief description of the meaning 
of the different sections without entering into detail:

    * queue: dictionary with the name of the queue and optional parameters 
      used to build/customize the header of the submission script

    * job: dictionary with the options used to prepare the environment before submitting the job

    * limits: dictionary with the constraints that must be fulfilled in order to run with this queueadapter

    * hardware: dictionary with information on the hardware available on this particular queue.


Use::

    abirun.py doc_manager

to get the complete documentation of the `manager.yml` file

In this simple case, we have one `QueueAdapter` named `gmac` that will submit `Tasks`
in a shell subprocess (`qtype: shell`) with mpirun. 
`env.sh` is the bash script I use to set the value of the environment variables 
(e.g. `PATH` and `LD_LIBRARY_PATH`) before running ABINIT.

Note that my laptop has 1 socket with 2 CPUs and 4 Gb of memory in total, hence I don't want
to run ABINIT tasks with more than 2 CPUs. This is the reason why `max_cores` is set to 2.
`Timelimit` is not used when you are using `qname=shell`, but it is very important when you 
are submitting jobs on a cluster because this value is used to generate the submission script.

At this point, you may wonder why we need to specify all these parameters in the configuration file.
The reason is that, before submitting a job to a resource manager, `AbiPy` will use the autoparal 
feature of ABINIT to get all the possible parallel configurations with `ncpus <= max_cores`. 
On the basis of these results, `AbiPy` selects the "optimal" one, and changes the ABINIT input file 
and the submission script accordingly .
(this is a very useful feature, especially for calculations done with `paral_kgb=1` that require 
the specification of `npkpt`, `npfft`, `npband`, etc).
If more than one `QueueAdapter` is specified, `AbiPy` will first compute all the possible 
configuration and then select the "optimal" `QueueAdapter` according to some kind of policy

Copy this example, change the entries in the `hardware` and the `limits` section according to
your machine, change `pre_run` so that the Abinit executables can be found in `$PATH`.
Save the file in the current working directory and run `abicheck.py`.
If everything is configured properly, you should see something like this in the terminal.

.. command-output:: abicheck.py --no-colors

This message tells us that everything is in place and we can finally run our first calculation with Abipy.
The directory `abipy/data/runs` contains python scripts to generate workflows for typical ab-initio calculations.
Here we focus on the configuration of the manager and the execution of the flow so we don't to explain how to 
generate input files and create Flow objects in python.

Let's start from the simplest example i.e. the `run_si_ebands.py` script that generates 
a flow to compute the band structure of silicon at the Kohn-Sham level 
(GS calculation to get the density followed by a NSCF run along a k-path in the first Brillouin zone).
Cd to ~abipy/data/runs and execute `run_si_ebands.py` to generate the flow::

    cd ~abipy/data/runs
    ./run_si_ebands.py

At this point, you should have a directory named `flow_si_ebands` with the following structure:

.. code-block:: shell

    $ tree flow_si_ebands/

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

`w0` is the directory containing the input files of the first workflow (well, we have only one workflow in our example).
`t0` and `t1` contain the input files need to run the SCF and the NSC run, respectively.

You might have noticed that each `Task` directory (w0/t0, w0/t1) presents the same structure:
    
   * run.abi: ABINIT input file
   * run.files: ABINIT files file
   * job.sh: Submission/shell script
   * outdata: Directory with output data files
   * indata: Directory with input data files 
   * tmpdata: Directory with temporary files

.. DANGER::
   `__AbinitFlow__.pickle` is the pickle file used to save the status of the `Flow`. Don't touch it! 

The `job.sh` has been generated using the information provided by `manager.yml`. 
In this case it's a simple shell script that executes the code but this is normal because we are using `qtype: shell`. 
The script will be more complicated when we start to submit jobs on a cluster with a resource manager.

We usually interact with the Abipy flow via the `abirun.py` script.
The script uses the syntax::

     `abirun.py FLOWDIR command [options]`

where `FLOWDIR` is the directory containing the flow and `command` defines the action to perform 
(use `--help` to get the list of possible commands).
`abirun.py` reconstruct the python Flow from the pickle file `__AbinitFlow__.pickle` located in `FLOWDIR`
and invokes the methods of the object depending on the options specified by the user on the command line.
Let's start to play with our flow.

Use::

    abirun.py flow_si_ebands status

to have a summary with the status of the different tasks and::

    abirun.py flow_si_ebands deps

to print the interconnection among the tasks in text format.

.. code-block:: shell

    <ScfTask, node_id=75244, workdir=flow_si_ebands/w0/t0>

    <NscfTask, node_id=75245, workdir=flow_si_ebands/w0/t1>
      +--<ScfTask, node_id=75244, workdir=flow_si_ebands/w0/t0>

.. TIP:: Alternatively one can use `abirun.py flow_si_ebands/ networkx`
	 to visualize the connections with the `networkx` package.

In this case, we have a flow with two tasks and the second task (w0/t1) 
depends on the ScfTask, more specifically on the density file produced by it.
This means that the second task cannot be executed/submitted until we have completed the first task. 
`abirun.py` knows the dependencies of our flow and will use this information to manage the submission/execution
of our tasks.

There are two commands that can be used to launch tasks: `single` and `rapid`.
The `single` command execute the first `Task` in the flow that is in the `READY` state that is a task
whose dependencies have been fulfilled while `rapid` submits all task of the flow that are in the `READY` state.
Let's try to run the flow with the `rapid` command and see what happens.

.. code-block:: shell

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
The `rapid` command tried to execute all tasks that are `READY` but since the second task depends on the first
one only the first task gets submitted.
Note that the SCF task (w0_t0) has been submitted with 2 MPI processors. 
Before submitting the task, indeed, AbiPy
invokes Abinit to get all the possible parallel configurations compatible with the constrains specified by the user,
select the "optimal" configuration according to some policy and submit the task with the optimized parameters.
At this point, there's no other task that can be executed, the script exits
and we have to wait for the SCF task before running the second part of the flow.

At each iteration, `abirun.py` prints a table with the status of the different tasks.
The meaning of the columns is as follows:

    * Queue: JobID @ QueueName (JobID == Process identifier if shell, job ID if we are submitting to QueueName)
    * MPI: Number of MPI processes used (computed automatically with autoparal, cannot exceed max_ncpus)
    * OMP: Number of OpenMP threads.
    * Gb: Memory requested in Gb (meaningless in this case because we're using the shell).
    * Warn: Number of warning messages found in the log file.
    * Com: Number of comments found in the log file.
    * Sub: Number of submissions (can be > 1 if Abipy encounters a problem and resubmit the task with different parameters
	without performing any operation that can change the physics of the system).
    * Rest: Number of restarts (Abipy can restart the job if convergence has not been reached)
    * Corr: Number of corrections performed. These operations can change the physics of the system.
    * Time: Time spent in the Queue (if ends with Q) or running time (if ends with R).
    * Node_ID: Node identifier used by Abipy to identify each node of the flow.

.. NOTE:: When the submission is done through the shell there's almost no difference between 
    job submission and job execution. The scenario is completely different if you are submitting 
    jobs to a resource manager because the task will get a priority value and will enter the queue.

If you execute `status` again, you should see that the first task is completed.
We can thus run `rapidfire` again to launch the `NscfTask`.
The second task won't take long and if you issue `status` again, you should see that the entire flow
completed successfully.

To understand what happened in more detail, use the `history` command to get the list of operations
performed by AbiPy on each task.

.. code-block:: shell

    $ abirun.py flow_si_ebands history

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


A closer inspection of the logs reveal that before submitting the first task, `abirun.py` has executed
Abinit in `autoparal` mode to get the list of possible parallel configuration and the calculation is then submitted.
At this point, `abirun.py` starts to look at the output files produced by the task to understand 
what's happening.
When the first task reaches completion, the second task is automatically changed to READY, 
the `irdden` input variable is added to the input file of the second task and a symbolic link to
the DEN file produced by the first task is created in the `indata` directory of the second task.
Another `autoparallel run` is now executed and the second task is finally submitted.

The command line interface is very flexible and sometimes it's the only tool available.
However, there are cases in which we would like to have a global view of what's happening 
The command::

    abirun.py flow_si_ebands notebook

generates a jupyter notebook with pre-defined calls that can be executed 
in order to get a graphical representation of the status of our flow inside a web browser
(requires `jupyter`, `nbformat` and, obviously, a web browser).
Expert users may want to use::

    abirun.py flow_si_ebands ipython

to open the `flow` in the `ipython` terminal so to have direct access to the API provided by the object.

------------------------------
How to configure the scheduler
------------------------------

In the previous example, we ran a simple band structure calculation for silicon in a few seconds 
on a laptop but one might have more complicated flows requiring hours or even days to complete.
For such cases, the `single` and `rapid` commands are not handy because we are supposed 
to monitor the evolution of the flow and re-run `abirun.py` when a new task is `READY`.
In these cases, it is much easier to delegate all the repetitive work to a python scheduler,
a sort of job that runs in the background and submits tasks automatically and perform the actions
required to complete the flow.

The parameters for the scheduler are declared in the YAML file `scheduler.yml`.
Also in this case, AbiPy will look first in the working directory and then inside `$HOME/.abinit/abipy`.
Crate a `scheduler.yml` in the working directory by copying the example below:

.. code-block:: yaml

    seconds: 10  # number of seconds to wait.
    #minutes: 0  # number of minutes to wait.
    #hours: 0    # number of hours to wait.
    #days: 0     # number of days to wait.

This file tells the scheduler to wake up every 10 seconds, inspect the status of the tasks
in the flow and perform the actions required for reach completion



.. IMPORTANT::

    Remember to set the time interval of the scheduler to a reasonable value.
    A small value leads to an increase of the submission rate but it also increases the CPU load 
    and the pressure on the hardware and on the resource manager.
    A too large time interval can have a detrimental effect on the throughput, especially 
    if you are submitting many small jobs.

At this point, we are ready to run our first calculation with the scheduler.
To make things more interesting, we execute a slightly more complicated flow that computes
the G0W0 corrections to the direct band gap of silicon at the Gamma point.
The flow consists of the following six tasks:

    1: ground state calculation to get the density
    2: NSCF calculation with several empty states. 
    3: calculation of the screening using the WFK produced by task 2
    4-5-6: Evaluation of the Self-Energy matrix elements with different values of nband 
     using the WFK produced by task 2 and the SCR file produced by task 3

Generate the flow with::

    ./run_si_g0w0.py

and let the scheduler manage the task submission with::

     abirun.py flow_si_g0w0 scheduler

You should see the following output on the terminal

.. code-block:: shell

    abirun.py flow_si_ebands scheduler

    Abipy Scheduler:
    PyFlowScheduler, Pid: 72038
    Scheduler options: {'seconds': 10, 'hours': 0, 'weeks': 0, 'minutes': 0, 'days': 0}

PID is the process identifier of the scheduler (also reported in the ... file)
We will see that the scheduler PID is extremely important when we start to run large flows on clusters. 

.. IMPORTANT:: Note that there must be only one scheduler associated to a given flow.

.. TIP:: Use `abirun.py . doc_scheduler` to get the full list of options supported by the scheduler.

.. command-output:: abirun.py doc_scheduler

------------------------------
Configuring AbiPy on a cluster
------------------------------

In this section we discuss how to configure the manager to run flows on a cluster.
The configuration depends on specific queue management system (Slurm, PBS, etc) so
we assume that you are already familiar with job submissions and you know the options 
that mush be specified in the job script in order to have your submission accepted 
by the management system (username, name of the queue ...)

Let's assume that your computing center uses Slurm and your jobs must be submitted to the `Oban` partition 
A `manager.yml` with a single `qadapter` will look like:

.. code-block:: yaml

    # Resource manager e.g slurm, pbs, shell
    qtype: slurm

    # Options passed to the resource manager 
    # (the syntax depends on qtype, consult the manual of your resource manager)
    qparams: 
      ntasks: 2
      time: 0:20:00
      partition: Oban
    
    # List of modules to import before running the calculation
    modules: 
	- intel/compilerpro/13.0.1.117
	- fftw3/intel/3.3

    mpi_runner: /path/to/mpirun
    
    # Shell environment
    shell_env: 
	 PATH: /home/user/local/bin/:$PATH
	 LD_LIBRARY_PATH: /home/user/local/lib:$LD_LIBRARY_PATH

    # Options for the automatic parallelization (Abinit autoparal feature)
    policy: 
	autoparal: 1
	max_ncpus: 2


Description:

`qtype` specifies the queue resource manager. This option tells Abipy how to generate submission
scripts, submit them, kill jobs in the queue and how to interpret the other options passed by the user. 

`qparams` is a dictionary with the parameters passed to the resource manager. 
 We use the *normalized* version of the options i.e dashes in the official name of the parameter 
 are replaced by underscores  (for the list of supported options see ...)

The complete list of options (`qparams`) supported by the `TaskManager` with Slurm  can be obtained with

.. command-output:: abirun.py . doc_manager slurm

`modules` is the list of modules to load, while `shell_env` allows the user 
to specify or to modify the values of the environment variables.

The `policy` section governs the automatic parallelization of the run: in this case abipy will use 
the `autoparal` features of abinit to determine an optimal configuration with **maximum** `max_ncpus` MPI nodes. 
Setting `autoparal` to 0 disables the automatic parallelization. **Other values of autoparal are not supported**.

In some cases, you may want to enforce some constraint on the "optimal" configuration. 
For example, you may want to select only those configurations whose parallel efficiency is greater than 0.7 
and whose number of MPI nodes is divisible by 4. 
One can easily enforce this constraint via the `condition` dictionary whose syntax is similar to the one used in `mongodb`

.. code-block:: yaml

    policy:
	autoparal: 1
	max_ncpus: 10
	condition: {$and: [ {"efficiency": {$gt: 0.7}}, {"tot_ncpus": {$divisible: 4}} ]}

The parallel efficiency is defined as $\epsilon = \dfrac{T_1}{T_N * N}$ where $N$ is the number 
of MPI processes and $T_j$ is the wall time needed to complete the calculation with $j$ MPI processes. 
For a perfect scaling implementation $\epsilon$ is equal to one.
The parallel speedup with N processors is given by $S = T_N / T_1$.
Note that `autoparal = 1` will automatically change your `job.sh` script as well as the input file 
so that we can run the job in parallel with the optimal configuration required by the user. 
For example, you can use `paral_kgb` in GS calculations and `abipy` will automatically set the values 
of `npband`, `npfft`, `npkpt` ... for you! 
Note that if no configuration fulfills the given condition, abipy will use the optimal configuration 
that leads to the highest parallel speedup (not necessarily the most efficient one).

Use::

    abirun.py FLOWDIR cancel

to cancel all tasks that have been submitted to the resource manager (the script asks for confirmation).
Abipy detects if there's a scheduler attached to the flow and it will also kill the scheduler

In the previous sections, we have discussed how to define, build and run a `Flow`, but there is a very 
important point that we haven't discussed yet.
It should be stressed, indeed, that `AbiPy` is only driving and monitoring the `Flow` while the actual calculation 
is delegated to ABINIT (a Fortran program that is usually executed in parallel on multiple CPUs that communicate 
via the network by means of the MPI protocol).
Besides CPUs and memory must be reserved in advance by sending a request to the resource manager 
installed on the clusters (SLURM, PBS, etc)

.. TIP:: nohup abirun.py FLOWDIR scheduler 2> sched.log

One can put this configuration file either in the configuration directory `$HOME/.abinit/abipy` 
or in the current working directory (the latter has precedence over the global configuration 
file located in `$HOME/.abinit/abipy`).

because it's possible to run the scheduler in the background with::

     nohup abirun.py FLOWDIR scheduler 2> sched.log

This shell command redirects the stdout/stderr of `abirun.py` to `sched.log` 
and kill the active session without killing the scheduler thanks to the `nohup` unix command.
In this case, the PID gives as a handle that can be used to check whether the scheduler
is still running or kill it when we login again.

---------------
Troubleshooting
---------------

There are two other `abirun` commands that are very handy, especially if something goes wrong:

Use::

    abirun.py FLOWDIR events

to print the events (Abinit Warnings/Errors/Comments) found in the log files and::

    abirun.py FLOWDIR debug

to analyze error files and log files for possible error messages.

To get information on the Abinit build, use

.. command-output:: abirun.py abibuild --verbose 

while::

    abirun.py flow_si_ebands handlers

show the so-called events handlers that have been installed in the flow 
(an event handler is an action that will be executed in response of a particular event

.. code-block:: shell

    $ abirun.py flow_si_ebands handlers --verbose

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

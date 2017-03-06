Introduction
============

Besides post-processing tools and a programmatic interface to generate input files,
AbiPy also provides a pythonic interface to execute Abinit or submit calculations on supercomputing centers.
This section discusses how to create the configuration files required to interface AbiPy with Abinit.

We assume that Abinit is already available on your machine and that you know how to setup
your environment (set the $PATH and $LD_LIBRARY_PATH environment variables, load modules with `module load`, etc.)
so that the operating system can load and execute code. 

.. IMPORTANT:: Please make sure that you can execute Abinit interactively with simple input files and 
                that the code works as expected before proceeding with the rest of the tutorial

How to configure the TaskManager and the Scheduler
==================================================

The `TaskManager` is responsible for the submission of the tasks 
(creation of the submission script, initialization of the shell environment) as well as for the 
optimization of the parameters used for parallel runs (number of MPI processes, number of OpenMP threads, 
automatic parallelization with ABINIT `autoparal` feature). 

AbiPy knows how to run/submit the code with the proper environment thanks to the options
specified in the `manager.yml` configuration file.

The configuration file for the `TaskManager` is written in `YAML <https://en.wikipedia.org/wiki/YAML>`_,
a human-readable data serialization language commonly used for configuration files
(a good introduction to the YAML syntax can be found `here <http://yaml.org/spec/1.1/#id857168>`_.
See also this `reference card <http://www.yaml.org/refcard.html>`_)

By default, AbiPy looks for a `manager.yml` file in the current working directory i.e.
the directory in which you execute your script and then inside `$HOME/.abinit/abipy`.
If no `manager.yml` is found, the code aborts immediately.
Configuration files for typical cases are available in ~abipy/data/managers
A typical example is reported below:

TaskManager for a personal computer
===================================

Let's start from the simplest case i.e. a personal computer in which we can 
execute Abinit directly within the shell.
In this case, the configuration is made simple by the fact that we don't have
to submit a script to allocate resources and we can execute the binary directly.

The configuration file I use on my laptop is

.. code-block:: yaml

    qadapters:
	# List of qadapters objects  (just one in this simplified example)
	- queue:
	    qtype: shell        # "Submit" jobs via the shell.
	    qname: localhost    # "Submit" to the `localhost` queue. (fake queue in this case)
	  priority: 1
	  job:
	    # source a script to setup the environment.
	    pre_run: "export PATH=$HOME/git_repos/abinit_eph/build_gcc/src/98_main:$PATH"
	    mpi_runner: mpirun
	  limits:
	    timelimit: 1:00:00
	    max_cores: 2
	  hardware:
	     num_nodes: 1
	     sockets_per_node: 1
	     cores_per_socket: 2
	     mem_per_node: 4 Gb

`AbiPy` gets all the information needed to submit the different `Tasks` from a 
configuration file, `manager.yml`, that is usually located in the directory `~/.abinit/abipy/`. 
`manager.yml` contains a list of `QueueAdapters` objects. 
Each `QueueAdapter` is responsible for all interactions with a specific queue management system (slurm, PBS, bash, etc).
This includes handling all details of queue script format as well as queue submission and management.

For the sake of brevity, we just try to give you a general overview of the meaning 
of the different sections without entering into detail.

    * queue: dictionary with the name of the queue and optional parameters 
      used to build/customize the header of the submission script

    * job: dictionary with the options used to prepare the enviroment before submitting the job

    * limits: dictionary with the constraints that must be fulfilled in order to run with this queueadapter

    * hardware: dictionary with information on the hardware available on this particular queue.

In this (simple) case, we have one `QueueAdapter` named `gmac` that will submit `Tasks`
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
your machine, change `pre_run` so that the abinit executables can be found in $PATH.
Save the file in the current working directory and run `abicheck.py`.
If everything is configured properly, you should see something like this in the terminal.

.. code-block:: shell

    $ abicheck.py
    AbiPy Manager:
    [Qadapter 0]
    ShellAdapter:localhost
    Hardware:
       num_nodes: 1, sockets_per_node: 1, cores_per_socket: 2, mem_per_node 4096,
    Qadapter selected: 0

    Abinitbuild:
    Abinit Build Information:
	Abinit version: 8.3.1
	MPI: True, MPI-IO: True, OpenMP: False
	Netcdf: True, ETSF-IO: False

    Abipy Scheduler:
    PyFlowScheduler, Pid: 71013
    Scheduler options: {'seconds': 10, 'hours': 0, 'weeks': 0, 'minutes': 0, 'days': 0}

    Installed packages:
    Package      Version
    -----------  -----------------------
    numpy        1.10.4
    scipy        0.17.0
    netCDF4      1.2.4
    apscheduler  2.1.0
    pydispatch   2.0.5
    yaml         3.11
    pymatgen     4.6.2
    matplotlib   1.5.1 (backend: Qt4Agg)


    Abipy requirements are properly configured

All went well and we can finally run our first calculation with Abipy.
The directory `abipy/data/runs` contains python scripts to generate workflows for typical ab-initio calculations.
Let's start from the simplest example i.e. the `run_si_ebands.py` script that generates 
a flow to compute the band structure of silicon at the Kohn-Sham level (GS calculation to get the
density followed by a Nscf run along a k-path in the first Brillouin zone).

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

`w0` is the directory containing the input files of the first workflow  (well, we have only one workflow in our example).
`t0` and `t1` contain the input files need to run the SCF and the NSC run, respectively.

You might have noticed that each `Task` directory (t0, t1) present the same structure:
    
   * run.abi: ABINIT input file
   * run.files: ABINIT files file
   * job.sh: Submission/shell script
   * outdata: Directory with output data files
   * indata: Directory with input data files 
   * tmpdata: Directory with temporary files

.. DANGER::
   `__AbinitFlow__.pickle` is the pickle file used to save the status of the `Flow`.
   Don't touch it! 

The `job.sh` has been generated using the information provided via `manager.yml`. 
In this case it's a simple shell script that executes the code but this is normal because we are using 
`qtype: shell`. The script will be more complicated when we start to submit jobs on a cluster with 
a resource manager.

We usually interact with the Abipy flow via the `abirun.py` script.
The script uses the syntax::

     `abirun.py FLOWDIR command [options]`

where `FLOWDIR` is the directory containing the flow generated by our script and
command defines the action to perform (use `--help` to get the list of possible commands).

Use::

	abirun.py flow_si_ebands/ status

to have a summary with the status of the different tasks and::

	abirun.py flow_si_ebands/ deps

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
abirun.py knows the dependencies of our flow and ...

.. code-block:: shell

    abirun.py flow_si_ebands/ rapid

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
The `rapid` command tried to executed all tasks that are READY. 
The SCF task has been submitted with 2 MPI processors

    * Queue: JobID @ QueueName  
    * MPI: Number of MPI processes used (computed automatically with autoparal, cannot exceed max_ncpus)
    * OMP: Number of OpenMP threads
    * Gb: Memory requested in Gb (meaningless in this case because we using the shell)
    * Warn: Number of warning messages found in the log file.
    * Com: Number of comments found in the log file.
    * Sub: Number of submissions (can be > 1 if Abipy encounters a problem and resubmit the task with different parameters
      without performing any operation that can change the physics of the system)
    * Rest: Number of restarts (Abipy can restart the job if convergence has not been reached)
    * Corr: Number of corrections performed. These operations can change the physics 
    * Time: Time spent in the Queue (if ends with Q) or running time (if ends with R)
    * Node_ID: Node identifier used by Abipy to identify each node of the flow. 

(when the submission is done through the shell there's no actual difference between submission and execution)

Now you see that one job is completed, run `rapidfire` again to execute the `NscfTask`,
wait a bit and then issue `abirun.py flow_si_ebands status` 

.. code-block:: shell

    $ abirun.py flow_si_ebands/ history

    Running on gmac2 -- system Darwin -- Python 2.7.12 -- abirun-0.1.0

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


.. TIP:: Use `abirun.py FLOWDIR events`  to print ABINIT events (Warning/Error/Comment).
         Use `abirun.py FLOWDIR debug`   to analyze error files and log files for possible error messages.

How to configure the scheduler
==============================

The `scheduler.yml` is another configuration file used to pass options to scheduler.
This file is much easier to understand and is needed only if you are running automatic workflows.
For this reason, we postpone the discussion of `scheduler.yml` and we focus on the
configuration of the task manager.

The other configuration file is named `scheduler.yml` and defines the parameters 
for the scheduler that will run/submit our jobs

Crate a `scheduler.yml` in the working directory with:

.. code-block:: yaml

    seconds: 10  # number of seconds to wait.
    #minutes: 0  # number of minutes to wait.
    #hours: 0    # number of hours to wait.
    #days: 0     # number of days to wait.

This file tells the scheduler to wake up every 10 seconds, inspect the status of the tasks
in the flow and perform the appropriate actions required for completion

.. TIP::

    Remember to set the time interval of the scheduler to a reasonable value.
    A small value leads to an increase of the submission rate but it also increases the CPU load 
    and the pressure on the hardware and on the resource manager.
    A too large time interval can have a detrimental effect on the throughput, especially 
    if you are submitting many small jobs.

.. TIP:: Use `abirun.py . doc_scheduler` to get the full list of options supported by the scheduler.

At this point, we are ready to run our first calculation with the scheduler.
Here we run a slight more complicated flow that computes the G0W0 corrections
to the direct band gap of silicon at the Gamma point.
The flow consists of the following tasks:

    1: ground state calculation to get the density
    2: NSCF calculation with several empty states. 
    3: calculation of the screening using the WFK produced by task 2
    4-5-6: Evaluation of the Self-Energy matrix elements with different values of nband 
     using the WFK produced by task 2 and the SCR file produced by task 3

Generate the flow with::

    ./run_si_g0w0.py

and let the scheduler manage the submission of the tasks with `abirun.py flow_si_g0w0 scheduler`
You should obtain the following output on the terminal

.. code-block:: shell

    abirun.py flow_si_ebands scheduler

    Abipy Scheduler:
    PyFlowScheduler, Pid: 72038
    Scheduler options: {'seconds': 2, 'hours': 0, 'weeks': 0, 'minutes': 0, 'days': 0}


Configuring AbiPy on a cluster
==============================

.. code-block:: yaml

    # Resource manager e.g slurm, pbs, shell
    qtype: slurm
    # Options passed to the resource manager (syntax depends on qtype, consult the manual of your resource manager)
    qparams:
	ntasks: 2
	time: 0:20:00
	partition: Oban
    # List of modules to import before running the calculation
    modules:
	- intel/compilerpro/13.0.1.117
	- fftw3/intel/3.3
    # Shell environment
    shell_env:
	 PATH: /home/user/local/bin/:$PATH
	 LD_LIBRARY_PATH: /home/user/local/lib:$LD_LIBRARY_PATH
    mpi_runner: /path/to/mpirun
    # Options for the automatic parallelization (Abinit autoparal feature)
    policy:
	autoparal: 1
	max_ncpus: 2


Description:

    * `qtype` specifies the queue resource manager (Slurm in this example). 

    * `qparams` is a dictionary with the parameters  passed to the resource manager. 
       We use the *normalized* version of the options i.e dashes in the official name of the parameter 
       are replaced by  underscores  (for the list of supported options see ...)

    * `modules` is the list of modules to load, while `shell_env` allows the user to specify or to modify 
       the values of the environment variables.

    * The `policy` section governs the automatic parallelization of the run: in this case abipy will use 
      the `autoparal` features of abinit to determine an optimal configuration with **maximum** `max_ncpus` MPI nodes. 
      Setting autoparal to 0 disables the automatic parallelization. **Other values of autoparal are not supported**.

One can put this configuration file either in the configuration directory `$HOME/.abinit/abipy` or in the current working directory (the latter has precedence over the global configuration file located in `$HOME/.abinit/abipy`).
The `TaskManager` can then be easily initialized by calling the class method `from_user_config`

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

The complete list of options supported by the `TaskManager` with slurm can be retrieved with the command::

    abirun.py . doc_manager slurm

In the previous sections, we have discussed how to define, build and run a `Flow`, but there is a very 
important point that we haven't discussed yet.
It should be stressed, indeed, that `AbiPy` is only driving and monitoring the `Flow` while the actual calculation 
is delegated to ABINIT (a Fortran program that is usually executed in parallel on multiple CPUs that communicate 
via the network by means of the MPI protocol).
Besides CPUs and memory must be reserved in advance by sending a request to the resource manager 
installed on the clusters (SLURM, PBS, etc)


.. TIP:: nohup abirun.py FLOWDIR scheduler 2> sched.log

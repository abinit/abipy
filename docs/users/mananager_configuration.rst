Introduction
============

Besides post-processing tools and a programmatic interface to generate Abinit input files,
AbiPy also provides a pythonic interface to invoke the code to retrieve
results and/or run calculations on supercomputing centers.
This section discusses how to create the configuration files required to interface AbiPy with Abinit.

We assume that Abinit is already available on your machine and that you know how to setup
your environment (set the PATH and LD_LIBRARY_PATH environment variables, load modules, etc.)
so that the operating system can load and execute the applications.
We will see that the configuration file requred by AbiPy

AbiPy knows how to run/submit the code with the proper environment thanks to the options
specified in the `manager.yml` configuration file.
`YAML <https://en.wikipedia.org/wiki/YAML>`_, a human-readable data serialization language commonly 
used for configuration files.

By default, AbiPy looks for a `manager.yml` file in the current working directory i.e.
the directory in which you execute your script and then inside `$HOME/.abinit/abipy`.
If no `manager.yml` is found, the code aborts immediately.

The `scheduler.yml` is another configuration file used to pass options to scheduler.
This file is much easier to understand and is needed only if you are running automatic workflows.
For this reason, we postpone the discussion of `scheduler.yml` and we focus on the
configuration of the task manager.

Configuration files for typical cases are available in ~abipy/data/managers

Configuring AbiPy on a personal computer
========================================

Let's start from the simplest case i.e. a personal computer in which we can 
execute Abinit directly within the shell.

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


Copy this example, change the entries in the `hardware` and the `limits` section according to
your machine, change `pre_run` so that the abinit executables can be found in $PATH.
Save the file in the current working directory and run `abicheck.py`.
If everything is configured properly, you should get something like this in the terminal.

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



Cd to abipy/data/runs, execute the `run_si_ebands.py` to generate a flow that 
computes the band structure of silicon.
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

Use `abirun flow_si_ebands status`


$ abirun.py flow_si_ebands/ status


 abirun.py flow_si_ebands/ rapid


.. code-block:: shell
		     `:-                                                               -:`
	     --`  .+/`                              `                                  `/+.  .-.
       `.  :+.   /s-                   `yy         .yo                                   -s/   :+. .`
     ./.  +o`   /s/           `-::-`   `yy.-::-`   `:-    .:::-`   -:`     .:`            /s/   :s- ./.
    .o.  /o:   .oo.         .oyo++syo. `yyys++oys. -ys  -syo++sy+` sy-     +y:            .oo-   oo` `o.
    ++   oo.   /oo          yy-    -yy `yy:    .yy`-ys .ys`    /yo sy-     +y:             oo/   /o:  ++
    +/   oo`   /oo         `yy.    .yy` yy.    `yy`-ys :ys     :yo oy/     oy:             +o/   :o:  /o
    -/   :+.   -++`         -sy+::+yyy` .sy+::+yy- -ys :yys/::oys. `oyo::/syy:            `++-   /+.  /:
     --  `//    /+-           -/++/-//    -/++/-   `+: :yo:/++/.     .:++/:oy:            -+/   `+-  --
      `.`  -:    :/`                                   :yo                 +y:           `/:`  `:. `.`
	    `..   .:.                                   .`                 `.           .:.  `..
		    ...                                                               ...

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

Now you see that one job is completed, run rapidfire again to execute the NscfTask.

At this point, we are ready to run our first calculation with the scheduler.
Crate a `scheduler.yml` in the working directory with:

.. code-block:: yaml

    # number of seconds to wait.
    seconds: 10

Remove the `flow_si_ebands` directory, regenerate the flow by re-running `run_si_ebands.py` and
execute the band structure calculation in an automatic way by issuing:

.. code-block:: shell

    abirun.py flow_si_ebands scheduler

.. code-block:: shell

    Abipy Scheduler:
    PyFlowScheduler, Pid: 72038
    Scheduler options: {'seconds': 2, 'hours': 0, 'weeks': 0, 'minutes': 0, 'days': 0}













.. _manager-examples:

****************
Manager Examples
****************

.. contents::
   :backlinks: top

Dragon1
-------

.. code-block:: yaml


	# dragon hardware: http://www.ceci-hpc.be/clusters.html#dragon
	hardware: &hardware
	   num_nodes: 26
	   sockets_per_node: 2
	   cores_per_socket: 8
	   mem_per_node: 112Gb
	
	job: &job
	    mpi_runner: mpirun
	    shell_env:
	        PATH: "$HOME/git_repos/abinit/_build_dragon1-intel-mpich-mkl.ac/src/98_main:$PATH"
	    modules:
	        - mpich/3.0.4/intel-13.0.0
	    pre_run: "ulimit -s unlimited"
	
	# queues
	qadapters:
	  - priority: 1
	    queue:
	       qtype: slurm
	       qname: Def
	    limits:
	       timelimit: 0-00:30:00
	       min_cores: 1
	       max_cores: 12
	       min_mem_per_proc: 1000
	       max_mem_per_proc: 2000
	       max_num_launches: 10
	    hardware: *hardware
	    job: *job


Gmac
----

.. code-block:: yaml


	qadapters:
	    - &batch
	      priority: 1
	      queue:
	        qname: gmac
	        qtype: shell
	      job:
	        mpi_runner: mpirun
	        pre_run:
	         - source ~/env.sh
	      limits:
	         min_cores: 1
	         max_cores: 1
	         timelimit: 0:10:0
	      hardware:
	         num_nodes: 1
	         sockets_per_node: 1
	         cores_per_socket: 2
	         mem_per_node: 4 Gb
	         # Optional
	         #condition: {"$eq": {omp_threads: 2}}
	
	batch_adapter: *batch


Hercules
--------

.. code-block:: yaml


	# hercules hardware: http://www.ceci-hpc.be/clusters.html#hercules
	hardware: &hardware
	   num_nodes: 65
	   sockets_per_node: 2
	   cores_per_socket: 8
	   mem_per_node: 54Gb
	
	job: &job
	    mpi_runner: mpirun
	    shell_env:
	        PATH: "$HOME/git_repos/abinit/_build_hercules.ac/src/98_main/:$PATH"
	    modules:
	        - impi/5.1.3.181-iccifort-2016.3.210-GCC-5.4.0-2.26
	        - imkl/11.3.3.210-iimpi-2016b
	    # here pre_run is a string in verbatim mode (note |)
	    pre_run: |
	        ulimit -s unlimited
	
	# queues
	qadapters:
	  - priority: 1
	    queue:
	       qtype: slurm
	       #qname: defq
	    limits:
	       timelimit: 0-00:30:00
	       min_cores: 1
	       max_cores: 12
	       min_mem_per_proc: 1000
	       max_mem_per_proc: 2000
	       max_num_launches: 10
	    hardware: *hardware
	    job: *job


Hmem
----

.. code-block:: yaml


	# hmem hardware: http://www.ceci-hpc.be/clusters.html#hmem
	# See also http://www.cism.ucl.ac.be/faq/index.php#hmem_specifics
	high: &high
	   num_nodes: 2
	   sockets_per_node: 4
	   cores_per_socket: 12
	   mem_per_node: 512Gb
	
	middle: &middle
	   num_nodes: 7
	   sockets_per_node: 4
	   cores_per_socket: 12
	   mem_per_node: 256Gb
	
	low: &low
	   num_nodes: 7
	   sockets_per_node: 4
	   cores_per_socket: 12
	   mem_per_node: 128Gb
	
	job: &job
	    mpi_runner: mpirun
	    shell_env:
	        PATH: "$HOME/git_repos/abinit/_build_hmem_intel_openmpi-mkl.ac/src/98_main/:$PATH"
	    modules:
	        - openmpi/1.5.3/intel-12.0.0.084
	    pre_run: "ulimit -s unlimited"
	
	# queues
	qadapters:
	  - priority: 3
	    #max_num_launches: 20
	    queue:
	       qname: High
	       qtype: slurm
	    limits:
	       timelimit: 10-0:0:0
	       min_cores: 1
	       max_cores: 48
	    hardware: *high
	    job: *job
	
	  - priority: 2
	    queue:
	       qname: Middle
	       qtype: slurm
	    limits:
	       timelimit: 5-0:0:0
	       min_cores: 1
	       max_cores: 48
	    hardware: *middle
	    job: *job
	
	  - priority: 1
	    queue:
	       qname: Low
	       qtype: slurm
	    limits:
	       timelimit: 5-0:0:0
	       min_cores: 1
	       max_cores: 48
	    hardware: *low
	    job: *job


Juqueen
-------

.. code-block:: yaml


	batch: &batch
	   num_nodes: 128
	   sockets_per_node: 1
	   cores_per_socket: 16
	   mem_per_node: 128Gb
	
	job: &job
	    mpi_runner: runjob
	    shell_env:
	        PATH: $HOME/abinit/801-private/bgq_xlf_legacy/src/98_main/:$PATH
	
	# List of qadapters
	# Note that on the BlueGeneQ we need at least two qadapters
	# One for submitting jobs to the computing nodes and another
	# one for executing small sequential ABINIT jobs on the frontend
	# The two qadapters have different shell environments, module files and binaries.
	qadapters:
	
	  # adapter for submitting jobs to the BlueGene.
	  - priority: 1
	    queue:
	       #qname: batch
	       qtype: bluegene
	       qparams:
	         # Mandatory on juqueen.
	         notification: error
	         mail_user: john@nowhere.com
	         environment: COPY_ALL
	    limits:
	       timelimit: 00:20:00
	       min_cores: 1
	       max_cores: 1024
	    hardware: *batch
	    job: *job
	
	  # shell adapter for small sequential jobs (e.g. autoparal tasks).
	  # Note that we need an Abinit executable that can be executed on the frontend
	  # TODO check priority
	  - priority: 10
	    queue:
	       qname: shell_adapter
	       qtype: shell
	    limits:
	       timelimit: 00:10:00
	       min_cores: 1
	       max_cores: 1
	    hardware:
	       num_nodes: 1
	       sockets_per_node: 1
	       cores_per_socket: 1
	       mem_per_node: 12Gb
	    job:
	        #mpi_runner: runjob
	        shell_env:
	            PATH: $HOME/abinit/801-private/bgq_frontend/src/98_main/:$PATH
	        modules:
	            gcc/4.8.3


Jureca
------

.. code-block:: yaml


	# See http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JURECA/Configuration/Configuration_node.html
	# and
	# http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JURECA/UserInfo/QuickIntroduction.html?nn=1803700#JURECABatchPart
	devel: &devel
	   num_nodes: 8
	   sockets_per_node: 2
	   cores_per_socket: 12
	   mem_per_node: 128Gb
	
	batch: &batch
	   num_nodes: 128
	   sockets_per_node: 2
	   cores_per_socket: 12
	   mem_per_node: 128Gb
	
	job: &job
	    # mpirun is not available on jureca.
	    # parallel applications must be executed with srun.
	    # shell_runner is used to run small sequential jobs on the frontend (e.g. autoparal jobs)
	    # None means that we should run the executable without prepending srun.
	    mpi_runner: srun
	    shell_runner: None
	    shell_env:
	        PATH: $HOME/abinit/801-private/jureca_mpi/src/98_main:$PATH
	    modules:
	        - intel-para/2015.07
	    pre_run: "ulimit -s unlimited"
	
	# queues
	qadapters:
	  - priority: 1
	    #max_num_launches: 20
	    queue:
	       qname: batch
	       qtype: slurm
	    limits:
	       timelimit: 0:10:0
	       min_cores: 1
	       max_cores: 12
	    hardware: *batch
	    job: *job


Lemaitre2
---------

.. code-block:: yaml


	# lemaitre2 hardware: http://www.ceci-hpc.be/clusters.html#lemaitre2
	hardware: &hardware
	   num_nodes: 112
	   sockets_per_node: 2
	   cores_per_socket: 6
	   mem_per_node: 48Gb
	
	job: &job
	    mpi_runner: mpirun
	    shell_env:  # Use your abinit exec
	        PATH: "$HOME/git_repos/abinit/_build_lemaitre2-intel-openmpi-mkl.ac/src/98_main/:$PATH"
	    modules: # Abinit compiled with abiconfig settings
	        - openmpi/1.6.5/intel-13.0.1.117
	    pre_run: "ulimit -s unlimited"
	
	# queues
	qadapters:
	  - priority: 1
	    queue:
	       qtype: slurm
	       qname: Def
	    limits:
	       timelimit: 0-0:30:00
	       min_cores: 1
	       max_cores: 12
	       min_mem_per_proc: 1000
	       max_mem_per_proc: 2000
	       max_num_launches: 10
	    hardware: *hardware
	    job: *job


Lemaitre3
---------

.. code-block:: yaml


	# lemaitre3 hardware: http://www.ceci-hpc.be/clusters.html#lemaitre3
	# For the configuration file see:
	#       https://github.com/abinit/abiconfig/blob/master/abiconfig/clusters/lemaitre3-intel-easybuild.ac
	hardware: &hardware
	   num_nodes: 80
	   sockets_per_node: 2
	   cores_per_socket: 12
	   mem_per_node: 95Gb
	
	job: &job
	    mpi_runner: mpirun
	    shell_env:  # Use your abinit exec
	        PATH: "$HOME/git_repos/abinit/_build_lemaitre3-intel-easybuild.ac/src/98_main/:$PATH"
	    modules: # Abinit compiled with abiconfig settings
	        - intel/2017b
	        - netCDF-Fortran/4.4.4-intel-2017b
	    pre_run: "ulimit -s unlimited"
	
	# queues
	qadapters:
	  - priority: 1
	    queue:
	       qtype: slurm
	       #qname: Def
	    limits:
	       timelimit: 0-0:30:00
	       min_cores: 1
	       max_cores: 12
	       min_mem_per_proc: 1000
	       max_mem_per_proc: 2000
	       max_num_launches: 10
	    hardware: *hardware
	    job: *job


Manneback
---------

.. code-block:: yaml


	# Hardware specification.
	Def: &Def
	   num_nodes: 672
	   sockets_per_node: 2
	   cores_per_socket: 4
	   mem_per_node: 24 Gb
	
	ObanAMD: &ObanAMD
	   num_nodes: 6
	   sockets_per_node: 4
	   cores_per_socket: 8
	   mem_per_node: 128 Gb
	
	ObanIntel: &ObanIntel
	   num_nodes: 3
	   sockets_per_node: 4
	   cores_per_socket: 8
	   mem_per_node: 256 Gb
	
	# Environment, modules, and parameters used to launch jobs.
	job: &job
	    mpi_runner: mpirun
	    shell_env:
	         PATH: "$HOME/git_repos/abinit/_build_manneback-gcc-openmpi.ac/src/98_main/:$PATH"
	    pre_run:
	        - "ulimit -s unlimited"
	        - "export OMP_NUM_THREADS=1"
	        - "unset SLURM_CPUS_PER_TASK"
	        - "module purge"
	        - "module load gompi/2016a FFTW/3.3.4-gompi-2016a"
	
	#policy:
	#   frozen_timeout: 0-12:0:0
	
	# List of qdapters.
	qadapters:
	  - priority: 1
	    queue:
	       qname: Def
	       qtype: slurm
	       qparams:
	                # This nodes must be excluded because they are not compatible with the Abinit build (SIGILL error).
	                exclude_nodes: mb-neh[070,201-212],mb-har[001-014],mb-har[101-116],mb-opt[111-116],mb-har[121-140],mb-sab[004,040,007,101-102],mb-wes[251-252],mb-ivy[205,206,208]
	    limits:
	       timelimit: 00:30:00
	       #timelimit_hard: 5-00:00:0
	       min_cores: 1
	       max_cores: 8
	       hint_cores: 4
	       min_mem_per_proc: 1000
	       max_mem_per_proc: 2000
	       max_num_launches: 5
	    job: *job
	    hardware: *Def


Nic4
----

.. code-block:: yaml


	# nic4 hardware. see http://www.ceci-hpc.be/clusters.html#nic4
	hardware: &hardware
	   num_nodes: 120
	   sockets_per_node: 2
	   cores_per_socket: 8
	   mem_per_node: 64Gb
	
	job: &job
	    mpi_runner: "mpirun"
	    mpi_runner_options: "--bind-to none"
	    shell_env:
	        PATH: "$HOME/git_repos/abinit/_build_nic4-intel-openmpi-mkl-hdf5.ac/src/98_main:$PATH"
	    pre_run: "ulimit -s unlimited"
	    modules:
	        - shared
	        - openmpi/1.7.5/intel2013_sp1.1.106
	        - intel/mkl/64/11.1/2013_sp1.1.106
	        - hdf5/1.8.13/openmpi-1.7.5-intel2013_sp1.1.106
	        - netcdf/4.3.2/openmpi-1.7.5-intel2013_sp1.1.106
	        - slurm/14.03.11
	
	# queues
	qadapters:
	  - priority: 1
	    queue:
	       qtype: slurm
	       qname: defq
	       qparams:
	          mail_type: FAIL
	          #mail_user: # Othere slurm options ...
	    limits:
	       timelimit: 0:30:0
	       min_cores: 1
	       max_cores: 16
	       min_mem_per_proc: 1000
	       max_mem_per_proc: 2000
	       max_num_launches: 5
	    hardware: *hardware
	    job: *job


Shell
-----

.. code-block:: yaml


	qadapters:
	    # List of qadapters objects
	    - priority: 1
	      queue:
	        qtype: shell
	        qname: localhost
	      job:
	        mpi_runner: mpirun
	        # source a script to setup the environment.
	        #pre_run: "source ~/env.sh"
	      limits:
	        timelimit: 1:00:00
	        max_cores: 2
	      hardware:
	         num_nodes: 1
	         sockets_per_node: 1
	         cores_per_socket: 2
	         mem_per_node: 4 Gb


Shell_nompi
-----------

.. code-block:: yaml


	qadapters:
	    # List of qadapters objects
	    - priority: 1
	      queue:
	        qtype: shell
	        qname: localhost
	      job:
	        mpi_runner: None
	        # source a script to setup the environment.
	        #pre_run: "source ~/env.sh"
	      limits:
	        timelimit: 1:00:00
	        max_cores: 1
	      hardware:
	         num_nodes: 1
	         sockets_per_node: 1
	         cores_per_socket: 2
	         mem_per_node: 4 Gb


Travis
------

.. code-block:: yaml


	qadapters:
	    -
	      priority: 1
	      queue:
	        qname: travis
	        qtype: shell
	      job:
	        mpi_runner: mpirun
	        pre_run:
	            - source activate test-environment
	      limits:
	         min_cores: 1
	         max_cores: 2
	         timelimit: 0:10:0
	      hardware:
	         num_nodes: 1
	         sockets_per_node: 1
	         cores_per_socket: 2
	         mem_per_node: 4 Gb


Ubu
---

.. code-block:: yaml


	qadapters:
	    # List of qadapters objects
	    - priority: 1
	      queue:
	        qtype: shell
	        qname: ubu
	      job:
	        modules:
	           - ubu_intel_16.0_mpich
	        mpi_runner: mpiexec
	        # source a script to setup the environment.
	        pre_run: "source ~/env.sh"
	      limits:
	        timelimit: 1:00:00
	        max_cores: 24
	      hardware:
	         num_nodes: 1
	         sockets_per_node: 1
	         cores_per_socket: 24
	         mem_per_node: 4 Gb


Vega
----

.. code-block:: yaml


	# vega hardware: http://www.ceci-hpc.be/clusters.html#vega
	hardware: &hardware
	   num_nodes: 44
	   sockets_per_node: 4
	   cores_per_socket: 16
	   mem_per_node: 256Gb
	
	job: &job
	    mpi_runner: mpirun
	    shell_env:
	        PATH: "$HOME/git_repos/abinit/_build_vega-intel-impi-mkl.ac/src/98_main/:$PATH"
	    modules:
	        - intel/2015a
	    #pre_run: "ulimit -s unlimited"
	
	# queues
	qadapters:
	  - priority: 1
	    queue:
	       qtype: slurm
	       qname: defq
	    limits:
	       timelimit: 0-0:30:0
	       min_cores: 1
	       max_cores: 16
	       min_mem_per_proc: 1000
	       max_mem_per_proc: 2000
	       max_num_launches: 5
	    hardware: *hardware
	    job: *job


Viper
-----

.. code-block:: yaml


	hardware: &hardware
	   num_nodes: 1
	   sockets_per_node: 2
	   cores_per_socket: 4
	   mem_per_node: 32Gb
	
	job: &job
	    mpi_runner: ~/bin/mpirun.openmpi
	    # pre_run is a string in verbatim mode (note |)
	    pre_run:
	        - "ulimit -s unlimited"
	        - "source ~/.bashrc"
	
	# queues
	qadapters:
	  - priority: 1
	    queue:
	       qname: euspec.q
	       qtype: sge
	       qparams:
	           parallel_environment: slots
	    limits:
	       timelimit: 0:10:0
	       min_cores: 1
	       max_cores: 8
	    hardware: *hardware
	    job: *job


Zenobe
------

.. code-block:: yaml


	# Hardware specification.
	westmere: &westmere
	   num_nodes: 274
	   sockets_per_node: 2
	   cores_per_socket: 6
	   mem_per_node: 24 Gb
	
	ivybridge: &ivybridge
	   num_nodes: 342
	   sockets_per_node: 2
	   cores_per_socket: 12
	   mem_per_node: 64 Gb
	
	# Environment, modules, and parameters used to launch jobs.
	job: &job
	    mpi_runner: mpirun
	    shell_env:
	         PATH: $HOME/git_repos/abinit_build_impi/src/98_main:$PATH
	    modules:
	        - compiler/intel/composerxe/2013_sp1.1.106
	        - intelmpi
	        - python/2.7
	    pre_run: "ulimit -s unlimited"
	
	# List of qdapters.
	qadapters:
	  # Westmere default.
	  - priority: 99
	    queue:
	       qname: main
	       qtype: pbspro
	       qparams:
	         group_list: napsimu
	         #qverbatim: |
	         #  #PBS -r y
	    limits:
	       timelimit: 15:0
	       min_cores: 1
	       max_cores: 24
	    job: *job
	    hardware: *westmere
	
	  # Ivybridge large.
	  - priority: 1
	    queue:
	       qname: large
	       qtype: pbspro
	       qparams:
	          group_list: napsimu
	          #qverbatim: |
	          #  #PBS -r y
	    limits:
	       timelimit: 1-0:0:0
	       min_cores: 96
	       max_cores: 3888
	    job: *job
	    hardware: *ivybridge


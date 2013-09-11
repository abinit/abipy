#!/usr/bin/env python
from __future__ import print_function, division

from abipy.abilab import qadapters 

"""
setup:
    String or list of commands to execute during the initial setup.
modules:
    String or list of modules to load before running the application.
shell_env:
    Dictionary with the shell environment variables to export
    before running the application.
omp_env:
    Dictionary with the OpenMP variables.
pre_run:
    String or list of commands to execute before launching the calculation.
post_run:
    String or list of commands to execute once the calculation is completed.
mpi_runner:
    Path to mpirun or `MpiRunner` instance. None if MPI is not used.
"""

# A simple shell script that defines two variables and run abinit in parallel with 2 CPUs
shell = qadapters.ShellAdapter(
    qparams=dict(MPI_NCPUS=2), 
    # Set the value of FOO, unset BAR, modify the value of PATH so that the OS can find our executable. 
    shell_env=dict(FOO=1, BAR=None, PATH="/home/user/bin:$PATH"),
    mpi_runner="mpirun",
    )

# Generate the script (we assume that executable and mpirun are located in $PATH.
script = shell.get_script_str(job_name="job.sh", launch_dir="/path/to/lauch_dir",
    executable="abinit", stdin="STDIN_FNAME", stdout="STDOUT_FNAME", stderr="STDERR_FNAME")

print(script)

# Here we use the SlurmAdapter to generate a submission script for Slurm.
slurm = qadapters.SlurmAdapter(
    qparams=dict(
        ntasks=12,
        partition="hmem",
        account='nobody@nowhere.org',
        time="119:59:59",
        #ntasks_per_node=None,
        #cpus_per_task=None,
        #ntasks=None,
        #time=None,
        #partition=None,
        #account=None,
    ),
    setup=["echo 'This is the list of commands executed during the initial setup'", "ssh user@node01"],
    modules=['intel-compilers/12.0.4.191', 'MPI/Intel/mvapich2/1.6', 'FFTW/3.3'],
    shell_env=dict(FOO=1, PATH="/home/user/bin:$PATH"),
    omp_env=dict(OMP_NUM_THREADS=1),
    pre_run=["echo 'List of command executed before launching the calculation'" ],
    post_run=["echo 'List of command executed once the calculation is completed'" ],
    mpi_runner="mpirun",
)

script = slurm.get_script_str(job_name="job.sh", launch_dir="/path/to/lauch_dir",
    executable="abinit", stdin="STDIN_FNAME", stdout="STDOUT_FNAME", stderr="STDERR_FNAME")

print(script)

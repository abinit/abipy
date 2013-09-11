#!/usr/bin/env python
from __future__ import print_function, division

from abipy.abilab import qadapters 


"""
setup:
    String or list of commands executed during the initial setup.
modules:
    String or list of modules to load before running the application.
shell_env:
    Dictionary with the shell environment variables to export
    before running the application.
omp_env:
    Dictionary with the OpenMP variables.
pre_run:
    String or list of commands executed before launching the calculation.
post_run:
    String or list of commands executed once the calculation is completed.
mpi_runner:
    Path to mpirun or `MpiRunner` instance. None if not used
"""


shell = qadapters.ShellAdapter(
    qparams={},
    setup=["echo 'List of commands executed during the initial setup'", "ssh user@node01"],
    # List of modules to load before running the application.
    modules= ['intel-compilers/12.0.4.191', 'MPI/Intel/mvapich2/1.6', 'FFTW/3.3'],
    # Dictionary with the shell environment variables to export before running the application.
    shell_env = dict(MPI_NCPUS=1, BAR=None, PATH="/home/user/bin:$PATH"),
    # OpenMP variables.
    omp_env = dict(OMP_NUM_THREADS=1),
    pre_run = ["echo 'List of command executed before launching the calculation'" ],
    post_run = ["echo 'List of command executed once the calculation is completed'" ],
    mpi_runner = "mpirun",
    )

script = shell.get_script_str(job_name="myjob", launch_dir="hello", 
    executable="abinit", stdin="STDIN", stdout="STDOUT", stderr="STDERR")

print(script)

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
    setup = ["echo 'This is the list of commands executed during the initial setup'", "ssh user@node01"],
    # List of modules to load before running the application.
    modules = ['intel-compilers/12.0.4.191', 'MPI/Intel/mvapich2/1.6', 'FFTW/3.3'],
    # Dictionary with the shell environment variables to export before running the application.
    shell_env = dict(FOO=1, PATH="/home/user/bin:$PATH"),
    # OpenMP variables.
    omp_env = dict(OMP_NUM_THREADS=1),
    pre_run = ["echo 'List of command executed before launching the calculation'" ],
    post_run = ["echo 'List of command executed once the calculation is completed'" ],
    mpi_runner= "mpirun",
)

script = slurm.get_script_str(job_name="myjob", launch_dir="hello", 
    executable="abinit", stdin="STDIN", stdout="STDOUT", stderr="STDERR")

print(script)

#for qad in [shell, slurm]:
#    for num in [2, 3]:
#        qad.set_mpi_ncpus(num)
#
#        script = qad.get_script_str(job_name="myjob", launch_dir="hello", 
#            executable="abinit", stdin="STDIN", stdout="STDOUT", stderr="STDERR")
#        #script_file = "job_file.sh"
#        #with open(script_file, "w") as fh:
#        #    fh.write(script)
#
#        #process = qad.submit_to_queue(script_file)
#        #import cPickle as pickle
#        #with open("test.pickle", "w") as fh:
#        #    #pickle.dump(qad, fh, protocol=0)
#        #    pickle.dump(cls, fh, protocol=0)

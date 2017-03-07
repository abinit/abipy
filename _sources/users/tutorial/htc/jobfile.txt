
############################################
Jobfile: Controling execution and submission
############################################

The job file serves both to execute abinit and to submit the task to a batch server.

Example 1: Setting execution parameters
=======================================

.. code-block:: python

    from abipy.htc import Launcher

    calc = Launcher('Silicon', pseudodir='Data/',
                                     bindir='/Users/Antonius/Work/Software/abinit/7.2.0-private/build1/src/98_main/')
    # Pseudos
    calc.set_pseudos('14si.pspnc')
    
    # Inputs
    calc.read('Data/Silicon_ground_state.in')
    
    
    # Run script
    calc.set_mpirun('openmpi -np 4')
    
    calc.set_modules('intel-compilers/12.0.4.191', 'FFTW/3.3')
    
    calc.set_lines_before('echo "Job starts at:"', 'date')
    calc.set_lines_after('echo "Job ends at:"', 'date')
    
    calc.execute()


Example 2: Making a submission script
=====================================

.. code-block:: python

    from abipy.htc import Launcher
    
    calc = Launcher('Silicon',
                    jobtype='PBS', # PBS, SGE, Slurm
                    pseudodir='Data/',
                    bindir='/Users/Antonius/Work/Software/abinit/7.2.0-private/build1/src/98_main/',
                    mpirun='mpiexec',
                    modules=['intel-compilers/12.0.4.191', 'MPI/Intel/mvapich2/1.6', 'FFTW/3.3'],
                    lines_before=['echo "Job starts at:"', 'date'],
                    lines_after=['echo "Job ends at:"', 'date'],
                    )
    
    # Pseudos
    calc.set_pseudos('14si.pspnc')
    
    # Inputs
    calc.read('Data/Silicon_ground_state.in')
    
    
    # Run script
    calc.set_jobname('MyJob')
    calc.set_nodes(4)
    calc.set_ppn(12)
    calc.set_runtime(24)
    
    calc.set_queue('hpcourte')
    calc.set_mail('gabriel.antonius@umontreal.ca')
    
    calc.set_submission_command('qsub')  # This is the default.
    
    calc.execute()


The available commands (set_nodes, set_ppn) will depend on the job type.
See the subclasses of :class:`~abipy.htc.JobFile`.

Submitting
----------

To submit with the command 'qsub', just use

.. code-block:: bash

    $ python myscript.py -v submit
    Submitting Silicon/calc

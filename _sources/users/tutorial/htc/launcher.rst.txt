
##############################################
Launcher: Creating and executing a calculation
##############################################

This tutorial aims to show how to use the :class:`~abipy.htc.Launcher` object
to write a calculations, and eventually...launch it.

First example
=============

Let’s start with our favourite example: silicon!

Writing the calculation
-----------------------

.. code-block:: python

    from abipy import Launcher

    # ============ Initialize calculation ============ #
    calc = Launcher('Silicon/rootname')

    # Set the executable
    calc.set_executable('/home/antonius/abinit/7.1.2/build/src/98_main/abinit')

    # Set the pseudopotentials files
    calc.set_pseudodir('/home/antonius/Atoms')
    calc.set_pseudos('14-Si.pspnc')


    # ============ Input variables ============ #
    calc.ndtset = 2

    # Ground state
    calc.iscf1 = 7
    calc.tolvrs1 = 1e-10

    # Non self-consistent run
    calc.iscf2 = -2
    calc.tolwfr2 = 1e-18
    calc.getden2 = -1

    # Basis set
    calc.ecut = 10.

    # K-points grid
    calc.kptopt = 1
    calc.ngkpt = [2, 2, 2]
    calc.nshiftk = 1
    calc.shiftk = [0.0, 0.0, 0.0]
    
    # Definition of the unit cell
    calc.acell = 3*[10.217]
    calc.rprim = [[.0, .5, .5],[.5, .0, .5],[.5, .5, .0]]
    calc.ntypat = 1
    calc.znucl = 14,
    calc.natom = 2
    calc.typat = [1, 1]
    calc.xred = [[.0, .0, .0], [.25,.25,.25]]
    calc.xcart = None


    # ============ Execution ============ #

    # Write the files.
    calc.make()

    # Run the calculation with abinit.
    calc.run()  

    # Print the calculation status.
    calc.report()


A few remarks:

1) All the calculation files are contained in a directory of its own,
   which is the name used to initialize Launcher. If not provided,
   the root name for the files is chosen automatically.

2) All pseudopotential files should be contained in the same directory.
   This directory is designated with the function :func:`~abipy.htc.FilesFile.set_pseudodir`.

3) The abinit executable is called through a bash script,
   written with the :class:`~abipy.htc.JobFile` object.


Launching the calculation from the shell
----------------------------------------

This script is quite complete. Running it should give:

.. code-block:: bash

    $ python myscript.py
    Writing Silicon/calc
    Running Silicon/calc
    Silicon/calc.out : Completed

Real life calculations however, don't always go smoothly,
and you might want to execute each step (writing, running, reporting)
at different times.

You can control the execution part from the shell.
To do so, replace all the last section of the file with the following:

.. code-block:: python

    # ============ Execution ============ #

    # Take action from the command line.
    calc.execute()

Now, when executing the script, you simply need to give the action at the end of the line.

.. code-block:: bash

    $ python myscript.py -v make
    Writing Silicon/calc
    $ python myscript.py -v run
    Running Silicon/calc
    $ python myscript.py -v report
    Silicon/calc.out : Completed

The '-v' option is for a verbose output (printing a line).

When the calculation is done, you can remove the data files and the log with

.. code-block:: bash

    $ python myscript.py -v clean
    Cleaning Silicon/calc
    About to remove the following files:
    Silicon/run/calc.log
    Silicon/run/out_data/odat_calc_DS1_DDB
    Silicon/run/out_data/odat_calc_DS1_DEN
    Silicon/run/out_data/odat_calc_DS1_EIG
    Silicon/run/out_data/odat_calc_DS1_WFK
    Silicon/run/out_data/odat_calc_DS2_DEN
    Silicon/run/out_data/odat_calc_DS2_EIG
    Silicon/run/out_data/odat_calc_DS2_WFK
    Do you want to proceed? (y/n)


===============================================

Using blocks of variables.
==========================

Rather than setting each variable individually, it is convenient to gather them in blocks,
using dictionaries.

.. code-block:: python

    from abipy.htc import Launcher
    
    calc = Launcher('ZnO', pseudodir='Data/',  # Directories set at initialization
                                bindir='/Users/Antonius/Work/Software/abinit/7.2.0-private/build1/src/98_main/')
    
    calc.set_pseudos('Zn.psp', 'O.psp')
    
    ZnO_unit_cell = {
        'acell' : 3*[8.6277],
        'rprim' : [[.0, .5, .5], [.5, .0, .5], [.5, .5, .0]],
        'ntypat' : 2,
        'znucl' : [30, 8],
        'natom' : 2,
        'typat' : [1, 2],
        'xred' : [[.0, .0, .0], [.25,.25,.25]]}
    
    shifted_kpt_grid = {
        'kptopt' : 1,
        'ngkpt' : [2,2,2],
        'nshiftk' : 4,
        'shiftk' : [[.5,.5,.5], [.5,.0,.0], [.0,.5,.0], [.0,.0,.5]]}
    
    unshifted_kpt_grid = {
        'kptopt' : 1,
        'ngkpt' : [2,2,2],
        'nshiftk' : 1,
        'shiftk' : [.0,.0,.0]}
    
    calc.ecut = 7.
    calc.nstep = 1
    
    calc.fuzzy = 1     # Watch out! Non existing variables are still written in the input file.
    calc.fuzzy = None  # Better unset it!
    
    calc.ndtset = 2
    calc.tolvrs1 = 1e-2
    calc.tolwfr2 = 1e-4
    calc.getden2 = -1
    calc.iscf2 = -2
    
    calc.set_variables(ZnO_unit_cell)                  # Add the blocks of input variables...
    calc.set_variables(shifted_kpt_grid)               #
    calc.set_variables(unshifted_kpt_grid, dataset=2)  # Add a dataset index to all variables.
    
    calc.execute()




Elaborated procedures: GW calculation
-------------------------------------

.. code-block:: python

    from abipy.htc import Launcher

    calc = Launcher('GW', pseudodir='Data/',
                                bindir='/Users/Antonius/Work/Software/abinit/7.2.0-private/build1/src/98_main/')
    
    calc.set_pseudos('Zn.psp', 'O.psp')
    
    nband_occ = 13
    nband_unocc = 20
    
    Ground_state = {
        'tolvrs' : 1e-2,
        'nstep' : 1,
        'iscf' : 7,
        'ixc' : 7,
        'nband' : nband_occ + 3,
        'kptopt' : 1,
        'ngkpt' : [2,2,2],
        'nshiftk' : 4,
        'shiftk' : [[.5,.5,.5], [.5,.0,.0], [.0,.5,.0], [.0,.0,.5]],
        }
    
    Wavefunctions = {
        'getden' : -1,
        'iscf' : -2,
        'tolwfr' : 1e-3,
        'nstep' : 1,
        'nband' : nband_occ + nband_unocc,
        }
    
    Screening = {
        'optdriver' : 3,
        }
    
    Sigma = {
        'optdriver' : 4,
        'nkptgw' : 1,
        'kptgw' : [0,0,0],
        'bdgw' : [nband_occ-2, nband_occ+3],
        }
    
    GW_parameters = {
        'gwcalctyp' : '00',
        'awtr' : 1,
        'symchi' : 1,
        'symsigma' : 0,
        'gwpara' : 2,
        'ppmodel' : 2,
        'ecutwfn' : 10.,
        'ecutsigx' : 20.,
        'ecuteps' : 2.,
        'nband' : nband_occ + nband_unocc,
        }
    
    Parameters = {
        'ecut' : 10.,
        'kptopt' : 1,
        'ngkpt' : [2,2,2],
        'nshiftk' : 1,
        'shiftk' : [.0,.0,.0],
        'istwfk' : '*1',
        }
    
    ZnO_unit_cell = {
        'acell' : 3*[8.6277],
        'rprim' : [[.0, .5, .5], [.5, .0, .5], [.5, .5, .0]],
        'ntypat' : 2,
        'znucl' : [30, 8],
        'natom' : 2,
        'typat' : [1, 2],
        'xred' : [[.0, .0, .0], [.25,.25,.25]]
        }
    
    
    # Common parameters
    calc.set_variables(ZnO_unit_cell)
    calc.set_variables(Parameters)
    calc.set_variables(GW_parameters)
    
    # Procedure
    calc.ndtset = 4
    calc.set_variables(Ground_state, 1)
    calc.set_variables(Wavefunctions, 2)
    calc.set_variables(Screening, 3)
    calc.set_variables(Sigma, 4)
    
    # Link files
    calc.link_io(idtset=3, odtset=2, datatype='WFK') # Here we link the output data files produced by a dataset
    calc.link_io(4, 2, 'WFK')                        # as an input data file for an other dataset.
    calc.link_io(4, 3, 'SCR')
    
    # Execution
    calc.execute()

The function :func:`~abipy.htc.AbinitInput.link_io` creates a symbolic link for the input data file
of a given dataset which points to the output data file of an other dataset.
If we want to link an external file as an input data file,
we use :func:`~abipy.htc.AbinitInput.link_idat` as like this:

.. code-block:: python

    # Procedure
    calc.ndtset = 2
    calc.jdtset = [3,4]
    calc.set_variables(Screening, 3)
    calc.set_variables(Sigma, 4)
    
    calc.link_idat(file='Data/odat_calc_DS2_WFK', dtset=3)
    calc.link_idat('Data/odat_calc_DS2_WFK', 4)
    calc.link_io(4, 3, 'SCR')


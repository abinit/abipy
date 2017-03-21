
####################################################
MassLauncher: Managing several calculations at once.
####################################################


Example 1: Multiple calculations
================================

.. code-block:: python

    from abipy.htc import MassLauncher
    
    calcs = MassLauncher(3, 'Example5/calc', pseudodir='Data/',
                                             bindir = '/Users/Antonius/Work/Software/abinit/7.2.0-private/build1/src/98_
    
    Si_unit_cell = {
        'acell' : 3*[10.261],
        'ntypat' : 1,
        'znucl' : [14],
        'natom' : 2,
        'typat' : [1, 1],
        'rprim' : [[.0, .5, .5], [.5, .0, .5], [.5, .5, .0]],
        'xred' : [[.0, .0, .0], [.25,.25,.25]],
        }
    
    Parameters = {
        'kptopt' : 1,
        'ngkpt' : [2,2,2],
        'nshiftk' : 4,
        'shiftk' : [[.5,.5,.5], [.5,.0,.0], [.0,.5,.0], [.0,.0,.5]],
        'nstep' : 1,
        'tolvrs' : 1e-2,
        'ecut' : 7.,
        }
    
    calcs.set_variables(Si_unit_cell)
    calcs.set_variables(Parameters)
    
    # Works as a list
    calcs[0].ixc = 7
    calcs[1].ixc = 7
    calcs[2].ixc = 11
    
    pseudos = ('14si.pspnc', '14-Si.LDA.fhi', '14-Si.GGA.fhi')
    for calc, psp in zip(calcs, pseudos):
        calc.set_pseudos(psp)
    
    
    calcs.execute()


Execution
---------

By default, an action given from the command line is execute for all calculations.

.. code-block:: bash

    $ python myscript.py -v make
    Writing Example5/calc1
    Writing Example5/calc2
    Writing Example5/calc3

But you can select a subset of calculations to perform an action.

.. code-block:: bash

    $ python myscript.py -v run -c 1 2
    Running Example5/calc1
    Running Example5/calc2

    $ python myscript.py -v report
    Example5/calc1.outC : Completed
    Example5/calc2.outG : Completed
    Example5/calc3 : Unstarted




Example 2: Convergence study of a GW calculation
================================================

.. code-block:: python

    from abipy.htc import MassLauncher
    
    calcs = MassLauncher(3, 'Example6/calc', pseudodir='Data/',
                                             bindir='/Users/Antonius/Work/Software/abinit/7.2.0-private/build1/src/98_ma
    
    calcs.set_pseudos('Zn.psp', 'O.psp')
    
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
        'xred' : [[.0, .0, .0], [.25,.25,.25]],
        }
    
    
    calcs.set_variables(ZnO_unit_cell)
    calcs.set_variables(Parameters)
    calcs.set_variables(GW_parameters)
    
    calcs.ndtset = 2
    calcs.jdtset = [3, 4]
    #calcs.set_variables(Ground_state, 1)
    #calcs.set_variables(Wavefunctions, 2)
    calcs.set_variables(Screening, 3)
    calcs.set_variables(Sigma, 4)
    
    calcs.link_idat('Data/odat_calc_DS2_WFK', 3)
    calcs.link_idat('Data/odat_calc_DS2_WFK', 4)
    calcs.link_io(4, 3, 'SCR')
    
    # Increasing ecuteps in each calculation
    ecuteps = 1.
    for calc in calcs:
        calc.ecuteps = ecuteps
        ecuteps += .5
    
    
    calcs.execute()


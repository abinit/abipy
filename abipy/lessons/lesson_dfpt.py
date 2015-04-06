#!/usr/bin/env python
"""Phonon band structure of AlAs."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import numpy as np
import abipy.abilab as abilab
import abipy.data as abidata  


def make_scf_input(ecut=2, paral_kgb=0):
    """
    This function constructs the input files for the phonon calculation: 
    GS input + the input files for the phonon calculation.
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    gs_inp = abilab.AbinitInput(structure=abidata.structure_from_ucell("AlAs"),
                                pseudos=abidata.pseudos("13al.981214.fhi", "33as.pspnc"))
    
    gs_inp.set_vars(
        nband=4,             
        ecut=ecut,         
        ngkpt=[4, 4, 4],
        nshiftk=4,
        shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
                0.0, 0.5, 0.0,
                0.5, 0.0, 0.0,
                0.5, 0.5, 0.5],
        #shiftk=[0, 0, 0],
        paral_kgb=paral_kgb,
        ixc=1,
        nstep=25,
        diemac=9.0,
        tolvrs=1.0e-10,
    )
    
    gs_inp.set_mnemonics(True)
    return gs_inp


#@abilab.flow_main
#def main(options):
#    qpoints = np.reshape([
#             0.00000000E+00,  0.00000000E+00,  0.00000000E+00, 
#             2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
#    ], (-1, 3))
#
#    #return make_phq_ecutconv_flow(workdir="flow_alas_qconv", qpt=qpoints[0])
#
#    scf_input = make_scf_input(ecut=2, paral_kgb=0)
#
#    from pymatgen.io.abinitio.flows import phonon_conv_flow
#    flow = phonon_conv_flow("flow_conv", scf_input, qpoints[0], 
#        params=["ecut", [1, 2, 3]]
#        #params=["tsmear", "ngkpt", [0.002, 0.003, 0.004], [[2, 2, 2], [4, 4, 4], [6, 6, 6]]],
#    )
#    from pymatgen.io.abinitio.flows import PhononFlow
#    flow = PhononFlow.from_scf_input("flow_alas_phonons_and_becs", scf_input, ph_ngqpt=(4,4,4))
#
#    flow.build_and_pickle_dump()
#    return flow
#
#
#if __name__ == "__main__":
#    sys.exit(main())

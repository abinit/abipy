#!/usr/bin/env python
"""
This example shows how to build multiple flows and use the BatchLauncher to execute them.
"""
from __future__ import division, print_function, unicode_literals

import sys
import os
import abipy.data as abidata  
import abipy.abilab as abilab


def make_scf_nscf_inputs(paral_kgb=1):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    #inp.set_mnemonics(True)
    structure = inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    ecut = 6
    global_vars = dict(ecut=ecut,
                       nband=8,
                       nstep=15,
                       paral_kgb=paral_kgb,
                    )

    if inp.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    inp.set_vars(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    inp[1].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    inp[2].set_kpath(ndivsm=6, kptbounds=kptbounds)
    inp[2].set_vars(tolwfr=1e-12)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()
    return scf_input, nscf_input

def build_flow(workdir, paral_kgb):
    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs(paral_kgb)

    # Build the flow.
    return abilab.bandstructure_flow(workdir, scf_input, nscf_input)


def main():
    batch_dir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    from pymatgen.io.abinitio.launcher import BatchLauncher
    batch = BatchLauncher(workdir=batch_dir)

    # Build multiple flows and add them to the BatchLauncher.
    # Each flow has a unique wordir inside batch.workdir
    for paral_kgb in range(2):
        flow_dir = os.path.join(batch.workdir, "flow_paral_kgb_%d" % paral_kgb)
        batch.add_flow(build_flow(workdir=flow_dir, paral_kgb=paral_kgb))

    # Submit to the queue in batch mode.
    # Use abibatch.py to inspect the status or resubmit.
    print("batch.submit() returned: ", batch.submit())

if __name__ == "__main__":
    sys.exit(main())

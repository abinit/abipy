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

    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"), pseudos=pseudos, ndtset=2)

    # Global variables
    ecut = 6
    global_vars = dict(
        ecut=ecut,
        nband=8,
        nstep=15,
        paral_kgb=paral_kgb,
    )

    if multi.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    multi[0].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
    multi[1].set_vars(tolwfr=1e-12)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = multi.split_datasets()
    return scf_input, nscf_input


def main():
    batch_dir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # intialize the BatchLauncher.
    from pymatgen.io.abinit.launcher import BatchLauncher
    batch = BatchLauncher(workdir=batch_dir)

    # Build multiple flows and add them to the BatchLauncher.

    for paral_kgb in range(2):
        #flow_dir = os.path.join(batch.workdir, "flow_paral_kgb_%d" % paral_kgb)

        # Get the SCF and the NSCF input and build the flow.
        scf_input, nscf_input = make_scf_nscf_inputs(paral_kgb)

        # Each flow will have a unique wordir inside batch.workdir. 
        # Note that we have to pass workdir=None and use set_name to specify the dir basename
        flow = abilab.bandstructure_flow(None, scf_input, nscf_input, allocate=False)
        flow.set_name("flow_paral_kgb_%d" % paral_kgb)

        batch.add_flow(flow)

    # Submit to the queue in batch mode.
    # Use abibatch.py to inspect the status or resubmit.
    job = batch.submit()
    print("batch.submit() returned: ", job.retcode)
    return job.retcode

if __name__ == "__main__":
    sys.exit(main())

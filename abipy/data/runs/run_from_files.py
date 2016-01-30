#!/usr/bin/env python
"""Flow for computing the band structure of silicon."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata  
import abipy.abilab as abilab


def make_scf_nscf_inputs(paral_kgb=1):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"), 
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)

    # Global variables
    ecut = 6
    global_vars = dict(ecut=ecut,
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


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs()

    # Build the flow.
    flow = abilab.Flow(workdir, manager=options.manager)

    # Create a Work, all tasks in work will read the file f
    # Note that the file must exist when the work is created
    # Use the standard approach based on tasks and works if
    # there's a node who needs a file produced in the future.
    work = abilab.Work()
    denfile = abidata.ref_file("si_DEN-etsf.nc")
    work.register(nscf_input, deps={denfile: "DEN"})
    flow.register_work(work)

    return flow
    

@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())

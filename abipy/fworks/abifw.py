#!/usr/bin/env python
from fireworks import FireTaskBase, FWAction, Firework, LaunchPad, ScriptTask
from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.core.rocket_launcher import launch_rocket

from abipy.fworks.tasks import *

def main():
    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)

    # Build the flow
    flow = build_flow()
    flow.build_and_pickle_dump()

    # create the Firework consisting of a single task
    firework = Firework(FireTaskWithFlow(flow=flow))

    # store workflow 
    launchpad.add_wf(firework)

    #launch it locally
    #launch_rocket(launchpad)

    return 0


import sys
import os
import abipy.data as abidata  
import abipy.abilab as abilab
from abipy.tools.devtools import profile


def make_scf_nscf_inputs():
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos("14si.pspnc")
    #pseudos = data.pseudos("Si.GGA_PBE-JTH-paw.xml")

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    print(inp.pseudos)
    structure = inp.set_structure_from_file(abidata.cif_file("si.cif"))

    # Global variables
    ecut = 6
    global_vars = dict(ecut=ecut,
                       nband=8,
                       timopt=-1,
                       istwfk="*1",
                       nstep=15,
                       paral_kgb=0,
                       #accesswff=3,
                    )

    if inp.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    inp[1].set_variables(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    inp[2].set_kpath(ndivsm=6, kptbounds=kptbounds)
    inp[2].set_variables(tolwfr=1e-12)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()
    return scf_input, nscf_input

def build_flow():
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = "flow_test"

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config()

    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs()

    # Build the flow.
    return abilab.bandstructure_flow(workdir, manager, scf_input, nscf_input)
    

if __name__ == "__main__":
    #profile("main()", local_vars=locals(), global_vars=globals())
    main()

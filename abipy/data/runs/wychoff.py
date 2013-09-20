#!/usr/bin/env python
from __future__ import division, print_function

import os

from abipy import abilab
import abipy.data as data

from abipy.abilab import FloatWithUnit

def special_positions(lattice, u):
    """Construct the crystalline `Structure` for given value of the internal parameter u."""
    frac_coords = {}

    frac_coords["Si"] = [ [1/2 + u, 1/2 - u, 0],
                          [u,       -u,      u],
                        ]

    frac_coords["O"] = [ [0, 0, u] ]
                            
    species, coords = [], []

    for symbol, positions in frac_coords.items():
        species += len(positions) * [symbol]
        coords += positions

    return abilab.Structure(lattice, species, coords, 
                            validate_proximity=True, coords_are_cartesian=False)

def build_workflows():
    pseudos = data.pseudos("14si.pspnc")

    base_structure = abilab.Structure.from_file(data.cif_file("si.cif"))

    news, uparams = [], [0.1, 0.2, 0.3]

    for u in uparams:
        new = special_positions(base_structure.lattice, u)
        print(new)
        news.append(new)
        #assert len(base_structure) == len(new)

    # Create the list of workflows. Each workflow defines a band structure calculation.
    works = []
    for new_structure, u in zip(news, uparams):

        s = new_structure.formula.replace(" ", "")
        workdir = os.path.join(os.path.dirname(__file__), s + "_u" + str(u))

        # Generate the workflow.
        works.append(build_workflow(workdir, new_structure, pseudos))

    return works


def build_workflow(workdir, structure, pseudos):

    gs_inp = abilab.AbiInput(pseudos=pseudos)
    gs_inp.set_structure(structure)

    # Variables global to the SCF and the NSCF run.
    global_vars = dict(ecut=FloatWithUnit(100, "eV").to("Ha"),
                       #nband=8,
                    )

    gs_inp.set_variables(**global_vars)

    # (GS run)
    gs_inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])

    gs_inp.set_variables(
        tolvrs=1e-6,
    )

    print(gs_inp)

    # (NSCF run)
    nscf_inp = abilab.AbiInput(pseudos=pseudos)

    nscf_inp.set_structure(structure)
    nscf_inp.set_variables(**global_vars)

    nscf_inp.set_variables(
        tolwfr=1e-12,
        kptopt=0,
        nkpt=1,
        kpt=[0, 0, 0],
    )

    print(nscf_inp)

    # Build the manager.
    #manager = abilab.TaskManager.simple_mpi(mpi_ncpus=2)

    manager = abilab.TaskManager(qtype="slurm",
       qparams=dict(
           ntasks=2,
           #partition="hmem",
           time="0:20:00",
           #account='nobody@nowhere.org',
           #ntasks_per_node=None,
           #cpus_per_task=None,
       ),
       #setup="SetEnv intel13_intel",
       modules = ["intel/compilerpro/13.0.1.117", "fftw3/intel/3.3"],
       shell_env=dict(
         PATH=("/home/naps/ygillet/NAPS/src/abinit-7.4.3-public/tmp_intel13/src/98_main/:" +
               "/home/naps/ygillet/NAPS/intel13/bin:$PATH"),
         LD_LIBRARY_PATH="/home/naps/ygillet/NAPS/intel13/lib:$LD_LIBRARY_PATH",
       ),
       mpi_runner="mpirun",
       #policy=dict(autoparal=1, max_ncpus=2),
    )

    # Initialize the workflow.
    work = abilab.Workflow(workdir, manager)

    # Register the input.
    gs_link = work.register(gs_inp)

    work.register(nscf_inp, links=gs_link.produces_exts("_DEN"))

    return work


if __name__ == "__main__":
    works = build_workflows()

    for work in works:
        work.build_and_pickle_dump()

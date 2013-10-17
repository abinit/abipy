#!/usr/bin/env python
"""
This example shows how to compute the band structure of a set of 
crystalline structures obtained by changin a set of internal paramaters
"""
from __future__ import division, print_function

import os
import sys 

from abipy import abilab
import abipy.data as data

from abipy.data.runs import enable_logging
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

def make_flow():
    pseudos = data.pseudos("14si.pspnc")

    base_structure = abilab.Structure.from_file(data.cif_file("si.cif"))

    news, uparams = [], [0.1, 0.2, 0.3]

    for u in uparams:
        new = special_positions(base_structure.lattice, u)
        news.append(new)

    workdir = os.path.join(os.path.dirname(__file__), base_structure.formula + "_WYCHOFF")
    manager = abilab.TaskManager.from_file("taskmanager.yaml")

    flow = abilab.AbinitFlow(workdir, manager)

    # Create the list of workflows. Each workflow defines a band structure calculation.
    for new_structure, u in zip(news, uparams):
        # Generate the workflow and register it.
        flow.register_work(make_workflow(new_structure, pseudos))

    return flow.allocate()


def make_workflow(structure, pseudos):
    """
    Return a `Workflow` object defining a band structure calculation
    for given `Structure`.
    """
    # Variables global to the SCF and the NSCF run.
    global_vars = dict(ecut=FloatWithUnit(100, "eV").to("Ha"),
                       #nband=8,
                    )

    # GS + NSCF run 
    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    inp.set_structure(structure)
    inp.set_variables(**global_vars)

    # (GS run)
    inp[1].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])

    inp[1].set_variables(
          tolvrs=1e-6)

    # (NSCF run)
    inp[2].set_variables(
        iscf=-2,
        tolwfr=1e-12,
        kptopt=0,
        nkpt=1,
        kpt=[0, 0, 0],
    )

    gs_inp, nscf_inp = inp.split_datasets()
    return abilab.BandStructureWorkflow(gs_inp, nscf_inp)


@enable_logging
def main():
    flow = make_flow()
    return flow.build_and_pickle_dump()

if __name__ == "__main__":
    sys.exit(main())

    
    
    

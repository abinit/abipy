#!/usr/bin/env python
"""
This example shows how to compute the band structure of a set of
crystalline structures obtained by changing a set of internal paramaters
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import os
import sys
import abipy.data as abidata

from abipy import abilab
from abipy import flowtk


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

    return abilab.Structure(lattice, species, coords, validate_proximity=True, coords_are_cartesian=False)


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_")

    pseudos = abidata.pseudos("14si.pspnc", "8o.pspnc")

    base_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    news, uparams = [], [0.2, 0.3]

    for u in uparams:
        new = special_positions(base_structure.lattice, u)
        news.append(new)

    flow = flowtk.Flow(workdir, manager=options.manager, remove=options.remove)

    # Create the list of workflows. Each workflow defines a band structure calculation.
    for new_structure, u in zip(news, uparams):
        # Generate the workflow and register it.
        flow.register_work(make_workflow(new_structure, pseudos))

    return flow


def make_workflow(structure, pseudos, paral_kgb=1):
    """
    Return a `Workflow` object defining a band structure calculation
    for given `Structure`.
    """
    # Variables global to the SCF and the NSCF run.
    global_vars = dict(
        ecut=abilab.FloatWithUnit(100, "eV").to("Ha"),
        paral_kgb=paral_kgb,
        #nband=8,
    )

    # GS + NSCF run
    multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=2)
    multi.set_vars(global_vars)

    # (GS run)
    multi[0].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    multi[0].set_vars(tolvrs=1e-6)

    # (NSCF run)
    multi[1].set_vars(
        iscf=-2,
        tolwfr=1e-12,
        kptopt=0,
        nkpt=1,
        kpt=[0, 0, 0],
    )

    gs_inp, nscf_inp = multi.split_datasets()

    return flowtk.BandStructureWork(gs_inp, nscf_inp)


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())

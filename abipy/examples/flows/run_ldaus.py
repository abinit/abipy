#!/usr/bin/env python
r"""
Flow for LDA+U calculations
===========================

LDA+U band structure of NiO for several values of U-J.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import numpy as np
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def make_scf_nscf_dos_inputs(structure, pseudos, luj_params, paral_kgb=1):
    # Input file taken from tldau_2.in
    multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=3)

    # Global variables
    global_vars = dict(
        #
        ecut=12,
        pawecutdg=30,
        nband=40,
        occopt=7,
        tsmear=0.015,
        nstep=50,
        paral_kgb=paral_kgb,
        #
        # Spin
        nsppol=1,
        nspden=2,
        nspinor=1,
        spinat=[0,  0,  1,
                0,  0, -1,
                0,  0,  0,
                0,  0,  0],
        # Kpoint Grid
        # The k point grid is not symmetric, but the calculations being
        # for the ground-state, this is not a problem.
    )

    multi.set_vars(global_vars)
    multi.set_vars(luj_params.to_abivars())

    # GS run.
    multi[0].set_vars(
        iscf=17,
        toldfe=1.0e-8,
        ngkpt=[2, 2, 2],
        chksymbreak=0,
    )

    # Band structure run.
    multi[1].set_kpath(ndivsm=6)
    multi[1].set_vars(tolwfr=1e-10)

    # Dos calculation.
    multi[2].set_vars(
        iscf=-3,   # NSCF calculation
        ngkpt=structure.calc_ngkpt(nksmall=8),
        shiftk=[0.0, 0.0, 0.0],
        nshiftk=1,
        tolwfr=1.e-8,
        #pawprtdos=1,
    )

    # Generate two input files for the GS and the NSCF run
    scf_input, nscf_input, dos_input = multi.split_datasets()

    return scf_input, nscf_input, dos_input


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(options.workdir, manager=options.manager)

    # Create the work for the band structure calculation.
    structure = abidata.structure_from_ucell("NiO")
    pseudos = abidata.pseudos("28ni.paw", "8o.2.paw")

    # The code below set up the parameters for the LDA+U calculation in NiO.
    #usepawu   1
    #lpawu   2 -1
    #upawu  8.0 0.0 eV
    #jpawu  0.8 0.0 eV
    usepawu = 1
    u_values = [5.0, 8.0]

    for u in u_values:
        # Apply U-J on Ni only.
        luj_params = abilab.LdauParams(usepawu, structure)
        luj_params.luj_for_symbol("Ni", l=2, u=u, j=0.1*u, unit="eV")

        scf_input, nscf_input, dos_input = make_scf_nscf_dos_inputs(structure, pseudos, luj_params)

        work = flowtk.BandStructureWork(scf_input, nscf_input, dos_inputs=dos_input)
        flow.register_work(work)

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("GENERATE_SPHINX_GALLERY", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).plot_networkx(tight_layout=True)


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())

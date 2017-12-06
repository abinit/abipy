#!/usr/bin/env python
r"""
Gruneisen parameters with DFPT phonons and finite difference
============================================================

The flow computes the dynamical matrix of AlAs for three different volumens
around the equilibrium configuration.
Gruneisen parameters are then computed in anaddb by approximating the derivative
of the dynamical matrix wrt volume with finite differences.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os
import numpy as np
import abipy.abilab as abilab
import abipy.data as abidata
from abipy import flowtk

def build_flow(options):
    """
    Create a `Flow` for phonon calculations:

        1) One workflow for the GS run.
        2) nqpt works for phonon calculations. Each work contains
           nirred tasks where nirred is the number of irreducible phonon perturbations
           for that particular q-point.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=options.workdir)

    # This function constructs the input files for the phonon calculation:
    # GS input + the input files for the phonon calculation.
    pseudos = abidata.pseudos("14si.pspnc")
    structure = abidata.structure_from_ucell("Si")

    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    # List of q-points for the phonon calculation.
    #qpoints = np.reshape(qpoints, (-1, 3))

    # Global variables used both for the GS and the DFPT run.
    gs_inp.set_vars(
        nband=4,
        ecut=2.0,
        ngkpt=[4, 4, 4],
        nshiftk=4,
        shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
                0.0, 0.5, 0.0,
                0.5, 0.0, 0.0,
                0.5, 0.5, 0.5],
        #shiftk=[0, 0, 0],
        paral_kgb=0,
        diemac=12.0,
        iomode=3,
        tolvrs=1.0e-18,
    )

    # Get the qpoints in the IBZ. Note that here we use a q-mesh with ngkpt=(4,4,4) and shiftk=(0,0,0)
    qpoints = gs_inp.abiget_ibz(ngkpt=(4, 4, 4), shiftk=(0, 0, 0)).points

    # Build input for GS calculation and register the first work.
    #work0 = flow.register_scf_task(scf_input)
    #ph_work = flowtk.PhononWork.from_scf_task(work0[0], qpoints=qpoints, is_ngqpt=True)
    #flow.register_work(ph_work)

    ph_inputs = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=len(qpoints))

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

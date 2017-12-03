#!/usr/bin/env python
r"""
Phonon Flow
===========

Phonon band structure of AlAs.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk

def make_scf_input(paral_kgb=0):
    """
    This function constructs the input files for the phonon calculation:
    GS input + the input files for the phonon calculation.
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    structure = abidata.structure_from_ucell("AlAs")
    pseudos = abidata.pseudos("13al.981214.fhi", "33as.pspnc")
    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

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
        paral_kgb=paral_kgb,
        tolvrs=1.0e-10,
        ixc=1,
        nstep=25,
        diemac=9.0,
        iomode=3,
    )

    return gs_inp


def build_flow(options):
    """
    Create a `Flow` for phonon calculations. The flow has two works.

    The first work contains a single GS task that produces the WFK file used in DFPT
    The second work contains multiple PhononTasks that are generated automatically
    in order to compute the dynamical matrix on a [4, 4, 4] mesh.
    Symmetries are taken into account: only q-points in the IBZ are generated and
    for each q-point only the independent atomic perturbations are computed.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=workdir)

    # Build input for GS calculation and register the first work.
    scf_input = make_scf_input()
    work0 = flow.register_scf_task(scf_input)

    # This call uses the information reported in the GS task (work0[0]) to
    # compute all the independent atomic perturbations corresponding to a [4, 4, 4] q-mesh.
    ph_work = flowtk.PhononWork.from_scf_task(work0[0], qpoints=[4, 4, 4], is_ngqpt=True)
    flow.register_work(ph_work)

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("GENERATE_SPHINX_GALLERY", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).plot_networkx()




@flowtk.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())

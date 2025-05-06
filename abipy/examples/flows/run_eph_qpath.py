#!/usr/bin/env python
r"""
e-ph matrix elements along q-path
=================================

This example demonstrates how to compute electron-phonon (e-ph) matrix elements
along a specified q-path, either by explicitly calculating the self-consistent potential perturbations
(Δq Vscf) for each q-point using ab initio methods, o
r by employing Fourier interpolation techniques starting from a coarse q-point mesh.
"""
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk
from abipy.flowtk.eph_flows import EphPotFlow


def build_flow(options):

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Build input for GS calculation
    # Fe BCC structure
    structure = abilab.Structure.from_abistring("""
natom 1
ntypat 1
typat  1
znucl 26

acell 3*5.201879033248
rprim  0.5  0.5  0.5   # BCC primitive vectors (to be scaled by acell)
      -0.5  0.5  0.5
      -0.5 -0.5  0.5
xred 3*0
""")

    pseudos = ["Fe.upf"]
    scf_input = abilab.AbinitInput(structure, pseudos=pseudos)
    num_ele = scf_input.num_valence_electrons

    scf_input.set_vars(
        ecut=50,
        nsppol=2,
        spinat=[0, 0, 6.0],
        nband=12,
        occopt=7,
        tsmear=0.0025,
        tolvrs=6.075e-16,
        nstep=300,
        paral_kgb=0,
        vloc_rcut=10,
        #toldfe1=5e-13,
        #ngfft=[30, 30, 30],
    )

    # K-point grid
    #scf_input.set_kmesh(ngkpt=[24, 24, 24], shiftk=[0, 0, 0])
    scf_input.set_kmesh(ngkpt=[2, 2, 2], shiftk=[0, 0, 0])

    # q-mesh for phonons.
    ngqpt = [2, 2, 2]
    #ngqpt = [1, 1, 1]

    # Lisg of q-points in reduced coordinates.
    qpath_list = [
        +0.10000,  +0.10000,  +0.10000,
        +0.00000,  +0.00000,  +0.00000,
        +0.10000,  +0.00000,  +0.10000,
    ]

    flow = EphPotFlow.from_scf_input(
        options.workdir,
        scf_input,
        ngqpt,
        qpath_list,
        ndivsm=0,                     # pass full list of q-points instead of boundaries.
        what_to_compute="gkq_qpath",  # Important!
        with_becs=False,              # Metal --> no becs, no dynamical quadrupoles.
        with_quad=False,
    )

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).graphviz_imshow()


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

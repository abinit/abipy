#!/usr/bin/env python
r"""
Flow to compute e-ph scattering potentials
==========================================

This example shows how to compute e-ph scattering potentials
along a q-path, merge the POT files in the DVDB file and use the
DVDB and the DDB file to analyze the average over the unit cell of the
periodic part as a function of q
"""
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk


def make_scf_input(ngkpt):
    """
    This function constructs the input file for the GS calculation:
    """
    structure = abilab.Structure.from_abistring(
"""
acell 1.0522E+01  1.0522E+01  1.0522E+01
rprim  0.0  0.5 0.5
       0.5  0.0 0.5
       0.5  0.5 0.0

ntypat 1
znucl  32
natom 2
typat 1 1
xred  0.0  0.0  0.0
      1/4  1/4  1/4
""")

    pseudos = abidata.pseudos("Ge.psp8")
    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    gs_inp.set_vars(
        nband=16,
        ecut=40.0,
        ngkpt=ngkpt,
        nshiftk=1,
        shiftk=[0, 0, 0],
        tolvrs=1.0e-8,
        nstep=150,
        paral_kgb=0,
    )

    return gs_inp


def build_flow(options):
    """
    Create a `Flow` for phonon calculations. The flow has two works.

    The first work contains a single GS task that produces the WFK file used in DFPT
    Then we have multiple Works that are generated automatically
    in order to compute the dynamical matrix on a [2, 2, 2] mesh.
    Symmetries are taken into account: only q-points in the IBZ are generated and
    for each q-point only the independent atomic perturbations are computed.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Use 2x2x2 both for k-mesh and q-mesh
    # Build input for GS calculation
    scf_input = make_scf_input(ngkpt=(8, 8, 8))

    # Create flow to compute all the independent atomic perturbations
    # corresponding to a [4, 4, 4] q-mesh.
    # Electric field and Born effective charges are also computed.
    from abipy.flowtk.eph_flows import EphPotFlow
    ngqpt = [4, 4, 4]

    qpath_list = [
        [+0.000, +0.000, +0.000],  # name: $\Gamma$, weight: 0.000
        [+0.500, +0.000, +0.500],  # name: X, weight: 0.000
        [+0.500, +0.250, +0.750],  # name: W, weight: 0.000
        [+0.375, +0.375, +0.750],  # name: K, weight: 0.000
        [+0.000, +0.000, +0.000],  # name: $\Gamma$, weight: 0.000
        [+0.500, +0.500, +0.500],  # name: L, weight: 0.000
        #[+0.625, +0.250, +0.625],  # name: U, weight: 0.000
        #[+0.500, +0.250, +0.750],  # name: W, weight: 0.000
        #[+0.500, +0.500, +0.500],  # name: L, weight: 0.000
        #[+0.375, +0.375, +0.750],  # name: K, weight: 0.000
        #[+0.625, +0.250, +0.625],  # name: U, weight: 0.000
        #[+0.500, +0.000, +0.500],  # name: X, weight: 0.000
    ]

    # Use small ndivsm to reduce computing time.
    flow = EphPotFlow.from_scf_input(options.workdir, scf_input,
                                     ngqpt, qpath_list, ndivsm=10,
                                     with_quads=True,
                                     with_becs=True)

    return flow


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

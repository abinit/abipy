#!/usr/bin/env python
r"""
Flexoelectric Tensor
=====================

This example shows how to compute the dynamical matrix of GaP on a user-defined q-mesh
including Born effective charges, the macroscopic dielectric matric and the dynamical quadrupoles Q*.
The final results (out_DDB, out_DVDB) will be produced automatically at the end of the run
and saved in ``flow_phonons_with_quad/outdata/``.

The Q* tensor may be needed to improve the accuracy of the Fourier interpolation of the phonon frequencies,
especially in the long-wavelength limit |q| --> 0.
This example is based on  <https://docs.abinit.org/tests/tutorespfn/Input/tlw_4.abi>
Note that only selected features are compatible with dynamical quadrupoles.
Please consult <https://docs.abinit.org/topics/longwave/>

https://docs.abinit.org/tests/tutorespfn/Input/tlw_1.abi
"""

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk


def build_flow(options):
    """
    Create a `Flow` for phonon calculations including
    eps_inf, Born effective charges and  dynamical quadrupoles
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    pseudos = abidata.pseudos("14si.fhi")

    # Initialize the structure from an abinit string
    # Other approaches (e.g. from cif files) are available as well.
    structure = abilab.Structure.from_abistring("""
#Definition of the unit cell
#***************************
acell 10.102 10.6 9.4
rprim   0.0000000000    0.5000000000   0.500000000
        0.5000000000    0.0000000000   0.500000000
        0.5000000000    0.5000000000   0.000000000


#Definition of the atom types and positions
#******************************************
ntypat 1
znucl 14
natom 2
typat 2*1
xred
        0.0500000     -0.0500000      0.0000000
        0.2000000      0.3000000      0.2000000
""")
    # Build input for GS calculation
    scf_input = abilab.AbinitInput(structure, pseudos=pseudos)

    # We use parameters similar to the ones in https://docs.abinit.org/tests/tutorespfn/Input/tlw_4.abi
    scf_input.set_vars(
        nband=4,
        ecut=4.0,
        ngkpt=[4, 4, 4],
        shiftk=[
            [0.5, 0.5, 0.5],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5],
        ],
        diemac=13.0,
        nstep=100,
        #tolvrs=1.0e-10,
        tolvrs=1.0e-18,  # This is the value used in tlw_4.abi
                          # but it is not always possible to reach this precision in more complex systems.
        useylm=1,
        ixc=7
        #iomode=3,
        #paral_kgb=1,
    )

    # At the time of writing, Flexoelectric calculations are implemented only for
    # NC LDA scalar-relativistic pseudos without non-linear core corrections.
    # This section shows how to use the Pseudo API to perform this kind of check
    # before runnnig the calculation.
    for pseudo in scf_input.pseudos:
        # print(pseudo)
        if not pseudo.isnc:
            raise RuntimeError("Only NC pseudos are compatible with Q*")
        if pseudo.has_nlcc:
            raise RuntimeError("NLCC is not compatible with Q*")
        if pseudo.xc.type != "LDA":
            raise RuntimeError("Only LDA is compatible with Q*")

    # Initialize the flow
    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    # Compute phonons on the ddb_ngqpt q-mesh.
    # Include Born effective charges and dynamical quadrupoles via `with_quad=True`.
    #
    ddb_ngqpt = [1, 1, 1]
    #ddb_ngqpt = [4, 4, 4]
    ph_work = flowtk.PhononWork.from_scf_input(scf_input, qpoints=ddb_ngqpt,
                                               is_ngqpt=True, with_becs=True, with_flexoe=True)

    # Add the phonon work to the flow
    flow.register_work(ph_work)

    return flow


# This block generates the thumbnails in the AbiPy gallery.
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

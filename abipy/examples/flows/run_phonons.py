#!/usr/bin/env python
r"""
Flow for phonons with DFPT
==========================

This example shows how to compute the phonon band structure of AlAs with AbiPy flows.
Symmetries are taken into account: only q-points in the IBZ are generated and
for each q-point only the independent atomic perturbations are computed.
The final results (out_DDB, out_DVDB) will be produced automatically at the end of the run
and saved in ``flow_phonons/outdata/``.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk

def make_scf_input(paral_kgb=0):
    """
    This function constructs the input file for the GS calculation:
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
        diemac=9.0,
        #iomode=3,
    )

    return gs_inp


def build_flow(options):
    """
    Create a `Flow` for phonon calculations. The flow has two works.

    The first work contains a single GS task that produces the WFK file used in DFPT
    Then we have multiple Works that are generated automatically
    in order to compute the dynamical matrix on a [4, 4, 4] mesh.
    Symmetries are taken into account: only q-points in the IBZ are generated and
    for each q-point only the independent atomic perturbations are computed.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Build input for GS calculation
    scf_input = make_scf_input()

    # Create flow to compute all the independent atomic perturbations
    # corresponding to a [4, 4, 4] q-mesh.
    # Electric field and Born effective charges are also computed.
    flow = flowtk.PhononFlow.from_scf_input(options.workdir, scf_input,
                                            ph_ngqpt=(4, 4, 4), with_becs=True)

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    #build_flow(options).plot_networkx(with_edge_labels=False, tight_layout=True)
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


############################################################################
#
# Run the script with:
#
#     run_phonons.py -s
#
# then use:
#
#    abirun.py flow_phonons history
#
# to get the list of actions perfomed by AbiPy to complete the flow.
# Note how the ``PhononWork`` has merged all the partial DDB files produced by the PhononTasks
#
# .. code-block:: bash
#
#    ===================================================================================================================================
#    ====================================== <PhononWork, node_id=241274, workdir=flow_phonons/w1> ======================================
#    ===================================================================================================================================
#    [Thu Dec  7 22:55:02 2017] Finalized set to True
#    [Thu Dec  7 22:55:02 2017] Will call mrgddb to merge [ .... ]
#
# Now open the final DDB file with:
#
#    abiopen.py flow_phonons/outdata/out_DDB
#
# and invoke anaddb to compute the phonon band structure and the phonon DOS with:
#
# .. code-block:: ipython
#
#     In [1]: phbst_file, phdos_file = abifile.anaget_phbst_and_phdos_files()
#     In [2]: %matplotlib
#     In [3]: phbst_file.plot_phbands()
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_phonons.png?raw=true
#    :alt: Phonon band structure of AlAs.
#

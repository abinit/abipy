#!/usr/bin/env python
r"""
Phonons with WFQ files (q-mesh denser than k-mesh)
==================================================

This example shows how to use WFQ files to compute phonons on a q-mesh
that is not necessarily commensurate with the k-mesh used for electrons.
Symmetries are taken into account: only q-points in the IBZ are generated.
Moreover WFQ files are computed only if k + q does not belong to the initial mesh and,
for each q-point, only the independent atomic perturbations are computed.
The final results (out_DDB, out_DVDB) will be produced automatically at the end of the run
and saved in the ``outdata/`` of the work.
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
        ngkpt=[2, 2, 2],
        nshiftk=1,
        shiftk=[0, 0, 0],
        #nshiftk=4,
        #shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
        #        0.0, 0.5, 0.0,
        #        0.5, 0.0, 0.0,
        #        0.5, 0.5, 0.5],
        paral_kgb=paral_kgb,
        tolvrs=1.0e-10,
        diemac=9.0,
    )

    return gs_inp


def build_flow(options):
    """
    Create a `Flow` for phonon calculations. The flow has two works.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=options.workdir)

    # Build input for GS calculation and create first work with 1 ScfTask.
    scf_input = make_scf_input()
    work = flow.register_scf_task(scf_input)
    scf_task = work[0]

    # Create work for phonon calculation with WFQ files with a [4, 4, 4] q-mesh.
    # Electric field and Born effective charges are also computed.
    wfkq_work = flowtk.PhononWfkqWork.from_scf_task(scf_task, ngqpt=[4, 4, 4], with_becs=True)
    flow.register_work(wfkq_work)

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
#     run_phonons_wfq.py -s
#
# then use:
#
#    abirun.py flow_phonons_wkq history
#
# to get the list of actions perfomed by AbiPy to complete the flow.
# Note how the ``PhononWfkqWork`` has merged all the partial DDB/DVDB files
# and removed the WFQ files at runtime to optimize the disk space.
#
# .. code-block:: bash
#
#    =========================================================================================================================
#    ============================= <PhononWfkqWork, node_id=360036, workdir=flow_phonons_wkq/w1> =============================
#    =========================================================================================================================
#    [Tue Sep 18 00:04:18 2018] Removing WFQ: flow_phonons_wkq/w1/t5/outdata/out_WFQ
#    [Tue Sep 18 00:04:54 2018] Removing WFQ: flow_phonons_wkq/w1/t14/outdata/out_WFQ
#
# Now open the final DDB file with:
#
#    abiopen.py flow_phonons_wkq/w1/outdata/out_DDB
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

#!/usr/bin/env python
r"""
Flow for e-Bands with frozen phonon
===================================

Electronic band structure of silicon in a distorted geometry (frozen phonon at q=0)
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import numpy as np
import abipy.data as data
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def make_scf_nscf_inputs(structure, paral_kgb=1):
    multi = abilab.MultiDataset(structure, pseudos=data.pseudos("14si.pspnc"), ndtset=2)

    # Global variables
    global_vars = dict(
        ecut=6,
        nband=8,
        timopt=-1,
        paral_kgb=0,
        #nstep=4, # This is not enough to converge. Used to test the automatic restart.
        nstep=10,
        iomode=3,
    )

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    multi[0].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
    multi[1].set_vars(tolwfr=1e-12)

    # Generate two input files for the GS and the NSCF run
    scf_input, nscf_input = multi.split_datasets()

    return scf_input, nscf_input


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # build the structures
    base_structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    modifier = abilab.StructureModifier(base_structure)

    etas = [-0.1, 0, +0.1]
    ph_displ = np.reshape(np.zeros(3*len(base_structure)), (-1,3))
    ph_displ[0,:] = [+1, 0, 0]
    ph_displ[1,:] = [-1, 0, 0]

    displaced_structures = modifier.displace(ph_displ, etas, frac_coords=False)

    flow = flowtk.Flow(options.workdir, manager=options.manager)

    for structure in displaced_structures:
        # Create the work for the band structure calculation.
        scf_input, nscf_input = make_scf_nscf_inputs(structure)

        work = flowtk.BandStructureWork(scf_input, nscf_input)
        flow.register_work(work)

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    #build_flow(options).plot_networkx(with_edge_labels=True, tight_layout=True)
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
#     run_phfrozen_ebands -s
#
# then use:
#
#    abirun.py flow_phfrozen_ebands/ structures -v
#
# to analyze the input/output structures including the atomic positions:
#
# .. code-block:: bash
#
#       Lattice parameters:
#                 formula  natom  angle0  angle1  angle2      a      b      c  volume  \
#       w0_t0_in      Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w0_t0_out     Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w0_t1_in      Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w1_t0_in      Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w1_t0_out     Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w1_t1_in      Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w2_t0_in      Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w2_t0_out     Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#       w2_t1_in      Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#
#                 abispg_num  P [GPa]  Max|F| eV/ang task_class              status
#       w0_t0_in        None      NaN            NaN    ScfTask  Completed
#       w0_t0_out         12   -3.638      3.180e+00    ScfTask  Completed
#       w0_t1_in        None      NaN            NaN   NscfTask  Completed
#       w1_t0_in        None      NaN            NaN    ScfTask  Completed
#       w1_t0_out        227   -5.212      7.430e-27    ScfTask  Completed
#       w1_t1_in        None      NaN            NaN   NscfTask  Completed
#       w2_t0_in        None      NaN            NaN    ScfTask  Completed
#       w2_t0_out         12   -4.192      2.095e+00    ScfTask  Completed
#       w2_t1_in        None      NaN            NaN   NscfTask  Completed
#
#       Atomic positions (columns give the site index):
#                                                    0  \
#       w0_t0_in   (Si, +0.970139 +0.000000 +0.014930)
#       w0_t0_out  (Si, +0.970139 +0.000000 +0.014930)
#       w0_t1_in   (Si, +0.970139 +0.000000 +0.014930)
#       w1_t0_in   (Si, +0.000000 +0.000000 +0.000000)
#       w1_t0_out  (Si, +0.000000 +0.000000 +0.000000)
#       w1_t1_in   (Si, +0.000000 +0.000000 +0.000000)
#       w2_t0_in   (Si, +0.029861 +0.000000 +0.985070)
#       w2_t0_out  (Si, +0.029861 +0.000000 +0.985070)
#       w2_t1_in   (Si, +0.029861 +0.000000 +0.985070)
#
#                                                    1              status
#       w0_t0_in   (Si, +0.279861 +0.250000 +0.235070)  Completed
#       w0_t0_out  (Si, +0.279861 +0.250000 +0.235070)  Completed
#       w0_t1_in   (Si, +0.279861 +0.250000 +0.235070)  Completed
#       w1_t0_in   (Si, +0.250000 +0.250000 +0.250000)  Completed
#       w1_t0_out  (Si, +0.250000 +0.250000 +0.250000)  Completed
#       w1_t1_in   (Si, +0.250000 +0.250000 +0.250000)  Completed
#       w2_t0_in   (Si, +0.220139 +0.250000 +0.264930)  Completed
#       w2_t0_out  (Si, +0.220139 +0.250000 +0.264930)  Completed
#       w2_t1_in   (Si, +0.220139 +0.250000 +0.264930)  Completed
#
#  Finally, we can plot the electronic bands with the command:
#
#    abirun.py flow_phfrozen_ebands ebands -p -t NscfTask
#
#  to select only the band structures produced by the NscfTask.
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_phfrozen_ebands.png?raw=true
#    :alt: Band structures of Si computed for different displacement amplitudes.
#

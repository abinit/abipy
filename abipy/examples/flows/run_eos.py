#!/usr/bin/env python
r"""
Flow for Equation of State
==========================

Flow to compute the equation of state by fitting E(V) at T = 0.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk

exclude_py_versions = ["2.7"]


def build_flow(options):
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Build GS input file.
    pseudos = abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    #silicon = abilab.Structure.zincblende(5.431, ["Si", "Si"], units="ang")
    silicon = abidata.cif_file("si.cif")

    scf_input = abilab.AbinitInput(silicon, pseudos)
    ecut = 12
    scf_input.set_vars(
        ecut=ecut,
        pawecutdg=40,
        nband=6,
        paral_kgb=0,
        iomode=3,
        toldfe=1e-9,
    )

    # K-point sampling (shifted)
    scf_input.set_autokmesh(nksmall=4)

    from abipy.flowtk.gs_works import EosWork
    flow = flowtk.Flow(options.workdir, manager=options.manager)

    # Si is cubic and atomic positions are fixed by symmetry so we
    # use move_atoms=False to compute E(V) with SCF-GS tasks instead of
    # performing a constant-volume optimization of the cell geometry.
    work = EosWork.from_scf_input(scf_input, move_atoms=False, ecutsm=0.5)
    flow.register_work(work)
    flow.allocate(use_smartio=True)

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
#     run_eos.py -s
#
# then use:
#
#    abirun.py flow_eos structures
#
# to get info about the initial and final structures.
#
# .. code-block:: bash
#
#        Lattice parameters:
#                  formula  natom  angle0  angle1  angle2      a      b      c  volume  \
#        w0_t0_in      Si2      2    60.0    60.0    60.0  3.854  3.854  3.854  40.479
#        w0_t0_out     Si2      2    60.0    60.0    60.0  3.854  3.854  3.854  40.479
#        w0_t1_in      Si2      2    60.0    60.0    60.0  3.857  3.857  3.857  40.582
#        w0_t1_out     Si2      2    60.0    60.0    60.0  3.857  3.857  3.857  40.582
#        w0_t2_in      Si2      2    60.0    60.0    60.0  3.861  3.861  3.861  40.684
#        w0_t2_out     Si2      2    60.0    60.0    60.0  3.861  3.861  3.861  40.684
#        w0_t3_in      Si2      2    60.0    60.0    60.0  3.864  3.864  3.864  40.786
#        w0_t3_out     Si2      2    60.0    60.0    60.0  3.864  3.864  3.864  40.786
#        w0_t4_in      Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#        w0_t4_out     Si2      2    60.0    60.0    60.0  3.867  3.867  3.867  40.888
#        w0_t5_in      Si2      2    60.0    60.0    60.0  3.870  3.870  3.870  40.991
#        w0_t5_out     Si2      2    60.0    60.0    60.0  3.870  3.870  3.870  40.991
#        w0_t6_in      Si2      2    60.0    60.0    60.0  3.873  3.873  3.873  41.093
#        w0_t6_out     Si2      2    60.0    60.0    60.0  3.873  3.873  3.873  41.093
#        w0_t7_in      Si2      2    60.0    60.0    60.0  3.877  3.877  3.877  41.195
#        w0_t7_out     Si2      2    60.0    60.0    60.0  3.877  3.877  3.877  41.195
#        w0_t8_in      Si2      2    60.0    60.0    60.0  3.880  3.880  3.880  41.297
#        w0_t8_out     Si2      2    60.0    60.0    60.0  3.880  3.880  3.880  41.297
#
#                  abispg_num  P [GPa]  Max|F| eV/ang task_class              status
#        w0_t0_in        None      NaN            NaN    ScfTask  Completed
#        w0_t0_out        227    1.327      9.506e-28    ScfTask  Completed
#        w0_t1_in        None      NaN            NaN    ScfTask  Completed
#        w0_t1_out        227    1.092      1.877e-27    ScfTask  Completed
#        w0_t2_in        None      NaN            NaN    ScfTask  Completed
#        w0_t2_out        227    0.859      4.649e-27    ScfTask  Completed
#        w0_t3_in        None      NaN            NaN    ScfTask  Completed
#        w0_t3_out        227    0.629      5.062e-27    ScfTask  Completed
#        w0_t4_in        None      NaN            NaN    ScfTask  Completed
#        w0_t4_out        227    0.403      5.945e-27    ScfTask  Completed
#        w0_t5_in        None      NaN            NaN    ScfTask  Completed
#        w0_t5_out        227    0.178      6.306e-27    ScfTask  Completed
#        w0_t6_in        None      NaN            NaN    ScfTask  Completed
#        w0_t6_out        227   -0.042      1.738e-27    ScfTask  Completed
#        w0_t7_in        None      NaN            NaN    ScfTask  Completed
#        w0_t7_out        227   -0.259      2.725e-27    ScfTask  Completed
#        w0_t8_in        None      NaN            NaN    ScfTask  Completed
#        w0_t8_out        227   -0.474      6.939e-27    ScfTask  Completed
#
#        Use `--verbose` to print atoms.
#
# Silicon is a cubic materials and the positions are fixed by symmeetry so we only need ScfTask to get E(V).
#
# The EOS results are saved in the JSON file in flow_eos/w0/outdata/eos_data.json.
#
# We can also use the GSR robot to compute and plot the EOS.
#
# Use
#
#    abirun.py flow_eos robot GSR
#
# to build a robot for GSR files and start an ipython shell.
# Then, inside ipython, type:
#
# .. code-block:: ipython
#
#        In [1]: %matplotlib
#        In [2]: robot.gridplot_eos()
#
# to produce
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_eos.png?raw=true
#    :alt: Equation of state of Si obtained with different models.
#

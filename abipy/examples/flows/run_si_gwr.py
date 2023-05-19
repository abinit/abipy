#!/usr/bin/env python
r"""
GWR flow with convergence studies
=================================

This script shows how to compute the G0W0 corrections in silicon.
More specifically, we build a flow to analyze the convergence of the QP corrections
wrt to the number of bands in the self-energy. More complicated convergence studies
can be implemented on the basis of this example.
"""

import os
import sys
import abipy.data as data
import abipy.abilab as abilab

from abipy import flowtk


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_","flow_")

    scf_input = abilab.AbinitInput(structure=data.cif_file("si.cif"),
                                   pseudos=data.pseudos("14si.pspnc"))

    # Global variables. gw_kmesh is used in all datasets except DATASET 1.
    scf_input.set_vars(
        ecut=6,
        tolvrs=1e-6,
        nband=4,
        paral_kgb=1,
    )

    # IMPORTANT: k-grid for GWR must be Gamma-centered.
    scf_input.set_kmesh(
        ngkpt=[2, 2, 2],
        shiftk=[0.0, 0.0, 0.0],
    )

    # GS-SCF to get the DEN, followed by direct diago to obtain green_nband bands.
    from abipy.flowtk.gwr_works import DirectDiagoWork

    flow = flowtk.Flow(workdir=options.workdir)

    green_nband = -1
    diago_work = DirectDiagoWork.from_scf_input(scf_input, green_nband)
    flow.register_work(diago_work)

    gwr_template = scf_input.make_gwr_qprange_input(gwr_ntau=6, nband=8, ecuteps=4)

    # Two possibilities:
    # 1) Change the value of one variable:

    varname_values = ("nband", [8, 12, 14])

    # or take the Cartesian product of two or more variables with e.g.:

    #varname_values = [
    #   ("gwr_ntau", [6, 8]),
    #   ("ecuteps", [2, 4]),
    #]

    from abipy.flowtk.gwr_works import GwrSigmaConvWork
    gwr_work = GwrSigmaConvWork.from_varname_values(
            varname_values, gwr_template, den_node=diago_work[0], wfk_node=diago_work[1])

    flow.register_work(gwr_work)

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
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())


############################################################################
#
# Run the script with:
#
#     run_si_gwr.py -s
#
# The last three tasks (``w0_t3``, ``w0_t4``, ``w0_t5``) are the SigmaTask who have produced
# a netcdf file with the GW results with different number of bands.
# We can check this with the command:
#
#    abirun.py flow_si_g0w0/ listext SIGRES
#
# .. code-block:: bash
#
#       Found 3 files with extension `SIGRES` produced by the flow
#       File                                        Size [Mb]    Node_ID  Node Class
#       ----------------------------------------  -----------  ---------  ------------
#       flow_si_g0w0/w0/t3/outdata/out_SIGRES.nc         0.05     241325  SigmaTask
#       flow_si_g0w0/w0/t4/outdata/out_SIGRES.nc         0.08     241326  SigmaTask
#       flow_si_g0w0/w0/t5/outdata/out_SIGRES.nc         0.13     241327  SigmaTask
#
# Let's use the SIGRES robot to collect and analyze the results:
#
#    abirun.py flow_si_g0w0/ robot SIGRES
#
# and then, inside the ipython terminal, type:
#
# .. code-block:: ipython
#
#       In [1]: df = robot.get_dataframe()
#       In [2]: df
#       Out[2]:
#                                                 nsppol     qpgap            ecutwfn  \
#       flow_si_g0w0/w0/t3/outdata/out_SIGRES.nc       1  3.627960  5.914381651684836
#       flow_si_g0w0/w0/t4/outdata/out_SIGRES.nc       1  3.531781  5.914381651684836
#       flow_si_g0w0/w0/t5/outdata/out_SIGRES.nc       1  3.512285  5.914381651684836
#
#                                                            ecuteps  \
#       flow_si_g0w0/w0/t3/outdata/out_SIGRES.nc  3.6964885323070074
#       flow_si_g0w0/w0/t4/outdata/out_SIGRES.nc  3.6964885323070074
#       flow_si_g0w0/w0/t5/outdata/out_SIGRES.nc  3.6964885323070074
#
#                                                          ecutsigx scr_nband  \
#       flow_si_g0w0/w0/t3/outdata/out_SIGRES.nc  5.914381651684846        25
#       flow_si_g0w0/w0/t4/outdata/out_SIGRES.nc  5.914381651684846        25
#       flow_si_g0w0/w0/t5/outdata/out_SIGRES.nc  5.914381651684846        25
#
#                                                sigma_nband gwcalctyp scissor_ene  \
#       flow_si_g0w0/w0/t3/outdata/out_SIGRES.nc          10         0         0.0
#       flow_si_g0w0/w0/t4/outdata/out_SIGRES.nc          20         0         0.0
#       flow_si_g0w0/w0/t5/outdata/out_SIGRES.nc          30         0         0.0
#
#                                                 nkibz
#       flow_si_g0w0/w0/t3/outdata/out_SIGRES.nc      6
#       flow_si_g0w0/w0/t4/outdata/out_SIGRES.nc      6
#       flow_si_g0w0/w0/t5/outdata/out_SIGRES.nc      6
#
#       In [3]: %matplotlib
#       In [4]: df.plot("sigma_nband", "qpgap", marker="o")
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_si_g0w0.png?raw=true
#    :alt: QP results in Si plotted vs the KS energy e0.
#

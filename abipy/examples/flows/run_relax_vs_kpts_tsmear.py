#!/usr/bin/env python
r"""
Relaxation of GaN with different K-meshes
=========================================

In this example, we employ the relaxation algorithms implemented in Abinit (``ionmov`` and ``optcell``)
to find the equilibrium configuration of GaN (atomic positions and lattice vectors).
The relaxation is done with different k-meshes to monitor the convergence of the results.
You will observe a change of the equilibrium parameters with respect to the k-point mesh.

Note the we are using pseudopotentials generated with the GGA which tends to
overestimate the lattice parameters and ecut is way too low.
If you replace GGA with LDA, you will observe that LDA tends to underestimate the parameters.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os

import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata


def relax_input(tsmear, nksmall):
    """
    Crystalline aluminum: optimization of the lattice parameter
    at fixed number of k points and broadening. Similar to tbase4_1.in with minor
    """
    #structure = abilab.Structure.fcc()
    inp = abilab.AbinitInput(structure=abidata.ucells.structure_from_ucell("Al"),
                             pseudos=abidata.pseudos("13al.981214.fhi"))

    # Define k-point sampling.
    # nshiftk and shift are automatically selected from the lattice and the number of divisions
    # for the smallest direction. nksmall 2 e.g. will automatically select
    #   ngkpt 2 2 2
    #   nshiftk 4
    #   shiftk
    #       0.5 0.5 0.5
    #       0.5 0.0 0.0
    #       0.0 0.5 0.0
    #       0.0 0.0 0.5
    inp.set_autokmesh(nksmall=nksmall)

    inp.set_vars(
        ecut=6,
        occopt=4,
        tsmear=tsmear,
        toldfe=1e-6,
        nstep=10,
        optcell=1,    # Optimization of the lattice parameters
        ionmov=2,
        ntime=10,
        dilatmx=1.05,
        ecutsm=0.5,
        ixc=1,
    )

    return inp

def build_flow(options):
    """
    Build and return a flow performing structural relaxations with different k-point samplings.
    """
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Let generate multiple inputs for different (tsmear, nksmall)
    # Product computes the Cartesian product of input iterables.
    # It's equivalent to nested for-loops
    tsmear_list = (0.01, 0.02, 0.03, 0.04)
    nksmall_list = (2, 4, 6)

    from itertools import product
    inputs = [relax_input(tsmear, nksmall) for tsmear, nksmall in product(tsmear_list, nksmall_list)]

    # Build flow form inputs.
    # As the calculations are independent, we can use Flow.from_inputs
    # Note the Flow.from_inputs is a simplified interface that, by default, builds tasks
    # for Ground-state calculation (GsTask).
    # Here we are performing a structural relaxation so we have to specify the task class explicitly.
    # AbiPy will use this piece of information to handle the restart of the RelaxTask that differs
    # from the one provided by GsTask.

    return flowtk.Flow.from_inputs(options.workdir, inputs=inputs, task_class=flowtk.RelaxTask)


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
# Run the script with:
#
#     run_relax_vs_kpts_tsmear.py -s
#
# then use:
#
#     abirun.py flow_relax_vs_kpts_tsmear hist -p
#
# to print (and plot) the structural relaxation for all the tasks:
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_relax_vs_kpts_tsmear_hist.png?raw=true
#    :alt: Commbiplot of the HIST files for the different (nkpt, tsmear) params
#
# To analyze the convergence of the relaxed lattice parameters, use:
#
# .. code-block:: bash
#
#	abirun.py flow_relax_vs_kpts_tsmear robot GSR
#
# to create a GSR robot for all the tasks in the flow and open an ipyton shell.
#
# Then, inside ipython, type:
#
# .. code-block:: ipython
#
#	 In [1]: %matplotlib
#	 In [2]: df = robot.get_dataframe()
#	 # Let's do some math with pandas to retrieve the Abinit acell from the a lattice parameter given in Ang.
#	 In [3]: import math
#	 In [4]: from abipy import abilab
#	 In [5]: df["acell"] = df["a"] * math.sqrt(2) * abilab.units.ang_to_bohr

#        In [7]: df["acell"]
#        Out[7]:
#        flow_relax_vs_kpts_tsmear/w0/t0/outdata/out_GSR.nc     7.558770
#        flow_relax_vs_kpts_tsmear/w0/t1/outdata/out_GSR.nc     7.505486
#        flow_relax_vs_kpts_tsmear/w0/t2/outdata/out_GSR.nc     7.496158
#        flow_relax_vs_kpts_tsmear/w0/t3/outdata/out_GSR.nc     7.558770
#        flow_relax_vs_kpts_tsmear/w0/t4/outdata/out_GSR.nc     7.505643
#        flow_relax_vs_kpts_tsmear/w0/t5/outdata/out_GSR.nc     7.495546
#        flow_relax_vs_kpts_tsmear/w0/t6/outdata/out_GSR.nc     7.558770
#        flow_relax_vs_kpts_tsmear/w0/t7/outdata/out_GSR.nc     7.501756
#        flow_relax_vs_kpts_tsmear/w0/t8/outdata/out_GSR.nc     7.496770
#        flow_relax_vs_kpts_tsmear/w0/t9/outdata/out_GSR.nc     7.558771
#        flow_relax_vs_kpts_tsmear/w0/t10/outdata/out_GSR.nc    7.504096
#        flow_relax_vs_kpts_tsmear/w0/t11/outdata/out_GSR.nc    7.499134
#        Name: acell, dtype: float64
#
#        # to plot the optimized acell vs nkpt for the different values of tsmear, use:
#	 In [6]: robot.plot_xy_with_hue(df, "nkpt", "acell", hue="tsmear")
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_relax_vs_kpts_tsmear.png?raw=true
#    :alt: optimized acell as function of nkpt and tsmear

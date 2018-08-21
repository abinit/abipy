#!/usr/bin/env python
r"""
Bethe-Salpeter Flow with factory functions
==========================================

Calculation of the BSE spectrum with the High-throuhput interface (factory functions).
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy import abilab


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize pseudos and Structure.
    pseudos = abidata.pseudos("14si.pspnc")
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    kppa = scf_kppa = 1
    nscf_nband = 6
    nscf_ngkpt = [4, 4, 4]
    nscf_shiftk = [0.1, 0.2, 0.3]
    bs_loband = 2
    bs_nband = nscf_nband
    mbpt_sciss = 0.7 * abilab.units.eV_to_Ha
    mdf_epsinf = 12
    ecuteps = 2
    ecut = 12

    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    # BSE calculation with model dielectric function.
    multi = abilab.bse_with_mdf_inputs(
        structure, pseudos,
        scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
        ecuteps, bs_loband, bs_nband, mbpt_sciss, mdf_epsinf,
        ecut=ecut,# pawecutdg=None,
        exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="unpolarized",
        smearing=None)
        #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

    work = flowtk.BseMdfWork(scf_input=multi[0], nscf_input=multi[1], bse_inputs=multi[2:])

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
#     run_ht_si_bsemdf.py -s
#
# then use:
#
#    abirun.py flow_ht_si_bsemdf/ status -v
#
# to get the status of the flow
#
# .. code-block:: bash
#
#        Work #0: <BseMdfWork, node_id=241906, workdir=flow_ht_si_bsemdf/w0>, Finalized=True
#        +--------+-----------+-----------------+--------------+------------+----------+-----------------+----------+-----------+
#        | Task   | Status    | Queue           | MPI|Omp|Gb   | Warn|Com   | Class    | Sub|Rest|Corr   | Time     |   Node_ID |
#        +========+===========+=================+==============+============+==========+=================+==========+===========+
#        | w0_t0  | Completed | 62430@localhost | 1|  1|2.0    | 2|  0      | ScfTask  | (1, 0, 0)       | 0:00:01R |    241907 |
#        +--------+-----------+-----------------+--------------+------------+----------+-----------------+----------+-----------+
#        | w0_t1  | Completed | 62438@localhost | 2|  1|2.0    | 3|  0      | NscfTask | (1, 0, 0)       | 0:00:05R |    241908 |
#        +--------+-----------+-----------------+--------------+------------+----------+-----------------+----------+-----------+
#        | w0_t2  | Completed | 62449@localhost | 2|  1|2.0    | 10|  6     | BseTask  | (1, 0, 0)       | 0:00:17R |    241909 |
#        +--------+-----------+-----------------+--------------+------------+----------+-----------------+----------+-----------+
#
#        all_ok reached
#
# The macroscopic dielectric function produced by the BseTask is stored in the ``out_MDF.nc`` file:
#
# .. code-block:: bash
#
#        abirun.py flow_ht_si_bsemdf/ listext MDF
#
#        Found 1 files with extension `MDF` produced by the flow
#        File                                          Size [Mb]    Node_ID  Node Class
#        ------------------------------------------  -----------  ---------  ------------
#        flow_ht_si_bsemdf/w0/t2/outdata/out_MDF.nc          0.4     241909  BseTask
#
# Open the file with:
#
#       abiopen.py flow_ht_si_bsemdf/w0/t2/outdata/out_MDF.nc
#
# and then inside the ipython terminal, issue:
#
# .. code-block:: ipython
#
#    In [1]: %matplotlib
#    In [2]: abifile.plot_mdfs()
#
# to produce:
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_ht_si_bsemdf.png?raw=true
#    :alt: Band structure of Si in the IBZ and along a k-path
#

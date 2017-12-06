#!/usr/bin/env python
r"""
Bethe-Salpeter Flow with factory functions
==========================================

Calculation of the BSE spectrum with the HT interface.
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
    nscf_ngkpt = [4,4,4]
    nscf_shiftk = [0.1, 0.2, 0.3]
    bs_loband = 2
    bs_nband = nscf_nband
    mbpt_sciss = 0.7
    mdf_epsinf = 12
    ecuteps = 2
    ecut = 12

    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    # BSE calculation with model dielectric function.
    multi = abilab.bse_with_mdf_inputs(
        structure, pseudos,
        scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
        ecuteps, bs_loband, bs_nband, mbpt_sciss, mdf_epsinf,
        ecut=ecut,#$ pawecutdg=None,
        exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="unpolarized",
        smearing=None)
        #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

    work = flowtk.BseMdfWork(scf_input=multi[0], nscf_input=multi[1], bse_inputs=multi[2:])

    #from pymatgen.io.abinit.calculations import bse_with_mdf_work
    #work = bse_with_mdf_work(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
    #                         ecuteps, bs_loband, bs_nband, mbpt_sciss, mdf_epsinf,
    #                         accuracy="normal", spin_mode="unpolarized", smearing=None,
    #                         charge=0.0, scf_solver=None, **extra_abivars)

    flow.register_work(work)
    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("GENERATE_SPHINX_GALLERY", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).plot_networkx(tight_layout=True)



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

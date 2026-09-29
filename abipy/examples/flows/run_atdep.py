#!/usr/bin/env python
r"""
Temperature-dependent phonons with aTDEP
========================================

This example shows how to build a flow that computes the temperature-dependent
interatomic force constants of MgO with aTDEP, starting from an AIMD trajectory (HIST.nc file).
The DDB file produced by aTDEP is then passed to anaddb to compute the phonon band structure and DOS.
"""

import os
import sys

import abipy.data as abidata
from abipy import flowtk
from abipy.abio.inputs import AnaddbInput, AtdepInput


def make_atdep_input(structure):
    """Build the input file for aTDEP."""
    inp = AtdepInput(structure)
    inp.set_vars(
        multiplicity=[[-2, 2, 2], [2, -2, 2], [2, 2, -2]],
        nstep_max=20,
        nstep_min=1,
        temperature=300,
        rcut=7.0,
        ngqpt1=[2, 2, 2],
        ngqpt2=[1, 1, 1],
    )

    return inp


def make_anaddb_input(structure):
    """Build the anaddb input file to compute phonon bands and DOS from the aTDEP DDB."""
    inp = AnaddbInput.phbands_and_dos(
        structure, ngqpt=[2, 2, 2], ndivsm=40, line_density=None,
        nqsmall=10, qppa=None, q1shft=(0, 0, 0), qptbounds=None,
        asr=2, chneut=1, dipdip=1, dipquad=0, quadquad=0,
        dos_method="tetra", lo_to_splitting="automatic",
        with_ifc=True, anaddb_kwargs={}, spell_check=False,
    )

    return inp


def build_flow(options):
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    structure = abidata.structure_from_mpid("mp-1265")  # MgO
    ddb_file = abidata.ref_file("MgO_dfpt_aimd/MgO_DDB.nc")
    hist_file = abidata.ref_file("MgO_dfpt_aimd/MgO_HIST.nc")

    atdep_inp = make_atdep_input(structure)
    anaddb_inp = make_anaddb_input(structure)
    anaddb_inp.set_vars(
        ngqpt=atdep_inp["ngqpt1"],
        rifcsph=atdep_inp["rcut"],
    )

    atdep_task = flowtk.AtdepTask(atdep_inp, hist_node=hist_file, ddb_node=ddb_file)
    anaddb_task = flowtk.AnaddbTask(anaddb_inp, ddb_node=atdep_task)

    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)
    flow.register_task(atdep_task)
    flow.register_task(anaddb_task, deps={atdep_task: "DDB"})

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

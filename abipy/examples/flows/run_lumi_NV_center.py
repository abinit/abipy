#!/usr/bin/env python
r"""
Delta SCF constrained occupation method calculation, to determine luminescent properties.
=========================================================================================

This example shows how to compute the luminescent properties of the NV- center in diamond.
It uses a 64 atoms supercell, where one C atom was replaced by one N atom and one vacancy was created.
See Fig.3 of https://doi.org/10.1103/PhysRevB.104.045303 for the setting of electron occupation
in the ground/excited state.
Steps:
1) Relaxation in the ground state
2) Relaxation in the excited state, starting from the relaxed ground state. Created at run-time
3) Scf computation in the relaxed/unrelaxed ground/excited state (4 computations).

Even if we use minimal settings, the workflow takes a few minutes to run on one core.
Filepaths of the 6 runs are stored in /w0/outdata/lumi.json
A quick post-processing is automatically done at the end of a LumiWork
and stored in /w0/outdata/Delta_SCF.json, with relevant luminescent properties
(ZPL energy, Stoke Shift, \Delta Q,...), see abipy/lumi/delta_scf.py
"""

import sys
import os
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata
from abipy.core.structure import Structure
from abipy.flowtk.lumi_works import LumiWork


def scf_inp(structure):
    pseudos = abidata.pseudos("N.psp8", "C.psp8")

    gs_scf_inp = abilab.AbinitInput(structure=structure, pseudos=pseudos)
    gs_scf_inp.set_vars(ecut=10, ### too low, just for example!
                        chksymbreak=0,
                        diemac=5,
                        prtwf=0,
                        nstep=300,
                        toldfe=1e-10,
                        chkprim=0,
                        cellcharge=-1 ### Negatively charged NV center
                    )

    ### Setting of the occupations, for spin up-dn in the ground/excited state
    ### Only valid for NV center in this particular cell.
    n_val = gs_scf_inp.num_valence_electrons
    n_cond = round(10)

    spin_up_gs = f"\n{int((n_val - 3) / 2)}*1 1 1   1 {n_cond}*0"
    spin_up_ex = f"\n{int((n_val - 3) / 2)}*1 1 1   1 {n_cond}*0"
    spin_dn_gs = f"\n{int((n_val - 3) / 2)}*1 1 0   0 {n_cond}*0"
    spin_dn_ex = f"\n{int((n_val - 3) / 2)}*1 0 0.5 0.5 {n_cond}*0"

    nsppol = 2

    #Dealing with supercell, Gamma only calculation
    shiftk = [0, 0, 0]
    ngkpt = [1, 1, 1]

    # Build SCF input for the ground state configuration.
    gs_scf_inp.set_kmesh_nband_and_occ(ngkpt, shiftk, nsppol, [spin_up_gs, spin_dn_gs])
    # Build SCF input for the excited configuration.
    exc_scf_inp = gs_scf_inp.deepcopy()
    exc_scf_inp.set_kmesh_nband_and_occ(ngkpt, shiftk, nsppol, [spin_up_ex, spin_dn_ex])

    return gs_scf_inp,exc_scf_inp


def relax_kwargs():

    # Dictionary with input variables to be added for performing structural relaxations.
    relax_kwargs = dict(
        ecutsm=0.5,
        toldff=1e-5, # TOO HIGH, just for testing purposes.
        tolmxf=1e-4, # TOO HIGH, just for testing purposes.
        ionmov=2,
        chkdilatmx=0,
    )
    # Relaxation settings could be different between excited and ground state...
    relax_kwargs_gs=relax_kwargs.copy()
    relax_kwargs_gs['optcell']=0 # in the ground state, no relaxation of the cell
    # Could be different!

    relax_kwargs_ex=relax_kwargs.copy()
    relax_kwargs_ex['optcell']=0 # in the excited state, no relaxation of the cell

    return relax_kwargs_gs, relax_kwargs_ex


def build_flow(options):

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(options.workdir, manager=options.manager)

    #Construct the structure

    stru=Structure.from_file(abidata.cif_file("NV_center_64_at_sc.cif"))

    ####### Delta SCF part of the flow #######

    gs_scf_inp,exc_scf_inp = scf_inp(stru)

    relax_kwargs_gs, relax_kwargs_ex = relax_kwargs()
    lumi_work=LumiWork.from_scf_inputs(gs_scf_inp, exc_scf_inp, relax_kwargs_gs, relax_kwargs_ex,four_points=True)

    flow.register_work(lumi_work)

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


if __name__ == '__main__':
    sys.exit(main())

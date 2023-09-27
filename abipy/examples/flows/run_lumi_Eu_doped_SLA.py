#!/usr/bin/env python
r"""
Delta SCF constrained occupation method calculation, to determine luminescent properties.
=========================================================================================

This example shows how to compute the luminescent properties of Eu doped phosphor.
It uses a 36 atoms cell of SrLiAl3N4. Two non-equivalent Sr sites are availabe for Eu, resulting
in two independent LumiWork. The creation of the supercells is done with make_doped_supercell().

Steps, for each structure:
1) Relaxation in the ground state
2) Relaxation in the excited state, starting from the relaxed ground state. Created at run-time
3) Scf computation in the relaxed/unrelaxed ground/excited state (4 computations).

Even if we use minimal settings, the workflow takes around one hour to run on one core.
Filepaths of the 6 runs are stored in outdata/lumi.json of each work
A quick post-processing is automatically done at the end of a LumiWork
and stored in outdata/Delta_SCF.json of each work, with relevant luminescent properties
(ZPL energy, Stoke Shift, \Delta Q,...), see abipy/lumi/delta_scf.py
"""

import sys
import os
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata
from abipy.core.structure import Structure
from abipy.flowtk.lumi_works import LumiWork

def get_non_eq_sites(structure,replaced_atom):
    ### return a list of positions of non-equivalent sites for the replaced atom. ###
    irred=structure.spget_equivalent_atoms().eqmap # mapping from inequivalent sites to atoms sites
    positions=structure.get_symbol2indices()[replaced_atom] # get indices of the replaced atom

    index_different_sites=[]

    for i in positions:
        if len(irred[i]) != 0:
            index_different_sites.append(irred[i][0])

    return(index_different_sites)


def make_doped_supercell(prim_structure,supercell_size,replaced_atom,dopant_atom):
    #return a list of doped supercell structure, one for each non-equivalent site of the replaced atom
    my_structure=prim_structure.copy()
    my_structure.make_supercell(supercell_size)

    list_ineq_pos=get_non_eq_sites(my_structure,replaced_atom)

    doped_structure_list=[]

    for pos in list_ineq_pos:
        final_structure=my_structure.copy()
        final_structure.replace(pos,dopant_atom)
        doped_structure_list.append(final_structure)

    return doped_structure_list


def scf_inp(structure):
    pseudos = abidata.pseudos("Eu.xml", "Sr.xml","Al.xml","N.xml","Li.xml")

    gs_scf_inp = abilab.AbinitInput(structure=structure, pseudos=pseudos)
    gs_scf_inp.set_vars(ecut=10,
                        pawecutdg=20,
                        chksymbreak=0,
                        diemac=5,
                        prtwf=0,
                        nstep=300,
                        toldfe=1e-10,
                        chkprim=0,
                        nbdbuf=5 # help convergence
                    )


    # Set DFT+U and spinat parameters according to chemical symbols.
    #symb2spinat = {"Eu": [0, 0, 7]}
    #symb2luj = {"Eu": {"lpawu": 3, "upawu": 7, "jpawu": 0.7}}

    #gs_scf_inp.set_usepawu(usepawu=1, symb2luj=symb2luj)
    #gs_scf_inp.set_spinat_from_symbols(symb2spinat, default=(0, 0, 0))

    ### Setting of the occupations, for spin up-dn in the ground/excited state
    ### Only valid for Eu doped
    n_val = gs_scf_inp.num_valence_electrons
    n_cond = round(15)

    spin_up_gs = f"\n{int((n_val - 7) / 2)}*1 7*1 {n_cond}*0"
    spin_up_ex = f"\n{int((n_val - 7) / 2)}*1 6*1 0 1 {n_cond - 1}*0"
    spin_dn = f"\n{int((n_val - 7) / 2)}*1 7*0 {n_cond}*0"

    nsppol = 2
    shiftk = [0, 0, 0]
    ngkpt = [1, 1, 1]

    # Build SCF input for the ground state configuration.
    gs_scf_inp.set_kmesh_nband_and_occ(ngkpt, shiftk, nsppol, [spin_up_gs, spin_dn])
    # Build SCF input for the excited configuration.
    exc_scf_inp = gs_scf_inp.deepcopy()
    exc_scf_inp.set_kmesh_nband_and_occ(ngkpt, shiftk, nsppol, [spin_up_ex, spin_dn])

    return gs_scf_inp,exc_scf_inp



def relax_kwargs():

    # Dictionary with input variables to be added for performing structural relaxations.
    relax_kwargs = dict(
        ecutsm=0.5,
        toldff=1e-4, # TOO HIGH, just for testing purposes.
        tolmxf=1e-3, # TOO HIGH, just for testing purposes.
        ionmov=2,
        chkdilatmx=0,
    )

    relax_kwargs_gs=relax_kwargs.copy()
    relax_kwargs_gs['optcell']=0 # in the ground state, allow relaxation of the cell

    relax_kwargs_ex=relax_kwargs.copy()
    relax_kwargs_ex['optcell']=0 # in the excited state, no relaxation of the cell

    return relax_kwargs_gs, relax_kwargs_ex


def build_flow(options):

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(options.workdir, manager=options.manager)

    #Construct the two structures (2 non-eq. sites for Sr) from the primitive cell of SLA (SrAlLi3N4)

    #prim_structure=structure.Structure.from_file('SLA_prim.cif')
    prim_structure=Structure.from_file(abidata.cif_file("SLA_prim.cif"))
    supercell_matrix=[1,1,1]  # Too small, just for test
    strus=prim_structure.make_doped_supercells(supercell_matrix,'Sr','Eu')


    ####### Delta SCF part of the flow #######

    # Create one LumiWork per structure
    for stru in strus:
        gs_scf_inp, exc_scf_inp = scf_inp(stru)
        relax_kwargs_gs, relax_kwargs_ex = relax_kwargs()
        lumi_work=LumiWork.from_scf_inputs(gs_scf_inp, exc_scf_inp, relax_kwargs_gs, relax_kwargs_ex)
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


#!/usr/bin/env python
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk

#def make_scf_input(scf_ngkpt, paral_kgb=0):
#    """
#    This function constructs the input file for the GS calculation:
#    on a scf_ngkpt Gamma-centered k-mesh.
#    """
#    structure = abilab.Structure.from_abivars(
#        acell=[7.7030079150, 7.7030079150, 7.7030079150],
#        rprim=[0.0000000000, 0.5000000000, 0.5000000000,
#               0.5000000000, 0.0000000000, 0.5000000000,
#               0.5000000000, 0.5000000000, 0.0000000000],
#        natom=2,
#        ntypat=2,
#        typat=[1, 2],
#        znucl=[3, 9],
#        xred=[0.0000000000, 0.0000000000, 0.0000000000,
#              0.5000000000, 0.5000000000, 0.5000000000]
#    )
#
#    pseudos = ["03-Li.LDA.TM.pspnc", "09-F.LDA.TM.pspnc"]
#
#    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)
#
#    gs_inp.set_vars(
#        ecut=50.0,
#        ngkpt=scf_ngkpt,
#        shiftk=[0, 0, 0],
#        nband=10,
#        paral_kgb=paral_kgb,
#        tolvrs=1.0e-10,
#        prtpot=1, # This is required for the Sternheimer method in the EPH part.
#    )
#
#    return gs_inp

def make_scf_input((scf_ngkpt, paral_kgb=0):
    """
    This function constructs the input file for the GS calculation:
    """

    # Initialize MgO structure from abinit variables.
    structure = abilab.Structure.from_abivars(
        acell=3 * [4.252718 * abilab.units.ang_to_bohr],
        rprim=[0.0000000000, 0.5000000000, 0.5000000000,
               0.5000000000, 0.0000000000, 0.5000000000,
               0.5000000000, 0.5000000000, 0.0000000000],
        natom=2,
        ntypat=2,
        typat=[1, 2],
        znucl=[12, 8],
        xred=[0.0000000000, 0.0000000000, 0.0000000000,
              0.5000000000, 0.5000000000, 0.5000000000]
    )

    # Input for GS part.
    # NC pseudos assumed in currect working directory.
    #pseudos = ["Mg-sp-gw.psp8", "O.psp8"]
    pseudos = abidata.pseudos("12mg.pspnc", "O.psp8")

    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    gs_inp.set_vars(
        nband=12,
        paral_kgb=paral_kgb,
        ecut=35.0,        # Too low. Should be ~50
        ngkpt=scf_ngkpt,  # Too coarse
        nshiftk=1,        # Gamma-centered mesh. Important to have the CBM/VBM!
        shiftk=[0, 0, 0],
        tolvrs=1.0e-10,
        diemac=9.0,
        nstep=150,
        nbdbuf=4,
        prtpot=1,        # Print potential for Sternheimer
        iomode=3,        # Produce output files in netcdf format.
    )

    return gs_inp

def build_flow(options):
    """
    Create a `Flow` for ZPR calculations.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=options.workdir)

    # Build input for the GS calculation and register the first GS ScfTask.
    gs_inp = make_scf_input(scf_ngkpt=[4, 4, 4])
    work = flow.register_scf_task(gs_inp)
    scf_task = work[0]

    # Build input for NSCF calculation along k-path (automatically selected by AbiPy)
    nscf_kpath_inp = gs_inp.make_ebands_input()
    work.register_nscf_task(nscf_kpath_inp, deps={scf_task: "DEN"})

    # This is the ab-initio q-mesh that should be COMMENSURATE with scf_ngkpt
    ddb_ngqpt = [2, 2, 2]

    # Now we build the NSCF tasks with different k-meshes and empty states.
    # Each NSCF task generates one of the WFK files used to build the EPH self-energy.
    # These WFK files are then used to erform convergece tests with respect to 
    # the interpolated fine q-mesh. 
    # NOTE that in the EPH we are gonna use k_mesh = q_mesh
    ngkpt_fine_list = [
        [8, 8, 8],
        #[12, 12, 12],
        #[32, 32, 32],
    ]

    nscf_empty_states_tasks = []
    for ngkpt_fine in ngkpt_fine_list:
        nscf_empty_kmesh_inp = gs_inp.new_with_vars(
            ngkpt=ngkpt_fine,
            #nband=630,      # Too low. ~300
            #nbdbuf=30,      # Reduces considerably the time needed to converge empty states!
            nband=40,      # Too low. ~300
            nbdbuf=5,      # Reduces considerably the time needed to converge empty states!
            tolwfr=1e-18,
            iscf=-2,
        )
        t = work.register_nscf_task(nscf_empty_kmesh_inp, deps={scf_task: "DEN"})
        nscf_empty_states_tasks.append(t)

    # Create work for phonon calculation on the coarse ddb_ngqpt q-mesh.
    # Electric field and Born effective charges are computed.
    ph_work = flowtk.PhononWfkqWork.from_scf_task(scf_task, ngqpt=ddb_ngqpt, with_becs=True)
    #for task in ph_work:
    #   task.input.set_vars(prtwf=-1)
    flow.register_work(ph_work)

    # Build template for e-ph self-energy calculation (real + imag part)
    # The k-points must be in the WFK file
    eph_template = gs_inp.new_with_vars(
        optdriver=7,             # Enter EPH driver.
        eph_task=4,              # Activate computation of EPH self-energy.
        ddb_ngqpt=ddb_ngqpt,     # Ab-initio q-mesh used to produce the DDB file.
        #nkptgw=2,
        #kptgw=[0, 0, 0,
        #       0.5, 5, 0],
        #bdgw=[1, 8, 1, 8],
        #gw_qprange=1,
        tmesh=[0, 200, 1],    # (start, step, num)
        zcut="0.01 eV",
        mixprec=1,
        boxcutmin=1.1,
    )

    # Set q-path for Fourier interpolation of phonons.
    eph_template.set_qpath(10)
    # Set q-mesh for phonons DOS.
    eph_template.set_phdos_qmesh(nqsmall=16, method="tetra")

    # Now we use the EPH template to perform a convergence study in which
    # we change the q-mesh used to integrate the self-energy and the number of bands.
    # The code will activate the Fourier interpolation of the DFPT potentials 
    # if eph_ngqpt_fine != ddb_ngqpt

    # Create empty work to contain EPH tasks with this value of eph_ngqpt_fine
    eph_work = flow.new_work()

    for ngkpt_fine, nscf_task in zip(ngkpt_fine_list, nscf_empty_states_tasks):
        new_inp = eph_template.new_with_vars(
            ngkpt=ngkpt_fine,
            eph_ngqpt_fine=eph_ngqpt_fine
        ) #, nband=nband)

        # The EPH code requires the GS WFK, the DDB file with all perturbations
        # and the DVDB file with the DFPT potentials (already merged by ph_work)
        deps = {nscf_task: "WFK", ph_work: ["DDB", "DVDB"]}
        eph_work.register_eph_task(new_inp, deps=deps)

    #flow.allocate()
    return flow


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
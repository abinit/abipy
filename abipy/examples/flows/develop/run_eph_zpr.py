#!/usr/bin/env python
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk

def make_scf_input(paral_kgb=0):
    """
    This function constructs the input file for the GS calculation:
    """
    structure = abilab.Structure.from_abivars(
        acell=[7.7030079150, 7.7030079150, 7.7030079150],
        rprim=[0.0000000000, 0.5000000000, 0.5000000000,
               0.5000000000, 0.0000000000, 0.5000000000,
               0.5000000000, 0.5000000000, 0.0000000000],
        natom=2,
        ntypat=2,
        typat=[1, 2],
        znucl=[3, 9],
        xred=[0.0000000000, 0.0000000000, 0.0000000000,
              0.5000000000, 0.5000000000, 0.5000000000]
    )

    pseudos = ["03-Li.LDA.TM.pspnc", "09-F.LDA.TM.pspnc"]

    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    gs_inp.set_vars(
        nband=10,
        ecut=50.0,
        istwfk="*1",
        ngkpt=[8, 8, 8],
        nshiftk=1,
        shiftk=[0, 0, 0],
        paral_kgb=paral_kgb,
        tolvrs=1.0e-12,
        diemac=9.0,
        nstep=150,
        nbdbuf=4,
    )

    return gs_inp


def build_flow(options):
    """
    Create a `Flow` for phonon calculations. The flow has two works.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=options.workdir)

    # Build input for GS calculation and create first work with 1 ScfTask.
    gs_inp = make_scf_input()
    work = flow.register_scf_task(gs_inp)
    scf_task = work[0]


    # Build new input for NSCF calculation with k-path (automatically selected by AbiPy)
    # Used to plot the KS band structure and interpolate the QP corrections.
    nscf_kpath_inp = gs_inp.new_with_vars(
        nband=12,
        tolwfr=1e-20,
        iscf=-2,
    )
    nscf_kpath_inp.set_kpath(ndivsm=10)
    work.register_nscf_task(nscf_kpath_inp, deps={scf_task: "DEN"})

    # Build another NSCF input with k-mesh and empty states.
    # This step generates the WFK file used to build the EPH self-energy.
    nscf_empty_kmesh_inp = gs_inp.new_with_vars(
        nband=630,     # Too low. ~300
        nbdbuf=30,     # Reduces considerably the time needed to converge empty states!
        tolwfr=1e-18,
        iscf=-2,
    )
    nscf_empty_task = work.register_nscf_task(nscf_empty_kmesh_inp, deps={scf_task: "DEN"})


    NGQPT = [32, 32, 32]
    nscf_empty_kmesh_inp.set_kmesh(
        ngkpt=NGQPT,
        shiftk=[0.0, 0.0, 0.0],
    )

    # Create work for phonon calculation with WFQ files with a [4, 4, 4] q-mesh.
    # Electric field and Born effective charges are also computed.
    ph_work = flowtk.PhononWfkqWork.from_scf_task(scf_task, ngqpt=NGQPT, with_becs=True)

    #for task in ph_work:
    #   task.input.set_vars(prtwf=-1)

    flow.register_work(ph_work)

    # Build template for self-energy calculation. See also v8/Input/t44.in
    # The k-points must be in the WFK file
    #
    eph_inp = gs_inp.new_with_vars(
        optdriver=7,             # Enter EPH driver.
        eph_task=4,              # Activate computation of EPH self-energy.
        ngkpt=NGQPT,
        ddb_ngqpt=NGQPT,         # q-mesh used to produce the DDB file (must be consistent with DDB data)
        symsigma=1,              # Use symmetries in self-energy integration (IBZ_k instead of BZ)
        # For more k-points...
        nkptgw=2,
        kptgw=[0, 0, 0,
               0.5, 5, 0],
        bdgw=[1, 8, 1, 8],
        #gw_qprange=-4,
        tmesh=[0, 200, 1],    # (start, step, num)
        zcut="0.1 eV",
    )

    # Set q-path for Fourier interpolation of phonons.
    eph_inp.set_qpath(10)

    # Set q-mesh for phonons DOS.
    eph_inp.set_phdos_qmesh(nqsmall=16, method="tetra")

    # EPH part requires the GS WFK, the DDB file with all perturbations
    # and the database of DFPT potentials (already merged by PhononWork)
    deps = {nscf_empty_task: "WFK", ph_work: ["DDB", "DVDB"]}

    # Now we use the EPH template to perform a convergence study in which
    # we change the q-mesh used to integrate the self-energy and the number of bands.
    # The code will activate the Fourier interpolation of the DFPT potentials if eph_ngqpt_fine != ddb_ngqpt

    #for eph_ngqpt_fine in [[4, 4, 4], [8, 8, 8]]:
    for eph_ngqpt_fine in [NGQPT]:
        # Create empty work to contain EPH tasks with this value of eph_ngqpt_fine
        eph_work = flow.register_work(flowtk.Work())
        for nband in [200, 300, 400, 500]:
            new_inp = eph_inp.new_with_vars(eph_ngqpt_fine=eph_ngqpt_fine, nband=nband)
            eph_work.register_eph_task(new_inp, deps=deps)

    flow.allocate()
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

#!/usr/bin/env python
r"""
Computation of Phonons, BECS and Eps_inf for MgO
==================================================

This example shows how to compute phonons with DFPT
Symmetries are taken into account: only q-points in the IBZ are generated.
The final results (out_DDB, out_DVDB) will be produced automatically at the end of the run
and saved in the ``outdata/`` directory of work[1].
"""

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk


def make_scf_input():
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

    # NC pseudos assumed in currect working directory.
    pseudos = ["Mg-sp-gw.psp8", "O.psp8"]

    # Input for GS part.
    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    gs_inp.set_vars(
        nband=16,
        paral_kgb=0,
        ecut=35.0,        # Too low. Should be ~50
        ngkpt=[4, 4, 4],  # Too coarse
        #nshiftk=1,
        shiftk=[0, 0, 0], # Gamma-centered mesh. Important to have the CBM/VBM!
        tolvrs=1.0e-10,
        diemac=9.0,
        nstep=150,
        nbdbuf=4,
        #prtpot=1,       # Print potential for Sternheimer
        iomode=3,        # Produce output files in netcdf format.
    )

    return gs_inp


def build_flow(options):
    """
    Create a `Flow` for phonon calculations. The flow has two works.

    - work[0]: GS + NSCF along a k-path
    - work[1]: DFPT work with phonons on a 4x4x4, BECS and eps_inf
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=options.workdir)

    # Build input for GS calculation and create first work with 1 ScfTask.
    gs_inp = make_scf_input()
    work = flow.register_scf_task(gs_inp)
    scf_task = work[0]

    # Build new input for NSCF calculation with k-path (automatically selected by AbiPy)
    # Used to plot the KS band structure
    nscf_kpath_inp = gs_inp.new_with_vars(
        nband=16,
        tolwfr=1e-18,
        iscf=-2,
    )
    nscf_kpath_inp.set_kpath(ndivsm=10)
    work.register_nscf_task(nscf_kpath_inp, deps={scf_task: "DEN"})

    # NSCF run with eph_ngkpt k-mesh to produce the WFK file for the EPH code.
    nscf_eph_inp = gs_inp.new_with_vars(
        tolwfr=1e-16,
        iscf=-2,
        prtpot=1,         # Needed for the Sternheimer.
                          # This POT file is actually a copy of the GS POT file.
    )

    eph_ngkpt = [32, 32, 32]
    nscf_eph_inp.set_kmesh(
        ngkpt=eph_ngkpt,
        shiftk=[0.0, 0.0, 0.0],
    )
    work.register_nscf_task(nscf_eph_inp, deps={scf_task: "DEN"})
    nscf_eph_work = work[-1]

    # Q-mesh for phonons. In this case k-mesh == q-mesh
    #ngkpt = [32, 32, 32]
    ngqpt = [4, 4, 4]

    # Create work for phonon calculationwith a [4, 4, 4] q-mesh.
    # Electric field and Born effective charges are also computed.
    ph_work = flowtk.PhononWork.from_scf_task(scf_task, ngqpt, is_ngqpt=True, with_becs=True, with_quad=False)

    flow.register_work(ph_work)

    # Build template for self-energy calculation. See also v8/Input/t44.in
    # The k-points must be in the WFK file
    #
    eph_inp = gs_inp.new_with_vars(
        optdriver=7,             # Enter EPH driver.
        eph_task=4,              # Activate computation of EPH self-energy.
        eph_stern=1,             # Use the sternheimer equations
        ngkpt=eph_ngkpt,
        ddb_ngqpt=ngqpt,         # q-mesh used to produce the DDB file (must be consistent with DDB data)
        symsigma=1,              # Use symmetries in self-energy integration (IBZ_k instead of BZ)
        # For more k-points...
        nkptgw=1,
        kptgw=[0, 0, 0],
        bdgw=[9, 9],
        #gw_qprange=-4,
        nfreqsp=8000,
        freqspmax="8.0 eV",
        tmesh=[0, 200, 1],    # (start, step, num)
        zcut="0.01 eV",
    )

    # Set q-path for Fourier interpolation of phonons.
    eph_inp.set_qpath(10)

    # Set q-mesh for phonons DOS.
    eph_inp.set_phdos_qmesh(nqsmall=16, method="tetra")

    # Convergence of the q-mesh grid
    # For each q-mesh check the frequency mesh and thus the time mesh density to obtain a stable spectral function
    for eph_ngqpt_fine in [ [8, 8, 8], [16, 16, 16], [32, 32, 32]]:
        # Create empty work to contain EPH tasks with this value of eph_ngqpt_fine
        eph_work = flow.register_work(flowtk.Work())

        for nfreq in [1000, 2000, 3000, 4000]:

            # EPH part requires the GS WFK, the DDB file with all perturbations
            # and the database of DFPT potentials (already merged by PhononWork)
            deps = {nscf_eph_work: ["WFK", "POT"], ph_work: ["DDB", "DVDB"]}
            new_inp = eph_inp.new_with_vars(eph_ngqpt_fine=eph_ngqpt_fine, nfreqsp=nfreq)
            eph_work.register_eph_task(new_inp, deps=deps)

            ce_inp = new_inp.new_with_vars(
                eph_task=9,
                tolcum=1e-3
            )
            # Cumulant expansion requires the SIGEPH output file obtained from the EPH calculation
            deps = {nscf_eph_work: ["WFK", "POT"], ph_work: ["DDB", "DVDB"], eph_work[-1]: "SIGEPH"}

            eph_work.register_eph_task(ce_inp, deps=deps)

    # Create workflow for the EPH task
    #eph_work.register_eph_task(eph_inp, deps=deps)

    #ce_inp = eph_inp.new_with_vars(
    #        eph_task=9,
    #        tolcum=1e-3
    #        )

    #deps = {nscf_eph_work: ["WFK","POT"], ph_work: ["DDB", "DVDB"], eph_work[0]: "SIGEPH"}

    #eph_work.register_eph_task(ce_inp, deps=deps)

    flow.allocate()
    flow.use_smartio()

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



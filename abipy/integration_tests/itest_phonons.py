"""Integration tests for phonon flows."""

import os
import numpy as np
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk
from abipy.core.testing import has_matplotlib

import logging
logger = logging.getLogger(__name__)


def scf_ph_inputs(tvars=None):
    """
    This function constructs the input files for the phonon calculation:
    GS input + the input files for the phonon calculation.
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    structure = abidata.structure_from_ucell("AlAs")

    # List of q-points for the phonon calculation (4,4,4) mesh.
    qpoints = [
             0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
             2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
             5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
             2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  2.50000000E-01,  0.00000000E+00,
            -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
            -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
            ]

    qpoints = np.reshape(qpoints, (-1, 3))

    # Global variables used both for the GS and the DFPT run.
    global_vars = dict(nband=4,
                       ecut=3.0,
                       ngkpt=[4, 4, 4],
                       shiftk=[0, 0, 0],
                       tolvrs=1.0e-6,
                       paral_kgb=0 if tvars is None else tvars.paral_kgb,
                    )

    multi = abilab.MultiDataset(structure=structure, pseudos=abidata.pseudos("13al.981214.fhi", "33as.pspnc"),
                               ndtset=1 + len(qpoints))

    multi.set_vars(global_vars)

    for i, qpt in enumerate(qpoints):
        # Response-function calculation for phonons.
        multi[i+1].set_vars(
            nstep=20,
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            kptopt=3,
            )

            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2      # Automatic generation of k points, taking

    # Split input into gs_inp and ph_inputs
    return multi.split_datasets()


def itest_phonon_flow(fwp, tvars):
    """
    Create a flow for phonon calculations:

        1) One work for the GS run.

        2) nqpt works for phonon calculations. Each work contains
           nirred tasks where nirred is the number of irreducible phonon perturbations
           for that particular q-point.
    """
    all_inps = scf_ph_inputs(tvars)
    scf_input, ph_inputs = all_inps[0], all_inps[1:]

    ph_ngqpt = (2, 2, 2)
    flow = flowtk.PhononFlow.from_scf_input(fwp.workdir, scf_input, ph_ngqpt=ph_ngqpt, with_becs=False)

    scheduler = flow.make_scheduler()
    assert scheduler.start() == 0
    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    with flow.open_final_ddb() as ddb:
        ddb_path = ddb.filepath
        ddb.to_string(verbose=2)
        assert len(ddb.structure) == 2
        #assert ddb.qpoints.frac_coords

    # Test PhononTask inspect method
    ph_task = flow[1][0]
    if has_matplotlib():
        assert ph_task.inspect(show=False)
    # Test get_results
    #ph_task.get_results()

    # Build new work with Anaddb tasks.
    # Construct a manager with mpi_procs==1 since anaddb do not support mpi_procs > 1 (except in elphon)
    shell_manager = fwp.manager.to_shell_manager(mpi_procs=1)
    awork = flowtk.Work(manager=shell_manager)

    # Phonons bands and DOS with gaussian method
    anaddb_input = abilab.AnaddbInput.phbands_and_dos(
        scf_input.structure, ngqpt=ph_ngqpt, ndivsm=5, nqsmall=10, dos_method="gaussian: 0.001 eV")

    atask = flowtk.AnaddbTask(anaddb_input, ddb_node=ddb_path, manager=shell_manager)
    awork.register(atask)

    # Phonons bands and DOS with tetrahedron method
    anaddb_input = abilab.AnaddbInput.phbands_and_dos(
        scf_input.structure, ngqpt=ph_ngqpt, ndivsm=5, nqsmall=10, dos_method="tetra")

    atask = flowtk.AnaddbTask(anaddb_input, ddb_node=ddb_path, manager=shell_manager)
    awork.register(atask)

    flow.register_work(awork)
    flow.allocate()
    flow.build()

    for i, atask in enumerate(awork):
        atask.history.info("About to run anaddb task: %d", i)
        atask.start_and_wait()
        assert atask.status == atask.S_DONE
        atask.check_status()
        assert atask.status == atask.S_OK

        # These output files should be produced in the task workdir.
        # Actually they should be in the outdir but anaddb uses different conventions.
        assert len(atask.wdir.list_filepaths(wildcard="*PHBST.nc")) == 1
        assert len(atask.wdir.list_filepaths(wildcard="*PHDOS.nc")) == 1


def itest_phonon_restart(fwp):
    """Test the restart of phonon calculations with the scheduler."""
    # Crystalline AlAs: computation of the second derivative of the total energy
    structure = abidata.structure_from_ucell("AlAs")

    # List of q-points for the phonon calculation (4,4,4) mesh.
    qpoints = np.reshape([
        0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
        2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
        #5.00000000E-01,  0.00000000E+00,  0.00000000E+00,  # XXX Uncomment this line to test restart from 1DEN
                                                            # Too long --> disabled
    ], (-1, 3))

    # Global variables used both for the GS and the DFPT run.
    global_vars = dict(
        nband=4,
        ecut=3.0,
        ngkpt=[4, 4, 4],
        shiftk=[0, 0, 0],
        tolvrs=1.0e-5,
    )

    multi = abilab.MultiDataset(structure=structure,
                                pseudos=abidata.pseudos("13al.981214.fhi", "33as.pspnc"),
                                ndtset=1 + len(qpoints))

    multi.set_vars(global_vars)

    for i, qpt in enumerate(qpoints):
        # Response-function calculation for phonons.
        multi[i+1].set_vars(
            rfphon=1,        # Will consider phonon-type perturbation.
            nqpt=1,          # One wavevector is to be considered.
            qpt=qpt,         # q-wavevector.
            kptopt=3,
            nstep=5,         # This is to trigger the phonon restart.
        )
        #rfatpol   1 1   # Only the first atom is displaced
        #rfdir   1 0 0   # Along the first reduced coordinate axis
        #kptopt   2      # Automatic generation of k points, taking

        # i == 0 --> restart from WFK
        if i == 1: multi[i+1].set_vars(prtwf=-1, nstep=5)  # Restart with WFK and smart-io.
        if i == 2: multi[i+1].set_vars(prtwf=0, nstep=8)   # Restart from 1DEN. Too long --> disabled.

    all_inps = multi.split_datasets()
    scf_input, ph_inputs = all_inps[0], all_inps[1:]

    flow = phonon_flow(fwp.workdir, scf_input, ph_inputs, manager=fwp.manager)
    flow.set_garbage_collector()

    for task in flow.iflat_tasks():
        task.manager.qadapter.max_num_launches = 20

    assert flow.make_scheduler().start() == 0

    flow.check_status(show=True, verbose=1)
    assert all(work.finalized for work in flow)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    assert sum(task.num_restarts for task in flow.iflat_tasks()) > 0


def phonon_flow(workdir, scf_input, ph_inputs, with_nscf=False, with_ddk=False, with_dde=False,
                manager=None, flow_class=flowtk.PhononFlow, allocate=True):
    """
    Build a :class:`PhononFlow` for phonon calculations.
    Args:
        workdir: Working directory.
        scf_input: Input for the GS SCF run.
        ph_inputs: List of Inputs for the phonon runs.
        with_nscf: add an nscf task in front of al phonon tasks to make sure the q point is covered
        with_ddk: add the ddk step
        with_dde: add the dde step it the dde is set ddk is switched on automatically
        manager: :class:`TaskManager` used to submit the jobs
                 Initialized from manager.yml if manager is None.
        flow_class: Flow class
    Returns:
        :class:`Flow` object
    """
    logger.critical("phonon_flow is deprecated and could give wrong results")
    if with_dde:
        with_ddk = True

    natom = len(scf_input.structure)

    # Create the container that will manage the different works.
    flow = flow_class(workdir, manager=manager)

    # Register the first work (GS calculation)
    # register_task creates a work for the task, registers it to the flow and returns the work
    # the 0the element of the work is the task
    scf_task = flow.register_task(scf_input, task_class=flowtk.ScfTask)[0]

    # Build a temporary work with a shell manager just to run
    # ABINIT to get the list of irreducible pertubations for this q-point.
    shell_manager = flow.manager.to_shell_manager(mpi_procs=1)

    if with_ddk:
        logger.info('add ddk')
        # TODO
        # MG Warning: be careful here because one should use tolde or tolwfr (tolvrs shall not be used!)
        ddk_input = ph_inputs[0].deepcopy()
        ddk_input.set_vars(qpt=[0, 0, 0], rfddk=1, rfelfd=2, rfdir=[1, 1, 1])
        ddk_task = flow.register_task(ddk_input, deps={scf_task: 'WFK'}, task_class=flowtk.DdkTask)[0]

    if with_dde:
        logger.info('add dde')
        dde_input = ph_inputs[0].deepcopy()
        dde_input.set_vars(qpt=[0, 0, 0], rfddk=1, rfelfd=2)
        dde_input_idir = dde_input.deepcopy()
        dde_input_idir.set_vars(rfdir=[1, 1, 1])
        dde_task = flow.register_task(dde_input, deps={scf_task: 'WFK', ddk_task: 'DDK'}, task_class=flowtk.DdeTask)[0]

    if not isinstance(ph_inputs, (list, tuple)):
        ph_inputs = [ph_inputs]

    for i, ph_input in enumerate(ph_inputs):
        fake_input = ph_input.deepcopy()

        # Run abinit on the front-end to get the list of irreducible pertubations.
        tmp_dir = os.path.join(workdir, "__ph_run" + str(i) + "__")
        w = flowtk.PhononWork(workdir=tmp_dir, manager=shell_manager)
        fake_task = w.register(fake_input)

        # Use the magic value paral_rf = -1 to get the list of irreducible perturbations for this q-point.
        abivars = dict(
            paral_rf=-1,
            rfatpol=[1, natom],  # Set of atoms to displace.
            rfdir=[1, 1, 1],     # Along this set of reduced coordinate axis.
        )

        fake_task.set_vars(abivars)
        w.allocate()
        w.start(wait=True)

        # Parse the file to get the perturbations.
        try:
            irred_perts = flowtk.yaml_read_irred_perts(fake_task.log_file.path)
        except Exception:
            print("Error in %s" % fake_task.log_file.path)
            raise

        logger.info(irred_perts)

        w.rmtree()

        # Now we can build the final list of works:
        # One work per q-point, each work computes all
        # the irreducible perturbations for a singe q-point.

        work_qpt = flowtk.PhononWork()

        if with_nscf:
            # MG: Warning this code assume 0 is Gamma!
            import copy
            nscf_input = copy.deepcopy(scf_input)
            nscf_input.set_vars(kptopt=3, iscf=-3, qpt=irred_perts[0]['qpt'], nqpt=1)
            nscf_task = work_qpt.register_nscf_task(nscf_input, deps={scf_task: "DEN"})
            deps = {nscf_task: "WFQ", scf_task: "WFK"}
        else:
            deps = {scf_task: "WFK"}

        if with_ddk:
            deps[ddk_task] = 'DDK'

        logger.info(irred_perts[0]['qpt'])

        for irred_pert in irred_perts:
            #print(irred_pert)
            new_input = ph_input.deepcopy()

            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            qpt = irred_pert["qpt"]
            idir = irred_pert["idir"]
            ipert = irred_pert["ipert"]

            # TODO this will work for phonons, but not for the other types of perturbations.
            rfdir = 3 * [0]
            rfdir[idir - 1] = 1
            rfatpol = [ipert, ipert]

            new_input.set_vars(
                #rfpert=1,
                qpt=qpt,
                rfdir=rfdir,
                rfatpol=rfatpol,
            )

            if with_ddk:
                new_input.set_vars(rfelfd=3)

            work_qpt.register_phonon_task(new_input, deps=deps)

        flow.register_work(work_qpt)

    if allocate: flow.allocate()

    return flow

"""Optical spectra with Optic."""

import pytest
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def make_inputs(tvars):
    """Construct input files."""
    structure = abidata.structure_from_ucell("GaAs")

    multi = abilab.MultiDataset(structure, pseudos=abidata.pseudos("31ga.pspnc", "33as.pspnc"), ndtset=5)

    # Global variables
    kmesh = dict(ngkpt=[4, 4, 4],
                 nshiftk=4,
                 shiftk=[[0.5, 0.5, 0.5],
                         [0.5, 0.0, 0.0],
                         [0.0, 0.5, 0.0],
                         [0.0, 0.0, 0.5]])

    global_vars = dict(ecut=2, paral_kgb=tvars.paral_kgb)
    global_vars.update(kmesh)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_vars(
        tolvrs=1e-6,
        nband=4)

    # NSCF run with large number of bands, and points in the the full BZ
    multi[1].set_vars(
        iscf=-2,
        nband=20,
        nstep=25,
        kptopt=1,
        tolwfr=1.e-8)
        #kptopt=3)

    # Fourth dataset: ddk response function along axis 1
    # Fifth dataset: ddk response function along axis 2
    # Sixth dataset: ddk response function along axis 3
    for idir in range(3):
        rfdir = 3 * [0]
        rfdir[idir] = 1

        multi[2+idir].set_vars(
            iscf=-3,
            nband=20,
            nstep=1,
            nline=0,
            prtwf=3,
            kptopt=3,
            nqpt=1,
            qpt=[0.0, 0.0, 0.0],
            rfdir=rfdir,
            rfelfd=2,
            tolwfr=1.e-9,
        )

    # scf_inp, nscf_inp, ddk1, ddk2, ddk3
    return multi.split_datasets()


def itest_optic_flow(fwp, tvars):
    """Test optic calculations."""
    if tvars.paral_kgb == 1:
        pytest.xfail("Optic flow with paral_kgb==1 is expected to fail (implementation problem)")

    """
    0.002         ! Value of the smearing factor, in Hartree
    0.0003  0.3   ! Difference between frequency values (in Hartree), and maximum frequency ( 1 Ha is about 27.211 eV)
    0.000         ! Scissor shift if needed, in Hartree
    0.002         ! Tolerance on closeness of singularities (in Hartree)
    1             ! Number of components of linear optic tensor to be computed
    11            ! Linear coefficients to be computed (x=1, y=2, z=3)
    2             ! Number of components of nonlinear optic tensor to be computed
    123 222       ! Non-linear coefficients to be computed
    """
    optic_input = abilab.OpticInput(
        broadening=0.002,
        domega=0.0003,
        maxomega=0.3,
        scissor=0.000,
        tolerance=0.002,
        num_lin_comp=1,
        lin_comp=11,
        num_nonlin_comp=2,
        nonlin_comp=(123, 222),
    )
    #print(optic_input)

    scf_inp, nscf_inp, ddk1, ddk2, ddk3 = make_inputs(tvars)

    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

    bands_work = flowtk.BandStructureWork(scf_inp, nscf_inp)
    flow.register_work(bands_work)

    # work with DDK tasks.
    ddk_work = flowtk.Work()
    for inp in [ddk1, ddk2, ddk3]:
        ddk_work.register_ddk_task(inp, deps={bands_work.nscf_task: "WFK"})

    flow.register_work(ddk_work)
    flow.allocate()
    flow.build_and_pickle_dump(abivalidate=True)

    # Run the tasks
    for task in flow.iflat_tasks():
        task.start_and_wait()
        assert task.status == task.S_DONE

    flow.check_status(show=True)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    # Optic does not support MPI with ncores > 1 hence we have to construct a manager with mpi_procs==1
    shell_manager = fwp.manager.to_shell_manager(mpi_procs=1)

    # Build optic task and register it
    optic_task1 = flowtk.OpticTask(optic_input, nscf_node=bands_work.nscf_task, ddk_nodes=ddk_work,
                                   manager=shell_manager)

    flow.register_task(optic_task1)
    flow.allocate()
    flow.build_and_pickle_dump(abivalidate=True)

    optic_task1.start_and_wait()
    assert optic_task1.status == optic_task1.S_DONE

    # Now we do a similar calculation but the dependencies are represented by
    # strings with the path to the input files instead of task objects.
    ddk_nodes = [task.outdir.has_abiext("1WF") for task in ddk_work]
    #ddk_nodes = [task.outdir.has_abiext("DDK") for task in ddk_work]
    #print("ddk_nodes:", ddk_nodes)
    assert all(ddk_nodes)

    #nscf_node = bands_work.nscf_task
    nscf_node = bands_work.nscf_task.outdir.has_abiext("WFK")
    assert nscf_node

    # This does not work yet
    optic_task2 = flowtk.OpticTask(optic_input, nscf_node=nscf_node, ddk_nodes=ddk_nodes)
    flow.register_task(optic_task2)
    flow.allocate()
    flow.build_and_pickle_dump(abivalidate=True)
    assert len(flow) == 4

    optic_task2.start_and_wait()
    assert optic_task2.status == optic_task2.S_DONE

    flow.check_status(show=True)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    assert all(work.finalized for work in flow)

    #assert flow.validate_json_schema()

    # Test get_results
    optic_task2.get_results()

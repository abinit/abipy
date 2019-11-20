"""Integration tests for GW flows."""

import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk

#from abipy.core.testing import has_abinit, has_matplotlib


def make_g0w0_inputs(ngkpt, tvars):
    """
    Input files for the calculation of the GW corrections.

    Returns: gs_input, nscf_input, scr_input, sigma_input
    """
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=4)

    # This grid is the most economical, but does not contain the Gamma point.
    scf_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap.
    gw_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.0, 0.0, 0.0,
                0.0, 0.5, 0.5,
                0.5, 0.0, 0.5,
                0.5, 0.5, 0.0]
    )

    # Global variables. gw_kmesh is used in all datasets except DATASET 1.
    ecut = 4

    multi.set_vars(
        ecut=ecut,
        pawecutdg=ecut*2 if multi.ispaw else None,
        istwfk="*1",
        paral_kgb=tvars.paral_kgb,
        gwpara=2,
    )
    multi.set_kmesh(**gw_kmesh)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(**scf_kmesh)
    multi[0].set_vars(tolvrs=1e-6, nband=4)

    # Dataset 2 (NSCF run)
    multi[1].set_vars(iscf=-2,
                      tolwfr=1e-10,
                      nband=10,
                      nbdbuf=2)

    # Dataset3: Calculation of the screening.
    multi[2].set_vars(
        optdriver=3,
        nband=8,
        ecutwfn=ecut,
        symchi=1,
        inclvkb=0,
        ecuteps=2.0,
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    kptgw = [
         -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
         -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
          5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
         -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
          5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
          0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
      ]

    multi[3].set_vars(
            optdriver=4,
            nband=10,
            ecutwfn=ecut,
            ecuteps=2.0,
            ecutsigx=2.0,
            symsigma=1,
            #gw_qprange=0,
    )

    bdgw = [4, 5]
    multi[3].set_kptgw(kptgw, bdgw)

    return multi.split_datasets()


def itest_g0w0_flow(fwp, tvars):
    """Test flow for G0W0 calculations."""
    scf, nscf, scr, sig = make_g0w0_inputs(ngkpt=[2, 2, 2], tvars=tvars)

    flow = flowtk.g0w0_flow(fwp.workdir, scf, nscf, scr, sig, manager=fwp.manager)
    # Will remove output files at run-time.
    flow.set_garbage_collector()
    flow.build_and_pickle_dump(abivalidate=True)

    for task in flow[0]:
        task.start_and_wait()

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    scf_task = flow[0][0]
    nscf_task = flow[0][1]
    scr_task = flow[0][2]
    sig_task = flow[0][3]

    # Test garbage_collector
    # The WFK|SCR file should have been removed because we call set_garbage_collector
    assert not scf_task.outdir.has_abiext("WFK")
    assert not nscf_task.outdir.has_abiext("WFK")
    assert not scr_task.outdir.has_abiext("SCR")
    assert not scr_task.outdir.has_abiext("SUS")

    # The sigma task should produce a SIGRES file.
    sigfile = sig_task.outdir.list_filepaths(wildcard="*SIGRES.nc")[0]
    assert sigfile
    with abilab.abiopen(sigfile) as sigres:
        sigres.to_string(verbose=2)
        assert sigres.nsppol == 1

    # Test SigmaTask inspect method
    #if has_matplotlib():
        #sig_task.inspect(show=False)

    # Test get_results for Sigma and Scr
    scr_task.get_results()
    sig_task.get_results()

    # Test SCR.nc file (this is optional)
    if scr_task.scr_path:
        with scr_task.open_scr() as scr:
            scr.to_string(verbose=2)
            assert len(scr.wpts) == 2
            assert scr.nwre == 1 and scr.nwim == 1
            for iq, qpoint in enumerate(scr.qpoints[:2]):
                #print(qpoint)
                qpt, iqcheck = scr.reader.find_qpoint_fileindex(qpoint)
                assert iqcheck == iq
                em1 = scr.get_em1(qpoint)
                #print(em1)

    # TODO Add more tests
    #assert flow.validate_json_schema()


def itest_g0w0qptdm_flow(fwp, tvars):
    """Integration test for G0W0WithQptdmFlow."""
    scf, nscf, scr, sig = make_g0w0_inputs(ngkpt=[2, 2, 2], tvars=tvars)

    flow = flowtk.G0W0WithQptdmFlow(fwp.workdir, scf, nscf, scr, sig, manager=fwp.manager)

    # Enable garbage collector at the flow level.
    # Note that here we have tp use this policy because tasks are created dynamically
    #flow.set_garbage_collector(policy="task")
    flow.set_garbage_collector(policy="flow")

    assert len(flow) == 3
    bands_work = flow[0]
    scr_work = flow[1]
    sigma_work = flow[2]

    assert scr_work.depends_on(bands_work.nscf_task)
    assert not scr_work.depends_on(bands_work.scf_task)

    for sigma_task in sigma_work:
        #print("sigma_task.deps", sigma_task.deps)
        assert sigma_task.depends_on(bands_work.nscf_task)
        assert not sigma_task.depends_on(bands_work.scf_task)
        assert sigma_task.depends_on(scr_work)

    flow.build_and_pickle_dump(abivalidate=True)
    flow.show_dependencies()
    # This call is needed to connect the node and enable
    # the callbacks, otherwise the scheduler enters a deadlock.
    flow.connect_signals()

    # Run the flow.
    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.show_status()
    assert all(work.finalized for work in flow)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    # Test set_garbage_collector
    # The WFK|SCR file should have been removed because we call set_garbage_collector
    #assert not scf_task.outdir.has_abiext("WFK")
    #assert not nscf_task.outdir.has_abiext("WFK")
    #assert not scr_task.outdir.has_abiext("SCR")
    #assert not scr_task.outdir.has_abiext("SUS")

    # The SCR file produced by scr_work should have been removed
    assert not scr_work.outdir.has_abiext("SCR")

    #assert flow.validate_json_schema()

    flow.finalize()


def itest_htc_g0w0(fwp, tvars):
    """Testing G0W0Work."""
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

    scf_kppa = 10
    nscf_nband = 10
    #nscf_ngkpt = [4,4,4]
    #nscf_shiftk = [0.0, 0.0, 0.0]
    ecut, ecuteps, ecutsigx = 4, 2, 3
    #scr_nband = 50
    #sigma_nband = 50

    extra_abivars = dict(
        ecut=ecut,
        istwfk="*1",
        paral_kgb=tvars.paral_kgb,
    )

    multi = abilab.g0w0_with_ppmodel_inputs(
        structure, pseudos,
        scf_kppa, nscf_nband, ecuteps, ecutsigx,
        ecut=ecut, pawecutdg=None,
        accuracy="normal", spin_mode="unpolarized", smearing=None,
        #ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
        #sigma_nband=None, gw_qprange=1):

    )
    multi.set_vars(paral_kgb=tvars.paral_kgb)

    scf_input, nscf_input, scr_input, sigma_input = multi.split_datasets()
    work = flowtk.G0W0Work(scf_input, nscf_input, scr_input, sigma_input)

    flow.register_work(work)
    flow.allocate()
    flow.connect_signals()

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions
    assert fwp.scheduler.nlaunch == 4

    # The sigma task should produce a SCR file.
    assert len(work[2].outdir.list_filepaths(wildcard="*SCR")) == 1

    flow.show_status()
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    assert all(work.finalized for work in flow)

    #assert flow.validate_json_schema()

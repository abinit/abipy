"""
Integration tests for flows/works/tasks that rely on external files e.g. DEN --> NscfTask.
"""
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk
from abipy.core.testing import AbipyTest


def make_scf_nscf_inputs(paral_kgb=1):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)

    # Global variables
    ecut = 6
    global_vars = dict(
        ecut=ecut,
        nband=8,
        nstep=15,
        paral_kgb=paral_kgb,
    )

    if multi.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(ngkpt=[8, 8, 8], shiftk=[0, 0, 0])
    multi[0].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
    multi[1].set_vars(tolwfr=1e-12)

    # Generate two input files for the GS and the NSCF run
    scf_input, nscf_input = multi.split_datasets()

    return scf_input, nscf_input


def itest_nscf_from_denfile(fwp, tvars):
    """Testing NscTask from pre-existent DEN file."""
    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs(paral_kgb=1)

    # Build the flow.
    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)

    # Create a Work, all tasks in work will start from the DEN file.
    # Note that the file must exist when the work is created
    # Use the standard approach based on tasks and works if
    # there's a node who needs a file produced in the future.
    # Need to copy DEN.nc to temp dir to avoid problem with multiple extensions.
    import shutil
    tmp_directory = AbipyTest.mkdtemp()
    den_filepath = os.path.join(tmp_directory, "si_DEN.nc")
    shutil.copyfile(abidata.ref_file("si_DEN.nc"), den_filepath)

    flow.register_nscf_task(nscf_input, deps={den_filepath: "DEN"}, append=True)

    flow.allocate(use_smartio=True)
    assert len(flow) == 1 and len(flow[0]) == 1

    task = flow[0][0]
    assert len(task.deps) == 1
    print(task.deps[0])
    filenode = task.deps[0].node
    assert filenode.is_file
    assert not filenode.is_task and not filenode.is_work and not filenode.is_flow
    assert not filenode.deps
    assert task.depends_on(filenode)
    assert not filenode.depends_on(task)
    assert filenode in task.get_parents()
    assert not filenode.get_parents()
    assert task in filenode.get_children()
    assert len(filenode.get_children()) == 1
    assert task.str_deps
    assert filenode.str_deps
    assert not task.get_children()
    #assert filenode.set_manager(fwp.manager)
    repr(filenode); str(filenode)
    assert filenode.filepath == den_filepath
    assert filenode.status == filenode.S_OK

    if AbipyTest.has_python_graphviz():
        assert callable(filenode.get_graphviz_dirtree().view)
        assert callable(task.get_graphviz_dirtree().view)
        assert callable(flow.get_graphviz().view)

    # Will remove output files (WFK)
    flow.set_garbage_collector()
    scheduler = flow.make_scheduler()
    assert scheduler.start() == 0
    assert not scheduler.exceptions
    assert scheduler.nlaunch == 1

    flow.check_status(show=True)
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()
    assert all(work.finalized for work in flow)

    # The WFK files should have been removed because we called set_garbage_collector
    # but the input DEN.nc should exist
    for task in flow[0]:
        assert not task.outdir.has_abiext("WFK")
    assert os.path.isfile(den_filepath)

    shutil.rmtree(tmp_directory)

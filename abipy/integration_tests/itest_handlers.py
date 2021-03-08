
import pytest
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def itest_tolsymerror_handler(fwp):
    """
    Test the handler of TolSymError. The test triggers:

        --- !TolSymError
        message: |
            Could not find the point group
        src_file: symptgroup.F90
        src_line: 236
        ...

    at the level of the symmetry finder and autoparal fails
    because it cannot find the parallel configurations.
    """
    pytest.xfail("tolsymerror_handler has been disabled because this problem has been fixed in v9.")
    structure = dict(
        acell=(1.0, 1.0, 1.0),
        xred=[
           1.0001907690, 1.0040151117, 0.0099335191,
           0.2501907744, 0.2540150788, 0.2599335332],
        rprim=[
          -6.2733366562, 0.0000000000, -3.6219126071,
          -6.2733366562, 0.0000000000,  3.6219126071,
          -4.1822244376, 5.9145585205,  0.0000000000],
        typat=(1, 1),
        ntypat=1,
        znucl=14,
        natom=2,
    )

    inp = abilab.AbinitInput(structure=structure, pseudos=abidata.pseudos("14si.pspnc"))

    inp.set_vars(
         ntime=5,
         tolrff=0.02,
         shiftk=[0, 0, 0],
         ngkpt=(4, 4, 4),
         chksymbreak=0,
         ecut=4,
         tolmxf=5e-05,
         nshiftk=1,
         #tolsym=1e-10,
    )

    flow = flowtk.Flow(workdir=fwp.workdir, manager=fwp.manager)
    flow.register_task(inp, task_class=flowtk.RelaxTask)

    flow.allocate()
    assert flow.make_scheduler().start() == 0

    flow.show_status()
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    task = flow[0][0]
    assert len(task.corrections) == 1
    assert task.corrections[0]["event"]["@class"] == "TolSymError"


def itest_dilatmxerror_handler(fwp):
    """Test the handler of DilatmxError. The test triggers:

        --- !DilatmxError
        message: |
            Dilatmx has been exceeded too many times (4)
            Restart your calculation from larger lattice vectors and/or a larger dilatmx
        src_file: mover.F90
        src_line: 840
        ...

    in variable cell structural optimizations.
    """
    #if fwp.on_travis:
    pytest.xfail("dilatmxerror_handler is not portable and it's been disabled!")

    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    structure.scale_lattice(structure.volume * 0.8)
    # Perturb the structure (random perturbation of 0.1 Angstrom)
    #structure.perturb(distance=0.1)

    inp = abilab.AbinitInput(structure=structure, pseudos=abidata.pseudos("14si.pspnc"))

    inp.set_vars(
        ecut=4,
        ngkpt=[4, 4, 4],
        shiftk=[0, 0, 0],
        nshiftk=1,
        chksymbreak=0,
        paral_kgb=1,
        optcell=1,
        ionmov=2,
        ecutsm=0.5,
        dilatmx=1.01,
        tolrff=0.02,
        tolmxf=5.0e-5,
        strfact=100,
        ntime=50,
        #ntime=5, To test the restart
        )

    # Create the flow
    flow = flowtk.Flow(fwp.workdir, manager=fwp.manager)
    flow.register_task(inp, task_class=flowtk.RelaxTask)

    flow.allocate()
    assert flow.make_scheduler().start() == 0

    flow.show_status()
    if not flow.all_ok:
        flow.debug()
        raise RuntimeError()

    task = flow[0][0]
    # Don't check the number of corrections as it's not portable.
    assert len(task.corrections)
    for i in range(task.num_corrections):
        assert task.corrections[i]["event"]["@class"] == "DilatmxError"
    #assert 0

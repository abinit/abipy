#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import abipy.abilab as abilab 
import abipy.flowtk as flowtk
import abipy.data as abidata


def gs_input(nsppol):
    # Fe normal bcc structure for test of a ferromagnetic calculation
    # The first dataset is without magnetization for comparison
    structure = abilab.Structure.from_abivars(dict(
        natom=1,
        ntypat=1,
        typat=1,
        znucl=26,
        acell=3*[5.42],
        rprim=[-0.5,  0.5,  0.5,
                0.5, -0.5,  0.5,
                0.5,  0.5, -0.5],
        xred=[0.0, 0.0, 0.0])
    )
    inp = abilab.AbinitInput(structure, pseudos=abidata.pseudos("26fe.pspnc"))

    inp.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0.5, 0.5, 0.5])

    # Optimization of the lattice parameters
    inp.set_vars(
        nsppol=nsppol,
        ecut=18,
        nband=8,
        occopt=3,
        tsmear=0.01,
        toldfe=1e-6, 
        nstep=50,
    )

    if nsppol == 2:
        inp.set_vars(spinat=[0.0, 0.0, 4.0])

    return inp


def afm_input():
    # Fe fcc structure with two atoms per unit cell for test of antiferromagnetic
    # This is the simplest fcc structure compatible with a X point spiral
    structure = abilab.Structure.from_abivars(dict(
        natom=2,
        ntypat=1,
        typat=[1, 1],
        znucl=26,
        acell=3*[6.60],
        rprim=[0.5, -0.5,  0.0,
               0.5,  0.5,  0.0,
               0.0,  0.0,  1.0],
        xred=[0.0, 0.0, 0.0,
              0.5, 0.0, 0.5],
        ))
    inp = abilab.AbinitInput(structure=structure, pseudos=abidata.pseudos("26fe.pspnc"))

    inp.set_kmesh(ngkpt=[6, 6, 4], shiftk=[0.5, 0.5, 0.5])

    # Antiferromagnetic order
    inp.set_vars(
        nsppol=1,
        nspden=2,
        spinat=[0.0, 0.0,  4.0,
                0.0, 0.0, -4.0],
        ecut=18, 
        nband=16,
        occopt=3,
        tsmear=0.01,
        tolwfr=1e-7, 
        nstep=70,
    )

    return inp


def gs_flow():
    inputs = [gs_input(nsppol) for nsppol in [1, 2]]
    flow = flowtk.Flow.from_inputs(workdir="flow_spin", inputs=inputs)
    flow.make_scheduler().start()

    with abilab.abirobot(flow, "GSR") as robot:
        data = robot.get_dataframe()
        print(data)
        robot.pairplot(x_vars="nsppol", y_vars=["energy", "a", "volume", "pressure"])

    #gstask_nospin, gstask_spin = flow[0][0], flow[0][1] 
    #data = abilab.PrettyTable(["property", "unpolarized", "polarized"])
    #with gstask_nospin.open_gsr() as gsr_nospin, gstask_spin.open_gsr() as gsr_spin:
    #    properties = ["energy", "pressure", "magnetization", "nelect_updown"]
    #    for p in properties:
    #        row = [p, getattr(gsr_nospin, p), getattr(gsr_spin, p)]
    #        data.add_row(row)

    #    plotter = abilab.ElectronDosPlotter()
    #    plotter.add_edos_from_file(gsr_spin.filepath, label="spin")
    #    plotter.add_edos_from_file(gsr_nospin.filepath, label="nospin")
    #    plotter.plot()

    #    #gsr_spin.plot_ebands()
    #    #gsr_nospin.plot_ebands_with_dos()
    #    #plotter = abilab.ElectronBandsPlotter()
    #    #plotter.add_ebands_from_file(gsr_spin.filepath, label="spin")
    #    #plotter.plot()
    #    #plotter = abilab.GSR_Plotter(gsr_nospin.filepath, gsr_spin.filepath)
    #    #plotter.add_ebands_from_file(gsr_nospin.filepath, label="nospin")
    #print(data)


def afm_flow():
    flow = flowtk.Flow.from_inputs(workdir="flow_afm", inputs=afm_input())

    flow.make_scheduler().start()

    with flow[0][0].open_gsr() as gsr:
        print("Energy: ", gsr.energy.to("Ha"))
        print("Magnetization: ",gsr.magnetization)


def tantalum_gsinput(nspinor=2):
    #  Single Ta atom in a big box (BCC), treated with spin-orbit coupling.
    structure = abilab.Structure.from_abivars(
        natom=1,
        ntypat=1,
        typat=[1],
        znucl=73,
        acell=3*[15.0],
        rprim=[ 0.5,  0.5, -0.5,
               -0.5,  0.5,  0.5,
                0.5, -0.5,  0.5],
        xred=[0.0, 0.0, 0.0]
    )
    inp = abilab.AbinitInput(structure=structure, pseudos=abidata.pseudos("73ta.hghsc"))

    inp.set_kmesh(ngkpt=[1, 1, 1], shiftk=[0.0, 0.0, 0.0])

    inp.set_vars(
        nspinor=nspinor,
        ecut=10, 
        ixc=2,
        istwfk=1,
        intxc=1,
        nband=26,
        occopt=7,
        tsmear=0.01,
        toldfe=1e-7, 
        nstep=70,
    )

    print(inp)
    return inp


def tantalum_flow():
    inputs = [tantalum_gsinput(nspinor) for nspinor in [1, 2]]
    flow = abilab.Flow.from_inputs(workdir="flow_tantalum", inputs=inputs)

    flow.make_scheduler().start()

    with abilab.GsrRobot.open(flow) as robot:
        data = robot.get_dataframe()
        print(data)
        robot.pairplot(x_vars="nspinor", y_vars=["energy", "magnetization", "pressure"])

    #for task in flow.iflat_tasks():
    #    with task.open_gsr() as gsr:
    #        print("Energy: ", gsr.energy.to("Ha"))
    #        print("Magnetization: ",gsr.magnetization)


if __name__ == "__main__":
    gs_flow()
    #afm_flow()
    #tantalum_flow()

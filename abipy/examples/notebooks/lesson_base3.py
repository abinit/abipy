#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab 
import abipy.data as abidata


def ngkpt_flow():
    ngkpt_list = [(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]

    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
    inp.set_structure_from_file(abidata.cif_file("si.cif"))

    # Global variables
    inp.set_variables(ecut=10, tolvrs=1e-9)

    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    flow = abilab.Flow.from_inputs("flow_base3_ngkpt", inputs=inp.split_datasets())
    flow.make_scheduler().start()

    with abilab.abirobot(flow, "GSR") as robot:
        #robot.ebands_plotter().plot()
        data = robot.get_dataframe()

    import matplotlib.pyplot as plt
    data.plot(x="nkpts", y="energy", title="Total energy vs nkpts", legend="Energy [eV]", style="b-o")
    plt.show()


def relax_flow():
    # Structural relaxation
    ngkpt_list = [(2, 2, 2), (4, 4, 4)]
    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
    inp.set_structure_from_file(abidata.cif_file("si.cif"))

    # Global variables
    inp.set_variables(
        ecut=10,
        tolvrs=1e-9,
        optcell=1,
        ionmov=3,
        ntime=10,
        dilatmx=1.05,
        ecutsm=0.5,
    )

    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    flow = abilab.Flow.from_inputs("flow_relax", inputs=inp.split_datasets(),
                                         task_class=abilab.RelaxTask)

    #flow.make_scheduler().start()
    flow.show_status()

    with abilab.GsrRobot.open(flow) as robot:
        data = robot.get_dataframe()

    import matplotlib.pyplot as plt
    import seaborn as sns
    #data.plot(x="nkpts", y="a", style="b-o")

    grid = sns.PairGrid(data, x_vars="nkpts", y_vars=["a", "volume"]) #, hue="tsmear")
    grid.map(plt.plot, marker="o")
    grid.add_legend()
    plt.show()


def bands_flow():
    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)
    inp.set_structure_from_file(abidata.cif_file("si.cif"))

    # Global variables
    inp.ecut = 10

    # Dataset 1
    inp[1].set_variables(tolvrs=1e-9)
    inp[1].set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])

    # Dataset 2
    inp[2].set_variables(tolwfr=1e-15)
    inp[2].set_kpath(ndivsm=5)

    scf_input, nscf_input = inp.split_datasets()

    flow = abilab.bandstructure_flow(workdir="flow_bands", scf_input=scf_input, nscf_input=nscf_input)
    flow.make_scheduler().start()
    
    nscf_task = flow[0][1]
    with nscf_task.open_gsr() as gsr:
        gsr.ebands.plot()


if __name__ == "__main__":
    ngkpt_flow()
    #relax_flow()
    #bands_flow()

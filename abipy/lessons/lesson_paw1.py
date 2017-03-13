#!/usr/bin/env python
from __future__ import division, print_function

import numpy as np
import abipy.abilab as abilab 
import abipy.flowapi as flowapi
import abipy.data as abidata


def gs_input(ecut, pawecutdg, acell_ang=3.567):
    # tpaw1_2.in
    # Input for PAW1 tutorial
    # Diamond at experimental volume
    structure = abilab.Structure.from_abivars(
        natom=2,
        ntypat=1,
        typat=2 * [1],
        znucl=6,
        acell=3*[abilab.Length(acell_ang, "ang").to("bohr")],
        rprim=[0.0,  0.5,  0.5,
               0.5,  0.0,  0.5,
               0.5,  0.5,  0.0],
        xred=[0.0, 0.0, 0.0,
              1/4, 1/4, 1/4],
    )
    inp = abilab.AbinitInput(structure=structure, pseudos=abidata.pseudos("6c.lda.atompaw"))

    # Optimization of the lattice parameters
    inp.set_vars(
        ecut=ecut, 
        pawecutdg=pawecutdg,
        ecutsm=0.5,
        nband=6,
        tolvrs=1e-10,
        nstep=20,
    )

    inp.set_autokmesh(nksmall=6) # ngkpt=[6, 6, 6], shiftk=[0.5, 0.5, 0.5])

    return inp


def flow_ecut_conv():
    inputs = [gs_input(ecut=ecut, pawecutdg=50) 
              for ecut in np.linspace(start=8, stop=24, num=9)]

    flow = flowapi.Flow.from_inputs("flow_ecut_conv", inputs)
    flow.make_scheduler().start()

    with abilab.abirobot(flow, "GSR") as robot:
        data = robot.get_dataframe()
        print(data)
        data.plot(x="ecut", y="energy", title="Energy vs ecut")


def flow_pawecutdg_conv_flow():
    inputs = [gs_input(ecut=12, pawecutdg=pawecutdg)
              for pawecutdg in np.linspace(start=12, stop=39, num=10)]

    flow = flowapi.Flow.from_inputs("flow_pawecutdg_conv", inputs)
    flow.make_scheduler().start()

    with abilab.abirobot(flow, "GSR") as robot:
        data = robot.get_dataframe()
        print(data)
        data.plot(x="pawecutdg", y="energy", title="Energy vs pawecutdg")


def flow_ecut_pawecutdg():
    import itertools
    ecut_list = np.linspace(start=8, stop=24, num=9)
    pawecutdg_list = [24, 30]
    inputs = [gs_input(ecut, pawecutdg) 
              for pawecutdg, ecut in itertools.product(pawecutdg_list, ecut_list)]

    flow = flowapi.Flow.from_inputs("flow_pawecutdg_ecut", inputs)
    flow.make_scheduler().start()

    with abilab.abirobot(flow, "GSR") as robot:
        data = robot.get_dataframe()
        print(data)
        robot.pairplot(x_vars="ecut", y_vars="energy", hue="pawecutdg")


def flow_eos():
    inputs = [gs_input(ecut=12, pawecutdg=24, acell_ang=acell_ang)
              for acell_ang in np.linspace(start=3.52, stop=3.55, num=7)]
    flow = flowapi.Flow.from_inputs("flow_eos", inputs)
    #flow.build()

    flow.make_scheduler().start()

    with abilab.abirobot(flow, "GSR") as robot:
        #fit = robot.eos_fit()
        fits, table = robot.eos_fit("all")
        print(table)

    print(fit)
    fit.plot()


if __name__ == "__main__":
    #flow_ecut_conv()
    #flow_pawecutdg_conv()
    #flow_eos()
    flow_ecut_pawecutdg()

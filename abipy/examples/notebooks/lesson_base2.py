#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab 
import abipy.data as abidata
import numpy as np

def h2_h_input(x=0.7, ecut=10, acell=(10, 10, 10)):
    """
    This file to optimize the H2 bond length, compute the associated total
    energy, then to compute the total energy of the isolated H atom.
    """
    inp = abilab.AbiInput(pseudos=abidata.pseudos("01h.pspgth"), ndtset=2)

    inp.set_variables(
        ecut=ecut, 
        nband=1,
        diemac=2.0,
        nstep=10,
        
    )

    inp.set_kmesh(ngkpt=(1,1,1), shiftk=(0,0,0))

    h2 = abilab.Structure.from_abivars(dict(
        natom=2,        
        ntypat=1,  
        typat=(1, 1),
        znucl=1,
        xcart=[-x, 0.0, 0.0,
               +x, 0.0, 0.0],
        acell=acell,
        rprim=[1, 0, 0, 0, 1, 0, 0, 0, 1],
    ))

    inp[1].set_structure(h2)
    inp[1].set_variables(
        ionmov=3,
        ntime=10,
        tolmxf=5e-4,
        toldff=5e-5,
    )

    h = abilab.Structure.from_abivars(dict(
        natom=1,        
        ntypat=1,  
        typat=1,
        znucl=1,
        xcart=[0.0, 0.0, 0.0],
        acell=acell,
        rprim=[1, 0, 0, 0, 1, 0, 0, 0, 1],
    ))

    inp[2].set_structure(h)
    inp[2].set_variables(
        nsppol=2,
        nband=(1, 1),
        occopt=2,
        occ=(1.0, 0.0),
        toldfe=1e-6,
        spinat=(0.0, 0.0, 1.0),
    )

    return inp.split_datasets()

def ecut_convergence_study():
    """
    H2 molecule in a big box
    Generate a flow to compute the total energy and forces as a function of the interatomic distance
    """
    inputs = []
    for ecut in range(10, 40, 5):
        inputs += h2_h_input(ecut=ecut)
    flow = abilab.AbinitFlow.from_inputs("flow_h2h_ecut", inputs)

    #flow.make_scheduler().start()

    import matplotlib.pyplot as plt
    import pandas as pd

    with abilab.GsrRobot.from_flow(flow) as robot:
        table = robot.get_dataframe()
        table = table[["formula", "ecut", "energy"]]
        print(table)

        grouped = table.groupby("ecut")
        rows = []
        for ecut, group in grouped:
            group = group.set_index("formula")
            atomization = 2 * group["energy"]["H1"] - group["energy"]["H2"]
            atomization = abilab.Energy(atomization, "eV").to("Ha")
            rows.append(dict(ecut=ecut, atomization=atomization))

        atable = pd.DataFrame(rows)
        print(atable)

        atable.plot("ecut", "atomization")
        plt.show()

def acell_convergence_study():
    """
    H2 molecule in a big box
    Generate a flow to compute the total energy and forces as a function of the interatomic distance
    """
    inputs = []
    for acell in range(8, 20, 2):
        inputs += h2_h_input(ecut=10, acell=3*[acell])
    flow = abilab.AbinitFlow.from_inputs("flow_h2h_acell", inputs)

    #flow.make_scheduler().start()

    def hh_dist(gsr):
        """This function receives a GSR file and computes the H-H distance"""
        if len(gsr.structure) == 1:
            l = None
        else:
            cart_coords = gsr.structure.cart_coords
            l = np.sqrt(np.linalg.norm(cart_coords[1] - cart_coords[0]))
        return "hh_dist", l

    import matplotlib.pyplot as plt
    import pandas as pd

    with abilab.GsrRobot.from_flow(flow) as robot:
        table = robot.get_dataframe(funcs=hh_dist)
        table = table[["formula", "a", "energy", "hh_dist"]]
        print(table)

        grouped = table.groupby("a")
        rows = []
        for a, group in grouped:
            group = group.set_index("formula")
            atomization = 2 * group["energy"]["H1"] - group["energy"]["H2"]
            atomization = abilab.Energy(atomization, "eV").to("Ha")
            # FIXME Values for hh_dist are wrong! why?
            rows.append(dict(a=a, atomization=atomization, hh_dist=group["hh_dist"]["H2"]))

        atable = pd.DataFrame(rows)
        print(atable)

        atable.plot("a", "atomization")
        plt.show()


if __name__ == "__main__":
    #ecut_convergence_study()
    acell_convergence_study()

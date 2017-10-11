#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import abipy.abilab as abilab 
import abipy.flowtk as flowtk
import abipy.data as abidata
import numpy as np


def h2_h_input(x=0.7, ecut=10, acell=(10, 10, 10)):
    """
    This file to optimize the H2 bond length, compute the associated total
    energy, then to compute the total energy of the isolated H atom.
    """
    h2 = abilab.Structure.from_abivars(
        natom=2,        
        ntypat=1,  
        typat=(1, 1),
        znucl=1,
        xcart=[-x, 0.0, 0.0,
               +x, 0.0, 0.0],
        acell=acell,
        rprim=[1, 0, 0, 0, 1, 0, 0, 0, 1],
    )

    h = abilab.Structure.from_abivars(
        natom=1,        
        ntypat=1,  
        typat=1,
        znucl=1,
        xcart=[0.0, 0.0, 0.0],
        acell=acell,
        rprim=[1, 0, 0, 0, 1, 0, 0, 0, 1],
    )

    global_vars = dict(
        ecut=ecut, 
        nband=1,
        diemac=2.0,
        nstep=10,
    )

    h2_inp = abilab.AbinitInput(structure=h2, pseudos=abidata.pseudos("01h.pspgth"))

    h2_inp.set_vars(global_vars)
    h2_inp.set_kmesh(ngkpt=(1,1,1), shiftk=(0,0,0))
    h2_inp.set_vars(
        ionmov=3,
        ntime=10,
        tolmxf=5e-4,
        toldff=5e-5,
    )

    h_inp = abilab.AbinitInput(structure=h, pseudos=abidata.pseudos("01h.pspgth"))
    h_inp.set_vars(global_vars)

    h_inp.set_vars(
        nsppol=2,
        nband=(1, 1),
        occopt=2,
        occ=(1.0, 0.0),
        toldfe=1e-6,
        spinat=(0.0, 0.0, 1.0),
    )

    return h2_inp, h_inp


def ecut_convergence_study(ecuts=range(10, 40, 5)):
    """
    H2 molecule in a big box
    Generate a flow to compute the total energy and forces as a function of the interatomic distance
    """
    inputs = []
    for ecut in ecuts:
        inputs += h2_h_input(ecut=ecut)

    flow = flowtk.Flow.from_inputs("flow_h2h_ecut", inputs)
    flow.make_scheduler().start()

    import matplotlib.pyplot as plt
    import pandas as pd

    with abilab.abirobot(flow, "GSR") as robot:
        frame = robot.get_dataframe()
        frame = frame[["formula", "ecut", "energy"]]
        print(frame)

        grouped = frame.groupby("ecut")
        rows = []
        for ecut, group in grouped:
            group = group.set_index("formula")
            atomization = 2 * group["energy"]["H1"] - group["energy"]["H2"]
            atomization = abilab.Energy(atomization, "eV").to("Ha")
            rows.append(dict(ecut=ecut, atomization=atomization))

        atomization_frame = pd.DataFrame(rows)
        print(atomization_frame)

        atomization_frame.plot("ecut", "atomization")
        plt.show()


def acell_convergence_study(acell_list=range(8, 20, 2), ecut=10):
    """
    H2 molecule in a big box
    Generate a flow to compute the total energy and the forces as function of the interatomic distance
    """
    inputs = []
    for acell in acell_list:
        inputs += h2_h_input(ecut=ecut, acell=3*[acell])
    flow = flowtk.Flow.from_inputs("flow_h2h_acell", inputs)
    flow.make_scheduler().start()

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

    with abilab.abirobot(flow, "GSR") as robot:
        frame = robot.get_dataframe(funcs=hh_dist)
        frame = frame[["formula", "a", "energy", "hh_dist"]]
        print(frame)

        grouped = frame.groupby("a")
        rows = []
        for a, group in grouped:
            group = group.set_index("formula")
            atomization = 2 * group["energy"]["H1"] - group["energy"]["H2"]
            atomization = abilab.Energy(atomization, "eV").to("Ha")
            # FIXME Values for hh_dist are wrong! why?
            rows.append(dict(a=a, atomization=atomization, hh_dist=group["hh_dist"]["H2"]))

        atomization_frame = pd.DataFrame(rows)
        print(atomization_frame)

        atomization_frame.plot("a", "atomization")
        plt.show()


if __name__ == "__main__":
    #ecut_convergence_study()
    acell_convergence_study()

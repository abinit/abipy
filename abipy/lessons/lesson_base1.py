#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab 
import abipy.data as abidata
import numpy as np


def gs_input(x=0.7, acell=(10, 10, 10)):
    """H2 molecule in a big box"""
    structure = abilab.Structure.from_abivars(
        ntypat=1,  
        znucl=1,
        natom=2,        
        typat=(1, 1),
        xcart=[-x, 0.0, 0.0,
               +x, 0.0, 0.0],
        acell=acell,
        rprim=[1, 0, 0, 0, 1, 0, 0, 0, 1]
    )

    inp = abilab.AbinitInput(structure=structure, pseudos=abidata.pseudos("01h.pspgth"))

    inp.set_vars(
        ecut=10, 
        nband=1,
        diemac=2.0,
        nstep=10,
        toldfe=1e-6,
    )

    inp.set_kmesh(ngkpt=(1,1,1), shiftk=(0,0,0))
    return inp


def analyze_flow(flow):
    def hh_dist(gsr):
        """This function receives a GSR file and computes the H-H distance"""
        cart_coords = gsr.structure.cart_coords
        l = np.sqrt(np.dot(cart_coords[1] - cart_coords[0]))
        return "hh_dist", l

    with abilab.abirobot(flow, "GSR") as robot:
        table = robot.get_dataframe(funcs=hh_dist)
        print(table)
        #robot.ebands_plotter().plot()

    return table


def build_flow():
    """
    H2 molecule in a big box:
    Generate a flow to compute the total energy and forces as a function of the interatomic distance
    """
    inputs = [gs_input(x) for x in np.linspace(0.5, 1.025, 21)]
    return abilab.Flow.from_inputs("flow_h", inputs)

    #table = abilab.PrettyTable(["length [Ang]", "energy [eV]"])
    #for task in flow.iflat_tasks():
    #    with task.open_gsr() as gsr:
    #        cart_coords = gsr.structure.cart_coords
    #        l = np.sqrt(np.linalg.norm(cart_coords[1] - cart_coords[0]))
    #        table.add_row([l, float(gsr.energy)])
    #print(table)
    #table.plot(title="Etotal vs interatomic distance")
    # Quadratic fit
    #fit = table.quadfit()
    #print(fit)
    #fit.plot()

if __name__ == "__main__":
    flow = build_flow()
    flow.make_scheduler().start()
    analyze_flow(flow)

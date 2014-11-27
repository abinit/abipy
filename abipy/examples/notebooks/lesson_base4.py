#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab 
import abipy.data as abidata


def relax_input(tsmear, nksmall):
    # Crystalline aluminum : optimization of the lattice parameter
    # at fixed number of k points and broadening.
    inp = abilab.AbiInput(pseudos=abidata.pseudos("13al.981214.fhi"))
    structure = abidata.ucells.structure_from_ucell("Al")
    inp.set_structure(structure)

    # nshiftk and shift are automatically selected from the lattice.
    #Definition of the k-point grid
    #ngkpt 2 2 2       # This is a 2x2x2 FCC grid, based on the primitive vectors
    #nshiftk 4         # of the reciprocal space. For a FCC real space lattice,
    #                  # like the present one, it actually corresponds to the
    #                  # so-called 4x4x4 Monkhorst-Pack grid, if the following shifts
    #                  # are used :
    #shiftk 0.5 0.5 0.5
    #       0.5 0.0 0.0
    #       0.0 0.5 0.0
    #       0.0 0.0 0.5
    inp.set_autokmesh(nksmall=nksmall)

    # Optimization of the lattice parameters
    inp.set_variables(
        ecut=6, 
        occopt=4,
        tsmear=tsmear,
        toldfe=1e-6, 
        nstep=10,
        optcell=1,
        ionmov=3,
        ntime=10,
        dilatmx=1.05,
        ecutsm=0.5,
        ixc=1,
    )

    #print(inp)
    return inp


def relax_flow():
    inp = relax_input(tsmear=0.05, nksmall=2)
    flow = abilab.AbinitFlow.from_inputs("flow_al_relax", inp)

    flow.make_scheduler().start()

    #table = abilab.PrettyTable(["nkpts", "etotal"])
    gs_task = flow[0][0]
    with gs_task.open_gsr() as gsr:
        print("input structure:\n", structure)
        print("relaxed structure:\n", gsr.structure)
        # TODO
        #print(gsr.energy_components)
        #return gsr


def convergence():
    # Cartesian product of input iterables. Equivalent to nested for-loops.
    from itertools import product
    tsmear_list = [0.01, 0.02, 0.03, 0.04]
    nksmall_list = [2, 4, 6]

    inputs = [relax_input(tsmear, nksmall) for tsmear, nksmall in product(tsmear_list, nksmall_list)]

    flow = abilab.AbinitFlow.from_inputs(workdir="flow_al_conv_relax", inputs=inputs)

    #flow.make_scheduler().start()

    with abilab.GsrRobot.open(flow) as robot:
        data = robot.get_dataframe()

    import matplotlib.pyplot as plt
    import seaborn as sns

    print(data)
    grid = sns.PairGrid(data, x_vars="nkpts", y_vars=["energy", "a", "volume"], hue="tsmear")
    grid.map(plt.plot, marker="o")
    grid.add_legend()

    plt.show()
    #grid = sns.FacetGrid(data, col="tsmear") 
    #grid.map(sns.pointplot, "nkpts", "a") 
    #sns.pairplot(data, x_vars="nkpts", y_vars=["energy", "a", "volume"], hue="tsmear")
    #grid.map(plt.scatter, s=50)


if __name__ == "__main__":
    #relax_flow()
    convergence()

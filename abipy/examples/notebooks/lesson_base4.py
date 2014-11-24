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

    print(inp)
    return inp

def relax_flow():
    inp = relax_input(tsmear=0.05, nksmall=2)
    flow = abilab.AbinitFlow.from_inputs("flow_al_relax", inp)
    flow.build()

    flow.make_scheduler().start()
    flow.show_status()

    #table = abilab.PrettyTable(["nkibz", "etotal"])
    gs_task = flow[0][0]
    with gs_task.open_gsr() as gsr:
        print("input structure:\n", structure)
        print("relaxed structure:\n", gsr.structure)
        # TODO
        #print(gsr.energy_components)
        #return gsr

def convergence():
    tsmear_list = [0.01, 0.02, 0.03, 0.04]
    nksmall_list = [2, 4, 6]

    work = abilab.Workflow()

    # Cartesian product of input iterables. Equivalent to nested for-loops.
    from itertools import product
    for tsmear, nksmall in product(tsmear_list, nksmall_list):
        inp = relax_input(tsmear, nksmall)
        task = work.register_relax_task(inp)
        task.set_user_info(tsmear=tsmear, nksmall=nksmall)

    flow = abilab.AbinitFlow(workdir="flow_al_conv_relax",)
    flow.register_work(work)
    flow.allocate()
    flow.build()

    #flow.make_scheduler().start()
    flow.show_status()

    #table = abilab.PrettyTable(["nkibz", "etotal"])
    #gs_task = flow[0][0]

    rows, row_names = [], []
    for task in flow.iflat_tasks():
        with task.open_gsr() as gsr:
            info = task.user_info
            structure = gsr.structure
            abc, angles = structure.lattice.abc, structure.lattice.angles
            nkibz = len(gsr.kpoints)
            row_names.append(task.pos_str)
            rows.append(dict(
                tsmear=info.tsmear, nksmall=info.nksmall, nkibz=nkibz, 
                nsppol=gsr.nsppol, nspinor=gsr.nspinor, nspden=gsr.nspden,
                energy=gsr.energy, magnetization=gsr.magnetization,
                a=abc[0], b=abc[1], c=abc[2], volume=structure.volume,
                angle0=angles[0], angle1=angles[1], angle2=angles[2],
                ))


    import pandas as pd
    data = pd.DataFrame(rows, index=row_names, columns=rows[0].keys())
    print(data)
    #data.plot(x="nkibz", y="energy")

    import matplotlib.pyplot as plt
    import seaborn as sns
    grid = sns.FacetGrid(data, col="tsmear", size=7) #, col_wrap=4, size=7)
    grid.map(sns.pointplot, "nkibz", "a", color=".3", ci=None)

    sns.pairplot(data, x_vars="nkibz", y_vars=["energy", "a", "volume"], hue="tsmear")
    plt.show()


from collections import Iterable, OrderedDict

class GsrRobot(Iterable):
    # TODO: Write mixin HasGsrFiles
    def __init__(self, *args):
        self._gsr_files = OrderedDict()
        for label, ncfile in args:
            self.add_file(label, ncfile)

    def add_file(self, label, ncfile)
        self._gsr_files[label] = ncfile

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement."""
        self.close()

    def close(self):
        """It automatically closes all the files that have been opened by self"""
        for i, gsr in enumerate(self):
            if self._doclose[i]: 
                gsr.close()

    def get_dataframe(self, **kwargs)
        rows, row_names = [], []
        for gsr in self:
            structure = gsr.structure
            abc, angles = structure.lattice.abc, structure.lattice.angles
            nkibz = len(gsr.kpoints)
            row_names.append(task.pos_str)
            rows.append(dict(
                #tsmear=info.tsmear, nksmall=info.nksmall, nkibz=nkibz, 
                nsppol=gsr.nsppol, nspinor=gsr.nspinor, nspden=gsr.nspden,
                energy=gsr.energy, magnetization=gsr.magnetization,
                a=abc[0], b=abc[1], c=abc[2], volume=structure.volume,
                angle0=angles[0], angle1=angles[1], angle2=angles[2],
                ))

        import pandas as pd
        return pd.DataFrame(rows, index=row_names, columns=rows[0].keys())


if __name__ == "__main__":
    #relax_flow()
    convergence()

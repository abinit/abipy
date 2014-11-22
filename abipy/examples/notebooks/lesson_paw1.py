#!/usr/bin/env python
from __future__ import division, print_function

import numpy as np
import abipy.abilab as abilab 
import abipy.data as abidata

def gs_input(ecut, pawecutdg, acell_ang=3.567):
    # tpaw1_2.in
    # Input for PAW1 tutorial
    # Diamond at experimental volume

    #inp = abilab.AbiInput(pseudos=abidata.pseudos("6c.pspnc"))
    inp = abilab.AbiInput(pseudos=abidata.pseudos("6c.lda.atompaw"))
    structure = abilab.Structure.from_abivars(dict(
        natom=2,
        ntypat=1,
        typat=2 * [1],
        znucl=6,
        acell=3*[abilab.Length(acell_ang, "ang").to("bohr")],
        rprim=[0.0,  0.5,  0.5,
               0.5,  0.0,  0.5,
               0.5,  0.5,  0.0],
        xred=[0.0, 0.0, 0.0,
              1/4, 1/4, 1/4])
    )
    inp.set_structure(structure)

    inp.set_autokmesh(nksmall=6) #ngkpt=[6, 6, 6], shiftk=[0.5, 0.5, 0.5])

    # Optimization of the lattice parameters
    inp.set_variables(
        ecut=ecut, 
        pawecutdg=pawecutdg,
        ecutsm=0.5,
        nband=6,
        tolvrs=1e-10, 
        nstep=20,
    )

    print(inp)
    return inp


def ecutconv_flow():
    inputs = [gs_input(ecut=ecut, pawecutdg=50) 
              for ecut in np.linspace(start=8, stop=24, num=9)]
    flow = abilab.AbinitFlow.from_inputs("flow_ecutconv", inputs)
    flow.build()
    flow.make_scheduler().start()

    table = abilab.PrettyTable(["ecut", "energy"])
    for task in flow.iflat_tasks():
        with task.open_gsr() as gsr:
            #table.add_row([gsr.ecut, gsr.energy])
            table.add_row([None, gsr.energy.to("Ha")])

    energies = table.get_column("energy")
    table.add_column("Delta_energy", [e - energies[-1] for e in energies])

    print(table)

def pawecutdgconv_flow():

    inputs = [gs_input(ecut=12, pawecutdg=pawecutdg)
              for pawecutdg in np.linspace(start=12, stop=39, num=10)]
    flow = abilab.AbinitFlow.from_inputs("flow_ecutconv", inputs)
    flow.build()

    flow.make_scheduler().start()

    table = abilab.PrettyTable(["pawecutdg", "energy"])
    for task in flow.iflat_tasks():
        with task.open_gsr() as gsr:
            #table.add_row([gsr.ecut, gsr.energy])
            table.add_row([None, gsr.energy.to("Ha")])

    energies = table.get_column("energy")
    table.add_column("Delta_energy", [e - energies[-1] for e in energies])

    print(table)


def eos_flow():
    inputs = [gs_input(ecut=12, pawecutdg=24, acell_ang=acell_ang)
              for acell_ang in np.linspace(start=3.52, stop=3.55, num=7)]
    flow = abilab.AbinitFlow.from_inputs("flow_eos", inputs)
    flow.build()

    flow.make_scheduler().start()

    energies, volumes = [], []
    for task in flow.iflat_tasks():
        with task.open_gsr() as gsr:
            energies.append(gsr.energy)
            volumes.append(gsr.structure.volume)

    eos = abilab.EOS(eos_name='birch_murnaghan')
    fit = eos.fit(volumes, energies)
    print(fit)
    fit.plot()

if __name__ == "__main__":
    #import pstats, cProfile
    #cProfile.runctx("ecutconv_flow()", globals(), locals(), "Profile.prof")
    #s = pstats.Stats("Profile.prof")
    #s.strip_dirs().sort_stats("time").print_stats()
    #ecutconv_flow()
    #pawecutdgconv_flow()
    eos_flow()

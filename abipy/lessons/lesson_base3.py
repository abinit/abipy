#!/usr/bin/env python
"""
Template for lessons: paragraph with a brief description of the lesson separated by a blank line from the text below.

Additional info go here. 
Users will import the module to have access to its public API
Example::

from abipy.tutorias import lesson_base1 as lesson

# To get help:
lesson.help()

# To build the flow:
flow = lesson.make_ngkpt_flow()

# To print the input files 
flow.show_inputs()

# start the flow with
flow.make_scheduler().start()

# Wait for completion and analyze the results.
flow.analyze()
"""
from __future__ import division, print_function

import sys
import abipy.abilab as abilab 
import abipy.data as abidata


def help(stream=sys.stdout):
    stream.write(__doc__)


class NgkptFlow(abilab.Flow):
    def analyze(self):
        with abilab.abirobot(self, "GSR") as robot:
            data = robot.get_dataframe()
            #robot.ebands_plotter().plot()

        import matplotlib.pyplot as plt
        data.plot(x="nkpts", y="energy", title="Total energy vs nkpts", legend="Energy [eV]", style="b-o")
        plt.show()


def make_ngkpt_flow():
    ngkpt_list = [(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]

    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
    inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    inp.set_variables(ecut=10, tolvrs=1e-9)

    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    return NgkptFlow.from_inputs(workdir="flow_base3_ngkpt", inputs=inp.split_datasets())


class RelaxFlow(abilab.Flow):
    def analyze(self):
        with abilab.GsrRobot.open(self) as robot:
            data = robot.get_dataframe()
            robot.pairplot(x_vars="nkpts", y_vars=["a", "volume"]) #, hue="tsmear")

            #grid = sns.PairGrid(data, x_vars="nkpts", y_vars=["a", "volume"]) #, hue="tsmear")
            #grid.map(plt.plot, marker="o")
            #grid.add_legend()
            #plt.show()


def make_relax_flow():
    # Structural relaxation for different k-point samplings.
    ngkpt_list = [(2, 2, 2), (4, 4, 4)]
    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
    inp.set_structure(abidata.cif_file("si.cif"))

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

    return RelaxFlow.from_inputs("flow_base3_relax", inputs=inp.split_datasets(), task_class=abilab.RelaxTask)

class EbandsFlow(abilab.Flow):
    def analyze(self):
        nscf_task = self[0][1]
        with nscf_task.open_gsr() as gsr:
            return gsr.ebands.plot()


def make_ebands_flow():
    """Band structure calculation."""
    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)
    inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    inp.ecut = 10

    # Dataset 1
    inp[1].set_variables(tolvrs=1e-9)
    inp[1].set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])

    # Dataset 2
    inp[2].set_variables(tolwfr=1e-15)
    inp[2].set_kpath(ndivsm=5)

    scf_input, nscf_input = inp.split_datasets()

    return abilab.bandstructure_flow(workdir="flow_base3_ebands", scf_input=scf_input, nscf_input=nscf_input, flow_class=EbandsFlow)


if __name__ == "__main__":
    flow = make_ngkpt_flow()
    #flow = make_relax_flow()
    #flow = make_ebands_flow()
    #flow.show_inputs()
    #flow.make_scheduler().start()
    #flow.analyze()

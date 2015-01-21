#!/usr/bin/env python
"""
Template for lessons: paragraph with a brief description of the lesson separated by a blank line from the text below.

Detailed description of the tutorial (goals, step-by-step description of the operations to be performed)
Users will import this module to access the public API. Each module should provide a help method
that prints this doc string. One of more factory functions (make_flow, make_convergence_flow ...)
that build and return a subclass of abilab.Flow. The Flow subclass provides a `analyze` method 
that performs the post-processing of the results and produces the final results (matplotlib plots, pandas dataframes, ...)
Users should be able to run the tutorial either via this script, or interactively inside ipython or ipython notebooks
The working directory of the flow should be named: flow_[name_of_the_lesson][_extra_info] so that each 
lesson will be done in different directories.

Example::

\033[91m Title\033[0m

\033[94m Background\033[0m

\033[94m The related abinit variables\033[0m

\033[1m ngkpt \033[0m

\033[94m The abipy flows in this lesson\033[0m

\033[94m The cource of this lesson\033[0m

Start this lessen by importing it in a new namespace:

\033[92m In []:\033[0m from abipy.tutorias import lesson_base1 as lesson

As always you can reread this lessons text using the command:

\033[92m In []:\033[0m lesson.help()

To build the flow:

\033[92m In []:\033[0m flow = lesson.make_flow()

To print the input files

\033[92m In []:\033[0m flow.show_inputs()

Start the flow with the scheduler and wait for completion.

\033[92m In []:\033[0m flow.make_scheduler().start()

To analyze the results.

\033[92m In []:\033[0m flow.analyze()
"""
from __future__ import division, print_function

import sys
import abipy.abilab as abilab 
import abipy.data as abidata
from abipy.lessons.lesson_helper_functions import help, abinit_help, get_local_copy


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
    flow.make_scheduler().start()
    flow.analyze()

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

from abipy.tutorias import lesson_base1 as lesson

# To get help:
lesson.help()

# To build the flow:
flow = lesson.make_ngkpt_flow()

# To print the input files 
flow.show_inputs()

# Start the flow with the scheduler and wait for completion.
flow.make_scheduler().start()

# To analyze the results.
flow.analyze()
"""
from __future__ import division, print_function

import sys
import abipy.abilab as abilab 
import abipy.flowtk as flowtk
import abipy.data as abidata


def help(stream=sys.stdout):
    stream.write(__doc__)


class NgkptFlow(flowtk.Flow):
    def analyze(self):
        with abilab.abirobot(self, "GSR") as robot:
            data = robot.get_dataframe()
            #robot.ebands_plotter().plot()

        import matplotlib.pyplot as plt
        data.plot(x="nkpts", y="energy", title="Total energy vs nkpts", legend="Energy [eV]", style="b-o")
        plt.show()


def make_ngkpt_flow():
    ngkpt_list = [(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]

    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
    # Global variables
    multi.set_vars(ecut=10, tolvrs=1e-9)

    for i, ngkpt in enumerate(ngkpt_list):
        multi[i].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    return NgkptFlow.from_inputs(workdir="flow_base3_ngkpt", inputs=multi.split_datasets())


class RelaxFlow(flowtk.Flow):
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
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))

    # Global variables
    multi.set_vars(
        ecut=10,
        tolvrs=1e-9,
        optcell=1,
        ionmov=3,
        ntime=10,
        dilatmx=1.05,
        ecutsm=0.5,
    )

    for i, ngkpt in enumerate(ngkpt_list):
        multi[i].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    return RelaxFlow.from_inputs("flow_base3_relax", inputs=multi.split_datasets(), task_class=flowtk.RelaxTask)


def make_ebands_flow():
    """Band structure calculation."""
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)
    # Global variables
    multi.set_vars(ecut=10)

    # Dataset 1
    multi[0].set_vars(tolvrs=1e-9)
    multi[0].set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])

    # Dataset 2
    multi[1].set_vars(tolwfr=1e-15)
    multi[1].set_kpath(ndivsm=5)

    scf_input, nscf_input = multi.split_datasets()

    return flowtk.bandstructure_flow(workdir="flow_base3_ebands", scf_input=scf_input, nscf_input=nscf_input)


if __name__ == "__main__":
    flow = make_ngkpt_flow()
    #flow = make_relax_flow()
    #flow = make_ebands_flow()

    flow.show_inputs()
    flow.make_scheduler().start()
    flow.analyze()

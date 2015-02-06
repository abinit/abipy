#!/usr/bin/env python
"""
\033[91m The calculation of the density of states and the bandstructure.\033[0m

\033[94m Background\033[0m

\033[94m The related abinit variables\033[0m

\033[1m ... \033[0m
\033[1m ... \033[0m
\033[1m ... \033[0m
\033[1m ... \033[0m

More info on the inputvariables and their use can be obtianed using the following function:

\033[92m In []:\033[0m lesson.abinit_help(inputvariable)

This will print the official abinit description of this inputvariable.

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
import os
import shutil
import abipy.abilab as abilab
import abipy.data as abidata

from abipy.lessons.lesson_helper_functions import abinit_help


def help(stream=sys.stdout):
    """
    Display the tutorial text.
    """
    stream.write(__doc__)


def get_local_copy():
    """
    Copy this script to the current working dir to explore and edit
    """
    dst = os.path.basename(__file__[:-1])
    if os.path.exists(dst):
        raise RuntimeError("file %s already exists. Remove it before calling get_local_copy" % dst)
    shutil.copyfile(__file__[:-1], dst)


class DosFlow(abilab.Flow):
    def analyze(self):
        #with abilab.abirobot(self, "GSR") as robot:
        #    data = robot.get_dataframe()
        #    #robot.ebands_plotter().plot()

        #import matplotlib.pyplot as plt
        #data.plot(x="energy", y="dos", title="Density Of States", legend="Energy [eV]", style="b-o")
        ##todo correct plot matteo could you please?
        #plt.show()

        plotter = abilab.ElectronDosPlotter()
        for task in self[0]:
            with task.open_gsr() as gsr:
                edos = gsr.ebands.get_edos(method="gaussian", step=0.01, width=0.1)
                ngkpt = task.strategy.abinit_input[0]["ngkpt"]
                plotter.add_edos("ngkpt %s" % str(ngkpt), edos)
        
        return plotter.plot()


def make_kptdos_flow():
    # Here we compute DOSes with SCF runs, a more standard and efficient approach would be 
    # to a SCF run and then multiple NSCF runs to get the DOS.
    ngkpt_list = [(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]

    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
    inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    inp.set_vars(ecut=10, tolvrs=1e-9)

    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    return DosFlow.from_inputs(workdir="flow_kptdos", inputs=inp.split_datasets())


class EbandsFlow(abilab.Flow):
    def analyze(self):
        nscf_task = self[0][1]
        with nscf_task.open_gsr() as gsr:
            return gsr.ebands.plot()

    def plot_edoses(self, method="gaussian", step=0.01, width=0.1, **kwargs):
        plotter = abilab.ElectronDosPlotter()
        for task in self.dos_tasks:
            with task.open_gsr() as gsr:
                edos = gsr.ebands.get_edos(method=method, step=step, width=width)
                ngkpt = task.get_inpvar("ngkpt")
                plotter.add_edos("ngkpt %s" % str(ngkpt), edos)

        return plotter.plot(**kwargs)

    def plot_ebands_and_dos(self, dos_idx=0, **kwargs):
        # plot dos
        with self.nscf_task.open_gsr() as gsr: 
            gs_ebands = gsr.ebands

        with self.dos_tasks[dos_idx].open_gsr() as gsr: 
            dos_ebands = gsr.ebands

        edos = dos_ebands.get_edos(method="gaussian", step=0.01, width=0.1)
        return gs_ebands.plot_with_edos(edos, **kwargs)


def make_electronic_structure_flow(ngkpts_for_dos = [(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]):
    """Band structure calculation."""
    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=2 + len(ngkpts_for_dos))
    inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    inp.ecut = 10

    # Dataset 1
    inp[1].set_vars(tolvrs=1e-9)
    inp[1].set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])

    # Dataset 2
    inp[2].set_vars(tolwfr=1e-15)
    inp[2].set_kpath(ndivsm=5)

    # Dataset 3
    for i, ngkpt in enumerate(ngkpts_for_dos):
        inp[3+i].set_vars(tolwfr=1e-15)
        inp[3+i].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    inputs = inp.split_datasets()
    scf_input, nscf_input, dos_input = inputs[0], inputs[1], inputs[2:]

    return abilab.bandstructure_flow(workdir="flow_base3_ebands", scf_input=scf_input, nscf_input=nscf_input,
                                     dos_inputs=dos_input, flow_class=EbandsFlow)


if __name__ == "__main__":
    #flow = make_electronic_structure_flow()
    #flow.build_and_pickle_dump()
    #flow = make_kptdos_flow()
    #flow.make_scheduler().start()
    #flow = abilab.Flow.pickle_load("flow_kptdos")
    flow = abilab.Flow.pickle_load("flow_base3_ebands")
    #flow.analyze()
    #flow.plot_ebands_and_dos()
    flow.plot_edoses()

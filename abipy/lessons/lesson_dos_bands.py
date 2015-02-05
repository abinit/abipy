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
        with abilab.abirobot(self, "GSR") as robot:
            data = robot.get_dataframe()
            #robot.ebands_plotter().plot()

        import matplotlib.pyplot as plt
        data.plot(x="energy", y="dos", title="Density Of States", legend="Energy [eV]", style="b-o")
        #todo correct plot matteo could you please?
        plt.show()


def make_kptdos_flow():
    ngkpt_list = [(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]

    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
    inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    inp.set_variables(ecut=10, tolvrs=1e-9)

    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    return DosFlow.from_inputs(workdir="flow_kptdos", inputs=inp.split_datasets())


class EbandsFlow(abilab.Flow):
    def analyze(self):
        nscf_task = self[0][1]
        with nscf_task.open_gsr() as gsr:
            return gsr.ebands.plot()

        # plot dos


def make_electronic_structure_flow():
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

    # Dataset 3
    inp[2].set_variables(tolwfr=1e-15)
    inp[1].set_kmesh(ngkpt=[6,6,6], shiftk=[0,0,0])

    scf_input, nscf_input, dos_input = inp.split_datasets()

    return abilab.bandstructure_flow(workdir="flow_base3_ebands", scf_input=scf_input, nscf_input=nscf_input,
                                     dos_inputs=dos_input, flow_class=EbandsFlow, manager=None)


if __name__ == "__main__":
    flow = make_electronic_structure_flow()
    flow.make_scheduler().start()
    flow.analyze()

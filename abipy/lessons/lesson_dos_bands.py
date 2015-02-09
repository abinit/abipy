#!/usr/bin/env python
"""
The calculation of the density of states and the bandstructure.
===============================================================

Background
----------

The related abinit variables
----------------------------

    * 1
    * 2

More info on the inputvariables and their use can be obtained using the following function:

    .. code-block :: python
        lesson.docvar("inputvariable")

This will print the official abinit description of this inputvariable.

The abipy flows in this lesson
------------------------------

The flow that we use in this lesson contains for the first time dependencies.
This means that some tasks in the flow can only be started if an other task is
ready. We will first perform one selfconsistend calculation to obtain a proper
density. Using this density we calculate

The cource of this lesson
-------------------------

Start this lessen by importing it in a new namespace:

    .. code-block :: python
        from abipy.lesson.lesson_dos_bands import Lesson()
        lesson = Lesson()

As always you can reread this lessons text using the command:

    .. code-block :: python
        lesson

To build the flow:

    .. code-block :: python
        flow = lesson.make_flow()

To print the input files

    .. code-block :: python
        flow.show_inputs()

Start the flow with the scheduler and wait for completion.

    .. code-block :: python
        flow.make_scheduler().start()

To analyze the results.

    .. code-block :: python
        flow.analyze()



Exercises
---------



Next
----

A logical next lesson would be lesson_g0w0


"""
from __future__ import division, print_function

import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.lessons.core import BaseLesson, get_pseudos


class EbandsDosFlow(abilab.Flow):
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

    def plot_ebands_with_edos(self, dos_idx=0, **kwargs):
        # plot dos
        with self.nscf_task.open_gsr() as gsr: 
            gs_ebands = gsr.ebands

        with self.dos_tasks[dos_idx].open_gsr() as gsr: 
            dos_ebands = gsr.ebands

        edos = dos_ebands.get_edos(method="gaussian", step=0.01, width=0.1)
        return gs_ebands.plot_with_edos(edos, **kwargs)


def make_electronic_structure_flow(ngkpts_for_dos=[(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]):
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
                                     dos_inputs=dos_input, flow_class=EbandsDosFlow)


class Lesson(BaseLesson):

    @property
    def doc_string(self):
        return __doc__

    @property
    def pyfile(self):
        return os.path.basename(__file__[:-1])

    @staticmethod
    def make_electronic_structure_flow():
        return make_electronic_structure_flow()


if __name__ == "__main__":
    l = Lesson()
    print(l.pyfile)

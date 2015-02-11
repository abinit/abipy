#!/usr/bin/env python
"""
The calculation of the density of states and the bandstructure.
===============================================================

Background
----------

This lesson focuses on calculating the density of states (DOS) and the band structure. On thing one should always
keep in mind is that these are the densities of states





The related abinit variables
----------------------------

    * 1
    * 2

"""
from __future__ import division, print_function


_ipython_lesson_ = """
More info on the inputvariables and their use can be obtained using the following function:

    .. code-block :: python

        lesson.docvar("inputvariable")

This will print the official abinit description of this inputvariable.

The abipy flows in this lesson
------------------------------

The flow that we use in this lesson contains for the first time dependencies.
This means that some tasks in the flow can only be started if an other task is
ready. We will first perform one self-consistent calculation to obtain a proper
density. Using this density we calculate in two more steps the DOS and the bandstructure.
For the DOS this not stricktly nessesary since the DOS will also be calculated on a regular grid.
In general the density will be converged already before the DOS is converged. For large systems it may become
nessesary to split. For the bandstructure we have a non-uniform grid so we do need to fix the density.

The course of this lesson
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

To visualize the dependencies in the flow:

    .. code-block :: python

        flow.show_dependencies()

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

_commandline_lesson_ = """
At this place they will not be discussed in detail. In stead you are
invited to read the abinit documentation on them. The full description,
directly from the abinit description is available via the following function:

    .. code-block :: shell

        abidocs.py man inputvariable

This will print the official abinit description of this inputvariable.


The course of this lesson
-------------------------
"""

import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.lessons.core import BaseLesson, get_pseudos


def make_electronic_structure_flow(ngkpts_for_dos=[(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]):
    """Band structure calculation."""
    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=2 + len(ngkpts_for_dos))
    inp.set_structure(abidata.cif_file("si.cif"))

    # Global variables
    inp.ecut = 10

    # Dataset 1
    inp[1].set_vars(tolvrs=1e-9)
    inp[1].set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0])

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
                                     dos_inputs=dos_input)


class Lesson(BaseLesson):

    @property
    def abipy_string(self):
        return __doc__+_ipython_lesson_

    @property
    def comline_string(self):
        return __doc__+_commandline_lesson_

    @property
    def pyfile(self):
        return os.path.basename(__file__)

    @staticmethod
    def make_electronic_structure_flow(**kwargs):
        return make_electronic_structure_flow(**kwargs)

    @staticmethod
    def analyze(my_flow):
        nscf_task = my_flow[0][1]
        with nscf_task.open_gsr() as gsr:
            return gsr.ebands.plot()

if __name__ == "__main__":
    l = Lesson()
    flow = l.make_electronic_structure_flow()
    flow.build_and_pickle_dump()
    l.manfile(l.comline_string)
    l.instruct()
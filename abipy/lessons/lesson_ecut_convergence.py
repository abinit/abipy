#!/usr/bin/env python
"""
Basis set convergence study and some more on flows, works, and tasks. 
=====================================================================

Background
----------

This lesson focuses on the convergence study of the completeness
of the basis set used. In our case the basis set consists of plane
waves. Plane waves are inherently well suited to capture the periodic
nature of a crystalline solid. In addition a plane wave basis set
has the advantage that it introduces only one convergence parameter,
the kinetic energy cutoff.

The sharp features of the wavefunctions near the nucleus are however
problematic for plane waves. Describing these features would require
very high frequency plane waves. In practice we will always use
pseudo-potentials in stead of the actual nuclear potential to improve
convergence. Effectively a pseudopotential replaces the sharp coulomb
potential of the nucleus and the core electrons by something more smooth
inside the pseudization region that connects smoothly to the real potential
outside the pseudization region.

Needless to say a different pseudo potential will require a different
cutoff for the calculation to be converged. In general norm-conserving
pseudos require a larger cut-off that ultra-soft pseudos and Projector
Augmented Wave 'pseudos' require even smaller cutoffs. Moreover two
pseudos of the same type for the same element may require different
cutoffs as well. 'Harder' (having a smaller pseudization radius) require
larger cutoffs than 'softer' pseudos. There are however many more
properties of a pseudo that determine the cutoff needed.

The related abinit variables
----------------------------

As said the most important parameter in the energy cutoff, in abinit ecut.
The most important input parameters concerning the basis set are:

    * ecut
    * dilatms
    * ecutsm
    * ecutdg


More info on the inputvariables and their use can be obtained using the
following function:

    .. code-block :: python

        lesson.docvar("inputvariable")

The abipy flows in this lesson
------------------------------

The course of this lesson
-------------------------

Start this lesson by importing it in a new namespace

    .. code-block :: python

        from abipy.lessons.lesson_ecut_convergence import Lesson()
        lesson = Lesson()

As always you can reread this lesson's text using the command:

    .. code-block :: python

        lesson

To build the flow:

    .. code-block :: python

        flow = lesson.make_ecut_flow()

To print the input files

    .. code-block :: python

        flow.show_inputs()

In this lesson we take a closer look at the structure of a Flow. In general
a flow is a container that contains 'works'. Works are (connected) series
of abinit executions we call tasks. To show the works contained in a flow
use the 'works()' method:

    .. code-block :: python

        flow.works()

to show the status of a flow:

    .. code-block :: python

        flow.show_status()

There are many more properties and methods of a flow than may also come in
handy. By typing [tab] in ipython after the period, you will be presented
with all the option. Feel free to experiment a bit at this point. By adding
a questionmark to the method or property ipython will show the information
and description of it:

    .. code-block :: python

        flow.open_files?

Will explain what this method is supposed to do.

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

A logical next lesson would be lesson_relaxation

"""
from __future__ import division, print_function

import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.lessons.core import BaseLesson


class EcutFlow(abilab.Flow):
    def analyze(self):
        with abilab.abirobot(self, "GSR") as robot:
            data = robot.get_dataframe()
            #robot.ebands_plotter().plot()

        import matplotlib.pyplot as plt
        data.plot(x="ecut", y="energy", title="Total energy vs ecut", legend="Energy [eV]", style="b-o")
        plt.show()


def make_ecut_flow(structure_file=None, ecut_list = (10, 12, 14, 16, 18)):
    #define the structure and add the necessary pseudos:
    if structure_file is None:
        inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ecut_list))
        inp.set_structure(abidata.cif_file("si.cif"))
        workdir = "lesson_Si_ecut_convergence"
    else:
        structure = abilab.Structure.from_file(structure_file)
        pseudos = abilab.PseudoTable()  ## todo fix this
        inp = abilab.AbiInput(pseudos=pseudos, ndtset=len(ecut_list))
        inp.set_structure(structure)
        workdir = "lesson_" + structure.composition.reduced_formula + "_ecut_convergence"

    # Global variables
    inp.set_vars(tolvrs=1e-9)
    inp.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0])

    for i, ecut in enumerate(ecut_list):
        inp[i+1].set_vars(ecut=ecut)

    return EcutFlow.from_inputs(workdir=workdir, inputs=inp.split_datasets())


class Lesson(BaseLesson):

    @property
    def doc_string(self):
        return __doc__

    @property
    def pyfile(self):
        return os.path.basename(__file__[:-1])

    @staticmethod
    def make_flow():
        return make_ecut_flow()


if __name__ == "__main__":
    l = Lesson()
    print(l.pyfile)

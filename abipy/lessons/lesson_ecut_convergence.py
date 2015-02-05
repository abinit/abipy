#!/usr/bin/env python
"""
\033[91m Basis set convergence study and some more on flows, works, and tasks. \033[0m

\033[94m Background\033[0m

This lesson focuses on the convergence study of the completeness
of the basis set used. In our case the basis consists of plane
waves. Plane waves are inherently well suited to capture the periodic
nature of a crystalline solid. In addition a plane wave basis set
has the advantage that it introduces only one convergence parameter,
the kinetic energy cutoff.

The sharp features of the wavefunctions near the nucleus are however
problematic for plane waves. Describing these features would require
very high frequency plane waves. In practice we will always use
pseudo-potentials in stead of the actual nuclear potential to improve
convergence. Effectively a speudopotential replaces the sharp coulomb
potential of the nucleus and the core electrons by someting more smooth
inside the pseudization region that connects smoothly to the real potential
outside the pseudization region.

Needless to say a different pseudo potential will require a different
cutoff for the calculation to be converged. In general normconserving
pseudos require a larger cut-off that ultra-soft pseudos and Projector
Augmented Wave 'pseudos' require even smaller cutoffs. Moreover two
pseudo's of the same type for the same element may requier different
cutoffs as well. 'Harder' (having a smaller pseudization radius) require
larger cutoffs than 'softer' pseudo's. There are however many more
properties of a pseudo tha determine the cutoff needed.

\033[94m The related abinit variables\033[0m

As said the most important parameter in the energy cutoff, in abinit ecut.

\033[1m ecut \033[0m
\033[1m dilatms \033[0m
\033[1m ecutsm \033[0m
\033[1m \033[0m

More info on the inputvariables and their use can be obtained using the
following function:

\033[92m In []:\033[0m lesson.abinit_help(inputvariable)

\033[94m The abipy flows in this lesson \033[0m

\033[94m The course of this lesson \033[0m

Start this lessen by importing it in a new namespace

\033[92m In []:\033[0m from abipy.lessons import lesson_ecut_convergence as lesson

As always you can reread this lesson's text using the command:

\033[92m In []:\033[0m lesson.help()

To build the flow:

\033[92m In []:\033[0m flow = lesson.make_ecut_flow()

To print the input files

\033[92m In []:\033[0m flow.show_inputs()

In this lesson we take a closer look at the structure of a Flow. In general
a flow is a container that contains 'works'. Works are (connected) series
of abinit executions we call tasks. To show the works contained in a flow
use the 'works()' method:

\033[92m In []:\033[0m flow.works()

to show the status of a flow:

\033[92m In []:\033[0m flow.show_status()

There are many more properties and methods of a flow than may also come in
handy. By typing [tab] in ipython after the period, you will be presented
with all the option. Feel free to experiment a bit at this point. By adding
a questionmark to the method or property ipython will show the information
and description of it:

\033[92m In []:\033[0m flow.open_files?

Will explain what this method is supposed to do.

Start the flow with the scheduler and wait for completion.

\033[92m In []:\033[0m flow.make_scheduler().start()

To analyze the results.

\033[92m In []:\033[0m flow.analyze()

\033[93m Exercises \033[0m



\033[93m Next \033[0m

A logical next lesson would be lesson_relaxation

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


class EcutFlow(abilab.Flow):
    def analyze(self):
        with abilab.abirobot(self, "GSR") as robot:
            data = robot.get_dataframe()
            #robot.ebands_plotter().plot()

        import matplotlib.pyplot as plt
        data.plot(x="ecut", y="energy", title="Total energy vs ecut", legend="Energy [eV]", style="b-o")
        plt.show()


def make_ecut_flow(structure_file=None):
    ecut_list = [10, 12, 14, 16, 18]

    #define the structure and add the nessesary pseudo's:
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
    inp.set_variables(tolvrs=1e-9)
    inp.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0])

    for i, ecut in enumerate(ecut_list):
        inp[i+1].set_variables(ecut=ecut)

    return EcutFlow.from_inputs(workdir=workdir, inputs=inp.split_datasets())


if __name__ == "__main__":
    flow = make_ecut_flow()
    flow.make_scheduler().start()
    flow.analyze()

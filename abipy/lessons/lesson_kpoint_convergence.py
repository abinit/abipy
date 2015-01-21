#!/usr/bin/env python
"""
\033[91m K-point convergence study for a semi-conductor and an introduction some of the basic concepts of the abipy library. \033[0m

\033[94m Background\033[0m

This lesson deals with the basic k-point convergence study that is needed in any DFT calculation of a solid. In a DFT
calculation of a solid the first Brillouin zone needs to be discritized to enable the integration of various quantities.
Effectively the integrals are turned in to sums over k-points. For any result to be converged we need a k-point mesh
that is dense enough, but at the same time as corse as possible to make for an efficient calculation. Various types
of materials require in general different densities of the k-point meshes. In general metals need denser meshes than
semiconductors. Your first investigation into a new compound will quit often be a k-point convergence study.

\033[94m The related abinit variables\033[0m

Yannick could you fill this? maybe even an interactive interface to you input parameter database work?

\033[1m ngkpt \033[0m
\033[1m kptshift \033[0m
\033[1m kptopt \033[0m


\033[94m The abipy flows of this lesson\033[0m

When performed manually, a k-point convergence study would require the preparation of a series of input-files, running
abinit for all the inputs and extracting and studying the quantity that is needed to be converged. This lesson shows
how this process can be greatly facilitated by using python scripts in the abipy framework. We fill construct a single
python object, a abipy flow, that contains all the information needed for the calculations but also provides methods
for acually running abinit, inspecting the input and output, and analyzing the results.

\033[94m The Course of this lesson\033[0m

This lesson can be started in ipython by importing it:

\033[92m In []:\033[0m from abipy.lessons import lesson_Si_kpoint_convergence as lesson

The lesson is now imported in your ipython session in its own namespace 'lesson'. This object now gives us all the
tools to follow this lesson. For instance the command:

\033[92m In []:\033[0m lesson.help()

displays this lessons information text, and can be recalled at any moment. The main object we use to pack
(connected series of) calculations is a flow. This lesson provides a method that returns a flow designed to perform
k-point convergence studies. This flow is made by the command:

\033[92m In []:\033[0m flow = lesson.make_ngkpt_flow()

'flow' is now an object that contains al the information needed to generate abinit input. In this case it is a special
flow for a k-point convergence study and since we did not specify anything when generating the flow the example case
of silicon is generated. Our flow however inherited from the abinit base flow so we have a lot of 'standard' methods
available. For instance:

\033[92m In []:\033[0m flow.show_inputs()

This will display all the inputs as they will be 'given' to abinit.

To start the the execution of calculations packed in this flow we an use the following command:

\033[92m In []:\033[0m flow.make_scheduler().start()

This starts the actual execution via a scheduler. The scheduler is a sort of daemon that starts to submit tasks that
are ready to run. In our case al the tasks in the flow are independent so the first cycle of the scheduler directly
submitted all of them. More complicated flows may have tasks that can only start using input from a previous task. We
will encounter some of those later.

The last step of analyzing the results can be done again in with a single command:

\033[92m In []:\033[0m flow.analyze()

This method of flow will open the necessary output files, retrieve the data, and produce a plot.

Finally, once you are through with this lesson and exited ipython:

\033[92m In []:\033[0m exit

You can see that in the directory that you were working there is now a subdir were the calculation have been performed.
Have a look at these folders and the files that are in them.


\033[93m Exercises \033[0m

As an exercise you can now start this lesson again but in stead of performing the convergence study for silicon study
the convergence for a metal. By using:

\033[92m In []:\033[0m flow = lesson.make_ngkpt_flow(structure_file=abidata.cif_file('al.cif'))

you will generate a flow for aluminum. Actually, you can pass the path to any cif file to perform a convergence study
on that material.

If you have time left it is also a good exercise to open the python file that contains this lesson and study the
implementations of the classes, methods and functions we used. You can get a copy of the file by using:

\033[92m In []:\033[0m lesson.get_local_copy()

Try to find what to change to change the set of k-point meshes that are used in the convergence study.
"""

from __future__ import division, print_function

import sys
import os
import shutil
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.lessons.lesson_helper_functions import abinit_help


# should n't we put these functions in a separate module and import them, we'll need them in any lesson...

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

class NgkptFlow(abilab.Flow):
    """
    A flow class for the study of k-point convergence studies. It inherits from the base class abilab.Flow, and in
    addion implements a method for analyzing specifically the k-point convergence data.
    """
    def analyze(self):
        with abilab.abirobot(self, "GSR") as robot:
            data = robot.get_dataframe()
            #robot.ebands_plotter().plot()

        import matplotlib.pyplot as plt
        data.plot(x="nkpts", y="energy", title="Total energy vs nkpts", legend="Energy [eV]", style="b-o")
        plt.show()


def make_ngkpt_flow(structure_file=None):
    """
    A 'factory function' (a function that returns an instance of the class defined above. If no specific system is
    specified, structure_file=None, an example flow for silicon in constructed and returned.
    """
    ngkpt_list = [(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]

    # Defining the structure and adding the appropriate pseudo potentials
    if structure_file is None:
        inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
        inp.set_structure(abidata.cif_file("si.cif"))
        workdir = "lesson_Si_kpoint_convergence"
    else:
        structure = abilab.Structure.from_file(structure_file)
        pseudos = abilab.p  ## todo fix this
        inp = abilab.AbiInput(pseudos=pseudos, ndtset=len(ngkpt_list))
        inp.set_structure(structure)
        workdir = "lesson_" + structure.composition.reduced_formula + "_kpoint_convergence"

    # Global variables
    inp.set_variables(ecut=10, tolvrs=1e-9)

    # Specific variables for the different calculations
    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])

    return NgkptFlow.from_inputs(workdir=workdir, inputs=inp.split_datasets())


if __name__ == "__main__":
    """
    this section of code will be executed when this script file would be executed directly. It can be seen as an
    example how the above method could be used in you onw scripts
    """
    flow = make_ngkpt_flow()
    flow.make_scheduler().start()
    flow.analyze()

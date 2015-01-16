#!/usr/bin/env python
"""
K-point convergence study for a semi-conductor and an introduction some of the basic concepts of the abipy library.

  Background.

This lesson deals with basic k-point convergence study that is needed in any DFT calculation of a solid. In a DFT
calculation of a solid the first Brillouin zone needs to be discritized to enable the integration of various quantities.
Effectively these integrals are turned in to sums over k-point. For any result to be converged we need a k-point mesh
that is dense enough, but at corse as possible to make for an efficient calculation. Various types of materials require
in general different densities of the k-point meshes. Metals need denser meshes than semiconductors. Your first
investigation into a new compound will quit often be a k-point convergence study.


  The related abinit variables.

Yannick could you fill this? maybe even an interactive interface to you input parameter database work?


  The abipy flows of this lesson.

When performed manually, a k-point convergence study would require the preparation of a series of input-files, running
abinit for all the inputs and extracting and studying the parameter that is needed to be converged. This lesson shows
how this process can be greatly facilitated by using python scripts in the abipy framework. We fill construct a single
python object, a abipy flow, that contains all the information needed for the calculations but also provides methods
for acually running abinit, inspecting the input and output, and analyzing the results.


  The Course of this lesson.

This lesson can be started in ipython by importing it:

[]: from abipy.lessons import lesson_Si_kpoint_convergence as lesson

The lesson is now imported in your ipython session in its own namespace 'lesson'. This object now gives us all the
tools to follow this lesson. For instance the command:

[]: lesson.help()

displays this lessons information text, and can be recalled at any moment. The main object we use to pack
(connected series of) calculations is a flow. This lesson provides a method that returns a flow designed to perform
k-point convergence studies. The flow is made by the command:

[]: flow = lesson.make_ngkpt_flow()

'flow' is now an object that contains al the information needed to generate abinit input. In this case it is a special
flow for a k-point convergence study and since we did not specify any thing in generating the flow the example case
of silicon is generated. Our flow however inherited from the base abinit flow so we have a lot of 'standard' methods
vailable. For instance:

[]: flow.show_inputs()

This will display all the inputs as they will be 'given' to abinit.

To start the the execution of calculations packed in this flow we an use the following command:

[]: flow.make_scheduler().start()

This starts the actual execution via a scheduler. The scheduler 'knows' about the computer we are on and how to run
calculations on this computer. In a later lesson we will revisit this topic in more detail. Now we just wait for the
calculations to complete.

The last step of analyzing the results can be done again in with a single command:

[]: flow.analyze()

This method of flow will open the nessesary output files, retrieve the data, and produce a plot.

Finally, once you are through with this lesson and exited ipython:

[]: exit

You can see that in the directory that you were working there is now a subdir were the calculation have been performed.
Have a look at these folders and the files that are in them.


Exercises

As an exercise you can now start this lesson again but in stead of performing the convergence study for silicon study
the convergence for a metal. By using:

[]: flow = lesson.make_ngkpt_flow(structure_file=abidata.cif_file('al.cif'))

you will generate a flow for aluminum. Actually, you can pass the path to any cif file to perform a

If you have time left it is also a good exercise to open the python file that contains this lesson and study the
implementations of the class, methods and functions we used. Try to find what to change to change the set of k-point
meshes that are used in the convergence study.
"""

from __future__ import division, print_function

import sys
import abipy.abilab as abilab 
import abipy.data as abidata


def help(stream=sys.stdout):
    """
    Display the tutorial text.
    """
    stream.write(__doc__)


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
        pseudos = abilab.PseudoTable()  ## todo fix this
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

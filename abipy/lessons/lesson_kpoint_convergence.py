#!/usr/bin/env python
"""
K-point convergence study for a semi-conductor and an introduction some of the basic concepts of the abipy library.
===================================================================================================================

Background
----------

This lesson deals with the basic k-point convergence study that
is needed in any DFT calculation of a solid. In a DFT calculation
of a solid the first Brillouin zone needs to be discretized to
enable the integration of various quantities. Effectively the
integrals are turned in to sums over k-points. For any result to
be converged we need a k-point mesh that is dense enough, but at
the same time as coarse as possible to make for an efficient
calculation. Various types of materials require in general different
densities of the k-point meshes. In general metals need denser meshes
than semiconductors. Your first investigation into a new compound
will quit often be a k-point convergence study.

The related abinit variables
----------------------------

The abinit parameters concerned with the k-point grid are:

    * ngkpt
    * shiftk
    * occopt (see exercises)
    * tsmear (see exercises)
    * kptopt (see exercises)

At this place they will not be discussed in detail. In stead you are
invited to read the abinit documentation on them. The full description,
directly from the abinit description is available via the following function:

    .. code-block :: python

        lesson.docvar("inputvariable")

This will print the official abinit description of this inputvariable.

The abipy flows of this lesson
------------------------------

When performed manually, a k-point convergence study would require
the preparation of a series of input-files, running abinit for all
the inputs and extracting and studying the quantity that is needed
to be converged. This lesson shows how this process can be greatly
facilitated by using python scripts in the abipy framework. We will
construct a single python object, a abipy flow, that contains all
the information needed for the calculations but also provides methods
for actually running abinit, inspecting the input and output, and
analyzing the results.

All calculations will however still be run in parallel.

The Course of this lesson
-------------------------

This lesson can be started in ipython by importing it:

    .. code-block :: python

        from abipy.lessons.lesson_kpoint_convergence import Lesson
        lesson = Lesson()

The lesson is now imported in your ipython session in its own
namespace 'lesson'. This object now gives us all the tools to
follow this lesson. For instance the command:

    .. code-block :: python

        lesson

displays this lessons information text, and can be recalled at
any moment. The main object we use to pack (connected series of)
calculations is a flow. This lesson provides a method that returns
a flow designed to perform k-point convergence studies. This flow
is made by the command:

    .. code-block :: python

        flow = lesson.make_ngkpt_flow()

'flow' is now an object that contains al the information needed
to generate abinit input. In this case it is a special flow for
a k-point convergence study and since we did not specify anything
when generating the flow the example case of silicon is generated.
Our flow however, inherited from the abinit base flow so we have
a lot of 'standard' methods available. For instance:

    .. code-block :: python

        flow.show_inputs()

This will display all the inputs as they will be 'given' to abinit.

To start the execution of calculations packed in this flow we
and use the following command:

    .. code-block :: python

        flow.make_scheduler().start()

This starts the actual execution via a scheduler. The scheduler is
a sort of daemon that starts to submit tasks that are ready to run.
In our case all the tasks in the flow are independent so the first
cycle of the scheduler directly submitted all of them. More
complicated flows may have tasks that can only start using input
from a previous task. We will encounter some of those later.

The last step of analyzing the results can be done again in with
a single command:

    .. code-block :: python

        flow.analyze()

This method of flow will open the necessary output files, retrieve
the data, and produce a plot.

Finally, once you are through with this lesson and exited ipython:

    .. code-block :: python

        exit

You can see that in the directory that you were working there is
now a subdir were the calculation have been performed. Have a look
at these folders and the files that are in them.


Exercises
---------

As an exercise you can now start this lesson again but in stead
of performing the convergence study for silicon study the
convergence for a metal. By using:

    .. code-block :: python

        flow = lesson.make_ngkpt_flow(structure_file=lesson.abidata.cif_file('al.cif'), metal=True)

you will generate a flow for aluminum. Actually, you can pass
the path to any cif file to perform a convergence study on that
material. Be careful however, aluminum is a metal and the default
parameters for occopt and tsmear are for semiconductors. The
keyword argument 'metal' fixes this. (you could also see what
happens if you don't put this flag :-) ) Look at the inputs to
see what has been changed and study the description of these
inputvariables using the abinit_help() method.

If you have time left it is also a good exercise to open the
python file that contains this lesson and study the implementations
of the classes, methods and functions we used. You can get a copy
of the file by using:

    .. code-block :: python

        lesson.get_local_copy()

Try to find what to change to change the set of k-point meshes that
are used in the convergence study.

Next
----

A logical next lesson would be lesson_ecut_convergence
"""

from __future__ import division, print_function

import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.lessons.core import BaseLesson, get_pseudos


class NgkptFlow(abilab.Flow):
    """
    A flow class for the study of k-point convergence studies. It inherits from the base class abilab.Flow, and in
    addition implements a method for analyzing specifically the k-point convergence data.
    """
    def analyze(self):
        with abilab.abirobot(self, "GSR") as robot:
            data = robot.get_dataframe()
            #robot.ebands_plotter().plot()

        import matplotlib.pyplot as plt
        data.plot(x="nkpts", y="energy", title="Total energy vs nkpts", legend=False, style="b-o")
        plt.show()


def make_ngkpt_flow(ngkpt_list=[(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)], structure_file=None, metal=False):
    """
    A `factory function` (a function that returns an instance of the class defined above. If no specific system is
    specified, structure_file=None, an example flow for silicon in constructed and returned.
    """
    # Defining the structure and adding the appropriate pseudo potentials
    if structure_file is None:
        inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
        inp.set_structure(abidata.cif_file("si.cif"))
        workdir = "flow_lesson_Si_kpoint_convergence"
    else:
        structure = abilab.Structure.from_file(structure_file)

        pseudos = get_pseudos(structure)
        inp = abilab.AbiInput(pseudos=pseudos, ndtset=len(ngkpt_list))
        inp.set_structure(structure)
        workdir = "flow_lesson_" + structure.composition.reduced_formula + "_kpoint_convergence"

    # Global variables
    inp.set_vars(ecut=10, tolvrs=1e-9)

    if metal:
        inp.set_vars(occopt=7, tsmear=0.04)

    # Specific variables for the different calculations
    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])

    return NgkptFlow.from_inputs(workdir=workdir, inputs=inp.split_datasets())


class Lesson(BaseLesson):

    @property
    def doc_string(self):
        return __doc__

    @property
    def pyfile(self):
        return os.path.basename(__file__[:-1])

    @staticmethod
    def make_flow():
        return make_ngkpt_flow()


if __name__ == "__main__":
    l = Lesson()
    print(l.pyfile)

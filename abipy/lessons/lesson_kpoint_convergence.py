#!/usr/bin/env python
"""
K-point convergence study for a semi-conductor. An introduction to some of the basic concepts of AbiPy
======================================================================================================

Background
----------

This lesson deals with the basic k-point convergence study that is needed in any DFT calculation for periodic systems. 
In such a systems, indeed, the first Brillouin zone (BZ) needs to be discretized when performing the 
integration of several important quantities e.g. the electronic density or the electronic energy.
Integrals over the BZ are therefore turned into sums over discrete k-points and the k-mesh should 
be dense enough, but at the same time as coarse as possible to make for an efficient calculation. 
Your first investigation into a new compound will quit often be a k-point convergence study.
It is worth stressing that the density of the k-point mesh needed to reach converged results is system-dependent.

Note moreover that metals need much denser meshes than semiconductors. 
The presence of the Fermi surface, indeed, induces discontinuities in the integrand functions and a 
fictitious broadening of the occupation factors (tsmear) should be introduced in order to accelerate 
the convergence of the integrals.

The related Abinit variables
----------------------------

The variables used to specify the k-point sampling are:

    * ngkpt
    * shiftk
    * kptopt (see exercises)

The variables used to specify the occupation scheme in metals are:

    * occopt (see exercises)
    * tsmear (see exercises)


"""
from __future__ import division, print_function

_ipython_lesson_ = """
At this place they will not be discussed in detail, you are invited to consult the abinit documentation. 
The full description, directly from the official abinit docs, is available via the following function:

    .. code-block:: python

        lesson.docvar("inputvariable")

that will print the official description of inputvariable.

The abipy flows of this lesson
------------------------------

When performed manually, a k-point convergence study would require
the preparation of several input files, running abinit for all
the inputs and then extracting and studying the quantity under investigation.
This lesson shows how this process can be greatly facilitated thanks to abipy. 
We will construct a single python object, an abipy flow, that contains all
the information needed for the calculations but also provides methods
for actually running abinit, inspecting the input and output, and analyzing the final results.

The Course of this lesson
-------------------------

This lesson can be started in ipython by importing it:

    .. code-block:: python

        from abipy.lessons.lesson_kpoint_convergence import Lesson
        lesson = Lesson()

The lesson is now imported in your ipython session in its own namespace 'lesson'. 
This object now gives us all the tools to follow this lesson. 
For instance the command:

    .. code-block:: python

        lesson

displays this lessons information text, and can be recalled at
any moment. The main object we use to connect different calculations is a flow. 
This lesson provides a method that returns a flow that performs k-point convergence studies. 
This flow is made by the command:

    .. code-block:: python

        flow = lesson.make_ngkpt_flow()

'flow' is now an object that contains all the information needed
to generate abinit input. In this case it is a special flow for
a k-point convergence study and since we did not specify anything
when generating the flow the example case of silicon is generated.
Our flow however, inherited from the abinit base flow so we have
a lot of 'standard' methods available. For instance:

    .. code-block:: python

        flow.show_inputs()

This will display all the inputs as they will be 'given' to abinit.

To start the execution of calculations packed in this flow we
and use the following command:

    .. code-block:: python

        flow.make_scheduler().start()

This starts the actual execution via a scheduler. The scheduler is
a sort of daemon that starts to submit tasks that are ready to run.
In our case all the tasks in the flow are independent so the first
cycle of the scheduler directly submitted all of them. More
complicated flows may have tasks that can only start using input
from a previous task. We will encounter some of those later.

The last step of analyzing the results can be done again in with a single command:

    .. code-block:: python

        lesson.analyze(flow)

This method of flow will open the necessary output files, retrieve
the data, and produce a plot.

Finally, once you are through with this lesson and exited ipython:

    .. code-block:: python

        exit

You can see that in the directory that you were working there is
now a subdir were the calculation have been performed. Have a look
at these folders and the files that are in them.


Exercises
---------

As an exercise you can now start this lesson again but instead
of performing the convergence study for silicon, we could perform a 
similar convergence analysis for a metal. 
By using:

    .. code-block:: python

        flow = lesson.make_ngkpt_flow(structure_file=lesson.abidata.cif_file('al.cif'), metal=True)

you will generate a flow for aluminum. Actually, you can pass
the path to any CIF file to perform a convergence study on that
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

    .. code-block:: python

        lesson.get_local_copy()

Try to find what to change to change the set of k-point meshes that
are used in the convergence study.

Next
----

A logical next lesson would be lesson_ecut_convergence
"""

_commandline_lesson_ = """
At this place they will not be discussed in detail. Instead you are
invited to read the abinit documentation on them. The full description,
directly from the abinit description is available via the following function:

    .. code-block:: shell

        abidoc.py man inputvariable

This will print the official abinit description of this inputvariable.


The course of this lesson
-------------------------

In the generation of this lesson by the python script all the input files have been generated automatically.
The input files have been organized in a workdir "flow_lesson_Si_kpoint_convergence". 
Inside this directory, you'll find a single work, w0, with four tasks, t0 - t3. 
Have a look at the input files, run.abi, of the four tasks to see what is different.

You'll see that also the files file and the jobs submission script are generated. 
In the job scripts you'll see that the jobs are prepared to run just on the front end.

You'll also see that the files file has been created as well.

To perform the k-point convergence study execute abinit with the four input sets.

Once the calculations are ready, you'll see three important output files.

    * run.out
    * run.log
    * run.err

The main summary of the calculation can be found in the .out file, we'll go there soon :-). The .err file should be
empty. If it's not something went wrong. If something went wrong read the .err. file. The .log file contains extensive
information on you calculation that could help to find out what went wrong in the case of errors. Especially there are
three types of messages that could help

    * COMMENT
    * WARNING
    * ERROR

In case of an error message abinit stopped the execution by itself, because of that error.

Now the .out file. Some interesting keywords to look for:

    * Symmetries
    * Citation for XC functional:
    * ETOT (the total energies during the electronic structure convergence)
    * Eigenvalues
    * Etotal (the total energy of an ionic step)

Obviously there is much more.

Collect the total energies of the four calculations and plot them as a function of the number of k-points in the
calculation.

Alternative to execution of the manual execution the calculations can also be executed using the abipy scheduler.

    .. code-block:: shell

        abirun.py flow_lesson_Si_kpoint_convergence scheduler

"""
import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.lessons.core import BaseLesson, get_pseudos


def make_ngkpt_flow(ngkpt_list=[(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)], structure_file=None, metal=False):
    """
    A `factory function` (a function that returns an instance of the class defined above. If no specific system is
    specified, structure_file=None, an flow for silicon in constructed and returned.
    """
    # Defining the structure and adding the appropriate pseudo potentials
    if structure_file is None:
        multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                    pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ngkpt_list))
        workdir = "flow_lesson_Si_kpoint_convergence"
    else:
        structure = abilab.Structure.from_file(structure_file)
        pseudos = get_pseudos(structure)
        multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=len(ngkpt_list))
        workdir = "flow_lesson_" + structure.composition.reduced_formula + "_kpoint_convergence"

    # Add mnemonics to input file.
    multi.set_mnemonics(True)

    # Global variables
    multi.set_vars(ecut=10, tolvrs=1e-9)

    if metal:
        multi.set_vars(occopt=7, tsmear=0.04)

    # Specific variables for the different calculations
    for i, ngkpt in enumerate(ngkpt_list):
        multi[i].set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])

    return abilab.Flow.from_inputs(workdir=workdir, inputs=multi.split_datasets())


class Lesson(BaseLesson):

    @property
    def abipy_string(self):
        return __doc__ + _ipython_lesson_

    @property
    def comline_string(self):
        return __doc__ + _commandline_lesson_

    @property
    def pyfile(self):
        return os.path.abspath(__file__).replace(".pyc", ".py")

    @staticmethod
    def make_ngkpt_flow(**kwargs):
        return make_ngkpt_flow(**kwargs)

    @staticmethod
    def analyze(my_flow, **kwargs):
        with abilab.abirobot(my_flow, "GSR") as robot:
            data = robot.get_dataframe()
        import matplotlib.pyplot as plt
        ax = data.plot(x="nkpts", y="energy", title="Total energy vs nkpts", legend=False, style="b-o")
        ax.set_xlabel('Number of k-points')
        ax.set_ylabel('Total Energy [eV]')
        return plt.show(**kwargs)


if __name__ == "__main__":
    l = Lesson()
    flow = l.make_ngkpt_flow()
    flow.build_and_pickle_dump()
    l.setup()

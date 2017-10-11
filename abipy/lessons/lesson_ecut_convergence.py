#!/usr/bin/env python
"""
Basis set convergence study and more info on flows, works, and tasks
====================================================================

Background
----------

This lesson focuses on the convergence study on the completeness of the plane wave (PW) basis set.
Plane waves are inherently well suited to capture the periodic nature of crystalline solids.
In addition, a PW basis set has the advantage that it introduces only one convergence parameter,
the kinetic energy cutoff (ecut).

The sharp features of the wavefunctions near the nucleus are however problematic for PWs.
Describing these features would require very high energy cutoff energies.
For this reason PW codes use pseudo-potentials in order to facilitate the convergence of the results.
A pseudopotential replaces the singular coulomb potential of the nucleus and the
core electrons by something smoother inside the so-called pseudization region.
The pseudopotential connects smoothly to the real all-electron potential outside the pseudization region.

Note that different pseudo potentials usually require a different cutoff energy to be converged.
In general norm-conserving pseudos require a larger cut-off than ultra-soft pseudos
or Projector Augmented Wave 'pseudos'.
Moreover two pseudos of the same type for the same element may require different  cutoff energies as well.
Pseudos with small pseudization radius usually require larger cutoffs than pseudos
with large pseudization radius.

The related abinit variables
----------------------------

    * ecut        (cutoff energy)
    * pawecutdg   (additional variable for the double-grid used in PAW)
    * ecutsm      (smoothing of the kinetic energy)

"""
from __future__ import division, print_function, unicode_literals, absolute_import

_ipython_lesson_ = """
More info on the input variables and their use can be obtained using:

    .. code-block:: python

        print(lesson.docvar("inputvariable"))

Description of the lesson
-------------------------
This lesson contains a factory function for a convergence study with respect to ecut.

Executing the lesson
--------------------

Start this lesson by importing it:

    .. code-block:: python

        from abipy.lessons.lesson_ecut_convergence import Lesson
        lesson = Lesson()

As usual, you can reread this text using the command:

    .. code-block:: python

        lesson

To build the flow:

    .. code-block:: python

        flow = lesson.make_ecut_flow()

To print the input files

    .. code-block:: python

        flow.show_inputs()

to show the status of a flow:

    .. code-block:: python

        flow.show_status()

There are many more properties and methods of the flow than may also come in handy.
By typing [tab] in ipython after the period, you will be presented
with all the options. Feel free to experiment a bit at this point.
In the ipython shell, one can get the description of the object by
adding a question mark at the end of the statement:

    .. code-block:: python

        flow.show_status?

Start the flow with the scheduler and wait for completion.

    .. code-block:: python

        flow.make_scheduler().start()

To analyze the results.

    .. code-block:: python

        lesson.analyze(flow)

Exercises
---------

Try to run the convergence study for Al.

Get a copy of the python script used in this lesson like before and look at the `analyze` method.
Use the code in `analyze` to build your Pandas dataframe and use its method to produce convergence plots:

Next
----

A logical next lesson would be lesson_relaxation

"""

_commandline_lesson_ = """
The full description, directly from the abinit documentation, is available via the shell command:

    .. code-block:: shell

        abidoc.py man inputvariable

that prints the official description of inputvariable.

The course of this lesson
-------------------------

As in the previous lesson, executing the python script creates the folder structure with the required input files.

One of the standard thing to look for to be converged in the total energy. We did that already in the previous lesson.
This time have a look at some of the other important properties. Look for instance at the convergence rate of the
forces, stress-tensor or the energies of the KS-orbitals.

Exercises
---------

Edit the input files to run the same convergence study for a different k-point mesh. Best to start small.
"""

import os
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk
from abipy.lessons.core import BaseLesson


def make_ecut_flow(structure_file=None, ecut_list = (10, 12, 14, 16, 18)):
    """
    Build and return a `Flow` to perform a convergence study wrt to ecut.

    Args:
        structure_file: (optional) file containing the crystalline structure.
            If None, crystalline silicon structure.
        ecut_list: List of cutoff energies to be investigated.
    """
    # Define the structure and add the necessary pseudos:
    if structure_file is None:
        multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                  pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(ecut_list))
        workdir = "flow_Si_ecut_convergence"
    else:
        structure = abilab.Structure.from_file(structure_file)
        pseudos = abilab.PseudoTable()  ## todo fix this
        multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=len(ecut_list))
        workdir = "flow_" + structure.composition.reduced_formula + "_ecut_convergence"

    # Add mnemonics to the input files.
    multi.set_mnemonics(True)

    # Global variables
    multi.set_vars(tolvrs=1e-9)
    multi.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0])

    # Here we set the value of ecut used by the i-th task.
    for i, ecut in enumerate(ecut_list):
        multi[i].set_vars(ecut=ecut)

    return flowtk.Flow.from_inputs(workdir=workdir, inputs=multi.split_datasets())


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
    def make_ecut_flow(**kwargs):
        return make_ecut_flow(**kwargs)

    @staticmethod
    def analyze(flow, **kwargs):
        with abilab.abirobot(flow, "GSR") as robot:
            data = robot.get_dataframe()

        import matplotlib.pyplot as plt
        ax = data.plot(x="ecut", y="energy", title="Total energy vs ecut", legend=False, style="b-o")
        ax.set_xlabel('Ecut [Ha]')
        ax.set_ylabel('Total Energy [eV]')
        return plt.show(**kwargs)


if __name__ == "__main__":
    lesson = Lesson()
    flow = lesson.make_ecut_flow()
    flow.build_and_pickle_dump()
    lesson.setup()

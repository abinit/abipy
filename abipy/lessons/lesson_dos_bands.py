#!/usr/bin/env python
"""
The calculation of the density of states and the bandstructure.
===============================================================

Background
----------

This lesson focuses on the calculation of the density of states (DOS) and 
the electronic band structure within the Kohn-Sham (KS) formalism. 

In contrast to the total energy and its derivatives, the energies of the KS-levels have no exact physical meaning,
except for the highest occupied state that actually would be the first ionization energy if the functional would be
exact. So why would we even want to calculate the KS-DOS and band structure? In most cases the KS spectrum is
qualitatively in agreement with the spectrum of ionization energies. Moreover in general we are able to make good
predictions on trends.

In lesson_g0w0.py, we discuss a more elaborated and accurate approach for the calculation of band energies and band gaps.

The related abinit variables
----------------------------

    * kptopt (negative values)
    * kptbounds (if you want to specify the bounds of the k-path)
    * ndivsm

"""
from __future__ import division, print_function


_ipython_lesson_ = """
More info on the inputvariables and their use can be obtained using the following function:

    .. code-block:: python

        lesson.docvar("inputvariable")

This will print the official abinit description of this inputvariable.

The abipy flows in this lesson
------------------------------

The flow that we use in this lesson contains for the first time dependencies.
This means that some tasks in the flow can only be started if an other task is
ready. We will first perform one self-consistent calculation to obtain a proper
density. Using this density we calculate in two more steps the DOS and the bandstructure.
For the DOS this not strictly necessary since the DOS will also be calculated on a regular grid.
In general the density will be converged already before the DOS is converged. For large systems it may become
nessesary to split. For the bandstructure, we have a non-uniform grid so we do need to fix the density.

The course of this lesson
-------------------------

Start this lesson by importing it in a new namespace:

    .. code-block:: python

        from abipy.lessons.lesson_dos_bands import Lesson
        lesson = Lesson()

As always you can reread this lessons text using the command:

    .. code-block:: python

        lesson

To build the flow:

    .. code-block:: python

        flow = lesson.make_flow()

To print the input files:

    .. code-block:: python

        flow.show_inputs()

To visualize the dependencies in the flow:

    .. code-block:: python

        flow.show_dependencies()

Start the flow with the scheduler and wait for completion.

    .. code-block:: python

        flow.make_scheduler().start()

To analyze the results.

    .. code-block:: python

        lesson.analyze(flow)

Exercises
---------

At this point, you may want to interact more with the underlying python objects
so that you can start to develop your own script or your post-processing tool.

Our flow consists of a BandStructureWork object that provides many tools for post-processing.
Use

    .. code-block:: python

            work = flow[0]

to have access to the band structure work and look at the `plot` methods that 
are available (hint: type work.plot in ipython and press TAB to get a list of methods)

1) Use the `plot_` methods to visualize the convergence of the DOS wrt to the number of k-points.
   Then change the value of the gaussian broadening (`width` parameter).

2) Plot bands and DOS on the same figure.

Rememeber that, in ipython, one can access the documentation of a method with `work.plot_edoses?`

Next
----

A logical next lesson would be lesson_g0w0
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
"""

import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.lessons.core import BaseLesson, get_pseudos


def make_electronic_structure_flow(ngkpts_for_dos=[(2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8)]):
    """Band structure calculation."""
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=2 + len(ngkpts_for_dos))
    # Global variables
    multi.set_vars(ecut=10)

    # Dataset 1
    multi[0].set_vars(tolvrs=1e-9)
    multi[0].set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0])

    # Dataset 2
    multi[1].set_vars(tolwfr=1e-15)
    multi[1].set_kpath(ndivsm=5)

    # Dataset 3
    for i, ngkpt in enumerate(ngkpts_for_dos):
        multi[2+i].set_vars(tolwfr=1e-15)
        multi[2+i].set_kmesh(ngkpt=ngkpt, shiftk=[0,0,0])

    inputs = multi.split_datasets()
    scf_input, nscf_input, dos_input = inputs[0], inputs[1], inputs[2:]

    return abilab.bandstructure_flow(workdir="flow_dos_bands", scf_input=scf_input, nscf_input=nscf_input,
                                     dos_inputs=dos_input)


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
    def make_flow(**kwargs):
        return make_electronic_structure_flow(**kwargs)

    @staticmethod
    def analyze(flow):
        nscf_task = flow[0][1]
        with nscf_task.open_gsr() as gsr:
            return gsr.ebands.plot()

if __name__ == "__main__":
    l = Lesson()
    flow = l.make_flow()
    flow.build_and_pickle_dump()
    l.setup()

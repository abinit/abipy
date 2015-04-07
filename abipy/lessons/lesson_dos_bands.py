#!/usr/bin/env python
"""
The calculation of the density of states and the bandstructure
==============================================================

Background
----------

This lesson focuses on the calculation of the density of states (DOS) and 
the electronic band structure within the Kohn-Sham (KS) formalism. 

In contrast to the total energy and its derivatives, the energies of the KS-levels have no exact physical meaning,
except for the highest occupied state that actually would be the first ionization energy if the DFT XC functional would be
exact. So why do we use the KS formalism to calculate electron DOSes and band structures? 

As a matter of fact, the KS energy spectrum is usually in qualitative agreement with experiments (let's ignore correlated systems).
Standard KS band structures with LDA or GGA are relatively cheap and KS calculations allow us to make good predictions on trends.
In lesson_g0w0.py, we discuss a more accurate and expensive approach for the calculation of band strtuctures and band gaps.

The related abinit variables
----------------------------

    * kptopt    (negative values if band structures are wanted.)
    * kptbounds (the boundaries of the k-path)
    * ndivsm    (number of points used to sample the smallest segment of the k-path)

"""
from __future__ import division, print_function


_ipython_lesson_ = """
More info on the inputvariables and their use can be obtained with:

    .. code-block:: python

        lesson.docvar("inputvariable")

This command will print the official description of `inputvariable`.

The abipy flow used in this lesson
----------------------------------

The flow used  in this lesson contains for the first time dependencies.
This means that some of the tasks in the flow can start only if its `parents` are completed.
We will first perform one self-consistent calculation to obtain a proper density. 
From this density we then calculate the DOS and the bandstructure in two independent tasks. 
The DOS calculation is done on a regular grid of k-points because the DOS is expressed in 
terms on a integral over the first Brillouin zone.
For the band structure task, we use a high symmetry path inside the BZ.

The course of this lesson
-------------------------

Start this lesson by importing it in a new namespace:

    .. code-block:: python

        from abipy.lessons.lesson_dos_bands import Lesson
        lesson = Lesson()

As usual, you can reread this text using the command:

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
so that you can start to develop your script or your post-processing tools.

Our flow consists of a `BandStructureWork` object that provides many tools for the post-processing.
Use

    .. code-block:: python

            work = flow[0]

to access the band structure work and look at the `plot` methods that 
are available (hint: type work.plot in ipython and press TAB to get a list of methods)

1) Use the `plot_` methods to visualize the convergence of the DOS wrt to the number of k-points.
   Then change the value of the gaussian broadening (`width` parameter).

2) Plot bands and DOS on the same figure.

Remember that, in ipython, one can access the documentation with `work.plot_edoses?`

Next
----

A logical next lesson would be lesson_g0w0
"""

_commandline_lesson_ = """
The full description, directly from the abinit documentation, is available via the following function:

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

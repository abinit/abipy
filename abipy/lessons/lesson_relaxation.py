#!/usr/bin/env python
"""
Relaxation of the unit cell with two different techniques
=========================================================

Background
----------

One of the tasks that is most performed using DFT is the relaxation of an
atomic structure. Effectively we search for that structure for which the
total energy is minimal. Since the total energy is in principal exact in
DFT the atomic position are in general rather good. 'In principal' means
if the exchange-correlation functional would be exact. However, since we
are comparing differences in total energies an certain amount of
error-cancellation can be expected.

In this lesson we focus on different types of structure relaxations.

The related abinit variables
----------------------------

    * ionmov
    * optcell
    * dilatmx
    * ecutsm
    * ntime

"""
from __future__ import division, print_function

_ipython_lesson_ = """

More info on the inputvariables and their use can be obtained using
the following function:

    .. code-block:: python

        lesson.docvar("inputvariable")


The abipy flows in this lesson
------------------------------

In this lesson we will use two different relaxation flows. One flow will
calculate the total energies of a compound at various volumes and fit an
equation of state to the energy v.s. volume data. Besides the optimal
volume, where the energy is minimal, this will also provide the bulk modulus,
the 'compressebility' of the systems. The other flow will automatically
optimize all degrees of freedom. In our first example Si, there is only
one degree of freedom, due to the symmetry of the crystal, the volume of
the unit cell, or the lattice parameter. In the second example, GaN, the
symmetry is lower and one internal degree of freedom appears.


The course of this lesson
-------------------------

Start this lesson by importing it in a new namespace:

    .. code-block:: python

        from abipy.lessons.lesson_relaxation import Lesson()
        lesson = Lesson()

As always you can reread this lessons text using the command:

    .. code-block:: python

        lesson

To build the flow for silicon

    .. code-block:: python

        flow = lesson.make_relax_eos_flow()

For Gallium Arsenide, use

    .. code-block:: python

        flow = lesson.make_relax_relax_flow()

To print the input files

    .. code-block:: python

        flow.show_inputs()

Start the flow with the scheduler and wait for completion.

    .. code-block:: python

        flow.make_scheduler().start()

To analyze the results.

    .. code-block:: python

        lesson.analyze(flow)

In the case of silicon, it will show a fit of the total energy vs the
volume of the unit cell. The minimum of this curve is the equilibrium
volume. From this fit, we can also obtain the bulk modulus.

Volume of the unit cell of silicon : XXX A^3 [ source ?]
Bulk modulus : 98.8 GPa [ source ? ]

In the case of gallium arsenide, you will see the change of equilibrium
volume and length of the box with respect to the k-point mesh.

Volume of the unit cell of GaN : XXX A^3 [ source ?]
Vertical distance between Ga and N : XXX A [ source ?]

Of course you will need to converge your results with respect to
the kpoint sampling and with respect with ecut...

Exercises
---------

As an exercise you can now try to get the equilibrium unit cell
of silicon automatically using abinit. You can inspire yourself
from the GaN relaxation. First download a local copy of the python
script.

    .. code-block:: python

        lesson.get_local_copy()

And have a look in make_relax_gan_flow(), try to do the same
with 'si.cif' file instead of 'gan.cif'

As a second exercice, you can try to converge the results obtained
here with respect to the k-point sampling and with respect to ecut
and compare the converged results with experimental data.

Next
----

A logical next lesson would be lesson_dos_bands
"""

_commandline_lesson_ = """
At this place they will not be discussed in detail. In stead you are
invited to read the abinit documentation on them. The full description,
directly from the abinit description is available via the following function:

    .. code-block:: shell

        abidocs.py man inputvariable

This will print the official abinit description of this inputvariable.


The course of this lesson
-------------------------
"""



import os
import numpy as np
import abipy.abilab as abilab
import abipy.data as abidata

from pymatgen.io.abinitio.eos import EOS
from abipy.core import Structure
from abipy.lessons.core import BaseLesson, get_pseudos


def get_dist(gsrfile):
    struct = gsrfile.structure
    red_dist = struct.frac_coords[2][2] - struct.frac_coords[0][2]
    return 'u',red_dist


class RelaxFlow(abilab.Flow):

    def analyze(self):
        with abilab.GsrRobot.open(self) as robot:
            data = robot.get_dataframe(funcs=get_dist)
            robot.pairplot(data, x_vars="nkpts", y_vars=["a", "c", "u"]) #, hue="tsmear")

            #grid = sns.PairGrid(data, x_vars="nkpts", y_vars=["a", "volume"]) #, hue="tsmear")
            #grid.map(plt.plot, marker="o")
            #grid.add_legend()
            #plt.show()


class EosFlow(abilab.Flow):
    def analyze(self):
        work = self.works[0]
        etotals = work.read_etotals(unit="eV")

        #eos_fit = EOS.DeltaFactor().fit(self.volumes, etotals)
        #eos_fit.plot()

        eos_fit = EOS.Birch_Murnaghan().fit(self.volumes, etotals)
        return eos_fit.plot()


def make_relax_flow(structure_file=None):
    # Structural relaxation for different k-point samplings.
    ngkpt_list = [[3, 3, 2], [6, 6, 4], [8, 8, 6]]

    if structure_file is None:
        structure = abilab.Structure.from_file(abidata.cif_file("gan2.cif"))
    else:
        structure = abilab.Structure.from_file(structure_file)

    inp = abilab.AbiInput(pseudos=get_pseudos(structure), ndtset=len(ngkpt_list))
    inp.set_structure(structure)

    # Add mnemonics to input file.
    inp.set_mnemonics(True)

    # Global variables
    inp.set_vars(
        ecut=20,
        tolrff=1.0e-2,
        nstep=100,
        optcell=2,
        ionmov=3,
        ntime=500,
        dilatmx=1.05,
        ecutsm=0.5,
    )

    for i, ngkpt in enumerate(ngkpt_list):
        inp[i+1].set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])

    return RelaxFlow.from_inputs("flow_gan_relax", inputs=inp.split_datasets(), task_class=abilab.RelaxTask)


def make_eos_flow(structure_file=None):
    # Structural relaxation for different k-point samplings.
    scale_volumes = np.arange(94, 108, 2) / 100.

    if structure_file is None:
        structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    else:
        structure = abilab.Structure.from_file(structure_file)

    inp = abilab.AbiInput(pseudos=get_pseudos(structure), ndtset=len(scale_volumes))
    inp.set_structure(structure)

    # Global variables
    inp.set_vars(
        ecut=16,
        tolvrs=1e-16
    )

    inp.set_kmesh(ngkpt=[4, 4, 4], shiftk=[[0.5, 0.5, 0.5], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]])

    for idt, scal_vol in enumerate(scale_volumes):
        new_lattice = structure.lattice.scale(structure.volume*scal_vol)

        new_structure = Structure(new_lattice, structure.species, structure.frac_coords)

        inp[idt+1].set_structure(new_structure)

    eos_flow = EosFlow.from_inputs("flow_si_relax", inputs=inp.split_datasets(), task_class=abilab.RelaxTask)
    eos_flow.volumes = structure.volume*scale_volumes
    return eos_flow


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
    def make_eos_flow(**kwargs):
        return make_eos_flow(**kwargs)

    @staticmethod
    def make_relax_flow(**kwargs):
        return make_relax_flow(**kwargs)

if __name__ == "__main__":
    l = Lesson()
    l = Lesson()
    flow = l.make_eos_flow()
    flow.build_and_pickle_dump()
    flow = l.make_relax_flow()
    flow.build_and_pickle_dump()
    l.manfile(l.comline_string)
    l.instruct()

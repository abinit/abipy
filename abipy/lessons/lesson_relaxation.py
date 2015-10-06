#!/usr/bin/env python
"""
Relaxation of the unit cell with two different techniques
=========================================================

Background
----------

In this lesson we discuss two different methods to find the equilibrium  structure of a system.
In the first method, we use the GS part of Abinit to calculate the total energy of silicon for 
different volumes and then we fit the energy vs the volume with a model for the equation of state (EOS).
The fit provides the optimal volume (i.e. the volume for which the total energy is minimal),
as well as the bulk modulus (the 'compressibility' of the system).
Note that this approach is only applicable to isotropic materials without any degree of freedom for the atomic positions. 
Indeed, the equation of state is obtained by performing a homogeneous compressions/dilatation of 
the initial Bravais lattice while keeping the atoms fixed in the initial high-symmetry positions.

In the second example, we find the equilibrium configuration of GaN.
In this case, the approach used for computing the EOS of silicon is not applicable because, 
one should optimize both the lattice parameters as well the distance between Ga and N.
For this reason, we employ the relaxation algorithms implemented in Abinit (`ionmov` and `optcell`)
in which the forces and the stresses obtained at the end of the SCF cycle are used to find the minimum energy configuration.

The related abinit variables
----------------------------

    * ionmov
    * optcell
    * dilatmx
    * ecutsm
    * ntime
    * tolmxf
    * tolrff

"""
from __future__ import division, print_function

_ipython_lesson_ = """
For a more detailed description of the variables, you are invited to consult the abinit documentation. 
The full description, directly from the official abinit docs, is available in ipython with the command:

    .. code-block:: python

        print(lesson.docvar("inputvariable"))


Description of the lesson
-------------------------

We will use two different AbiPy flows to find the equilibrium configuration.
The first flow, si_flow, calculates the total energy of silicon at different volumes 
and computes the equation of state E(V).
The other flow, gan_flow, uses Abinit to optimize all degrees of freedom (atomic positions
and lattice vectors).

Executing the lesson
--------------------

Start this lesson by importing it with the commands:

    .. code-block:: python

        from abipy.lessons.lesson_relaxation import Lesson
        lesson = Lesson()

As usual, you can reread this text using the command:

    .. code-block:: python

        lesson

To build the flow for silicon, use

    .. code-block:: python

        si_flow = lesson.make_eos_flow()

For Gallium Nitride, use

    .. code-block:: python

        gan_flow = lesson.make_relax_flow()

To print the input files

    .. code-block:: python

        si_flow.show_inputs()

Start the flow with the scheduler and wait for completion.

    .. code-block:: python

        si_flow.make_scheduler().start()

To analyze the results.

    .. code-block:: python

        # For Silicon
        lesson.analyze_eos_flow(si_flow)

        # For Gallium Nitride, use
        lesson.analyze_eos_flow(gan_flow)

In the case of silicon, python will show a fit of the total energy vs the
volume of the unit cell. The minimum of this curve is the equilibrium
volume. From this fit, we can also obtain the bulk modulus.
Note that this approach is only applicable to isotropic materials since the
equation of state has been obtained by performing 
a homogeneous compressions/dilatation of the initial Bravais lattice.

Try to compare the results with these experimental results:
Volume of the unit cell of silicon: 40.05 A^3 [NSM]
Bulk modulus: 98 GPa [NSM]

In the case of gallium nitride, you will observe a change of the equilibrium
parameters with respect to the k-point mesh.

Try to compare the results with these experimental results:

    * Volume of the unit cell of GaN: 45.73 A^3 [Schulz & Thiemann 1977]
    * Lattice parameters of GaN: a = 3.190 A, c = 5.189 A [Schulz & Thiemann 1977]
    * Vertical distance between Ga and N : about 0.377 * c [ Schulz & Thiemann, 1977]

Of course you will need to converge your results with respect to the k-point sampling and the 
cutoff energy ecut.

Note the we are using pseudopotentials generated with the GGA which tends to 
overestimate the lattice parameters. 
If you use LDA-type pseudopotentials, you will observe that LDA tends to underestimate the parameters.

Exercises
---------

As an exercise you can now try to get the equilibrium unit cell of silicon automatically using abinit. 
You can use the code for the relaxation of GaN as template.
First download a local copy of the python script.

    .. code-block:: python

        lesson.get_local_copy()

and have a look at the code in make_relax_gan_flow().
Try to do the same with 'si.cif' file instead of 'gan.cif'

Pay attention to the fact that for silicon, you cannot use tolrff
to stop your self-consistent cycle. 
Silicon has no internal degree of freedom, the forces are zero by symmetry. 
and hence the tolrff criterion makes no sense.

As a second exercise, you can try to converge the results for silicon with respect 
to the k-point sampling and ecut.
Compare the converged results with experimental data.

Next
----

A logical next lesson would be lesson_dos_bands
"""


_commandline_lesson_ = """
The full description of the variables, directly from the abinit documentation is available via the shell command:

    .. code-block:: shell

        abidoc.py man inputvariable

that prints the official description of `inputvariable`.

As in the previous lessons, executing the python script creates the folder structure with the input files for this
lesson.

For the flow_si_relax folder, look in particular to the changes in the unit cell (rprim) in the input files and the
corresponding change in the unit cell volume (ucvol), total energy (etotal) and stresses (strten) in the output file.
For the flow_gan_relax, observe in the output files how the automatic relaxation takes place.
At each step of the relaxation, a full SCF-cycle is performed and forces and the stress are computed. 
The ions are then moved according to the forces and a new SCF-cycle is started.
The procedure is interated until convergence is achieved. 
This is the reason why there are two stopping criteria for structural relaxation:
tolrff or tolvrs are used for the SCF cycle whereas tolmxf govers the relaxation algorithm.

Exercises
---------

Edit the input files to run the same jobs with different values of ecut. 

You can also try to change the stopping criterion to see if this affects the final results.

Finally, try to generate the input file for silicon, and try to guess why setting the stopping 
criterion on the forces won't work in this case!
"""

import os
import numpy as np
import abipy.abilab as abilab
import abipy.data as abidata

from pymatgen.io.abinit.eos import EOS
from abipy.core import Structure
from abipy.lessons.core import BaseLesson, get_pseudos


def make_relax_flow(structure_file=None):
    """
    Build and return a flow that perform a structural relaxation for different k-point samplings.
    """
    ngkpt_list = [[3, 3, 2], [6, 6, 4], [8, 8, 6]]

    if structure_file is None:
        structure = abilab.Structure.from_file(abidata.cif_file("gan2.cif"))
    else:
        structure = abilab.Structure.from_file(structure_file)

    multi = abilab.MultiDataset(structure=structure, pseudos=get_pseudos(structure), ndtset=len(ngkpt_list))

    # Add mnemonics to input file.
    multi.set_mnemonics(True)

    # Global variables
    multi.set_vars(
        ecut=20,
        tolrff=5.0e-2,
        nstep=30,
        optcell=2,
        ionmov=3,
        ntime=50,
        dilatmx=1.05,
        ecutsm=0.5,
        tolmxf=5.0e-5,
    )

    for i, ngkpt in enumerate(ngkpt_list):
        multi[i].set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])

    return abilab.Flow.from_inputs("flow_gan_relax", inputs=multi.split_datasets(), task_class=abilab.RelaxTask)


def make_eos_flow(structure_file=None):
    """
    Build and return a Flow to compute the equation of state 
    of an isotropic material for different k-point samplings.
    """
    scale_volumes = np.arange(94, 108, 2) / 100.

    if structure_file is None:
        structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    else:
        structure = abilab.Structure.from_file(structure_file)

    multi = abilab.MultiDataset(structure=structure, pseudos=get_pseudos(structure), ndtset=len(scale_volumes))

    # Global variables
    multi.set_vars(
        ecut=16,
        tolvrs=1e-16
    )

    multi.set_kmesh(ngkpt=[4, 4, 4], shiftk=[[0.5, 0.5, 0.5], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]])

    for idt, scal_vol in enumerate(scale_volumes):
        new_lattice = structure.lattice.scale(structure.volume*scal_vol)

        new_structure = Structure(new_lattice, structure.species, structure.frac_coords)

        multi[idt].set_structure(new_structure)

    eos_flow = abilab.Flow.from_inputs("flow_si_relax", inputs=multi.split_datasets())
    eos_flow.volumes = structure.volume * scale_volumes
    return eos_flow


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
    def make_eos_flow(**kwargs):
        return make_eos_flow(**kwargs)

    @staticmethod
    def analyze_eos_flow(flow, **kwargs):
        work = flow[0]
        etotals = work.read_etotals(unit="eV")
        eos_fit = EOS.Birch_Murnaghan().fit(flow.volumes, etotals)
        return eos_fit.plot(**kwargs)

    @staticmethod
    def make_relax_flow(**kwargs):
        return make_relax_flow(**kwargs)

    @staticmethod
    def analyze_relax_flow(flow, **kwargs):
        def get_dist(gsrfile):
            struct = gsrfile.structure
            red_dist = struct.frac_coords[2][2] - struct.frac_coords[0][2]
            return 'u', red_dist

        with abilab.GsrRobot.open(flow) as robot:
            data = robot.get_dataframe(funcs=get_dist)
            return robot.pairplot(data, x_vars="nkpts", y_vars=["a", "c", "volume", "u"]) 


if __name__ == "__main__":
    l = Lesson()
    flow = l.make_eos_flow()
    flow.build_and_pickle_dump()
    flow = l.make_relax_flow()
    flow.build_and_pickle_dump()
    l.setup()

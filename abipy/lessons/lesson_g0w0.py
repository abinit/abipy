#!/usr/bin/env python
"""
$G_0W_0$ band structure with an energy-dependent scissors operator
==================================================================

Background
----------

Standard functionals (LDA and GGA), systematically underestimate band gaps, giving values 
that are about 30-40% smaller than experimental data.
The inability of standard Kohn-Sham (KS) theory to give band gaps close to experiment is often referred to as the **band-gap problem**. 
From a theoretical point of view this is not surprising since KS eigenvalues are not supposed to give the correct band energies.
The band structure of a crystal is rigorously defined as the energies needed to add or subtract electrons from the many-body system
which, in turn, are related to the difference between total energies of many-body states differing by one electron. 

An alternative, more traditional, approach to the study of exchange-correlation effects in 
many-body systems is provided by Many-Body Perturbation Theory (MBPT) which defines a rigorous approach to the description of excited-state properties, 
based on the Green's function formalism.
In this lesson, we discuss how to use the MBPT part of ABINIT to compute the band-structure of silicon within the so-called $G_0W_0$ approximation.

For a very brief introduction see MBPT_NOTES_. 

.. _MBPT_NOTES: http://www.abinit.org/documentation/helpfiles/for-v7.10/tutorial/theory_mbt.html 

Related ABINIT variables
------------------------

    * optdriver
    * ecuteps
    * ecutsigx
    * nband
    * gwcalctyp
    * gw_qprange
    * all gw** variables

"""
from __future__ import division, print_function, unicode_literals

_ipython_lesson_ = """
More info on the input variables and their use can be obtained
using the following function:

    .. code-block:: python

        lesson.abinit_help(inputvariable)

This will print the official abinit description of this variables.

To open the python script in ipython use:

    .. code-block:: python

        %load $lesson.pyfile

The abipy flows of this lesson
------------------------------

In this lesson, we will construct an `abipy` flow made of two works.
The first work is a standard KS band-structure calculation that consists of 
an initial GS calculation to get the density followed by two NSCF calculations.
The first NSCF task computes the KS eigenvalues on a high-symmetry path in the BZ,
whereas the second NSCF task is done on a homogeneous k-mesh so that one can calculate 
the DOS from the KS eigenvalues. 

The second work represents the real GW workflow in which we read the density computed in the first task of 
the previous work to compute the KS bands for many empty states. 
The WFK file produced in this step is then used to compute the screened interaction $W$. 
Finally we do a self-energy calculation in which we use the $W$ produced
in the previous step and the WFK file to compute the matrix elements of the self-energy and 
the $G_0W_0$ corrections for all the k-points in the IBZ and 8 bands (4 occupied + 4 empty)

Once the flow is completed, we can interpolate the $G_0W_0$ corrections as function of the initial KS energy 
to obtain an energy-dependent scissors operator. 
At this point, we can apply the scissors operator onto the KS band structure to obtain an approximated $G_0W_0$
band dispersion.

The course of this lesson
-------------------------

This lesson can be started in ipython by importing it:

    .. code-block:: python

        from abipy.lessons.lesson_g0w0 import Lesson()
        lesson = Lesson()

The lesson is now imported in your ipython session in its own namespace 'lesson'. 
This object now gives us all the tools to follow this lesson. As before:

    .. code-block:: python

        lesson

displays this lessons information text. This lesson provides a
factory function that returns a flow designed to perform a standard
G0W0 calculation.

In the previous lesson we have actually been running job directly on
the front-end. These calculations were so small that this was not a
problem. GW calculations, however, (even the under converged examples
we are using here) are much more involved. To run submit calculations
to the actual nodes of the cluster we only need to provide abipy
with different manager settings. First have a look at the current
manager.yml file. This one tells abipy what it needs to know to run
shell jobs. Next copy the file we prepared for this cluster:

cp /data/euspec/doc/abinit-templates/manager_viper.yml.

Have a look at this file as well. It may look complicated but if fact it
is just a translation of the user manual of the cluster. For a new cluster
one person has to create it once. Also note the it only mentions which queueing
systems is installed how to use this systems is programmed in abipy. To
use this manager move it to manager.yml. (abipy will first look for a manager
file in you current folder and secondly in ~/.abinit/abipy, so you can put
one there an don't bother about it for every calculation)

displays this lessons information text, and can be recalled at any moment.
The main object we use to pack (connected series of) calculations is a flow. 
This lesson provides a method that returns a flow designed to perform k-point convergence studies. 

This flow is made by the command:

    .. code-block:: python

        flow = lesson.make_g0w0_scissors_flow()

`flow` is now an object that contains al the information needed to generate abinit inputs. 

    .. code-block:: python

        flow.show_inputs()

will display all the inputs as they will be 'given' to abinit. In
previous lessons we ran the flows each time directly inside ipython.
For relatively small calculations this is very practical. There are
however other ways more suited for large calculations.

To start the execution of calculations packed in this flow we use the following command:

    .. code-block:: python

        flow.make_scheduler().start()

This starts the actual execution via a scheduler. 

The last step of analyzing the results can be done again in with a single command:

    .. code-block:: python

        lesson.analyze(flow)

This method of flow will open the necessary output files, retrieve the data, and produce a plot.

Finally, once you are through with this lesson and exited ipython:

    .. code-block:: python

        exit

You can see that in the directory that you were working there is
now a subdir were the calculation have been performed. 
Have a look at these folders and the files that are in them.

Exercises
---------

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
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.lessons.core import BaseLesson

abinit_help = abilab.abinit_help


def make_inputs(ngkpt, paral_kgb=0):
    # Crystalline silicon
    # Calculation of the GW band structure with the scissors operator.
    # Dataset 1: ground state run to get the density.
    # Dataset 2: NSCF run to get the KS band structure on a high-symmetry k-path.
    # Dataset 3: NSCF run with a homogeneous sampling of the BZ to compute the KS DOS.
    # Dataset 4: NSCF run with empty states to prepare the GW steps.
    # Dataset 5: calculation of the screening from the WFK file computed in dataset 4.
    # Dataset 6: Use the SCR file computed at step 5 and the WFK file computed in dataset 4 to get the GW corrections.

    inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=6)
    inp.set_structure(abidata.cif_file("si.cif"))

    # Add mnemonics to input file.
    inp.set_mnemonics(True)

    # This grid is the most economical, but does not contain the Gamma point.
    scf_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    dos_kmesh = dict(
        ngkpt=(6, 6, 6),
        shiftk=[0.0, 0.0, 0.0])

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap. 
    gw_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.0, 0.0, 0.0,  
                0.0, 0.5, 0.5,  
                0.5, 0.0, 0.5,  
                0.5, 0.5, 0.0]
    )
       
    # Global variables
    ecut = 6
    inp.set_vars(
        ecut=ecut,
        istwfk="*1",
        paral_kgb=paral_kgb,
        gwpara=2,
    )

    # Dataset 1 (GS run to get the density)
    inp[1].set_kmesh(**scf_kmesh)
    inp[1].set_vars(
        tolvrs=1e-6,
        nband=4,
    )
    inp[1].set_kmesh(**scf_kmesh)

    # Dataset 2 (NSCF run)
    inp[2].set_vars(iscf=-2,
                    tolwfr=1e-12,
                    nband=8,
                   )
    inp[2].set_kpath(ndivsm=8)

    # Dataset 3 (DOS NSCF)
    inp[3].set_vars(iscf=-2,
                    tolwfr=1e-12,
                    nband=35,
                    #nband=10,
                   )
    inp[3].set_kmesh(**dos_kmesh)

    # Dataset 4 (NSCF run for GW)
    inp[4].set_vars(iscf=-2,
                    tolwfr=1e-12,
                    nband=35,
                   )
    inp[4].set_kmesh(**gw_kmesh)

    # Dataset3: Calculation of the screening.
    inp[5].set_vars(
        optdriver=3,   
        nband=25,    
        ecutwfn=ecut,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
    )
    inp[5].set_kmesh(**gw_kmesh)

    inp[6].set_vars(
            optdriver=4,
            nband=10,      
            ecutwfn=ecut,
            ecuteps=4.0,
            ecutsigx=6.0,
            symsigma=1,
            gw_qprange=-4,  # Compute GW corrections for all kpts in IBZ, all occupied states and 4 empty states,
        )
    inp[6].set_kmesh(**gw_kmesh)

    return inp.split_datasets()


def make_g0w0_scissors_flow(workdir="flow_lesson_g0w0"):
    # Change the value of ngkpt below to perform a GW calculation with a different k-mesh.
    scf, bands_nscf, dos_nscf, gw_nscf, scr, sig = make_inputs(ngkpt=[2,2,2])

    flow = abilab.Flow(workdir=workdir)
    work0 = abilab.BandStructureWork(scf, bands_nscf, dos_inputs=dos_nscf)
    flow.register_work(work0)

    work1 = abilab.Work()
    gw_nscf_task = work1.register_nscf_task(gw_nscf, deps={work0[0]: "DEN"})
    scr_task = work1.register_scr_task(scr, deps={gw_nscf_task: "WFK"})
    sigma_task = work1.register_sigma_task(sig, deps={gw_nscf_task: "WFK", scr_task: "SCR"})
    flow.register_work(work1)

    return flow.allocate()

class Lesson(BaseLesson):

    @property
    def abipy_string(self):
        return __doc__ + _ipython_lesson_

    @property
    def comline_string(self):
        return __doc__ + _commandline_lesson_

    @property
    def pyfile(self):
        return __file__

    @staticmethod
    def make_flow(**kwargs):
        return make_g0w0_scissors_flow(**kwargs)

    @staticmethod
    def analyze_flow(flow, domains_spin=[[-10, 6.02], [6.1, 20]]):
        # Read the G0W0 correction form the output file of the sigma_task
        # and construct the scissors_builder object.
        sigma_task = flow[1][2]
        builder = sigma_task.get_scissors_builder()

        # Plot G0W0 results as function of the initial KS energy.
        builder.plot_qpe_vs_e0()

        # Build the scissors operator with a polyfit in the regions specified by domains_spin
        builder.build(domains_spin=domains_spin)

        # Plot the fit.
        builder.plot_fit()

        # Here we pass the KS bands to the scissors-builder.
        # plot_qpbands with apply the scissors onto the input band structure and plot the final results. 
        bands_task = flow[0][1]
        bands_filepath = bands_task.gsr_path
        builder.plot_qpbands(bands_filepath, bands_label="KS Bands", title="Silicon Bands (KS and KS+scissors)")

        # TODO: Fix problems with boundaries!
        #dos_task = flow[0][2]
        #dos_filepath = dos_task.outdir.has_abiext("GSR")
        #builder.plot_qpbands(bands_filepath, dos_filepath=dos_filepath,
        #                     title="Silicon Bands and DOS (KS and KS+scissors)")


if __name__ == "__main__":
    l = Lesson()
    flow = l.make_flow()
    flow.build_and_pickle_dump()
    l.manfile(l.comline_string)
    l.instruct()

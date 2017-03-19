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
many-body systems is provided by Many-Body Perturbation Theory (MBPT) which defines a rigorous approach to the description
of excited-state properties, based on the Green's function formalism.
In this lesson, we discuss how to use the MBPT part of ABINIT to compute the band-structure of silicon
within the so-called $G_0W_0$ approximation.

For a very brief introduction to the many-body formalism, see MBPT_NOTES_.

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
More info on the input variables and it's usage can be obtained using the command:

    .. code-block:: python

        lesson.abinit_help(inputvariable)

that prints the official description of the variables.

Description of the lesson
-------------------------

In this lesson, we will construct an `abipy` flow made of two works.
The first work is a standard KS band-structure calculation that consists of
an initial GS calculation to get the density followed by two NSCF calculations.
The first NSCF task computes the KS eigenvalues on a high-symmetry path in the BZ,
whereas the second NSCF task employs a homogeneous k-mesh so that one can compute
the DOS from the KS eigenvalues.
This work is similar to the one we have already encountered in lesson_dos_bands.

The second work represents the real GW workflow that uses the density computed in the first task of
the previous work  to compute the KS bands for many empty states.
The WFK file produced in this step is then used to compute the screened interaction $W$.
Finally, we perform a self-energy calculation that uses the $W$ produced
in the previous step and the WFK file to compute the matrix elements of the self-energy and
the $G_0W_0$ corrections for all the k-points in the IBZ and 8 bands (4 occupied + 4 empty)

Once the flow is completed, we can interpolate the $G_0W_0$ corrections as function of the initial KS energy
to obtain an energy-dependent scissors operator.
At this point, we can apply the scissors operator onto the KS band structure to obtain
an approximated $G_0W_0$ band dispersion.

Don't worry if there are steps of the entire procedure that are not clear to you.
GW calculations are much more complicated than standard KS band structures and
the main goal of this lesson is to give you an overview of the Abipy capabilities.

Executing the lesson
--------------------

This lesson can be started in ipython by importing it with:

    .. code-block:: python

        from abipy.lessons.lesson_g0w0 import Lesson
        lesson = Lesson()

The `lesson` object gives us all the tools needed to execute this tutorial. As usual:

    .. code-block:: python

        lesson

displays this text.

The lesson module provides a factory function that returns a flow designed to perform standard G0W0 calculations.
To build the flow, use

    .. code-block:: python

        flow = lesson.make_flow()

`flow` is the object containing al the information needed to generate abinit inputs.

    .. code-block:: python

        flow.show_inputs()

displays all the inputs that will be 'passed' to abinit.

To start the execution of calculations packed in this flow we use the following command:

    .. code-block:: python

        flow.make_scheduler().start()

This starts the actual execution via a scheduler.

The last step of analyzing the results can be done again in with a single command:

    .. code-block:: python

        lesson.analyze(flow)

This method of flow will open the necessary output files, retrieve the data, and produce a plot.

Finally, once you have completed this lesson you can exit ipython with:

    .. code-block:: python

        exit

You can see that in the working directory, here is
now a subdir were the calculation have been performed.
Have a look at these folders and the files that are in them.

Next
----

A logical next lesson would be lesson_bse.
Please consult the ipython notebook available on the abipy website
"""

_commandline_lesson_ = """
The full description, directly from the abinit documentation, is available via the following shell command:

    .. code-block:: shell

        abidoc.py man inputvariable

This command will print the official description of inputvariable.

The course of this lesson
-------------------------

"""
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowapi as flowapi

from abipy.lessons.core import BaseLesson


def make_inputs(ngkpt, paral_kgb=0):
    """
    Crystalline silicon: calculation of the G0W0 band structure with the scissors operator.

    Args:
        ngkpt: list of 3 integers. Abinit variable defining the k-point sampling.
        paral_kgb: Option used to select the eigensolver in the GS part.

    Return:
        Six AbinitInput objects:

        0: ground state run to get the density.
        1: NSCF run to get the KS band structure on a high-symmetry k-path.
        2: NSCF run with a homogeneous sampling of the BZ to compute the KS DOS.
        3: NSCF run with empty states to prepare the GW steps.
        4: calculation of the screening from the WFK file computed in dataset 4.
        5: Use the SCR file computed at step 5 and the WFK file computed in dataset 4 to get the GW corrections.
    """

    multi = abilab.MultiDataset(abidata.cif_file("si.cif"),
                                pseudos=abidata.pseudos("14si.pspnc"), ndtset=6)
    # Add mnemonics to input file.
    multi.set_mnemonics(True)

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
    multi.set_vars(
        ecut=ecut,
        istwfk="*1",
        paral_kgb=paral_kgb,
        gwpara=2,
    )

    # Dataset 1 (GS run to get the density)
    multi[0].set_kmesh(**scf_kmesh)
    multi[0].set_vars(
        tolvrs=1e-6,
        nband=4,
    )
    multi[0].set_kmesh(**scf_kmesh)

    # Dataset 2 (NSCF run)
    multi[1].set_vars(iscf=-2,
                      tolwfr=1e-12,
                      nband=8,
                      )
    multi[1].set_kpath(ndivsm=8)

    # Dataset 3 (DOS NSCF)
    multi[2].set_vars(iscf=-2,
                      tolwfr=1e-12,
                      nband=35,
                      #nband=10,
                      )
    multi[2].set_kmesh(**dos_kmesh)

    # Dataset 4 (NSCF run for GW)
    multi[3].set_vars(iscf=-2,
                      tolwfr=1e-12,
                      nband=35,
                     )
    multi[3].set_kmesh(**gw_kmesh)

    # Dataset3: Calculation of the screening.
    multi[4].set_vars(
        optdriver=3,
        nband=25,
        ecutwfn=ecut,
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,
    )
    multi[4].set_kmesh(**gw_kmesh)

    multi[5].set_vars(
            optdriver=4,
            nband=10,
            ecutwfn=ecut,
            ecuteps=4.0,
            ecutsigx=6.0,
            symsigma=1,
            gw_qprange=-4,  # Compute GW corrections for all kpts in IBZ,
                            # all occupied states and 4 empty states,
        )
    multi[5].set_kmesh(**gw_kmesh)

    return multi.split_datasets()


def make_g0w0_scissors_flow(workdir="flow_lesson_g0w0", ngkpt=(2,2,2)):
    # Change the value of ngkpt below to perform a GW calculation with a different k-mesh.
    scf, bands_nscf, dos_nscf, gw_nscf, scr, sig = make_inputs(ngkpt=ngkpt)

    flow = flowapi.Flow(workdir=workdir)
    work0 = flowapi.BandStructureWork(scf, bands_nscf, dos_inputs=dos_nscf)
    flow.register_work(work0)

    work1 = flowapi.Work()
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
        return os.path.abspath(__file__).replace(".pyc", ".py")

    @staticmethod
    def make_flow(**kwargs):
        return make_g0w0_scissors_flow(**kwargs)

    @staticmethod
    def analyze(flow, domains_spin=((-10, 6.02), (6.1, 20))):
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
    lesson = Lesson()
    flow = lesson.make_flow()
    flow.build_and_pickle_dump()
    lesson.setup()

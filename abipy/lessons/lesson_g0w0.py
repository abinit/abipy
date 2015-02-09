#!/usr/bin/env python
"""
\033[91m Computation of the G0W0 band structure of silicon with a energy-dependent scissors operator.\033[0m

\033[94m Background\033[0m

Standard approximations for XC functionals (LDA, GGA) severely underestimate band-gaps and BLABLABLA

\033[94m The related abinit variables\033[0m

\033[1m optdriver \033[0m
\033[1m ecuteps \033[0m
\033[1m ecutsigx \033[0m

More info on the input variables and their use can be obtained
using the following function:

\033[92m In []:\033[0m lesson.abinit_help(inputvariable)

This will print the official abinit description of this variables.

\033[94m The abipy flows of this lesson\033[0m

In this lesson, we will construct an abipy flow made of two works.
The first work is a standard KS band-structure calculation that
consists of an initial GS calculation to get the density followed
by two NSCF calculations. The first NSCF task computes the KS
eigenvalues of a high-symmetry path in the BZ, whereas the second
NSCF task is done on a homogeneous k-mesh so that one can calculate
the DOS from the KS eigenvalues. These two NSCF tasks 

The second work represents a typical GW workflow in which we
read the density computed in the first task of the previous work
to obtain and compute the KS eigenvalues and eigenvectors for many
empty states. The WFK file produced in this step is then used to
compute the screened interaction W. Finally we do a self-energy
calculation in which we use the W produced in the previous step
and the WFK file to compute the matrix elements of the self-energy
and the G0W0 corrections for all the k-points in the IBZ and 8
bands (4 occupied + 4 empty).

Once the flow is completed, we can post-process the results and
compute the G0W0 band-structure of silicon with a scissors operator...

\033[94m The Course of this lesson\033[0m

This lesson can be started in ipython by importing it:

\033[92m In []:\033[0m from abipy.lessons import lesson_g0w0 as lesson

The lesson is now imported in your ipython session in its own
namespace 'lesson'. This object now gives us all the tools to
follow this lesson. As before:

\033[92m In []:\033[0m lesson.help()

displays this lessons information text. This lesson provides a
factory function that returns a flow designed to perform a standard
G0W0 calculation.

In the previous lesson we have actually been running job directly on
the frontend. These calculations were so small that this was not a
problem. GW calculations, however, (even the underconverged examples
we are using here) are much more involved. To run submit calculations
to the actual worknodes of the cluster we only need to provide abipy
with different manager settings. First have a look at the current
manager.yml file. This one tells abipy what it needs to know to run
shell jobs. Next copy the file we prepared for this cluster:

\033[92m In []:\033[0m cp /data/euspec/doc/abinit-templates/manager_viper.yml .

Have a look at this file as well. It may look complicated but if fact it
is just a translation of the user manual of the cluster. For a new cluster
one person has to create it once. Also note the it only mentions which queueing
systems is installed how to use this systems is programmed in abipy. To
use this manager move it to manager.yml. (abipy will first look for a manager
file in you current folder and secondly in ~/.abinit/abipy, so you can put
one there an don't bother about it for every calculation)

This flow is made by the command:

\033[92m In []:\033[0m flow = lesson.make_g0w0_scissors_flow()

'flow' is now an object that contains al the information needed
to generate abinit inputs.

\033[92m In []:\033[0m flow.show_inputs()

will display all the inputs as they will be 'given' to abinit. In
previous lessons we ran the flows each time directly inside ipython.
For relatively small calculations this is very practical. There are
however other ways more suited for large calculations.

To start the execution of calculations packed in this flow we use the following command:

\033[92m In []:\033[0m flow.make_scheduler().start()

This starts the actual execution via a scheduler. 

The last step of analyzing the results can be done again in with a single command:

\033[92m In []:\033[0m flow.analyze()

This method of flow will open the necessary output files, retrieve
the data, and produce a plot.

Finally, once you are through with this lesson and exited ipython:

\033[92m In []:\033[0m exit

You can see that in the directory that you were working there is
now a subdir were the calculation have been performed. Have a look
at these folders and the files that are in them.

\033[93m Exercises \033[0m

As an exercise you can now start this lesson again but in stead
of performing the convergence study for silicon study the
convergence for a metal. By using:

\033[92m In []:\033[0m flow = lesson.make_ngkpt_flow(structure_file=lesson.abidata.cif_file('al.cif'), metal=True)

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
of the classes, methods and functions we used. 
You can get a copy of the file by using:

\033[92m In []:\033[0m lesson.get_local_copy()

"""
from __future__ import division, print_function, unicode_literals

import os
import sys
import shutil
import abipy.data as abidata  
import abipy.abilab as abilab


def help(stream=sys.stdout):
    """
    Display the tutorial text.
    """
    stream.write(__doc__)


def get_local_copy():
    """
    Copy this script to the current working dir to explore and edit
    """
    dst = os.path.basename(__file__[:-1])
    if os.path.exists(dst):
        raise RuntimeError("file %s already exists. Remove it before calling get_local_copy" % dst)
    shutil.copyfile(__file__[:-1], dst)


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


def make_g0w0_scissors_flow():
    # Change the value of ngkpt below to perform a GW calculation with a different k-mesh.
    scf, bands_nscf, dos_nscf, gw_nscf, scr, sig = make_inputs(ngkpt=[2,2,2])

    flow = abilab.Flow(workdir="flow_lesson_g0w0")
    work0 = abilab.BandStructureWork(scf, bands_nscf, dos_inputs=dos_nscf)
    flow.register_work(work0)

    work1 = abilab.Work()
    gw_nscf_task = work1.register_nscf_task(gw_nscf, deps={work0[0]: "DEN"})
    scr_task = work1.register_scr_task(scr, deps={gw_nscf_task: "WFK"})
    sigma_task = work1.register_sigma_task(sig, deps={gw_nscf_task: "WFK", scr_task: "SCR"})
    flow.register_work(work1)

    return flow.allocate()


def analyze_flow(flow, domains_spin=[[-10, 6.02], [6.1, 20]]):
    sigma_task = flow[1][2]
    builder = sigma_task.get_scissors_builder()

    #builder.plot_qpe_vs_e0()
    builder.build(domains_spin=domains_spin)
    builder.plot_fit()

    bands_task = flow[0][1]
    bands_filepath = bands_task.outdir.has_abiext("GSR")
    builder.plot_qpbands(bands_filepath, title="Silicon Bands (KS and KS+scissors)")

    # TODO: Fix problems with boundaries!
    #dos_task = flow[0][2]
    #dos_filepath = dos_task.outdir.has_abiext("GSR")
    #builder.plot_qpbands(bands_filepath, dos_filepath=dos_filepath,
    #                     title="Silicon Bands and DOS (KS and KS+scissors)")


if __name__ == "__main__":
    #flow = make_g0w0_scissors_flow()
    #flow.show_inputs()
    #flow.make_scheduler().start()

    flow = abilab.Flow.pickle_load("flow_lesson_g0w0")
    analyze_flow(flow)


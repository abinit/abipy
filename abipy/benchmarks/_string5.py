#!/usr/bin/env python
"""
Hydronium ion + NH3 molecule. String method.
Moving the proton from H2O to NH3 keeping O and H atoms fixed.
Based on tutoparal/tstring_04.in
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import operator
import numpy as np
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input(): # pragma: no cover
    """
    """
    pseudos = abidata.pseudos("8o_hard.paw", "7n.paw", "1h.paw")

    # Crystal structure and definition of the path
    xangst = np.fromstring("""
        0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
       -3.7593832509E-01 -2.8581911534E-01  8.7109635973E-01
       -3.8439081179E-01  8.6764073738E-01 -2.8530130333E-01
        4.0000000000E+00  0.0000000000E+00  0.0000000000E+00
        4.3461703447E+00 -9.9808458269E-02 -9.5466143436E-01
        4.3190273240E+00 -7.8675247603E-01  5.6699786920E-01
        4.3411410402E+00  8.7383785043E-01  4.0224838603E-01
        1.0280313162E+00  2.2598784215E-02  1.5561763093E-02""", sep=" ").reshape((-1,3))

    xangst_lastimg = np.fromstring("""
        0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
       -3.0400286349E-01 -1.9039526061E-01  9.0873550186E-01
       -3.2251946581E-01  9.0284480687E-01 -1.8824324581E-01
        4.0000000000E+00  0.0000000000E+00  0.0000000000E+00
        4.4876385468E+00 -1.4925704575E-01 -8.9716581956E-01
        4.2142401901E+00 -7.8694929117E-01  6.3097154506E-01
        4.3498225718E+00  8.7106686509E-01  4.2709343135E-01
        2.9570301511E+00  5.5992672027E-02 -1.3560839453E-01""", sep=" ").reshape((-1,3))

    structure = abilab.Structure.from_abivars(
        acell = abilab.ArrayWithUnit([10, 5, 5], "ang").to("bohr"),
        rprim=np.eye(3),
        typat=[1, 3, 3, 2, 3, 3, 3, 3], # Type of atoms (H2O + NH3 + H)
        znucl=[8.0, 7.0, 1.0],
        xangst=xangst,
    )

    inp = abilab.AbinitInput(structure, pseudos)

    inp.set_vars(
        # Options for parallelism
        paral_kgb=1,

        # Input/output options
        prtwf=0,
        prtden=0,
        prteig=0,
        prtdensph=0,
        timopt=-1,

        # Convergence parameters
        ecut=20.,
        pawecutdg=40.,

        # Control of the relaxation    # TO BE SUPPRESSED WHEN USING IMGMOV keyword
        #ionmov 3                      # BFGS (Broyden) algorithm for ions relaxation
        #optcell 0                     # No cell optimization
        #ntime 20                      # Max. number of "time" steps
        #tolmxf 5.0d-5                 # Stopping criterion of relaxation cycle

        # Control of the SCF cycle
        toldff=5.0e-7,
        nstep=50,

        # Electronic configuration
        nband=10,
        occopt=1,
        ixc="-001009",              # Select LDA XC functional (LDA PZ from LibXC)

        # BZ sampling
        kptopt=0,                 # Scheme for k-points generation
        nkpt=1,
        kpt=3*[0.],               # Explicit k-point (gamma point)
        nsym=1,                   # No symmetry

        natfix=2,
        iatfix=[1, 4],            # Keep O and N atoms fixed
        charge=1.0,               # Charge of the simulation cell

        # IMAGE section.
        nimage=12,                # Number of points along the string
        imgmov=2,                 # Selection of "String Method" algo
        ntimimage=50,             # Max. number of relaxation steps of the string
        tolimg=0.0001,            # Tol. criterion (will stop when average energy of cells < tolimg)
        dynimage="0 10*1 0",      # Keep first and last images fixed
        fxcartfactor=1.0,         # Time step for evolution step of string method.
        prtvolimg=2,              # Printing volume (0=full, 1=intermediate, 2=minimal)
        xangst_lastimg=xangst_lastimg,
    )

    return inp


def build_flow(options): # pragma: no cover
    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    template = make_input()

    # Processor distribution.
    #pconfs = [
    #  dict
    #npimage=10,  # CPU distribution over images
    #npband=10,
    #npfft=2,
    #bandpp=1,    # CPU distribution for 20 CPU cores per image
    #]

    # Get the list of possible parallel configurations from abinit autoparal.
    #max_ncpus, min_eff = options.max_ncpus, options.min_eff
    #print("Getting all autoparal configurations up to max_ncpus: ",max_ncpus," with efficiency >= ",min_eff)
    #pconfs = template.abiget_autoparal_pconfs(max_ncpus, autoparal=1, verbose=options.verbose)
    #if options.verbose:
    #print(pconfs)

    #%% nprocs_to_test = 200

    #work = flowtk.Work()
    #for d, omp_threads in product(pconfs, options.omp_list):
    #    mpi_procs = reduce(operator.mul, d.values(), 1)
    #    if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
    #    manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
    #    print("wfoptalg:", wfoptalg, "done with MPI_PROCS:", mpi_procs, "and:", d)
    #    inp = template.new_with_vars(d)
    #    work.register_scf_task(inp, manager=manager)

	#flow.register_work(work)

    return flow.allocate()


@bench_main
def main(options):  # pragma: no cover
    if options.info:
        # print doc string and exit.
        print(__doc__)
        return

    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())

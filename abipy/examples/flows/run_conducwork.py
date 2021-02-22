#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import sys
import abipy.data as abidata
import abipy.abilab as abilab
from abipy import flowtk
import abipy.abio.factories as factory


def make_scf_input(structure, pseudos, ngkpt=(2,2,2), shiftk=(0,0,0),
                   **variables):
    """Build and return SCF input given the structure and pseudopotentials"""

    scf_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    # Global variables
    scf_inp.set_vars(**variables)

    # Dataset 1 (GS run)
    scf_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
    #scf_inp.set_vars(toldfe=1e-10)

    return scf_inp


def make_nscf_input(structure, pseudos, ngkpt=(2,2,2), shiftk=(0,0,0),
                    **variables):
    """Build and return NSCF input given the structure and pseudopotentials"""

    scf_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    # Global variables
    scf_inp.set_vars(**variables)

    # Dataset 1 (GS run)
    scf_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
    scf_inp.set_vars(iscf=-2)

    return scf_inp


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Get structure and pseudos from the abipy database
    structure = abidata.structure_from_ucell("Al")
    pseudos = abidata.pseudos("13al.pspnc")

    # Variables
    variables = dict(
            ecut=20,
            tsmear=0.05,
            nband=12,
            nbdbuf=2,
            occopt=3,
            iomode=1,
            nstep=20)

    ngkpt = [4, 4, 4]
    ngkpt_fine = [8, 8, 8]
    shiftk = [0.0, 0.0, 0.0]
    ngqpt = [2, 2, 2]
    tmesh = [0, 30, 11] # Conductivity at temp from 0K to 300K by increment of 30

    #Kerange Variables
    nbr_proc = 4
    ngqpt_fine = [16, 16, 16] # The sigma_ngkpt grid must be divisible by the qpt grid
    sigma_ngkpt = [16, 16, 16]
    einterp = [1, 5, 0, 0] # Star functions Interpolation
    sigma_erange = [-0.3, -0.3, "eV"] # Negative value for metals

    # Nom de mon flow
    flow = flowtk.Flow(workdir=options.workdir)

    # Create inputs Object
    scf_input = make_scf_input(structure, pseudos,
                               tolvrs=1e-12,
                               ngkpt=ngkpt,
                               shiftk=shiftk,
                               **variables)

    nscf_input = make_nscf_input(structure, pseudos,
                                 tolwfr=1e-18,
                                 ngkpt=ngkpt_fine,
                                 shiftk=shiftk,
                                 **variables)

    # Create Work Object
    # Work 0 : Calcul SCF
    gs_work = flowtk.Work()
    gs_work.register_scf_task(scf_input)
    flow.register_work(gs_work)

    # Work 1 : Calcul DDB et DVDB
    ph_work = flowtk.PhononWork.from_scf_task(gs_work[0],
                                              qpoints=ngqpt, is_ngqpt=True,
                                              tolerance={"tolvrs": 1e-8})
    flow.register_work(ph_work)

    # Work 2 : Conduc with Kerange
    multi = factory.conduc_kerange_from_inputs(scf_input=scf_input,
                               nscf_input=nscf_input,
                               tmesh=tmesh,
                               ddb_ngqpt=ngqpt,
                               eph_ngqpt_fine=ngqpt_fine,
                               sigma_ngkpt=sigma_ngkpt,
                               sigma_erange=sigma_erange,
                               einterp=einterp)

    # Here we can change multi to change the variable of a particular dataset

    conduc_work = flowtk.ConducWork.from_phwork(phwork=ph_work, # Linking the DDB and DVDB via a PhononWork
                                                multi=multi, # The multidataset object
                                                nbr_proc=nbr_proc, # Needed to parallelize the calculation
                                                flow=flow,
                                                withKerange=True, # Using Kerange
                                                skipInter=True, # Doing DVDB interpolation during the conductivity task
                                                omp_nbr_thread=1) # The default value, no need to specify it in this case
    # If you already have the DDB and DVDB, use from_filepath(DDB, DVDB, multi, ...) instead of from_phwork

    flow.register_work(conduc_work)

    return flow.allocate(use_smartio=True)


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())

############################################

#!/usr/bin/env python
"""Band structure of silicon with the HT interface."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import abipy.data as abidata  

from abipy import abilab
from pymatgen.io.abinitio.calculations import bandstructure_work


def build_flow(options):
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    scf_kppa = 40
    nscf_nband = 6
    ndivsm = 5
    #dos_ngkpt = [4,4,4]
    #dos_shiftk = [0.1, 0.2, 0.3]

    extra_abivars = dict(
        ecut=6, 
        timopt=-1,
        accesswff=3, 
        istwfk="*1",
    )

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

    # Initialize the flow.
    # FIXME  Abistructure is not pickleable with protocol -1
    flow = abilab.Flow(workdir=workdir, manager=manager, pickle_protocol=0)

    work = bandstructure_work(structure, abidata.pseudos("14si.pspnc"), scf_kppa, nscf_nband, ndivsm, 
                              spin_mode="unpolarized", smearing=None, **extra_abivars)

    flow.register_work(work)
    return flow.allocate()

    #dos_kppa = 10
    #bands = bandstructure_work("hello_dos", runmode, structure, pseudos, scf_kppa, nscf_nband,
    #                      ndivsm, accuracy="normal", spin_mode="polarized",
    #                      smearing="fermi_dirac:0.1 eV", charge=0.0, scf_solver=None,
    #                      dos_kppa=dos_kppa)


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())

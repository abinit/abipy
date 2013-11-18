#!/usr/bin/env python
"""Flow for computing the band structure of silicon."""
from __future__ import division, print_function

import sys
import os
import abipy.data as abidata  
import abipy.abilab as abilab

from abipy.data.runs import enable_logging, AbipyTest, MixinTest

class EbandsFlowTest(AbipyTest, MixinTest):
    """
    Unit test for the flow defined in this module.  
    Users who just want to learn how to use this flow can ignore this section.
    """
    def setUp(self):
        super(FlowTest, self).setUp()
        self.init_dirs()
        self.flow = build_bands_flow(self.workdir)

    def move_files(self):
        pass
        # (wi, ti) --> (out_file, ref_file)
        files_to_move = {
            (0, 0): {"out_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc"}
        }
        for (wi, ti), d in files_to_move.items():
            task = flow[wi][ti]
            for out_file, ref_file in d.items():
                print("Will move out_file %s to ref_file %s" % (out_file, ref_file))
                #task.rename(out_file, ref_file)

        #Remove all files except those matching these regular expression.
        #work.rmtree(exclude_wildcard="*.abi|*.abo|*_WFK*|*_GSR.nc|*DEN-etsf.nc")
        #work = flow[0]
        #                                                                                   
        #work[0].rename("out_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
        #work[0].rename("out_DEN-etsf.nc", "si_DEN-etsf.nc")
        #work[0].rename("out_GSR.nc", "si_scf_GSR.nc")
        #
        #work[1].rename("out_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
        #work[1].rename("out_GSR.nc", "si_nscf_GSR.nc")
        #work[1].remove_files("out_DS2_DEN-etsf.nc")
        #flow.rmtree()


def make_scf_nscf_inputs():
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos("14si.pspnc")
    #pseudos = data.pseudos("Si.GGA_PBE-JTH-paw.xml")

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    print(inp.pseudos)
    structure = inp.set_structure_from_file(abidata.cif_file("si.cif"))

    # Global variables
    ecut = 6
    global_vars = dict(ecut=ecut,
                       nband=8,
                       timopt=-1,
                       accesswff=3,
                       istwfk="*1",
                    )

    if inp.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    inp[1].set_variables(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    inp[2].set_kpath(ndivsm=6, kptbounds=kptbounds)
    inp[2].set_variables(tolwfr=1e-12)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()
    return scf_input, nscf_input

def build_bands_flow(workdir):

    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs()

    # Instantiate the TaskManager from `taskmanager.yml`.
    manager = abilab.TaskManager.from_user_config()
                                                               
    # Build the flow.
    flow = abilab.bandstructure_flow(workdir, manager, scf_input, nscf_input)
    return flow.allocate()
    

@enable_logging
def main():
    flow = build_bands_flow(workdir="tmp_si_ebands")
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())

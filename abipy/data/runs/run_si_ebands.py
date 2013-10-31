#!/usr/bin/env python
"""Flow for computing the band structure of silicon."""
from __future__ import division, print_function

import sys
import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester, enable_logging
from abipy.core.testing import AbipyTest

###########################################################
# Begin unit test for the flow defined in this module. 
# Users who just want to learn how to use this flow can ignore this section.
############################################################
import unittest

class MixinTest(object):
    #DORUN = False
    DORUN = True

    def init_dirs(self):
        import inspect
        import os
        ## Get the filename of the calling module.
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        #
        apath = os.path.abspath(mod.__file__)
        print(apath)
        #assert 0
                                                                             
        #base = os.path.basename(apath).replace(".py","").replace("run_","")
        # Results will be produced in workdir. 
        # refdir contains the reference data (might not be present if this 
        # is the first time we execute the AbinitFlow.
        #self.workdir = os.path.join(os.path.dirname(apath), "tmp_" + base)
        #self.refdir = os.path.join(os.path.dirname(apath), "data_" + base)

    def test_pickle(self):
        """Testing whether the flow object is pickleable."""
        self.serialize_with_pickle(self.flow, protocols=None, test_eq=True)

    def test_run(self):
        """Running the flow with PyFlowsScheduler."""
        if not self.DORUN:
            print("Skipping test_run")
            return 

        from pymatgen.io.abinitio.launcher import PyFlowsScheduler
        self.flow.build_and_pickle_dump()

        sched_options = dict(
            weeks=0,
            days=0,
            hours=0,
            minutes=0,
            seconds=5,
        )

        sched = PyFlowsScheduler(**sched_options)
        sched.add_flow(self.flow)
        sched.start()

        all_ok = self.flow.all_ok
        self.assertTrue(all_ok)

    def tearDown(self):
        if not self.DO_RUN:
            return 

        if not self.flow.all_ok:
            return 


class FlowTest(AbipyTest, MixinTest):
    def setUp(self):
        super(FlowTest, self).setUp()
        self.flow = build_bands_flow("__hello_test__")

        self.init_dirs()

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
#########################
# End unit test section #
#########################


def make_scf_nscf_inputs():
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=2)
    structure = inp.set_structure_from_file(data.cif_file("si.cif"))

    # Global variables
    global_vars = dict(ecut=6,
                       nband=8,
                       timopt=-1,
                       accesswff=3,
                       istwfk="*1",
                    )

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
    tester = Tester()
    flow = build_bands_flow(tester.workdir)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())

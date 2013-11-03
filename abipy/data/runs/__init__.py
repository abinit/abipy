from __future__ import division, print_function

import os
import inspect
import functools

from abipy.core.testing import AbipyTest
from pymatgen.util.decorators import enable_logging
from pymatgen.io.abinitio.launcher import PyFlowsScheduler

__all__ = [

]

class MixinTest(object):
    #DORUN = False
    DORUN = True

    def init_dirs(self):
        ## Get the filename of the calling module.
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        #
        apath = os.path.abspath(mod.__file__)
        print(apath)
        #assert 0
                                                                             
        base = os.path.basename(apath).replace(".py","").replace("run_","")
        # Results will be produced in workdir. 
        # refdir contains the reference data (might not be present if this 
        # is the first time we execute the AbinitFlow.
        self.workdir = os.path.join(os.path.dirname(apath), "tmp_" + base)
        #self.refdir = os.path.join(os.path.dirname(apath), "data_" + base)

    def test_pickle(self):
        """Testing whether the flow object is pickleable."""
        self.serialize_with_pickle(self.flow, protocols=None, test_eq=True)

    def test_run(self):
        """Running the flow with PyFlowsScheduler..."""
        if not self.DORUN:
            print("Skipping test_run")
            return 

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

    #def move_files(self):
    #    pass
    #    # (wi, ti) --> (out_file, ref_file)
    #    files_to_move = {
    #        (0, 0): {"out_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc"}
    #    }
    #    for (wi, ti), d in files_to_move.items():
    #        task = flow[wi][ti]
    #        for out_file, ref_file in d.items():
    #            print("Will move out_file %s to ref_file %s" % (out_file, ref_file))
    #            #task.rename(out_file, ref_file)

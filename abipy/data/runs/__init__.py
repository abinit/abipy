import os
import inspect
import functools
#import abipy.abilab as abilab

#from pymatgen.io.abinitio.launcher import PyFlowsScheduler
from pymatgen.util.decorators import enable_logging

__all__ = [
    "Tester",
]

class Tester(object):
    def __init__(self):
        # Get the filename of the calling module.
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        
        apath = os.path.abspath(mod.__file__)

        base = os.path.basename(apath).replace(".py","").replace("run_","")

        # Results will be produced in workdir. 
        # refdir contains the reference data (might not be present if this 
        # is the first time we execute the AbinitFlow.
        self.workdir = os.path.join(os.path.dirname(apath), "tmp_" + base)
        self.refdir = os.path.join(os.path.dirname(apath), "data_" + base)

    #def __str__(self):
    #    return str(self.__dict__)

    #@property
    #def retcode(self):
    #    """Return code."""
    #    try:
    #        return self._retcode
    #    except AttributeError:
    #        return 66

    #def set_flow_and_run(self, flow):
    #    self.set_flow(flow)
    #    self.run()

    #def run(self):
    #    raise NotImplementedError("")
    #    #retcodes = PyResourceManager(self.flow, max_ncpus=2, sleep_time=2.5).run()
    #    #self._retcode = max(retcodes)

    #def finalize(self):
        #pass
        #self.work.rm_indatadir()
        #self.work.rm_tmpdatadir()

        #refdir, workdir = self.refdir, self.workdir
        #if not os.path.exists(refdir):
        #    self.work.move(refdir)

        #else:
        #    diffs = {}
        #    for dirpath, dirnames, filenames in os.walk(refdir):
        #        for fname in filenames:
        #            ref_path = os.path.join(dirpath, fname)
        #            new_path = os.path.join(workdir, os.path.relpath(ref_path, start=refdir))
        #            diffs[ref_path] = self.compare(ref_path, new_path)
        #    print(diffs)

    #@staticmethod
    #def compare(ref_path, new_path):
    #    if ref_path.endswith(".abo"):
    #        with open(ref_path, "r") as ref, open(new_path,"r") as new:
    #            ref_lines = ref.readlines()
    #            new_lines = new.readlines()
    #            return ref_lines != new_lines
    #    else:
    #        return "no comparison for file %s" % ref_path


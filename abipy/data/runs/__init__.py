import os
import inspect

from pymatgen.io.abinitio.launcher import PyResourceManager
import abipy.abilab as abilab


class Tester(object):
    def __init__(self):
        
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])

        apath = os.path.abspath(mod.__file__)

        base = os.path.basename(apath).replace(".py","").replace("run_","")

        self.refdir = os.path.join(os.path.dirname(apath), "data_" + base)
        self.workdir = os.path.join(os.path.dirname(apath), "tmp_" + base)

    def __str__(self):
        return str(self.__dict__)

    @property
    def retcode(self):
        """Return code the workflow."""
        try:
            return self._retcode
        except AttributeError:
            return 66

    def set_work(self, work):
        self.work = work

    def make_manager(self):
        return abilab.TaskManager.simple_mpi(mpi_ncpus=1)

    def run(self):
        retcodes = PyResourceManager(self.work, max_ncpus=2, sleep_time=5).run()
        self._retcode = max(retcodes)

    def set_work_and_run(self, work):
        self.set_work(work)
        self.run()
                                                                                
    def finalize(self):
        self.work.rm_indatadir()
        self.work.rm_tmpdatadir()

        refdir, workdir = self.refdir, self.workdir

        if not os.path.exists(refdir):
            self.work.move(refdir)

        else:
            diffs = {}
            for dirpath, dirnames, filenames in os.walk(refdir):
                for fname in filenames:
                    ref_path = os.path.join(dirpath, fname)
                    new_path = os.path.join(workdir, os.path.relpath(ref_path, start=refdir))
                    diffs[ref_path] = self.compare(ref_path, new_path)
            print(diffs)

    @staticmethod
    def compare(ref_path, new_path):

        if ref_path.endswith(".about"):
            with open(ref_path, "r") as ref, open(new_path,"r") as new:
                ref_lines = ref.readlines()
                new_lines = new.readlines()
                return ref_lines != new_lines
        else:
            return "no comparison for file %s" % ref_path

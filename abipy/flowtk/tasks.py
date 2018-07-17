from __future__ import print_function, division, unicode_literals, absolute_import

from pymatgen.io.abinit.tasks import *


class BoxcuttedPhononTask(PhononTask):
    """
    This task compute phonons with a two-step algorithm.
    The first DFPT run is done with low-accuracy settings for boxcutmin and ecut
    The second DFPT run uses boxcutmin 2.0 and normal ecut and restarts from
    the 1WFK file generated previously.
    """
    @classmethod
    def patch_flow(cls, flow):
        for task in flow.iflat_tasks():
            if isinstance(task, PhononTask): task.__class__ = cls

    def setup(self):
        super(BoxcuttedPhononTask, self).setup()
        self.final_dfp_done = False if not hasattr(self, "final_dfp_done") else self.final_dfp_done
        if not self.final_dfp_done:
            # First run: use boxcutmin 1.5 and low-accuracy hints (assume pseudos with hints).
            pseudos = self.input.pseudos
            ecut = max(p.hint_for_accuracy("low").ecut for p in pseudos)
            pawecutdg = max(p.hint_for_accuracy("low").pawecutdg for p in pseudos) if self.input.ispaw else None
            self.set_vars(boxcutmin=1.5, ecut=ecut, pawecutdg=pawecutdg, prtwf=1)

    def _on_ok(self):
        results = super(BoxcuttedPhononTask, self)._on_ok()
        if not self.final_dfp_done:
            # Second run: use exact box and normal-accuracy hints (assume pseudos with hints).
            pseudos = self.input.pseudos
            ecut = max(p.hint_for_accuracy("normal").ecut for p in pseudos)
            pawecutdg = max(p.hint_for_accuracy("normal").pawecutdg for p in pseudos) if self.input.ispaw else None
            self.set_vars(boxcutmin=2.0, ecut=ecut, pawecutdg=pawecutdg, prtwf=-1)
            self.finalized = True
            self.final_dfp_done = True
            self.restart()

        return results


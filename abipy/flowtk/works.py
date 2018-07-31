# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from pymatgen.io.abinit.works import Work, MergeDdb
from abipy.abio import input_tags as tags


class ElasticWork(Work, MergeDdb):
    """
    Work for the calculation of elastic constants and (optionally piezoelectric tensor
    """

    @classmethod
    def from_scf_input(cls, scf_input, tolerance=None, with_internal_strain=True, with_piezoelectric=False, manager=None):
        """
        Similar to `from_scf_task`, the difference is that this method requires
        an input for SCF calculation instead of a ScfTask. All the tasks (Scf + Phonon)
        are packed in a single Work whereas in the previous case we usually have multiple works.
        """
        new = cls(manager=manager)

        # Register task for SCF calculation.
        scf_task = new.register_scf_task(scf_input)

        multi = scf_task.input.make_strain_perts_inputs(tolerance=tolerance, manager=manager,
                phonon_pert=with_internal_strain, kptopt=2)

        ddk_tasks = []
        if with_piezoelectric:
            ddk_multi = scf_task.input.make_ddk_inputs(tolerance=None, manager=manager)
            #for inp in multi.filter_by_tags(tags=tags.DDK):
            for inp in ddk_multi:
                ddk_task = new.register_ddk_task(inp, deps={scf_task: "WFK"})
                ddk_tasks.append(ddk_task)
            assert len(ddk_tasks) == 3
            ddk_deps = {ddk_task: "DDK" for ddk_task in ddk_tasks}

        if with_internal_strain:
            ph_deps = {scf_task: "WFK"}
            if with_piezoelectric: ph_deps.update(ddk_deps)
            #for inp in multi.filter_by_tags(tags=tags.PHONON):
            for inp in multi:
                if inp.get("rfphon", 0) == 1:
                    #new.register_phonon_task(inp, deps=ph_deps)
                    new.register_bec_task(inp, deps=ph_deps)

        elast_deps = {scf_task: "WFK"}
        if with_piezoelectric: elast_deps.update(ddk_deps)
        #for inp in multi.filter_by_tags(tags=tags.STRAIN)
        for inp in multi:
            # FIXME
            #new.register_elastic_task(inp, deps=elast_deps)
            new.register_bec_task(inp, deps=elast_deps)

        return new

    def on_all_ok(self):
        """
        This method is called when all the tasks reach S_OK.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files(delete_source_ddbs=False)

        results = self.Results(node=self, returncode=0, message="DDB merge done")
        return results

# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from pymatgen.io.abinit.works import Work, MergeDdb


class ElasticWork(Work, MergeDdb):
    """
    This Work computes the elastic constants and (optionally) the piezoelectric tensor.
    It consists of Response function calculations for:

        * rigid-atom elastic tensor
        * rigid-atom piezoelectric tensor
        * interatomic force constants at gamma
        * Born effective charges

    The structure is assumed to be already relaxed

    Create a `Flow` for phonon calculations. The flow has one works with:

	- 1 GS Task
	- 3 DDK Task
	- 4 Phonon Tasks (Gamma point)
	- 6 Elastic tasks (3 uniaxial + 3 shear strain)

    The Phonon tasks and the elastic task will read the DDK produced at the beginning
    """

    @classmethod
    def from_scf_input(cls, scf_input, tolerance=None, with_internal_strain=True,
                       with_piezoelectric=False, manager=None):
        """
        Similar to `from_scf_task`, the difference is that this method requires
        an input for SCF calculation instead of a ScfTask. All the tasks (Scf + Phonon)
        are packed in a single Work whereas in the previous case we usually have multiple works.
        """
        new = cls(manager=manager)
        #scf_input.deepcopy().pop_relax_vars()

        # Register task for SCF calculation.
        scf_task = new.register_scf_task(scf_input)

        if with_piezoelectric:
            # Calculate the ddk wf's needed for piezoelectric tensor and Born effective charges.
            ddk_tolerance = {"tolwfr": 1.0e-20}
            ddk_multi = scf_task.input.make_ddk_inputs(tolerance=ddk_tolerance, manager=manager)
            ddk_tasks = []
            for inp in ddk_multi:
                ddk_task = new.register_ddk_task(inp, deps={scf_task: "WFK"})
                ddk_tasks.append(ddk_task)
            ddk_deps = {ddk_task: "DDK" for ddk_task in ddk_tasks}

        # Build input files for strain and (optionally) phonons.
        tolerance = {"tolvrs": 1e-10}
        strain_multi = scf_task.input.make_strain_perts_inputs(tolerance=tolerance, manager=manager,
            phonon_pert=with_internal_strain, kptopt=2)

        if with_internal_strain:
            # Phonon perturbation (read DDK if piezo).
            ph_deps = {scf_task: "WFK"}
            if with_piezoelectric: ph_deps.update(ddk_deps)
            for inp in strain_multi:
                if inp.get("rfphon", 0) == 1:
                    new.register_phonon_task(inp, deps=ph_deps)

        # Finally compute strain pertubations (read DDK if piezo)
        elast_deps = {scf_task: "WFK"}
        if with_piezoelectric: elast_deps.update(ddk_deps)
        for inp in strain_multi:
            if inp.get("rfstrs", 0) != 0:
                new.register_elastic_task(inp, deps=elast_deps)

        return new

    def on_all_ok(self):
        """
        This method is called when all the tasks of the Work reach S_OK.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files(delete_source_ddbs=False, only_dfpt_tasks=False)

        results = self.Results(node=self, returncode=0, message="DDB merge done")
        return results

# coding: utf-8
"""Work subclasses related to DFTP."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from .works import Work, MergeDdb


class NscfDdksWork(Work):
    """
    This work requires a DEN file and computes the KS energies with a non self-consistent task
    with a dense k-mesh and empty states.
    This task is then followed by the computation of the DDK matrix elements with nstep = 1
    (the first order change of the wavefunctions is not converged but we only need the matrix elements)
    Mainly used to prepare optic calculations or other post-processing steps requiring the DDKs.
    """

    @classmethod
    def from_scf_task(cls, scf_task, ddk_ngkpt, ddk_shiftk, ddk_nband, manager=None):
        """
        Build NscfDdksWork from a scf_task.

        Args:
            scf_task: GS task. Must produce the DEN file required for the NSCF run.
            ddk_ngkpt: k-mesh used for the NSCF run and the non self-consistent DDK tasks.
            ddk_shiftk: k-mesh shifts
            ddk_nband: Number of bands (occupied + empty) used in the NSCF task and the DDKs tasks.
            manager: TaskManager instance. Use default if None.

        Return: NscfDdksWork instance
        """
        new = cls(manager=manager)

        # NSCF task with nband states and points in the IBZ (note kptopt = 1)
        nscf_inp0 = scf_task.input.deepcopy()
        nscf_inp0.set_vars(nband=ddk_nband, prtwf=1)
        nscf_inp0.set_kmesh(ddk_ngkpt, ddk_shiftk, kptopt=1)
        nscf_task0 = new.register_nscf_task(nscf_inp0, deps={scf_task: "DEN"})

        # NSCF run with nband states and points in the IBZ defined by time-reversal only (as required by DDK)
        # This is gonna be quick because Abinit will symmetrize states from the previous WFK file.
        # Time-reversal symmetry can be used in optic.
        #nscf_inp1 = nscf_inp0.deepcopy()
        #nscf_inp0.set_kmesh(ddk_ngkpt, ddk_shiftk, kptopt=2)
        #nscf_task1 = new.register_nscf_task(nscf_inp1)

        # This is the task producing the KS energies for optic
        new.task_with_ks_energies = nscf_task0

        # Build task for one-shot DDKs (note kptopt 2)
        ddk_inputs = nscf_inp0.make_ddk_inputs(kptopt=2)
        new.ddk_tasks = []
        for ddk_inp in ddk_inputs:
            # FIXME: prtwfk should be set to 0 but need to replace DDK.nc
            ddk_inp.set_vars(nstep=1, nline=0, prtwf=1)
            #new.register_ddk_task(ddk_inp, deps={nscf_task0: "WFK"})
            # FIXME: Here I have a conflict with DDK.nc and DDK
            t = new.register_task(ddk_inp, deps={nscf_task0: "WFK"})
            new.ddk_tasks.append(t)

        return new


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
    def from_scf_input(cls, scf_input, with_relaxed_ion=True, with_piezo=False, with_dde=False,
                       tolerances=None, den_deps=None, manager=None):
        """
        Args:
            scf_input:
            with_relaxed_ion:
            with_piezo:
            with_dde: Compute electric field perturbations.
            tolerances: Dict of tolerances
            den_deps:
            manager:

        Similar to `from_scf_task`, the difference is that this method requires
        an input for SCF calculation instead of a ScfTask. All the tasks (Scf + Phonon)
        are packed in a single Work whereas in the previous case we usually have multiple works.
        """
        if tolerances is None: tolerances = {}
        new = cls(manager=manager)

        # Register task for WFK0 calculation (either SCF or NCSCF if den_deps is given)
        if den_deps is None:
            wfk_task = new.register_scf_task(scf_input)
        else:
            tolwfr = 1.0e-20
            if "nscf" in tolerances:
                tolwfr = tolerances["nscf"]["tolwfr"]
            nscf_input = scf_input.new_with_vars(iscf=-2, tolwfr=tolwfr)
            wfk_task = new.register_nscf_task(nscf_input, deps=den_deps)

        if with_piezo or with_dde:
            # Calculate the ddk wf's needed for piezoelectric tensor and Born effective charges.
            #ddk_tolerance = {"tolwfr": 1.0e-20}
            ddk_tolerance = tolerances.get("ddk", None)
            ddk_multi = scf_input.make_ddk_inputs(tolerance=ddk_tolerance, manager=manager)
            ddk_tasks = []
            for inp in ddk_multi:
                ddk_task = new.register_ddk_task(inp, deps={wfk_task: "WFK"})
                ddk_tasks.append(ddk_task)
            ddk_deps = {ddk_task: "DDK" for ddk_task in ddk_tasks}

        if with_dde:
            # Add tasks for electric field perturbation.
            #dde_tolerance = None
            dde_tolerance = tolerances.get("dde", None)
            dde_multi = scf_input.make_dde_inputs(tolerance=dde_tolerance, use_symmetries=True, manager=manager)
            dde_deps = {wfk_task: "WFK"}
            dde_deps.update(ddk_deps)
            for inp in dde_multi:
                new.register_dde_task(inp, deps=dde_deps)

        # Build input files for strain and (optionally) phonons.
        #strain_tolerance = {"tolvrs": 1e-10}
        strain_tolerance = tolerances.get("strain", None)
        strain_multi = scf_input.make_strain_perts_inputs(tolerance=strain_tolerance, manager=manager,
            phonon_pert=with_relaxed_ion, kptopt=2)

        if with_relaxed_ion:
            # Phonon perturbation (read DDK if piezo).
            ph_deps = {wfk_task: "WFK"}
            if with_piezo: ph_deps.update(ddk_deps)
            for inp in strain_multi:
                if inp.get("rfphon", 0) == 1:
                    new.register_phonon_task(inp, deps=ph_deps)

        # Finally compute strain pertubations (read DDK if piezo).
        elast_deps = {wfk_task: "WFK"}
        if with_piezo: elast_deps.update(ddk_deps)
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

# coding: utf-8
"""
Workflows for calculations within the quasi-harmonic approximation.

See [Phys. Rev. B 110, 014103](https://doi.org/10.1103/PhysRevB.110.014103)
"""
from __future__ import annotations

import numpy as np

from abipy.tools.serialization import mjson_write
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.flows import Flow


class VzsisaFlow(Flow):
    """
    Flow for QHA calculations with the VZSISA approach.
    Main entry point for client code.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VzsisaFlow
    """

    @classmethod
    def from_scf_input(cls,
                       workdir,
                       scf_input,
                       bo_vol_scales,
                       ph_vol_scales,
                       ngqpt,
                       with_becs: bool,
                       with_quad: bool,
                       ndivsm=-20,
                       edos_ngkpt=None,
                       manager=None) -> VzsisaFlow:
        """
        Build a flow for QHA calculations from an |AbinitInput| for GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate other inputs.
            bo_vol_scales: List of scaling factors for the BO volumes.
            ph_vol_scales: List of scaling factors for the BO volumes.
            ngqpt: Three integers defining the q-mesh for phonon calculation.
            with_becs: Activate calculation of Electric field and Born effective charges.
            with_quad: Activate calculation of dynamical quadrupoles. Require `with_becs`
                Note that only selected features are compatible with dynamical quadrupoles.
                Please consult <https://docs.abinit.org/topics/longwave/>
            ndivsm: if > 0, it's the number of divisions for the smallest segment of the path (Abinit variable).
                if < 0, it's interpreted as the pymatgen `line_density` parameter in which the number of points
                in the segment is proportional to its length. Typical value: -20.
                This option is the recommended one if the k-path contains two consecutive high symmetry k-points
                that are very close as ndivsm > 0 may produce a very large number of wavevectors.
                if 0, deactivate band structure calculation.
            edos_ngkpt: Three integers defining the the k-sampling for the computation of the
                electron DOS with the relaxed structures. Useful for metals or small gap semiconductors
                in which the electronic contribution should be included.
                None disables the computation of the e-DOS.
            manager: |TaskManager| instance. Use default if None.
        """
        flow = cls(workdir=workdir, manager=manager)

        flow.register_work(VzsisaWork.from_scf_input(scf_input, bo_vol_scales, ph_vol_scales, ngqpt, with_becs, with_quad,
                                                     ndivsm, ionmov=2, edos_ngkpt=edos_ngkpt))
        return flow

    def finalize(self):
        """
        This method is called by the scheduler when the flow is completed.
        Here we write the vszisa.json file in the outdata directory of the flow.
        This json file can be used to build the Vzsisa object used to post-process the data.
        """
        work = self[0]
        data = {"bo_vol_scales": work.bo_vol_scales, "ph_vol_scales": work.ph_vol_scales}
        #data["initial_structure"] = work.initial_scf_input.structure

        # Build list of strings with path to the relevant output files ordered by V.
        data["gsr_relax_paths"] = [task.gsr_path for task in work.relax_tasks_vol]

        gsr_relax_entries, gsr_relax_volumes = [], []
        for task in work.relax_tasks_vol:
            with task.open_gsr() as gsr:
                gsr_relax_entries.append(dict(
                    volume=gsr.structure.volume,
                    energy_eV=float(gsr.energy),
                    pressure_GPa=float(gsr.pressure),
                    #structure=gsr.structure,
                ))
                gsr_relax_volumes.append(gsr.structure.volume)

        data["gsr_relax_entries"] = gsr_relax_entries
        data["gsr_relax_volumes_ang3"] = gsr_relax_volumes

        data["ddb_relax_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in work.ph_works]
        data["ddb_relax_volumes_ang3"] = [ph_work[0].input.structure.volume for ph_work in work.ph_works]

        data["gsr_relax_edos_paths"] = [] if not work.edos_work else [task.gsr_path for task in work.edos_work]

        # Write json file.
        mjson_write(data, self.outdir.path_in("vzsisa.json"), indent=4)

        return super().finalize()


class VzsisaWork(Work):
    """
    This work performs a structural relaxation of the initial structure, then a set of distorted
    structures is genenerated and the relaxed structures are used
    to compute phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VzsisaWork
    """

    @classmethod
    def from_scf_input(cls, scf_input, bo_vol_scales, ph_vol_scales, ngqpt,
                       with_becs: bool, with_quad: bool, ndivsm: int,
                       ionmov: int, edos_ngkpt=None) -> VzsisaWork:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.
        See VzsisaFlow for rhe meaning of the arguments.
        """
        work = cls()

        # Save attributes in work
        work.initial_scf_input = scf_input
        work.bo_vol_scales = np.array(bo_vol_scales)
        work.ph_vol_scales = np.array(ph_vol_scales)
        work.ngqpt = ngqpt
        work.with_becs = with_becs
        work.with_quad = with_quad
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))
        work.ndivsm = ndivsm

        # Create input for relaxation and register the relaxation task.
        work.relax_template = relax_template = scf_input.deepcopy()

        # optcell = 3: constant-volume optimization of cell geometry
        relax_template.pop_tolerances()
        relax_template.set_vars(optcell=3, ionmov=ionmov, tolvrs=1e-8, tolmxf=1e-6)
        relax_template.set_vars_ifnotin(ecutsm=1.0, dilatmx=1.05)

        work.initial_relax_task = work.register_relax_task(relax_template)

        return work

    def on_ok(self, sender):
        """
        This method is called when one task reaches status `S_OK`.
        It executes on_all_ok when all tasks in self have reached `S_OK`.
        """
        if sender == self.initial_relax_task:
            # Get relaxed structure
            relaxed_structure = sender.get_final_structure()
            v0 = relaxed_structure.volume

            # build new tasks for structural relaxation at fixed volume.
            relax_template = self.relax_template
            self.relax_tasks_vol = []
            for bo_scale in self.bo_vol_scales:
                new_structure = relax_template.structure.scale_lattice(v0 * bo_scale)
                task = self.register_relax_task(relax_template.new_with_structure(new_structure))
                self.relax_tasks_vol.append(task)

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """
        This callback is called when all tasks in the Work reach status `S_OK`.
        Here we add a new PhononWork for each volume using the relaxed structure.
        """
        self.edos_work = None if self.edos_ngkpt is None else Work()
        self.ph_works = []

        # Build phonon works for the different relaxed structures.
        for task, bo_scale in zip(self.relax_tasks_vol, self.bo_vol_scales):
            if all(abs(bo_scale - self.ph_vol_scales)) > 1e-3: continue
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)

            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, with_quad=self.with_quad,
                                                ndivsm=0 if bo_scale != 0.0 else self.ndivsm)

            # Reduce the number of files produced in the DFPT tasks to avoid possible disk quota issues.
            for task in ph_work[1:]:
                task.input.set_vars(prtden=0, prtpot=0)

            self.flow.register_work(ph_work)
            self.ph_works.append(ph_work)

            # Add task for electron DOS calculation to edos_work
            if self.edos_work is not None:
                edos_input = scf_input.make_edos_input(self.edos_ngkpt, prtwf=-1)
                self.edos_work.register_nscf_task(edos_input, deps={ph_work[0]: "DEN"})

        if self.edos_ngkpt is not None:
            self.flow.register_work(self.edos_work)

        self.flow.allocate(build=True)

        return super().on_all_ok()

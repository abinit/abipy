# coding: utf-8
"""
Workflows for calculations within the quasi-harmonic approximation.
"""
from __future__ import annotations

import numpy as np

from abipy.tools.serialization import mjson_write
from abipy.dfpt.deformation_utils import generate_deformations
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.flows import Flow


class vZSISAFlow(Flow):
    """
    Flow for QHA calculations with the vZSISA approach

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: vZSISAFlow
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, bo_scales, ph_scales, ngqpt, with_becs, with_quad,
                       edos_ngkpt=None, manager=None) -> vZSISAFlow:
        """
        Build a flow for QHA calculations from an |AbinitInput| for GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            bo_scales:
            ph_scales:
            ngqpt: Three integers defining the q-mesh for phonon calculation.
            with_becs: Activate calculation of Electric field and Born effective charges.
            with_quad: Activate calculation of dynamical quadrupoles. Require `with_becs`
                Note that only selected features are compatible with dynamical quadrupoles.
                Please consult <https://docs.abinit.org/topics/longwave/>
            edos_ngkpt: Three integers defining the the k-sampling for the computation of the
                electron DOS with the relaxed structures. Useful for metals or small gap semiconductors
                in which the electronic contribution should be included.
                None disables the computation of the e-DOS.
            manager: |TaskManager| instance. Use default if None.
        """
        flow = cls(workdir=workdir, manager=manager)

        work = vZSISAWork.from_scf_input(scf_input, bo_scales, ph_scales, ngqpt, with_becs, with_quad,
                                         optcell=2, ionmov=2, edos_ngkpt=edos_ngkpt)

        flow.register_work(work)
        return flow

    def finalize(self):
        """
        This method is called when the flow is completed.
        It performs some basic post-processing of the results to facilitate further analysis.
        """
        work = self[0]
        data = {}

        # Build list of strings with path to the relevant output files ordered by V.
        data["gsr_relax_paths"] = [task.gsr_path for task in work.relax_tasks_vol]
        data["gsr_relax_volumes"] = [task.gsr_path for task in work.relax_tasks_vol]
        data["ddb_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in work.ph_works]
        data["gsr_edos_path"] = []
        if work.edos_work:
            data["gsr_edos_paths"] = [task.gsr_path for task in work.edos_work]

        #entries = []
        #items = zip(self.relax_and_phonon_work.relax_tasks_vol, self.relax_and_phonon_work.ph_works_vol)
        #for ivol, (relax_task, ph_work) in enumerate(items):
        #    ddb_path = ph_work.outdir.has_abiext("DDB")
        #    with relax_task.open_gsr() as gsr:
        #        entries.append(dict(
        #            volume=gsr.structure.volume,
        #            energy_eV=float(gsr.energy),
        #            pressure_GPa=float(gsr.pressure),
        #            structure=gsr.structure,
        #            gsr_edos_path=None,
        #            ddb_path=ddb_path,
        #        )

        #data["entries_for_each_volume"] = entries

        mjson_write(data, self.outdir.path_in("vzsisa.json"), indent=4)

        return super().finalize()


class vZSISAWork(Work):
    """
    This work performs a structural relaxation of the initial structure, then a set of distorted
    structures is genenerated and the relaxed structures are used
    to compute phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: vZSISAWork
    """

    @classmethod
    def from_scf_input(cls, scf_input, bo_scales, ph_scales, ngqpt,
                       with_becs: bool, with_quad: bool,
                       optcell: int, ionmov: int, edos_ngkpt=None) -> vZSISAWork:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            bo_scales:
            ph_scales:
            with_becs: Activate calculation of Electric field and Born effective charges.
            optcell, ionmov: Abinit input variables.
            edos_ngkpt: Three integers defining the the k-sampling for the computation of the
                electron DOS with the relaxed structures. Useful for metals
                in which the electronic contribution should be included.
                None disables the computation of the e-DOS.
        """
        work = cls()

        # Save attributes in work
        work.initial_scf_input = scf_input
        work.bo_scales = np.array(bo_scales)
        work.ph_scales = np.array(ph_scales)
        work.ngqpt = ngqpt
        work.with_becs = with_becs
        work.with_quad = with_quad
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))

        # Create input for relaxation and register the relaxation task.
        work.relax_template = relax_template = scf_input.deepcopy()
        relax_template.pop_tolerances()
        relax_template.set_vars(tolvrs=1e-8, toldff=1.e-6, optcell=optcell, ionmov=ionmov)
        if optcell is not None and optcell != 0:
            relax_template.set_vars_ifnotin(ecutsm=0.5, dilatmx=1.05)

        work.initial_relax_task = work.register_relax_task(relax_template)

        return work

    def on_ok(self, sender):
        """
        This method is called when one task reaches status `S_OK`.
        It executes on_all_ok when all tasks in self have reached `S_OK`.
        """
        if sender == self.initial_relax_task:
            # Get relaxed structure and build new task for structural relaxation at fixed volume.
            relaxed_structure = sender.get_final_structure()
            v0 = relaxed_structure.volume

            relax_template = self.relax_template
            self.relax_tasks_vol = []
            for bo_scale in self.bo_scales:
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
        for task, bo_scale in zip(self.relax_tasks_vol, self.bo_scales):
            if all(abs(bo_scale - self.ph_scales)) > 1e-3: continue
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)
            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, with_quad=self.with_quad,
                                                ddk_tolerance=None)

            self.flow.register_work(ph_work)
            self.ph_works.append(ph_work)

            # Add electron DOS calculation.
            if self.edos_work is not None:
                edos_input = scf_input.make_edos_input(self.edos_ngkpt)
                self.edos_work.register_nscf_task(edos_input, deps={ph_work[0]: "DEN"})

        if self.edos_ngkpt is not None: self.flow.register_work(self.edos_work)
        self.flow.allocate(build=True)

        return super().on_all_ok()


class ThermalRelaxTask(RelaxTask):

    def _on_ok(self):
        results = super()._on_ok()
        # Check for convergence.
        #if not self.collinear_done:
        #    self.input.set_vars(strtarget=strtarget)
        #    self.finalized = False
        #    self.restart()

        return results


class ThermalRelaxWork(Work):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ThermalRelaxWork
    """

    @classmethod
    def from_relax_input(cls, relax_input, qha, temperatures, pressures):
        """
        """
        work = cls()

        work.temperatures = temperatures
        work.pressures = pressures
        #work.qha = qha

        for pressure in pressures:
            for temperature in temperatures:
                strtarget = qha.get_strtarget(temperature, pressure)
                new_input = relax_input.new_with_vars(strtarget=strtarget)
                work.register_relax_task(new_input)

        return work

    def on_all_ok(self):
        """
        Implement the post-processing step at the end of the Work.
        """
        return super().on_all_ok()

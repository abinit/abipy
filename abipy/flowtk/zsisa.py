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


class ZsisaFlow(Flow):
    """
    Flow for QHA calculations with the VZSISA approach.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaFlow
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, eps, ngqpt, with_becs, with_quad,
                       edos_ngkpt=None, manager=None) -> ZsisaFlow:
        """
        Build a flow for QHA calculations from an |AbinitInput| for GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            eps:
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
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

        # optcell = 2: full optimization of cell geometry (
        work = ZsisaWork.from_scf_input(scf_input, eps, ngqpt, with_becs, with_quad,
                                        ionmov=2, edos_ngkpt=edos_ngkpt)

        flow.register_work(work)
        return flow

    def finalize(self):
        """
        This method is called when the flow is completed.
        It performs some basic post-processing of the results to facilitate further analysis.
        """
        work = self[0]
        data = {"eps": work.eps}

        ## Build list of strings with path to the relevant output files ordered by V.
        #data["gsr_relax_paths"] = [task.gsr_path for task in work.relax_tasks_vol]

        #entries, gsr_relax_volumes = [], []
        #for task in work.relax_tasks_vol:
        #    with task.open_gsr() as gsr:
        #        entries.append(dict(
        #            volume=gsr.structure.volume,
        #            energy_eV=float(gsr.energy),
        #            pressure_GPa=float(gsr.pressure),
        #            #structure=gsr.structure,
        #        ))
        #        gsr_relax_volumes.append(gsr.structure.volume)

        #data["gsr_entries"] = entries
        #data["gsr_relax_volumes"] = gsr_relax_volumes

        #data["ddb_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in work.ph_works]
        #data["ddb_relax_volumes"] = [ph_work[0].input.structure.volume for ph_work in work.ph_works]

        data["gsr_edos_path"] = [] if not work.edos_work else [task.gsr_path for task in work.edos_work]

        # Write json file
        mjson_write(data, self.outdir.path_in("zsisa.json"), indent=4)

        return super().finalize()


class ZsisaWork(Work):
    """
    This work performs a structural relaxation of the initial structure, then a set of distorted
    structures is genenerated and the relaxed structures are used
    to compute phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaWork
    """

    @classmethod
    def from_scf_input(cls, scf_input, eps, ngqpt,
                       with_becs: bool, with_quad: bool,
                       ionmov: int, edos_ngkpt=None) -> ZsisaWork:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            with_becs: Activate calculation of Electric field and Born effective charges.
            ionmov: Abinit input variables.
            edos_ngkpt: Three integers defining the the k-sampling for the computation of the
                electron DOS with the relaxed structures. Useful for metals
                in which the electronic contribution should be included.
                None disables the computation of the e-DOS.
        """
        work = cls()

        # Save attributes in work
        work.initial_scf_input = scf_input
        work.eps = float(eps)
        work.ngqpt = ngqpt
        work.with_becs = with_becs
        work.with_quad = with_quad
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))

        # Create input for relaxation and register the relaxation task.
        work.relax_template = relax_template = scf_input.deepcopy()
        relax_template.pop_tolerances()
        # optcell = 2: full optimization of cell geometry
        relax_template.set_vars(optcell=2, ionmov=ionmov, tolvrs=1e-8, toldff=1.e-6)
        #if optcell is not None and optcell != 0:
        #relax_template.set_vars_ifnotin(ecutsm=1.0, dilatmx=1.05)

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
            self.deformed_structures_dict = generate_deformations(relaxed_structure, eps=self.eps)

            relax_template = self.relax_template
            self.relax_tasks_vol = []
            for structure in self.deformed_structures_dict.values():
                # Relax deformed structure at fixed unit cell.
                new_input = relax_template.new_with_structure(structure, optcell=0)
                task = self.register_relax_task(new_input)
                self.relax_tasks_vol.append(task)

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """
        This callback is called when all tasks in the Work reach status `S_OK`.
        Here we add a new PhononWork for each volume using the relaxed structure.
        """
        self.edos_work = Work()

        # Build phonon works for the different relaxed structures.
        self.ph_works = []

        start = 1
        for task, deform_name in zip(self[start:], self.deformed_structures_dict.keys(), strict=True):
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)
            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, ddk_tolerance=None)

            # Reduce the number of files produced in the DFPT tasks to avoid possible disk quota issues.
            prtvars = dict(prtden=0, prtpot=0)
            for task in ph_work[1:]:
                task.input.set_vars(**prtvars)

            ph_work.set_name(deform_name)
            self.ph_works.append(ph_work)
            self.flow.register_work(ph_work)

            # Add task for electron DOS calculation to edos_work
            if self.edos_ngkpt is not None:
                edos_input = scf_input.make_edos_input(self.edos_ngkpt)
                t = self.edos_work.register_nscf_task(edos_input, deps={ph_work[0]: "DEN"})
                t.set_name(deform_name)

        if self.edos_ngkpt is not None:
            self.flow.register_work(self.edos_work)

        self.flow.allocate(build=True)

        return super().on_all_ok()



class ThermalRelaxWork(Work):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ThermalRelaxWork
    """

    @classmethod
    def from_relax_input(cls, relax_input, zsisa, temperatures, pressures):
        """
        """
        work = cls()

        work.temperatures = np.array(temperatures)
        work.pressures_gpa = np.array(pressures_gpa)
        work.zsisa = zsisa

        for pressure_gpa in work.pressures_gpa:
            for temperature in work.temperatures:
                strtarget = zsisa.get_strtarget(temperature, pressure)
                new_input = relax_input.new_with_vars(strtarget=strtarget)
                task = work.register_task(new_input, task_class=ThermalRelaxTask)
                # Attach pressure and temperature
                task.pressure_gpa = pressure_gpa
                task.temperature = temperature

        return work

    def on_all_ok(self):
        """
        Implement the post-processing step at the end of the Work.
        """
        return super().on_all_ok()


class ThermalRelaxTask(RelaxTask):

    def _on_ok(self):
        results = super()._on_ok()

        zsisa = self.work.zsisa

        # Check for convergence.
        #if not self.collinear_done:
        #    self.input.set_vars(strtarget=strtarget)
        #    self.finalized = False
        #    self.restart()

        return results



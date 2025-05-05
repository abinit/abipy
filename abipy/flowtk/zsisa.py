# coding: utf-8
"""
Workflows for calculations within the quasi-harmonic approximation.
"""
from __future__ import annotations

import numpy as np
import abipy.core.abinit_units as abu

from abipy.tools.serialization import mjson_write, mjson_load
from abipy.dfpt.deformation_utils import generate_deformations
from abipy.abio.inputs import AbinitInput
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.flows import Flow


class ZsisaFlow(Flow):
    """
    Flow for QHA calculations with the ZSISA approach.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaFlow
    """

    @classmethod
    def from_scf_input(cls,
                       workdir: PathLike,
                       scf_input: AbinitInput,
                       eps: float,
                       ngqpt,
                       with_becs: bool,
                       with_quad: bool,
                       ndivsm=-20,
                       edos_ngkpt=None,
                       manager=None) -> ZsisaFlow:
        """
        Build a flow for QHA calculations from an |AbinitInput| for GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            eps:
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
        flow.register_work(ZsisaWork.from_scf_input(scf_input, eps, ngqpt, with_becs, with_quad,
                                                    ndivsm, ionmov=2, edos_ngkpt=edos_ngkpt))
        return flow

    def finalize(self):
        """
        This method is called when the flow is completed.
        It performs some basic post-processing of the results to facilitate further analysis.
        """
        work = self[0]
        data = {"eps": work.eps, "spgrp_number": work.spgrp_number}

        # Build list of strings with path to the relevant output files ordered by V.
        data["gsr_relax_paths"] = [task.gsr_path for task in work.relax_tasks_strained]
        data["strain_inds"] = work.strain_inds

        gsr_relax_entries, gsr_relax_volumes = [], []
        for task in work.relax_tasks_strained:
            with task.open_gsr() as gsr:
                gsr_relax_entries.append(dict(
                    volume=gsr.structure.volume,
                    energy_eV=float(gsr.energy),
                    pressure_GPa=float(gsr.pressure),
                    #structure=gsr.structure,
                ))
                gsr_relax_volumes.append(gsr.structure.volume)

        data["gsr_relax_entries"] = gsr_relax_entries

        data["ddb_relax_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in work.ph_works]
        data["gsr_relax_edos_paths"] = [] if not work.edos_work else [task.gsr_path for task in work.edos_work]
        data["gsr_relax_ebands_paths"] = [] if work.ndivsm == 0 else [ph_work.ebands_task.gsr_path for ph_work in work.ph_works]

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
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       eps: float,
                       ngqpt,
                       with_becs: bool,
                       with_quad: bool,
                       ndivsm: int,
                       ionmov: int,
                       edos_ngkpt=None) -> ZsisaWork:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.
        See ZsisaFlow for the meaning of the arguments.
        """
        work = cls()

        # Save attributes in work
        work.initial_scf_input = scf_input
        work.eps = float(eps)
        work.ngqpt = ngqpt
        work.with_becs = with_becs
        work.with_quad = with_quad
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))
        work.ndivsm = ndivsm

        # Create input for relaxation and register the relaxation task.
        work.relax_template = relax_template = scf_input.deepcopy()

        # optcell = 2: full optimization of cell geometry
        relax_template.pop_tolerances()
        relax_template.set_vars(optcell=2, ionmov=ionmov, tolvrs=1e-8, tolmxf=1e-6)
        relax_template.set_vars_ifnotin(ecutsm=1.0, dilatmx=1.05)

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

            self.strained_structures_dict, self.strain_inds, self.spgrp_number = generate_deformations(relaxed_structure, self.eps)

            self.relax_tasks_strained = []
            for structure in self.strained_structures_dict.values():
                # Relax deformed structure with fixed unit cell.
                task = self.register_relax_task(self.relax_template.new_with_structure(structure, optcell=0))
                self.relax_tasks_strained.append(task)

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

        for task, strain_name, strain_ind in zip(self[1:], self.strained_structures_dict.keys(), self.strain_inds, strict=True):
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)
            #scf_input.pop_vars(["dilatmx"])
            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, with_quad=self.with_quad,
                                                ndivsm=0 if np.any(strain_ind != 0) else self.ndivsm)

            # Reduce the number of files produced in the DFPT tasks to avoid possible disk quota issues.
            for task in ph_work[1:]:
                task.input.set_vars(prtden=0, prtpot=0)

            ph_work.set_name(strain_name)
            self.ph_works.append(ph_work)
            self.flow.register_work(ph_work)

            # Add task for electron DOS calculation to edos_work
            if self.edos_ngkpt is not None:
                edos_input = scf_input.make_edos_input(self.edos_ngkpt)
                self.edos_work.register_nscf_task(edos_input, deps={ph_work[0]: "DEN"}).set_name(strain_name)

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
    def from_relax_input(cls,
                         relax_input: AbinitInput,
                         zsisa,
                         temperatures: list,
                         pressures: list) -> ThermalRelaxWork:
        """
        Args:
            relax_input:
            zsisa:
            temperatures:
            pressures
        """
        work = cls()

        work.temperatures = np.array(temperatures)
        work.pressures_gpa = np.array(pressures_gpa)
        work.zsisa = zsisa

        for pressure_gpa in work.pressures_gpa:
            for temperature in work.temperatures:
                new_input = relax_input.new_with_vars(strtarget=strtarget)
                task = work.register_task(new_input, task_class=ThermalRelaxTask)
                # Attach pressure and temperature to task.
                task.pressure_gpa = pressure_gpa
                task.temperature = temperature

        return work

    def on_all_ok(self):
        """
        Implement the post-processing step at the end of the Work.
        """
        # Build work for elastic properties (clamped-ions)
        # activate internal strain and piezoelectric part.
        #from abipy.flowtk.dfpt import ElasticWork
        #elastic_work = ElasticWork.from_scf_input(scf_input, with_relaxed_ion=True, with_piezo=True)
        return super().on_all_ok()


class ThermalRelaxTask(RelaxTask):

    def _on_ok(self):
        results = super()._on_ok()
        zsisa = self.work.zsisa
        tol_gpa = 0.01

        json_path = self.outdir.path_in("thermal_history.json")
        if os.path.exists(json_path):
            thermal_hist = mjson_load(json_path)
        else:
            thermal_hist = {"initial_structure": self.input.structure, "tol_gpa": tol_gpa, "history": []}

        with self.open_gsr() as gsr:
            relaxed_structure = gsr.structure
            # Stress tensor is in GPa units
            cart_therm_stress = zsisa.get_cart_thermal_stress(relaxed_structure, self.temperature, self.pressure_gpa)
            converged = np.all(np.abs(cart_therm_stress - gsr.cart_stress_tensor)) < tol_gpa

            thermal_hist["history"].append(dict(
                structure=relaxed_structure,
                cart_therm_stress=cart_therm_stress,
                cart_bo_stress=gsr.cart_stress_tensor,
                converged=converged,
            ))

        mjson_write(thermal_hist, json_path, indent=4)

        # Check for convergence.
        if not converged:
            # In fortran notation The components of the stress tensor must be stored according to:
            # (1,1) → 1; (2,2) → 2; (3,3) → 3; (2,3) → 4; (3,1) → 5; (1,2) → 6
            # TODO: strtarget refers to Cartesian coords I suppose! Also, check sign!
            strtarget = np.empty(6)
            strtarget[0] = cart_therm_stress[0,0]
            strtarget[1] = cart_therm_stress[1,1]
            strtarget[2] = cart_therm_stress[2,2]
            strtarget[3] = cart_therm_stress[1,2]
            strtarget[4] = cart_therm_stress[2,0]
            strtarget[5] = cart_therm_stress[0,1]
            strtarget /= abu.HaBohr3_GPa
            self.input.set_vars(strtarget=strtarget)
            self.finalized = False
            self.restart()

        return results

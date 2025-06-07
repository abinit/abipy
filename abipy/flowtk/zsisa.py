# coding: utf-8
"""
Workflows for calculations within the ZSISA approximation to the QHA.
"""
from __future__ import annotations

import itertools
import numpy as np
import abipy.core.abinit_units as abu

from abipy.tools.serialization import mjson_write, mjson_load
from abipy.tools.typing import PathLike
from abipy.abio.inputs import AbinitInput
from abipy.dfpt.deformation_utils import generate_deformations
from abipy.dfpt.qha_general_stress import QHA_ZSISA
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.flows import Flow
from abipy.flowtk.dfpt_works import ElasticWork


class ZsisaFlow(Flow):
    """
    Flow for QHA calculations with the ZSISA approximation.
    This is the main entry point for client code.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaFlow
    """

    @classmethod
    def from_scf_input(cls,
                       workdir: PathLike,
                       scf_input: AbinitInput,
                       eps: float,
                       mode: str,
                       ngqpt,
                       with_becs: bool,
                       with_quad: bool,
                       ndivsm: int = -20,
                       edos_ngkpt=None,
                       ionmov: int = 2,
                       tolmxf=1e-5,
                       manager=None) -> ZsisaFlow:
        """
        Build a flow for QHA calculations from an |AbinitInput| for GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            eps: Strain magnitude to be applied to the lattice. Strains will be applied to the reference lattice.
            mode: 'TEC' for thermal expansion coefficients.
                'ECs" for thermal expansion coefficients + elastic constants.
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
            ionmov: Algorithm for ion optimization.
            tolmxf: Tolerance of Max force.
            manager: |TaskManager| instance. Use default if None.
        """
        flow = cls(workdir=workdir, manager=manager)
        flow.register_work(ZsisaWork.from_scf_input(scf_input, eps, mode, ngqpt, with_becs, with_quad,
                                                    ndivsm, ionmov, tolmxf,
                                                    edos_ngkpt=edos_ngkpt))
        return flow

    def finalize(self):
        """
        This method is called when the flow is completed.
        It performs some basic post-processing of the results to facilitate further analysis.
        """
        work = self[0]
        data = {"eps": work.eps, "mode": work.mode, "spgrp_number": work.spgrp_number}
        data["gsr_bo_path"] = work.initial_relax_task.gsr_path

        # Build list of strings with paths to the relevant output files.
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
        data["gsr_relax_ebands_paths"] = [] if work.ndivsm == 0 else \
            [ph_work.ebands_task.gsr_path for ph_work in work.ph_works if ph_work.ebands_task is not None]

        # Write json file
        mjson_write(data, self.outdir.path_in("zsisa.json"), indent=4)

        return super().finalize()


class ZsisaWork(Work):
    """
    This work performs a structural relaxation of the initial structure,
    then a set of distorted structures is genenerated and the relaxed structures are used to compute
    phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaWork
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       eps: float,
                       mode: str,
                       ngqpt,
                       with_becs: bool,
                       with_quad: bool,
                       ndivsm: int = -20,
                       ionmov: int = 2,
                       tolmxf=1e-5,
                       edos_ngkpt=None) -> ZsisaWork:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.
        See ZsisaFlow for the meaning of the arguments.
        """
        work = cls()

        # Save attributes in work.
        work.initial_scf_input = scf_input
        work.eps = float(eps)
        work.mode = mode
        work.ngqpt = ngqpt
        work.with_becs = with_becs
        work.with_quad = with_quad
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))
        work.ndivsm = ndivsm

        # Create input for relaxation and register the initial relaxation task.
        if "tolvrs" not in scf_input:
            raise ValueError("tolvrs should be specified in scf_input.")

        work.relax_template = relax_template = scf_input.deepcopy()

        # optcell = 2 --> full optimization of cell geometry.
        relax_kwargs = dict(optcell=2, ionmov=ionmov, tolmxf=tolmxf)
        relax_template.set_vars(**relax_kwargs)
        relax_template.set_vars_ifnotin(ecutsm=1.0, dilatmx=1.05)

        work.initial_relax_task = work.register_relax_task(relax_template)

        return work

    def on_ok(self, sender):
        """
        This method is called when one task reaches status `S_OK`.
        """
        if sender == self.initial_relax_task:
            # Get relaxed structure and build new task for structural relaxation at fixed volume.
            relaxed_structure = sender.get_final_structure()

            self.strained_structures_dict, self.strain_inds, self.spgrp_number = generate_deformations(
                relaxed_structure, self.eps, mode=self.mode)

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
        # Build phonon works for the different relaxed structures.
        self.ph_works = []
        self.edos_work = Work()

        for task, strain_name, strain_ind in zip(self[1:], self.strained_structures_dict.keys(), self.strain_inds, strict=True):
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)
            #scf_input.pop_vars(["dilatmx"])
            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, with_quad=self.with_quad,
                                                ndivsm=0 if np.any(strain_ind != 0) else self.ndivsm)

            # Reduce the number of files produced by the DFPT tasks to avoid possible disk quota issues.
            for ph_task in ph_work[1:]:
                ph_task.input.set_vars(prtden=0, prtpot=0)

            ph_work.set_name(strain_name)
            self.ph_works.append(ph_work)
            self.flow.register_work(ph_work)

            # Add task for electron DOS calculation to edos_work.
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

    #@classmethod
    #def from_zsisa_flow(self, zsisa_flow,
    #                    temperatures: list,
    #                    pressures: list) -> ThermalRelaxWork:

    @classmethod
    def from_relax_input(cls,
                         relax_input: AbinitInput,
                         mode: str,
                         phdos_paths,
                         gsr_bo_path
                         temperatures: list,
                         pressures: list) -> ThermalRelaxWork:
        """
        Args:
            relax_input:
            mode:
            phdos_paths:
            gsr_bo_path:
            temperatures: List of temperatures in K.
            pressures: List of pressures in GPa.
        """
        work = cls()

        work.mode = mode
        work.phdos_paths = phdos_paths
        work.gsr_bo_path = gsr_bo_path

        work.temperatures = np.array(temperatures)
        work.pressures_gpa = np.array(pressures_gpa)

        #zsisa = self.get_zsisa(self, guess_path)
        work.thermal_relax_tasks = []
        for pressure_gpa, temperature in itertools.product(work.pressures_gpa, work.temperatures):
            #new_structure = zsisa.structure_guess
            #new_input = relax_input.new_with_structure(new_structure)
            new_input = relax_input.new_with_vars(strtarget=strtarget)

            # Register ThermalRelaxTask and attach pressure and temperature to the task.
            task = work.register_task(new_input, task_class=ThermalRelaxTask)
            task.mode = mode
            task.pressure_gpa = pressure_gpa
            task.temperature = temperature

            work.thermal_relax_tasks.append(task)

        return work

    def get_zsisa(self, guess_path) -> QHA_ZSISA:
        return QHA_ZSISA.from_files(self.phdos_paths, guess_path, self.gsr_bo_path)

    #def on_ok(self, sender):
    #    """
    #    This method is called when one task reaches status `S_OK`.
    #    """
    #    # Find sender in self.thermal_relax_tasks.
    #    for relax_task in self.thermal_relax_tasks:
    #        if sender.node_id == relax_task.node_id:
    #            break
    #    else:
    #        raise RuntimeError(f"Cannot find {sender=} in {self}")

    #    #if work.mode == "elastic":
    #    # Build work for elastic properties
    #    #relaxed_structure = relax_task.get_final_structure()
    #    #scf_input = relax_task.input.new_with_structure(relaxed_structure)
    #    #scf_input.pop_relax_vars()
    #    #elastic_work = ElasticWork.from_scf_input(scf_input, with_relaxed_ion=True, with_piezo=True)
    #    #self.flow.register_work(elastic_work)
    #    #self.flow.allocate(build=True)

    #    return super().on_all_ok()


class ThermalRelaxTask(RelaxTask):

    def _on_ok(self):
        results = super()._on_ok()

        # === Main loop to calculate stress and trigger elastic constants ===
        guess_path = self.gsr_path
        zsisa = self.work.get_zsisa(guess_path)

        elastic_path = self.outdir.path_in("elastic_constant.txt")
        chk_converge, stress = zsisa.cal_stress(self.temperature, self.pressure_gpa,
                                                mode=self.mode, elastic_path=elastic_path)
        print(stress)
        if not chk_converge:
            # Write new input for relaxation.
            extra_vars = {
                "strtarget": stress,
                #"ionmov": 2,
                #"ntime": 60,
                #"optcell": 2,
                #"dilatmx": 1.04,
                #"tolmxf": 1.0e-5,
                #"strfact": 1000.,
                #"prtden": 0,
                #"prtwf": 0,
                #"prteig": 0
            }
            self.input.set_vars(**extra_vars)
            self.finalized = False
            # Restart will take care of using the output structure in the input.
            self.restart(submit=False)

        #else:
        #    if self.mode == "ECs":
        #        # Build work for elastic properties
        #        relaxed_structure = self.get_final_structure()
        #        scf_input = self.input.new_with_structure(relaxed_structure)
        #        scf_input.pop_relax_vars()
        #        elastic_work = ElasticWork.from_scf_input(scf_input, with_relaxed_ion=True, with_piezo=True)
        #        self.flow.register_work(elastic_work)
        #        self.flow.allocate(build=True)

        return results

        # TODO Write files with results
        #tol_gpa = 0.01
        #json_path = self.outdir.path_in("thermal_history.json")
        #if os.path.exists(json_path):
        #    thermal_hist = mjson_load(json_path)
        #else:
        #    thermal_hist = {"initial_structure": self.input.structure, "tol_gpa": tol_gpa, "history": []}

        #with self.open_gsr() as gsr:
        #    relaxed_structure = gsr.structure
        #    # Stress tensor is in GPa units
        #    cart_therm_stress = zsisa.get_cart_thermal_stress(relaxed_structure, self.temperature, self.pressure_gpa)
        #    converged = np.all(np.abs(cart_therm_stress - gsr.cart_stress_tensor)) < tol_gpa

        #    #def cal_stress(self, temp, pressure = 0, mode = "TEC" , elastic_path = "elastic_constant.txt" ):

        #    thermal_hist["history"].append(dict(
        #        structure=relaxed_structure,
        #        cart_therm_stress=cart_therm_stress,
        #        cart_bo_stress=gsr.cart_stress_tensor,
        #        converged=converged,
        #    ))

        #mjson_write(thermal_hist, json_path, indent=4)

        ## Check for convergence.
        #if not converged:
        #    # In fortran notation The components of the stress tensor must be stored according to:
        #    # (1,1) → 1; (2,2) → 2; (3,3) → 3; (2,3) → 4; (3,1) → 5; (1,2) → 6
        #    # TODO: strtarget refers to Cartesian coords I suppose! Also, check sign!
        #    strtarget = np.empty(6)
        #    strtarget[0] = cart_therm_stress[0,0]
        #    strtarget[1] = cart_therm_stress[1,1]
        #    strtarget[2] = cart_therm_stress[2,2]
        #    strtarget[3] = cart_therm_stress[1,2]
        #    strtarget[4] = cart_therm_stress[2,0]
        #    strtarget[5] = cart_therm_stress[0,1]
        #    strtarget /= abu.HaBohr3_GPa
        #    self.input.set_vars(strtarget=strtarget)
        #    self.finalized = False
        #    self.restart()


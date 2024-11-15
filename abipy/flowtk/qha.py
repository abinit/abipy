# coding: utf-8
"""
Workflows for calculations within the quasi-harmonic approximation.
"""
from __future__ import annotations

import numpy as np

from abipy.tools.serialization import mjson_write # mjson_load, mjson_loads,
from abipy.dfpt.deformation_utils import generate_deformations
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.flows import Flow


class RelaxAndAddPhononWorks(Work):
    """
    This work performs a structural relaxation of the initial structure, then a set of distorted
    structures is genenerated and the relaxed structures are used
    to compute phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: RelaxAndAddPhononWorks
    """

    @classmethod
    def from_scf_input(cls, scf_input, eps, ngqpt, with_becs: bool, optcell: int, ionmov: int,
                       edos_ngkpt=None) -> RelaxAndAddPhononWorks:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
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
        work.eps = eps
        work.ngqpt = ngqpt
        work.with_becs = with_becs
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
            # Get relaxed structure and build deformed structures.
            relaxed_structure = sender.get_final_structure()
            self.deformed_structures_dict = generate_deformations(relaxed_structure, eps=self.eps)

            # Add new tasks to relax the deformed structures at fixed cell.
            for structure in self.deformed_structures_dict.values():
                new_input = self.relax_template.new_with_structure(structure, optcell=0)
                self.register_relax_task(new_input)

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """
        This callback is called when all tasks in the Work reach status `S_OK`.
        Here we add a new PhononWork for each volume using the relaxed structure.
        """
        if self.edos_ngkpt is not None: self.edos_work = Work()

        # Build phonon works for the different relaxed structures.
        self.ph_works = []
        start = 1
        for task, def_name in zip(self[start:], self.deformed_structures_dict.keys(), strict=True):
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)
            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, ddk_tolerance=None)

            ph_work.set_name(def_name)
            self.ph_works.append(ph_work)
            self.flow.register_work(ph_work)

            # Add electron DOS calculation.
            if self.edos_ngkpt is not None:
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


class QhaFlow(Flow):
    """
    Flow for QHA calculations.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: QhaFlow
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, ngqpt, with_becs, eps=0.005,
                       edos_ngkpt=None, manager=None) -> QhaFlow:
        """
        Build a flow for QHA calculations from an |AbinitInput| for GS-SCF calculations.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            ngqpt: Three integers defining the q-mesh for phonon calculation.
            with_becs: Activate calculation of Electric field and Born effective charges.
            eps:
            edos_ngkpt: Three integers defining the the k-sampling for the computation of the
                electron DOS with the relaxed structures. Useful for metals or small gap semiconductors
                in which the electronic contribution should be included.
                None disables the computation of the e-DOS.
            manager: |TaskManager| instance. Use default if None.
        """
        flow = cls(workdir=workdir, manager=manager)

        new_work = RelaxAndAddPhononWorks.from_scf_input(scf_input, eps, ngqpt,
                                                         with_becs=with_becs, optcell=2, ionmov=2,
                                                         edos_ngkpt=edos_ngkpt)
        flow.register_work(new_work)
        flow.relax_and_add_phonon_works = new_work

        return flow

    #def on_all_ok(self):
    #    if self.on_all_ok_num_calls > 0: return True
    #    self.on_all_ok_num_calls += 1

    #    print(f"In on_all_ok with {self.on_all_ok_num_calls =}")

    #    # implement_logic_to_create_new_work
    #    #for work in self.relax_and_add_phonon_works.ph_works

    #    """
    #    self.register_work(work)
    #    self.allocate()
    #    self.build_and_pickle_dump()
    #    """

    #    #from abipy.dfpt.qha_approximation import QHA_App
    #    #qha = QHA_App.from_files_app_ddb(cls, ddb_paths, phdos_paths)

    #    # The scheduler will keep on running the flow.
    #    return False

    def finalize(self):
        """
        This method is called when the flow is completed.
        It performs some basic post-processing of the results to facilitate further analysis.
        """
        data = {}
        #work = self[0]
        #data = dict(
        #    initial_structure=work.initial_structure,
        #    relaxed_structure=work.final_structure,
        #    relaxed_structure_energy_eV=work.initial_structure,
        #    initial_deformed_structures=work.initial_structure,
        #    relaxed_deformed_structures_energy_ev=work.initial_structure,
        #)

        # QHA.from_files(cls, gsr_paths, phdos_paths)
        # Build list of strings with path to the relevant output files ordered by V.
        #data["gsr_relax_paths"] = [task.gsr_path for task in self.relax_and_phonon_work.relax_tasks_vol]
        data["ddb_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in self.relax_and_add_phonon_works.ph_works]
        #data["gsr_edos_path"] = []
        #if self.relax_and_phonon_work.edos_ngkpt is not None:
        #    data["gsr_edos_paths"] = [task.gsr_path for task in self.relax_and_phonon_work.edos_work]

        #entries = []
        #items = zip(self.relax_and_phonon_work.relax_tasks_vol, self.relax_and_phonon_work.ph_works_vol)
        #for ivol, (relax_task, ph_work) in enumerate(items):
        #    ddb_path = ph_work.outdir.has_abiext("DDB")
        #    with relax_task.open_gsr() as gsr:
        #        entry = dict(
        #            volume=gsr.structure.volume,
        #            energy_eV=gsr.energy,
        #            pressure_GPa=gsr.pressure,
        #            structure=gsr.structure,
        #            gsr_edos_path=None,
        #            ddb_path=ddb_path,
        #        )

        #    if self.relax_and_phonon_work.edos_ngkpt is not None:
        #        task = self.relax_and_phonon_work.edos_work[ivol]
        #        entry["gsr_edos_path"] = task.gsr_path

        #    entries.append(entry)
        #data["entries_for_each_volume"] = entries

        mjson_write(data, self.outdir.path_in("qha_data.json"), indent=4)

        return super().finalize()

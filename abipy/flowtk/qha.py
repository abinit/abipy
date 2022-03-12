# coding: utf-8
"""
Workflows for calculations within the quasi-harmonic approximation.
"""
import numpy as np

from monty.json import jsanitize
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.flows import Flow


class RelaxAndPhononWork(Work):
    """
    This work performs a structural relaxation for different volumes, then it uses
    the relaxed structures to compute phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: RelaxAndPhononWork
    """

    @classmethod
    def from_scf_input(cls, scf_input, volumes, ngqpt, with_becs, optcell, ionmov, edos_ngkpt=None):
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            volumes: List of volumes in Ang**3
            with_becs: Activate calculation of Electric field and Born effective charges.
            optcell, ionmov: Abinit input variables.
            edos_ngkpt: Three integers defining the the k-sampling for the computation of the
                electron DOS with the relaxed structures. Useful for metals or small gap semiconductors
                in which the electronic contribution should be included.
                None disables the computation of the e-DOS.
        """
        work = cls()
        work.initial_scf_input = scf_input
        work.ngqpt = np.reshape(ngqpt, (3,))
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))
        work.volumes = np.array(volumes)
        work.with_becs = with_becs

        # Create input for relaxation and register the relaxation tasks.
        relax_template = scf_input.deepcopy()
        relax_template.pop_tolerances()
        relax_template.set_vars(tolvrs=1e-10, toldff=1.e-6, optcell=optcell, ionmov=ionmov)
        if optcell is not None and optcell != 0:
            relax_template.set_vars_ifnotin(ecutsm=0.5, dilatmx=1.05)

        # Construct len(volumes) works. Each work performs a structural relaxation
        # at fixed volume followed by a DFPT calculation with the relaxed structure.
        work.relax_tasks_vol = []
        for new_volume in work.volumes:
            new_structure = relax_template.structure.scale_lattice(new_volume)
            new_input = relax_template.new_with_structure(new_structure)
            task = work.register_relax_task(new_input)
            work.relax_tasks_vol.append(task)

        work.ph_works_vol = len(volumes) * [None]
        work.edos_work = None
        return work

    def on_all_ok(self):
        """
        This callback is called when all tasks in the Work reach status `S_OK`.
        Here are add a new PhononWork for each volume using the relaxed structure.
        """
        if self.edos_ngkpt is not None: self.edos_work = Work()

        # Build phonon works for the different volumes.
        for ivol, task in enumerate(self.relax_tasks_vol):
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)
            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, ddk_tolerance=None)
            self.ph_works_vol[ivol] = ph_work

            # Add e-DOS calculation.
            if self.edos_ngkpt is not None:
                edos_input = scf_input.make_edos_input(self.edos_ngkpt)
                self.edos_work.register_nscf_task(edos_input, deps={ph_work[0]: "DEN"})

            self.flow.register_work(ph_work)

        if self.edos_ngkpt is not None: self.flow.register_work(self.edos_work)
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

        return super().on_all_ok()


class QhaFlow(Flow):
    """

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: QhaFlow
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, volumes, ngqpt, with_becs,
                       edos_ngkpt=None, metadata=None, manager=None):
        """
        Build a |Flow| for QHA calculations from an |AbinitInput| for GS-SCF calculations.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            volumes: List of volumes in Ang**3.
            ngqpt: Three integers defining the q-mesh for phonon calculation.
            with_becs: Activate calculation of Electric field and Born effective charges.
            edos_ngkpt: Three integers defining the the k-sampling for the computation of the
                electron DOS with the relaxed structures. Useful for metals or small gap semiconductors
                in which the electronic contribution should be included.
                None disables the computation of the e-DOS.
            metadata: Dictionary with metadata to be be added to the final JSON file.
            manager: |TaskManager| instance. Use default if None.
        """
        ngqpt = np.reshape(ngqpt, 3)
        flow = cls(workdir=workdir, manager=manager)
        flow.metadata = jsanitize(metadata) if metadata is not None else None

        work = RelaxAndPhononWork.from_scf_input(scf_input, volumes, ngqpt,
                                                 with_becs=with_becs, optcell=3, ionmov=3, edos_ngkpt=edos_ngkpt)
        flow.relax_and_phonon_work = work
        flow.register_work(work)

        return flow

    def finalize(self):
        """
        This method is called when the flow is completed.
        It performs some basic post-processing of the results to facilitate further analysis.
        """
        d = {}
        if self.metadata is not None: d.update({"metadata": self.metadata})

        # QHA.from_files(cls, gsr_paths, phdos_paths)
        # Build list of strings with path to the relevant output files ordered by V.
        d["gsr_relax_paths"] = [task.gsr_path for task in self.relax_and_phonon_work.relax_tasks_vol]
        d["ddb_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in self.relax_and_phonon_work.ph_works_vol]
        d["gsr_edos_path"] = []
        if self.relax_and_phonon_work.edos_ngkpt is not None:
            d["gsr_edos_paths"] = [task.gsr_path for task in self.relax_and_phonon_work.edos_work]

        from abipy import abilab
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
        #d["entries_for_each_volume"] = entries

        abilab.mjson_write(d, self.outdir.path_in("qha.json"), indent=4)

        #pyscript = """
        #from abipy import abilab
        #from abipy.qha import Qha
        #data = abilab.mjson_read("qha.json")
        #data
        #qha = from_files(gsr_files_paths, phdos_files_paths):
        #"""

        return super().finalize()

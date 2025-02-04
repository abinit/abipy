# coding: utf-8
"""
Workflows for calculations within the quasi-harmonic approximation
with two degrees of freedom. The main entry point is Qha2dFlow.
"""
from __future__ import annotations

import itertools
import numpy as np

from abipy.tools.serialization import mjson_write
from abipy.dfpt.deformation_utils import generate_deformations
from abipy.abio.inputs import AbinitInput
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.flows import Flow


class Qha2dFlow(Flow):
    """
    Flow for QHA calculations with two degrees of freedom.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Qha2dFlow
    """

    @classmethod
    def from_scf_input(cls,
                       workdir: PathLike,
                       scf_input: AbinitInput,
                       bo_strains_ac: list[list],
                       phdos_strains_ac: list[list],
                       ngqpt,
                       with_becs: bool,
                       with_quad: bool,
                       ndivsm=-20,
                       edos_ngkpt=None,
                       manager=None) -> Qha2dFlow:
        """
        Build a flow for QHA calculations from an |AbinitInput| for GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            bo_strains_ac
            phdos_strains_ac
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
        flow.register_work(Qha2dWork.from_scf_input(scf_input, bo_strains_ac, phdos_strains_ac,
                                                    ngqpt, with_becs, with_quad,
                                                    ndivsm, ionmov=2, edos_ngkpt=edos_ngkpt))
        return flow

    def finalize(self):
        """
        This method is called when the flow is completed.
        It performs some basic post-processing of the results to facilitate further analysis.
        """
        work = self[0]
        data = {"bo_strains_ac": work.bo_strains_ac, "phdos_strains_ac": work.phdos_strains_ac}

        # Build list of strings with path to the relevant output files ordered by V.
        data["gsr_relax_paths"] = [task.gsr_path for task in work.relax_tasks_strained]

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
        mjson_write(data, self.outdir.path_in("qha_2d.json"), indent=4)

        return super().finalize()


class Qha2dWork(Work):
    """
    This work performs a structural relaxation of the initial structure, then a set of distorted
    structures is genenerated and the relaxed structures are used
    to compute phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Qha2dWork
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       bo_strains_ac,
                       phdos_strains_ac,
                       ngqpt,
                       with_becs: bool,
                       with_quad: bool,
                       ndivsm: int,
                       ionmov: int,
                       edos_ngkpt=None) -> Qha2dWork:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.
        See Qha2dFlow for the meaning of the arguments.
        """
        work = cls()

        # Save attributes in work
        work.initial_scf_input = scf_input

        # Make sure ench row is a numpy array.
        work.bo_strains_ac = bo_strains_ac
        work.phdos_strains_ac = phdos_strains_ac
        for i in range(2):
            work.bo_strains_ac[i] = np.array(bo_strains_ac[i])
            work.phdos_strains_ac[i] = np.array(phdos_strains_ac[i])

        work.ngqpt = ngqpt
        work.with_becs = with_becs
        work.with_quad = with_quad
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))
        work.ndivsm = ndivsm

        # Create input for relaxation and register the relaxation task.
        work.relax_template = relax_template = scf_input.deepcopy()

        # optcell = 3 --> constant-volume optimization of cell geometry.
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
            # Get relaxed structure and build new task for structural relaxation at fixed volume.
            relaxed_structure = sender.get_final_structure()

            self.relax_tasks_strained = []

            for s1, s3 in itertools.product(self.bo_strains_ac[0], self.bo_strains_ac[1]):
                strain_name = f"{s1=}, {s3=}"
                # Apply strain to the structure
                strain_tensor = np.diag([s1, s1, s3])
                strained_structure = relaxed_structure.apply_strain(strain_tensor, inplace=False)
                #print("strained_structure:", strained_structure)

                # Relax deformed structure with fixed unit cell.
                task = self.register_relax_task(self.relax_template.new_with_structure(strained_structure, optcell=0))

                task.bo_strain = np.array((s1, s3))
                task.in_phdos_strains = np.any(np.abs(s1 - self.phdos_strains_ac[0]) < 1e-3) and \
                                        np.any(np.abs(s3 - self.phdos_strains_ac[1]) < 1e-3)
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

        for task in self.relax_tasks_strained:
            s1, s3 = task.bo_strain
            strain_name = f"{s1=}, {s3=}"

            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)

            if task.in_phdos_strains:
                ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                    with_becs=self.with_becs, with_quad=self.with_quad,
                                                    ndivsm=0 if np.any(task.bo_strain != 0) else self.ndivsm)

                # Reduce the number of files produced in the DFPT tasks to avoid possible disk quota issues.
                for task in ph_work[1:]:
                    task.input.set_vars(prtden=0, prtpot=0)

                ph_work.set_name(strain_name)
                self.ph_works.append(ph_work)
                self.flow.register_work(ph_work)

            # Add task for electron DOS calculation to edos_work.
            if self.edos_ngkpt is not None:
                edos_input = scf_input.make_edos_input(self.edos_ngkpt)
                self.edos_work.register_nscf_task(edos_input, deps={ph_work[0]: "DEN"}).set_name(deformationname)

        if self.edos_ngkpt is not None:
            self.flow.register_work(self.edos_work)

        self.flow.allocate(build=True)

        return super().on_all_ok()

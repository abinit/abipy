# coding: utf-8
"""
Workflows for calculations within the ZSISA approximation to the QHA.
"""
from __future__ import annotations

import itertools
import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from abipy.tools.serialization import mjson_load, Serializable
from abipy.tools.typing import PathLike, Figure
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.abio.inputs import AbinitInput
from abipy.electrons import GsrFile
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.deformation_utils import generate_deformations
from abipy.dfpt.qha_general_stress import QHA_ZSISA
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.flows import Flow
from abipy.flowtk.dfpt_works import ElasticWork


@dataclasses.dataclass(kw_only=True)
class ZsisaResults(Serializable):
    eps: float
    mode: str
    spgrp_number: int
    gsr_bo_path: str
    strain_inds: np.ndarray
    gsr_relax_paths: list[str]
    gsr_relax_entries: list[dict]
    ddb_relax_paths: list[str]
    gsr_relax_edos_paths: list[str]
    gsr_relax_ebands_paths: list[str]


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
        Build a flow for ZSISA calculations from an |AbinitInput| representing a GS-SCF calculation.

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
        # Consistency check.
        if "ngkpt" in scf_input:
            ngkpt = np.array(scf_input["ngkpt"], dtype=int)
            ngqpt = np.array(ngqpt, dtype=int)

            if np.any(ngkpt % ngqpt != 0):
                raise ValueError(f"ngqpt should be a divisor of ngkpt but got {ngqpt=} and {ngkpt=}")

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
        data = {
            "eps": work.eps,
            "mode": work.mode,
            "spgrp_number": work.spgrp_number,
            "gsr_bo_path": work.initial_relax_task.gsr_path,
            "strain_inds": work.strain_inds,
        }

        # Build list of strings with paths to the relevant output files.
        data["gsr_relax_paths"] = [task.gsr_path for task in work.relax_tasks_strained]

        gsr_relax_entries = []
        for task in work.relax_tasks_strained:
            with task.open_gsr() as gsr:
                gsr_relax_entries.append(dict(
                    volume=gsr.structure.volume,
                    energy_eV=float(gsr.energy),
                    pressure_GPa=float(gsr.pressure),
                    #structure=gsr.structure,
                ))

        data["gsr_relax_entries"] = gsr_relax_entries

        data["ddb_relax_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in work.ph_works]
        data["gsr_relax_edos_paths"] = [] if not work.edos_work else [task.gsr_path for task in work.edos_work]
        data["gsr_relax_ebands_paths"] = [] if work.ndivsm == 0 else \
            [ph_work.ebands_task.gsr_path for ph_work in work.ph_works if ph_work.ebands_task is not None]

        # Write json file.
        ZsisaResults(**data).json_write(self.outdir.path_in("ZsisaResults.json"), indent=4)

        return super().finalize()


class ZsisaWork(Work):
    """
    This work performs a structural relaxation of the initial structure,
    then a set of distorted structures is generated and the relaxed structures are used to compute
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
                       #mode: str
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
        Here we take the relaxed structure from initial_relax_task and generate deformed
        structures that will be relaxed at fixed unit cell before running phonon calculations.
        """
        if sender == self.initial_relax_task:
            # Get relaxed structure and build new task for structural relaxation at fixed volume.
            relaxed_structure = sender.get_final_structure()

            # Generate deformed structures with the associated indices in the 6d matrix.
            self.strained_structures_dict, self.strain_inds, self.spgrp_number = generate_deformations(
                relaxed_structure, self.eps, mode=self.mode)

            # Relax each deformed structure with fixed unit cell (optcell 0).
            self.relax_tasks_strained = []
            for structure in self.strained_structures_dict.values():
                task = self.register_relax_task(self.relax_template.new_with_structure(structure, optcell=0))
                self.relax_tasks_strained.append(task)

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """
        This callback is called when all tasks in the Work reach status `S_OK`.
        Here we add a new PhononWork for each deformed structure that has been relaxed at fixed cell.
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
            #for ph_task in ph_work[1:]:
            #    ph_task.input.set_vars(prtpot=0)

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
    A work made of ThermalRelaxTask tasks.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ThermalRelaxWork
    """

    @classmethod
    def from_zsisa_flow(cls,
                        zsisa_flow_dir,
                        temperatures: list,
                        pressures_gpa: list,
                        verbose: int = 0) -> ThermalRelaxWork:
        """
        Args:
            relax_input:
            mode:
            phdos_paths:
            gsr_bo_path:
            temperatures: List of temperatures in K.
            pressures_gpa: List of pressures in GPa.
        """
        flow = Flow.as_flow(zsisa_flow_dir)

        json_filepath = flow.outdir.path_in("ZsisaResults.json")
        data = mjson_load(json_filepath)

        work = cls()

        work.temperatures = np.array(temperatures, dtype=float)
        work.pressures_gpa = np.array(pressures_gpa, dtype=float)
        work.mode = data.mode
        #work.phdos_paths = phdos_paths
        #work.gsr_bo_path = gsr_bo_path

        # Build zsisa object.
        nqsmall_or_qppa = 2
        work.zsisa = QHA_ZSISA.from_json_file(json_filepath, nqsmall_or_qppa, verbose=verbose)

        scf_input = flow[0].initial_scf_input
        relax_template = flow[0].relax_template

        # Generate initial ThermalRelaxTask tasks.
        work.thermal_relax_tasks = []
        for pressure_gpa, temperature in itertools.product(work.pressures_gpa, work.temperatures):
            converged, stress = work.zsisa.cal_stress(temperature, pressure_gpa,
                                                      mode=work.mode, elastic_path=None)

            #print(f"{temperature=}, {pressure_gpa=}, {stress=}")
            # TODO: Relax options with ecutsm and strfact?
            new_relax_input = relax_template.new_with_vars(strtarget=stress)

            # Attach pressure and temperature to the task.
            task = work.register_task(new_relax_input, task_class=ThermalRelaxTask)
            task.mode = work.mode
            task.pressure_gpa = pressure_gpa
            task.temperature = temperature
            task.elastic_work = None

            work.thermal_relax_tasks.append(task)

        return work

    #def on_ok(self, sender):
    #    """
    #    This method is called when one task reaches status `S_OK`.
    #    """
    #    # Find sender in self.thermal_relax_tasks.
    #    for relax_task in self.thermal_relax_tasks:
    #        if sender == relax_task:
    #            break
    #    else:
    #        raise RuntimeError(f"Cannot find {sender=} in {self}")

    def on_all_ok(self):
        """
        Callback triggered when all tasks in the ThermalRelaxWork reach the `S_OK` status.

        This method writes a JSON file containing the paths to the GSR files associated with each
        temperature (T) and pressure (P) point. If elastic constants have also been computed,
        a similar JSON entry is generated for the corresponding DDB files.
        """

        data = dict(
            #initial_structure=self.input.structure,
            mode=self.mode,
            #phdos_paths=self.phdos_paths,
            #gsr_bo_path=self.gsr_bo_path,
            #temperatures=self.temperatures,
            #pressures_gpa=self.pressures_gpa,
        )

        data["task_entries"] = []
        for task in self.thermal_relax_tasks:
            entry = dict(
                pressure_gpa=task.pressure_gpa,
                temperature=task.temperature,
                gsr_path=task.gsr_path,
                elastic_ddb_path=None,
            )
            # Add path to the DDB file with the 2nd order derivatives wrt strain.
            if task.elastic_work is not None:
                entry["elastic_ddb_path"] = task.elastic_work.outdir.path_in("out_DDB")

            data["task_entries"].append(entry)

        ThermalRelaxResults(**data).json_write(self.outdir.path_in("ThermalRelaxResults.json"), indent=4)

        return super().on_all_ok()


class ThermalRelaxTask(RelaxTask):
    """
    This task implements an iterative method to find the optimal lattice parameters
    at a given temperature T and external pressure P_ext.
    Starting from an initial lattice guess, the thermal and Born-Oppenheimer (BO) stresses are computed.
    A target stress is defined using the thermal stress and P_ext.
    Then, the lattice and atomic positions are relaxed repeatedly until the BO stress matches the target stress

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ThermalRelaxTask
    """

    def _on_ok(self):
        print("In _on_ok")
        results = super()._on_ok()

        # Get relaxed structure and stress_guess from GSR file.
        with self.open_gsr() as gsr:
            relaxed_structure = gsr.structure
            stress_guess = gsr.cart_stress_tensor / 29421.02648438959

        guess_path = self.gsr_path
        zsisa = self.work.zsisa
        zsisa.set_structure_stress_guess(relaxed_structure, stress_guess)

        #elastic_path = self.outdir.path_in("elastic_constant.txt")
        converged, stress = zsisa.cal_stress(self.temperature, self.pressure_gpa,
                                             mode=self.mode, elastic_path=None)
        print(f"{stress_guess=}")
        print(f"{stress=}")
        print(f"{converged=}")

        if not converged:
            # Change strtarget and restart the task but without submitting it.
            extra_vars = {
                "strtarget": stress,
                "ionmov": 2,
                "ntime": 100,
                "optcell": 2,
                "dilatmx": 1.04,
                "tolmxf": 1.0e-5,
                "strfact": 1000.,
            }
            self.input.set_vars(**extra_vars)
            self.finalized = False
            # Restart will take care of using the output structure in the input.
            self.restart()

        else:
            if self.mode == "ECs":
                # Build work for elastic properties and attach it to the task.
                scf_input = self.input.new_with_structure(relaxed_structure, ionmov=0, optcell=0)
                self.elastic_work = ElasticWork.from_scf_input(scf_input, with_relaxed_ion=True, with_piezo=True)
                self.flow.register_work(self.elastic_work)
                self.flow.allocate(build=True)

        return results


@dataclasses.dataclass(kw_only=True)
class ThermalRelaxResults(Serializable):
    """
    """
    mode: str
    #phdos_paths: list[str]
    #gsr_bo_path: str
    #temperatures: np.ndarray
    #pressures_gpa: np.ndarray
    task_entries: list[dict]

    def get_dataframe(self) -> pd.DataFrame:
        """
        Build pandas DataFrame with temperature, pressure, lattice parameters and spacegroup info.
        """
        rows = []
        for entry in self.task_entries:
            row = {k: entry[k] for k in ("temperature", "pressure_gpa")}
            with GsrFile(entry["gsr_path"]) as gsr:
                row.update(gsr.structure.get_dict4pandas())
            rows.append(row)

        return pd.DataFrame(rows).sort_values(by="temperature")

    #def get_elastic_dataframe(self, **kwargs) -> pd.DataFrame:
    #    """Build pandas DataFrame with temperatures, pressure and elastic constants."""
    #    rows = []
    #    for entry in self.task_entries:
    #        ddb_path = entry["elastic_ddb_path"]
    #        if ddb_path is None:
    #            raise ValueError("elastic_ddb_path is None --> elastic constants are not available!")

    #        raise NotImplementedError()
    #        row = {k: entry[k] for k in ("temperature", "pressure_gpa")}
    #        with DdbFile(ddb_path) as ddb:
    #            elastic_data = ddb.anaget_elastic(**kwargs)
    #            df = elastic_data.get_elastic_dataframe(tensor_name)
    #            #row.update(gsr.structure.get_dict4pandas())
    #        rows.append(row)

    #    return pd.DataFrame(rows).sort_values(by="temperature")

    @add_fig_kwargs
    def plot_lattice_vs_temp(self, **kwargs) -> Figure:
        angles = ["alpha", "beta", "gamma"]
        lengths = ["a", "b", "c"]
        volume = "volume"
        data = self.get_dataframe()
        #hue = "pressure"
        hue = None

        #ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
        #                                        sharex=False, sharey=True, squeeze=False)
        #ax_list = ax_list.ravel()

        ax, fig, plt = get_ax_fig_plt(ax=None)

        import seaborn as sns
        sns.lineplot(data=data, x="temperature", y="a", ax=ax, hue=hue,
                     markers=True, dashes=False)

        #for angle in angles:

        return fig

    #@add_fig_kwargs
    #def plot_elastic_vs_temp(self, **kwargs) -> Figure:
    #     df = self.get_elastic_dataframe()

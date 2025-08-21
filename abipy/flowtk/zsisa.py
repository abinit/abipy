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

from functools import cached_property
from abipy.tools.serialization import Serializable
from abipy.tools.typing import PathLike, VectorLike, Figure
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, plot_xy_with_hue, set_visible, set_grid_legend
from abipy.abio.inputs import AbinitInput
from abipy.electrons import GsrFile
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.deformation_utils import generate_deformations
from abipy.dfpt.qha_general_stress import QHA_ZSISA, spgnum_to_crystal_system, cmat_inds_names
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.works import Work, PhononWork
from abipy.flowtk.dfpt_works import ElasticWork
from abipy.flowtk.flows import Flow


class ZsisaFlow(Flow):
    """
    Flow for QHA calculations with the ZSISA approximation.
    This is the main entry point for client code.
    See <https://arxiv.org/pdf/2503.09531>.

    This AbiPy flow performs the following operations:

    1) The initial structure is relaxed and a minimum set of deformed configurations is generated.

    2) The atomic positions in the deformed structures are relaxed at fixed cell,
       and the relaxed configurations are then used to compute phonons, Born-effective charges,
       the electronic dielectric tensor (eps_inf), and dynamical, quadrupoles with DFPT.

    3) Phonon DOSes are computed for all the deformed structures by interpolating the interatomic
       forces constants on a much denser q-mesh, and the results are used to compute the thermal stresse.

    4) An iterative process is used for determining lattice parameters at temperature T and external pressure P.
       The process begins with an initial guess for the lattice configuration,
       followed by the computation of thermal and BO stresses.
       A target stress is defined based on thermal stress and P, and the lattice and atomic positions
       are relaxed iteratively until the target stress matches the BO stress and the BO forces are zeroed,
       ensuring convergence to the optimal lattice configuration.

    5) Finally, relaxed-ions elastic constants for the thermal-relaxed configurations
       are optionally computed with DFPT.

    The final results are stored in the ``outdata`` directory of the flow.
    Two files are produced at the end of the run:

    zsisa_data.csv:

        lattice parameters, thermal expansion, and elastic tensor components for each T and P in CSV format.

    ZsisaResults.json:

        JSON file with the location of the different files (GSR.nc, DDB) on the filesystem.
        To reconstruct a python object from file, use:


    .. code-block:: python

        from abipy.zsisa import ZsisaResults
        data = ZsisaResults.json_load("ABSPATH_TO_FILE")

        # Post-processing methods.
        data.get_dataframe()
        data.plot_lattice_vs_temp()
        data.plot_thermal_expansion()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaFlow
    """

    @classmethod
    def from_scf_input(cls,
                       workdir: PathLike,
                       scf_input: AbinitInput,
                       eps: float,
                       mode: str,
                       ngqpt: VectorLike,
                       with_becs: bool,
                       with_quad: bool,
                       temperatures: VectorLike,
                       pressures_gpa: VectorLike,
                       nqsmall_or_qppa: int,
                       with_elastic: bool = True,
                       ndivsm: int = 0,
                       edos_ngkpt: VectorLike | None = None,
                       ionmov: int = 2,
                       tolmxf: float = 1e-5,
                       qha_model: str = "zsisa",
                       with_piezo: bool = False,
                       manager=None) -> ZsisaFlow:
        """
        Build a flow for ZSISA calculations from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            eps: Strain magnitude to be applied to the lattice. Strains will be applied to the reference lattice.
            mode: "TEC" for thermal expansion coefficients.
                "ECs" for thermal expansion coefficients + elastic constants.
            ngqpt: Three integers defining the q-mesh for phonon calculation.
            with_becs: Activate calculation of Electric field and Born effective charges.
            with_quad: Activate calculation of dynamical quadrupoles. Require `with_becs`
                Note that only selected features are compatible with dynamical quadrupoles.
                Please consult <https://docs.abinit.org/topics/longwave/>
            temperatures: List of temperatures in K.
            pressures_gpa: List of pressures in GPa.
            nqsmall_or_qppa: Defines the q-mesh for the computation of the phonon DOS.
                if > 0, it is interpreted as nqsmall
                if < 0, it is interpreted as qppa.
            with_elastic: False to disable to computation of elastic constants. Only available if mode == "TEC".
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
            qha_model: Specifies the QHA model type. Possible values are:
                'zsisa': Standard ZSISA model.
                'v_zsisa': v-ZSISA model.
                'zsisa_slab': ZSISA model adapted for slab geometries.
            with_piezo: True to compute piezoelectric tensor.
            manager: |TaskManager| instance. Use default if None.
        """
        # Consistency check.
        if "ngkpt" in scf_input:
            ngkpt = np.array(scf_input["ngkpt"], dtype=int)
            ngqpt = np.array(ngqpt, dtype=int)

            if np.any(ngkpt % ngqpt != 0):
                raise ValueError(f"ngqpt should be a divisor of ngkpt but got {ngqpt=} and {ngkpt=}")

        flow = cls(workdir=workdir, manager=manager)

        # Store temperatures and pressures in flow.
        flow.temperatures = np.array(temperatures, dtype=float)
        flow.pressures_gpa = np.array(pressures_gpa, dtype=float)
        flow.qha_model = qha_model
        flow.nqsmall_or_qppa = nqsmall_or_qppa
        flow.with_piezo = with_piezo
        flow.with_elastic = with_elastic

        if not with_elastic and mode == "ECs":
            raise ValueError(f"with_elastic = False is not compatible with {mode=}")

        flow.register_work(ZsisaWork.from_scf_input(scf_input, eps, mode, ngqpt, with_becs, with_quad,
                                                    nqsmall_or_qppa, ndivsm, ionmov, tolmxf,
                                                    edos_ngkpt=edos_ngkpt))
        return flow

    def on_all_ok(self):
        """
        This method is called when all the works in the flow have reached S_OK.
        Here we write a json file with the paths to the DDB/GSR files produced so far
        then we create a new work for the self-consistent relaxation under thermal stress
        for each temperature and pressure.
        """
        json_filepath = self.outdir.path_in("ZsisaResults.json")
        def _path(basename: str):
            return self.outdir.path_in(basename)

        self.on_all_ok_num_calls += 1

        if self.on_all_ok_num_calls == 1:
            # Write json file with metadata and empty list of results.
            self._write_zsisa_results_step1(json_filepath)
            # Relaxation with thermal stress for the different (T, P) values.
            self.thermal_relax_work = ThermalRelaxWork.from_zsisa_flow(self, self.temperatures, self.pressures_gpa)
            self.register_work(self.thermal_relax_work)
            self.allocate(build=True)
            return False

        if self.on_all_ok_num_calls == 2:
            # Add results to the json file.
            self.thermal_relax_work.update_zsisa_results(json_filepath)

            # Produce csv file in the outdir of the flow.
            results = ZsisaResults.json_load(json_filepath)
            results.get_dataframe().to_csv(_path("zsisa_data.csv"))

        return True

    def _write_zsisa_results_step1(self, json_filepath: str) -> None:
        """
        Write json file with initial data.
        """
        work = self[0]
        data = {
            "spgrp_number": work.spgrp_number,
            "qha_model": work.flow.qha_model,
            "eps": work.eps,
            "mode": work.mode,
            "gsr_bo_path": work.initial_relax_task.gsr_path,
            "inds_6d": work.inds_6d,
        }

        # Build list of strings with paths to the relevant output files.
        data["gsr_relax_paths"] = [task.gsr_path for task in work.relax_tasks_strained]

        data["ddb_relax_paths"] = [ph_work.outdir.has_abiext("DDB") for ph_work in work.ph_works]
        data["gsr_relax_edos_paths"] = [] if not work.edos_work else [task.gsr_path for task in work.edos_work]
        data["gsr_relax_ebands_paths"] = [] if work.ndivsm == 0 else \
            [ph_work.ebands_task.gsr_path for ph_work in work.ph_works if ph_work.ebands_task is not None]

        # Init with empty list.
        data["thermal_relax_entries"] = []

        # Write json file.
        ZsisaResults(**data).json_write(json_filepath, indent=4)


class ZsisaWork(Work):
    """
    This work performs a structural relaxation of the initial structure,
    then a set of distorted structures is generated and the relaxed structures
    are used to compute phonons, BECS and the dielectric tensor with DFPT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaWork
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       eps: float,
                       mode: str,
                       ngqpt: VectorLike,
                       with_becs: bool,
                       with_quad: bool,
                       nqsmall_or_qppa: int,
                       ndivsm: int = 0,
                       ionmov: int = 2,
                       tolmxf: float = 1e-5,
                       edos_ngkpt: VectorLike | None = None,
                       manager=None) -> ZsisaWork:
        """
        Build the work from an |AbinitInput| representing a GS-SCF calculation.
        See ZsisaFlow for the meaning of the arguments.
        """
        work = cls(manager=manager)

        # Save attributes in work.
        work.initial_scf_input = scf_input
        work.eps = float(eps)
        work.mode = mode
        work.ngqpt = ngqpt
        work.with_becs = with_becs
        work.with_quad = with_quad
        work.edos_ngkpt = edos_ngkpt if edos_ngkpt is None else np.reshape(edos_ngkpt, (3,))
        work.nqsmall_or_qppa = nqsmall_or_qppa
        work.ndivsm = ndivsm

        # Create input for relaxation and register the initial relaxation task.
        if "tolvrs" not in scf_input:
            raise ValueError("tolvrs must be specified in scf_input.")

        work.relax_template = relax_template = scf_input.deepcopy()

        # optcell = 2 --> full optimization of cell geometry.
        relax_kwargs = dict(optcell=2, ionmov=ionmov, tolmxf=tolmxf)
        relax_template.set_vars(**relax_kwargs)
        relax_template.set_vars_ifnotin(ecutsm=1.0, dilatmx=1.05)

        work.initial_relax_task = work.register_multi_relax_task(relax_template)

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

            # Generate deformed structures with the associated indices in the phdos_6d matrix.
            self.strained_structures_dict, self.inds_6d, self.spgrp_number = generate_deformations(
                relaxed_structure, self.eps, mode=self.mode)

            # Relax each deformed structure with fixed unit cell (optcell 0).
            self.relax_tasks_strained = []
            for structure in self.strained_structures_dict.values():
                task = self.register_relax_task(self.relax_template.new_with_structure(structure, optcell=0, dilatmx=1.0))
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

        for task, strain_name, strain_ind in zip(self[1:], self.strained_structures_dict.keys(), self.inds_6d, strict=True):
            relaxed_structure = task.get_final_structure()
            scf_input = self.initial_scf_input.new_with_structure(relaxed_structure)
            ph_work = PhononWork.from_scf_input(scf_input, self.ngqpt, is_ngqpt=True, tolerance=None,
                                                with_becs=self.with_becs, with_quad=self.with_quad,
                                                ndivsm=0 if np.any(strain_ind != 0) else self.ndivsm)

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
                        zsisa_flow: Flow | PathLike,
                        temperatures: VectorLike,
                        pressures_gpa: VectorLike,
                        nqsmall_or_qppa: int | None = None,
                        verbose: int = 0) -> ThermalRelaxWork:
        """
        Args:
            zsisa_flow: instance of ZsisaFlow or directory hosting the Flow.
            temperatures: List of temperatures in K.
            pressures_gpa: List of pressures in GPa.
            nqsmall_or_qppa: Define the q-mesh for the computation of the PHDOS.
                if > 0, it is interpreted as nqsmall
                if < 0, it is interpreted as qppa.
            verbose: Verbosity level.
        """
        flow = Flow.as_flow(zsisa_flow)

        json_filepath = flow.outdir.path_in("ZsisaResults.json")
        data = ZsisaResults.json_load(json_filepath)

        work = cls()
        work.temperatures = np.array(temperatures, dtype=float)
        work.pressures_gpa = np.array(pressures_gpa, dtype=float)
        work.mode = data.mode

        # Build zsisa object.
        if nqsmall_or_qppa is None:
            nqsmall_or_qppa = flow[0].nqsmall_or_qppa

        print(f"Computing Phonon DOS using {nqsmall_or_qppa=} ...")
        work.zsisa = QHA_ZSISA.from_json_file(json_filepath, nqsmall_or_qppa, verbose=verbose)

        work0 = flow[0]

        scf_input = work0.initial_scf_input
        relax_template = work0.relax_template

        # Get relaxed structure and stress_guess and energy_guess from the GSR file.
        with work0.initial_relax_task.open_gsr() as gsr:
            structure_guess = gsr.structure
            stress_guess = gsr.cart_stress_tensor * abu.GPa_to_au
            energy_guess = gsr.energy

        # Generate initial ThermalRelaxTask tasks.
        work.thermal_relax_tasks = []
        for pressure_gpa, temperature in itertools.product(work.pressures_gpa, work.temperatures):
            tdata = work.zsisa.get_tstress(temperature, pressure_gpa, work.mode,
                                           structure_guess, stress_guess, energy_guess,
                                           bo_elastic_voigt=None)

            # TODO: Relax options with ecutsm and strfact?
            extra_vars = {
                "strtarget": tdata.stress_au,
                "ionmov": 2,
                "ntime": 100,
                "optcell": 2,
                "dilatmx": 1.04,
                "tolmxf": 1.0e-5,
                "strfact": 1000.,  # Give more weight to the stress in the relaxation.
            }

            new_relax_input = relax_template.new_with_vars(**extra_vars)

            # Attach pressure and temperature to the task.
            task = work.register_task(new_relax_input, task_class=ThermalRelaxTask)

            task.mode = work.mode
            task.pressure_gpa = pressure_gpa
            task.temperature = temperature
            task.elastic_work = None

            work.thermal_relax_tasks.append(task)

        return work

    def update_zsisa_results(self, json_filepath: str) -> None:
        """
        Compute elastic constants C at the end of the calculation.
        Use C to obtain to reduce noise in thermal expansion coefficients.
        Also, compute C(T, P) if self.mode == "ECs".

        Args:
            json_filepath: Path to json file.
        """
        results = ZsisaResults.json_load(json_filepath)
        zsisa = self.zsisa

        for task in self.thermal_relax_tasks:
            # Call anaddb to get elastic tensor from the DDB file.
            ddb_filepath = task.elastic_work.outdir.path_in("out_DDB")
            with DdbFile(ddb_filepath) as ddb:
                edata = ddb.anaget_elastic()
                elastic_relaxed = edata.elastic_relaxed.zeroed(1e-3)

            # Get relaxed structure and stress_guess and energy_guess from the GSR file.
            with task.open_gsr() as gsr:
                relaxed_structure = gsr.structure
                stress_guess = gsr.cart_stress_tensor * abu.GPa_to_au
                energy_guess = gsr.energy

            tdata = zsisa.get_tstress(task.temperature, task.pressure_gpa, self.mode,
                                      relaxed_structure, energy_guess, stress_guess,
                                      bo_elastic_voigt=elastic_relaxed.voigt)
            #print(tdata)

            # Init entry and add it to list.
            entry = dict(
                nqsmall_or_qppa=task.flow.nqsmall_or_qppa,
                pressure_gpa=task.pressure_gpa,
                temperature=task.temperature,
                gsr_path=task.gsr_path,
                # Add path to the DDB file with the 2nd order derivatives wrt strain.
                elastic_ddb_path=ddb_filepath,
                therm=tdata.therm,
                elastic=tdata.elastic,
                gibbs_atom=tdata.gibbs / len(relaxed_structure),
            )
            entry = ThermalRelaxEntry(**entry)

            results.thermal_relax_entries.append(entry)

        results.json_write(json_filepath, indent=4)


class ThermalRelaxTask(RelaxTask):
    """
    This task implements an iterative method to find the optimal lattice parameters
    at a given temperature T and external pressure P_ext.
    Starting from an initial lattice guess, the thermal and Born-Oppenheimer (BO) stresses are computed.
    A target stress is defined using the thermal stress and P_ext.
    Then, the lattice and atomic positions are relaxed repeatedly until the BO stress matches the target stress.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ThermalRelaxTask
    """
    def _post_init_(self):
        self.num_converged = 0

    def _on_ok(self):
        results = super()._on_ok()

        # Get relaxed structure, stress_guess and energy_guess from the GSR file.
        with self.open_gsr() as gsr:
            relaxed_structure = gsr.structure
            stress_guess = gsr.cart_stress_tensor * abu.GPa_to_au
            energy_guess = gsr.energy

        zsisa = self.work.zsisa
        tdata = zsisa.get_tstress(self.temperature, self.pressure_gpa, self.mode,
                                  relaxed_structure, energy_guess, stress_guess,
                                  bo_elastic_voigt=None)
        if tdata.converged:
            self.num_converged += 1
        else:
            self.num_converged -= 1
            self.num_converged = max(self.num_converged, 0)

        if not tdata.converged or (tdata.converged and self.num_converged == 1):
            # Change strtarget and restart the relaxation task.
            self.input.set_vars(strtarget=tdata.stress_au)
            if tdata.converged:
                self.input.set_vars(dilatmx=1.0)

            self.finalized = False
            # NB: Restart will take care of using the output structure as input.
            self.restart()

        else:
            if self.flow.with_elastic:
                # Build work for elastic constants and attach it to the task.
                scf_input = self.input.new_with_structure(relaxed_structure, ionmov=0, optcell=0)
                # Remove all irdvars before running. Important!
                scf_input.pop_irdvars()
                self.elastic_work = ElasticWork.from_scf_input(
                    scf_input, with_relaxed_ion=True, with_piezo=self.flow.with_piezo
                )
                self.flow.register_work(self.elastic_work)
                self.flow.allocate(build=True)

        return results


# TODO: Why Voigt notation for alpha? The matrix is not necessarily symmetric
_ALPHA_COMPS = (
    'alpha_xx',
    'alpha_yy',
    'alpha_zz',
    'alpha_yz',
    'alpha_xz',
    'alpha_xy',
)


@dataclasses.dataclass(kw_only=True)
class ThermalRelaxEntry:

    nqsmall_or_qppa: int          # Define the q-mesh for the computation of the PHDOS.
    pressure_gpa: float           # Pressure in GPa.
    temperature: float            # Temperature in K.
    gsr_path: str                 # Path to the GSR file.
    elastic_ddb_path: str         # Path to the DDB file with the 2nd order derivatives wrt strain.
    therm: np.ndarray | None      # Thermal_expansion (Voigt notation).
    elastic: np.ndarray | None    # Elastic constants (6,6) matrix in GPa (Voigt notation, relaxed-ions).
    gibbs_atom: float             # Gibbs free energy per atom in eV.

    def get_dict4pandas(self) -> dict:
        """
        Return dictionary to build pandas dataframes.
        """
        dct = {k: getattr(self, k) for k in ("temperature", "pressure_gpa", "gibbs_atom", "nqsmall_or_qppa")}
        with GsrFile(self.gsr_path) as gsr:
            dct.update(gsr.structure.get_dict4pandas())

        if self.therm is not None:
            # Add thermal expansion coefficients.
            dct.update({name: self.therm[i] for i, name in enumerate(_ALPHA_COMPS)})

        if self.elastic is not None:
            # Add elastic constants.
            for inds, value in np.ndenumerate(self.elastic):
                inds = np.array(inds, dtype=int)
                # Start to count from 1.
                key = f"C_{inds[0]+1}{inds[1]+1}"
                dct[key] = value

        return dct


@dataclasses.dataclass(kw_only=True)
class ZsisaResults(Serializable):
    """
    Main entry point for post-processing and visualizing the results of a Zsisa calculation.

    This object stores the locations of the GSR/DDB files produced by a ZSISA calculation
    so that we can easily redo a thermal relaxation run if needed.
    For instance, one may want to increase the list of temperatures, pressures or increase
    the q-mesh used to compute the phonon DOS via Fourier interpolation.

    To read the object from file use:

    .. code-block:: python

        from abipy.flowtk.zsisa import ZsisaResults
        data = ZsisaResults.json_load("json_filepath")
        data.get_dataframe()
        data.plot_lattice_vs_temp()
        data.plot_thermal_expansion()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ZsisaResults
    """

    spgrp_number: int                  # Space group number.
    eps: float                         # Strain magnitude to be applied to the lattice.
    mode: str                          # "TEC" or "ECs".
    qha_model: str

    # TODO: Should have inds_6d as well as strain_inds
    inds_6d: np.ndarray                # List of indices in the phdos 6D grid used for finite differences.

    gsr_bo_path: str                   # Path to GSR file with the relaxed BO configuration.
    gsr_relax_paths: list[str]         # Paths to GSR files for the deformed structures after relaxation.
    ddb_relax_paths: list[str]         # Paths to DDB files with phonons for the deformed structures after relaxation.
    gsr_relax_edos_paths: list[str]    # Paths to GSR files with electron DOS for the deformed structures after relaxation.
    gsr_relax_ebands_paths: list[str]  # Paths to GSR files with electron bands for the deformed structures after relaxation.

    thermal_relax_entries: list[ThermalRelaxEntry]

    def __post_init__(self):
        """
        It seems MontyDecoder does not support nested dataclasses.
        Here we convert dict to ThermalRelaxEntry.
        """
        for ie, entry in enumerate(self.thermal_relax_entries):
            if isinstance(entry, ThermalRelaxEntry): continue
            self.thermal_relax_entries[ie] = ThermalRelaxEntry(**entry)

    @property
    def has_thermal_expansion(self) -> bool:
        """True if thermal_expansion has been computed with elastic constants."""
        return all(entry.therm is not None for entry in self.thermal_relax_entries)

    @property
    def has_elastic(self) -> bool:
        """True if T-dep elastic constants are available."""
        return all(entry.elastic is not None for entry in self.thermal_relax_entries)

    @cached_property
    def cycle_markers(self):
        """Create a list of markers you want to cycle through."""
        from itertools import cycle
        return cycle(('o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h', '+', 'x'))

    def get_dataframe(self) -> pd.DataFrame:
        """
        Build pandas DataFrame with temperature, pressure_gpa, lattice parameters, thermal expansion
        and spacegroup info.
        """
        rows = [entry.get_dict4pandas() for entry in self.thermal_relax_entries]
        return pd.DataFrame(rows).sort_values(by="temperature")

    @cached_property
    def col2label(self) -> dict:
        """Dict mapping pandas column names to latex labels."""
        col2label = {
            "temperature": "T (K)",
            "a": r"a ($\AA$)",
            "b": r"b ($\AA$)",
            "c": r"c ($\AA$)",
            "alpha": r"$\alpha$ (degrees)",
            "beta": r"$\beta$ (degrees)",
            "gamma": r"$\gamma$ (degrees)",
            "gibbs_atom": r"$G$ (eV/atom)",
        }

        # Labels for thermal expansion coefficients.
        for alpha_comp in _ALPHA_COMPS:
            alpha, comp = alpha_comp.split("_")
            col2label[alpha_comp] = r"${\%s}_{%s}$ (K$^{-1}$)" % (alpha, comp)

        # Labels for elastic constants.
        for i, j in itertools.product(range(1, 7), range(1, 7)):
            key = f"C_{i}{j}"
            col2label[key] = "$C_{%s%s}$ (GPa)" % (i, j)

        return col2label

    def _get_cnames_from_c_select(self, c_select: str) -> list[str]:
        """
        Return list of strings with the names of the elastic tensor components from `c_select`.
        """
        if c_select == "symmetry":
            # Get names of columns associated to the independent elastic constants.
            sym = spgnum_to_crystal_system(self.spgrp_number)
            _, c_names = cmat_inds_names(sym, self.mode)
        elif c_select == "all":
            c_names = [k for k in df.keys() if k.startswith("C_")]
        else:
            raise ValueError(f"Invalid {c_select=}")

        return c_names

    @staticmethod
    def _add_pga_col(df: pd.DataFrame) -> tuple[str, pd.DataFrame]:
        """Add new column to df with p_key to have nicer labels."""
        p_key = "P (GPa)"
        df = df.copy()
        df[p_key] = df["pressure_gpa"]
        return p_key, df

    @add_fig_kwargs
    def plot_lattice_vs_temp(self,
                             df: pd.DataFrame | None = None,
                             fontsize: int = 8,
                             **kwargs) -> Figure:
        """
        Plot lattice parameters and angles as a function of T grouped by pressure P.

        Args:
            df: dataframe with data. None to compute it inside the function.
            fontsize: fontsize for legends and titles.
        """
        angles = ["alpha", "beta", "gamma"]
        lengths = ["a", "b", "c"]

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=3, ncols=2,
                                               sharex=True, sharey=False, squeeze=False)
        if df is None:
            df = self.get_dataframe()
        one_pressure = df["pressure_gpa"].nunique() == 1

        # Add new column with p_key to have nicer labels.
        p_key, df = self._add_pga_col(df)

        plt_kwargs = dict(fontsize=fontsize,
                          hue=p_key if not one_pressure else None,
                          col2label=self.col2label,
                          marker="o",
                          show=False
                          )

        for ix, length in enumerate(lengths):
            ax = ax_mat[ix, 0]
            plot_xy_with_hue(df, "temperature", length, ax=ax, **plt_kwargs)
            if ix != len(lengths) - 1:
                set_visible(ax, False, *["xlabel"])
            if ix != 0:
                set_visible(ax, False, *["legend"])

        for ix, angle in enumerate(angles):
            ax = ax_mat[ix, 1]
            plot_xy_with_hue(df, "temperature", angle, ax=ax, **plt_kwargs)
            if ix != len(angles) - 1:
                set_visible(ax, False, *["xlabel"])
            set_visible(ax, False, *["legend"])

        if "title" not in kwargs:
            fig.suptitle("Lattice parameters and angles")

        return fig

    @add_fig_kwargs
    def plot_thermal_expansion(self,
                               df: pd.DataFrame | None = None,
                               fontsize: int = 8,
                               **kwargs) -> Figure:
        """
        Plot thermal expansion alpha as a function of T grouped by pressure P.

        Args:
            df: dataframe with data. None to compute it inside the function.
            fontsize: fontsize for legends and titles.
        """
        if not self.has_thermal_expansion:
            raise ValueError("Thermal expansion coefficients are not available!")

        nrows, ncols = 3, 2
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        if df is None:
            df = self.get_dataframe()
        one_pressure = df["pressure_gpa"].nunique() == 1

        # Add new column with p_key to have nicer labels.
        p_key, df = self._add_pga_col(df)

        plt_kwargs = dict(fontsize=fontsize,
                          hue=p_key if not one_pressure else None,
                          col2label=self.col2label,
                          marker="o",
                          show=False
                          )

        alpha_comps_mat = np.reshape(_ALPHA_COMPS, (2, 3)).T

        for ii, jj in itertools.product(range(nrows), range(ncols)):
            ax = ax_mat[ii, jj]
            alpha_name = alpha_comps_mat[ii, jj]
            plot_xy_with_hue(df, "temperature", alpha_name, ax=ax, **plt_kwargs)
            if ii != nrows - 1:
                set_visible(ax, False, *["xlabel"])
            if jj == 1:
                set_visible(ax, False, *["ylabel"])
            if (ii, jj) != (0, 0):
                set_visible(ax, False, *["legend"])

        if "title" not in kwargs:
            fig.suptitle("Thermal expansion coefficients")

        return fig

    @add_fig_kwargs
    def plot_elastic_vs_t(self,
                          pressure_gpa: float | None = None,
                          c_select: str = "symmetry",
                          df: pd.DataFrame | None = None,
                          ax=None,
                          colormap: str = "jet",
                          fontsize: int = 8,
                          **kwargs) -> Figure:
        """
        Plot elastic constants as a function of T for fixed pressure on a single figure.

        Args:
            pressure_gpa: Pressure to select. If None, the minimum pressure is used.
            c_select: "symmetry" if only the non-zero components should be plotted.
                "all" to plot all components.
            df: dataframe with data. None to compute it inside the function.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: Color map. Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            fontsize: fontsize for legends and titles.
        """
        if not self.has_elastic:
            raise ValueError("T-dependent elastic constants are not available!")

        if df is None:
            df = self.get_dataframe()

        if pressure_gpa is None:
            pressure_gpa = df["pressure_gpa"].values.min()

        # Add new column with p_key to have nicer labels.
        p_key, df = self._add_pga_col(df)
        df = df[df[p_key] == pressure_gpa]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        c_names = self._get_cnames_from_c_select(c_select)
        for i, c_name in enumerate(c_names):
            plt_kwargs = dict(
                marker=next(self.cycle_markers),
                color=cmap(float(i) / len(c_names)),
                label=c_name,
            )
            ax.plot(df["temperature"], df[c_name], **plt_kwargs)

        set_grid_legend(ax, fontsize, xlabel="T (K)", ylabel="Elastic constant (GPa)")

        if "title" not in kwargs:
            fig.suptitle(f"Temperature-dependent elastic constants at P={pressure_gpa} (GPa)")

        return fig

    @add_fig_kwargs
    def plot_elastic_vs_t_and_p(self,
                                c_select: str = "symmetry",
                                df: pd.DataFrame | None = None,
                                ax=None,
                                colormap: str = "jet",
                                fontsize: int = 8,
                                **kwargs) -> Figure:
        """
        Plot elastic constants as a function of T grouped by pressure.
        One subplot for each element of the C tensor.

        Args:
            c_select: "symmetry" if only the non-zero components should be plotted.
                "all" to plot all components.
            df: dataframe with data. None to compute it inside the function.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: Color map. Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            fontsize: fontsize for legends and titles.
        """
        if not self.has_elastic:
            raise ValueError("T-dependent elastic constants are not available!")

        if df is None:
            df = self.get_dataframe()

        # Add new column with p_key to have nicer labels.
        p_key, df = self._add_pga_col(df)

        c_names = self._get_cnames_from_c_select(c_select)
        nrows, ncols = len(c_names), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=True)
        cmap = plt.get_cmap(colormap)

        for ix, (c_name, ax) in enumerate(zip(c_names, ax_list, strict=True)):
            plt_kwargs = dict(
                fontsize=fontsize,
                marker=next(self.cycle_markers),
                hue=p_key,
                col2label=self.col2label,
                show=False,
            )

            plot_xy_with_hue(df, "temperature", c_name, ax=ax, **plt_kwargs)

            if ix != 0:
                set_visible(ax, False, *["legend"])

            if ix != len(c_names) - 1:
                set_visible(ax, False, *["xlabel"])

            #set_grid_legend(ax, fontsize, xlabel="T (K)", ylabel="Elastic constant (GPa)")

        #if "title" not in kwargs:
        #    fig.suptitle(f"Temperature-dependent elastic constants at P={pressure_gpa} (GPa)")

        return fig

    @add_fig_kwargs
    def plot_gibbs(self,
                   df: pd.DataFrame | None = None,
                   ax=None,
                   colormap: str = "jet",
                   fontsize: int = 8,
                   **kwargs) -> Figure:
        """
        Plot Gibbs energy per atom in eV as a function of T grouped by pressure in Gpa.

        Args:
            df: dataframe with data. None to compute it inside the function.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: Color map. Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            fontsize: fontsize for legends and titles.
        """
        if df is None:
            df = self.get_dataframe()

        # Add new column with p_key to have nicer labels.
        p_key, df = self._add_pga_col(df)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        plt_kwargs = dict(
            fontsize=fontsize,
            hue=p_key,
            col2label=self.col2label,
            marker="o",
            show=False,
        )

        plot_xy_with_hue(df, "temperature", "gibbs_atom", ax=ax, **plt_kwargs)

        #if "title" not in kwargs:
        #    fig.suptitle(f"Temperature-dependent elastic constants at P={pressure_gpa} (GPa)")

        return fig

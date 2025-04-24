# coding: utf-8
"""Work subclasses related to GS calculations."""
from __future__ import annotations

# TODO: Should we allow for relax and relax_opts?

import sys
#import json
import itertools
#import pickle
import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from dataclasses import field
from typing import Optional
#from monty.json import MSONable
from monty.string import list_strings #, marquee
from monty.functools import lazy_property
from pymatgen.analysis.elasticity.strain import Strain
from abipy.core.structure import Structure
from abipy.tools.numtools import build_mesh
from abipy.tools.derivatives import central_fdiff_weights, check_num_points_for_order # finite_diff
from abipy.tools.typing import Figure
from abipy.abio.inputs import AbinitInput
from abipy.abio.outputs import BerryPhasePolarization
from abipy.abio.enums import StrEnum
from abipy.tools.serialization import HasPickleIO
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_grid_legend,
    rotate_ticklabels, set_visible, set_ax_xylabels)
#from abipy.tools.serialization import mjson_write #, pmg_serialize
from .works import Work

#def centered_indices(n):
#    half = n // 2
#    if n % 2 == 0:
#        return list(range(-half, half)), half
#    #else:
#    #return list(range(-half, half + 1)),


def idir2s(idir: int):
    """Convert direction index to string."""
    return {0: "x", 1: "y", 2: "z"}[idir]


def dir2str(coeffs, variables='xyz'):
    """
    >>> dir2str((1, 2, 3)))
    "x + 2y + 3z"
    >>> dir2str((0, -1, 4)))
    "-y + 4z"
    >>> dir2str((1, 0, 0)))
    "x"
    >>> dir2str((0, 0, 0)))   # Output: (empty string)
    ""
    """
    terms = []
    for coeff, var in zip(coeffs, variables, strict=True):
        if coeff == 0:
            continue  # Skip terms with a coefficient of 0
        elif coeff == 1:
            terms.append(f"{var}")
        elif coeff == -1:
            terms.append(f"-{var}")
        else:
            terms.append(f"{coeff}{var}")

    return ' + '.join(terms).replace('+ -', '- ')


def mat33_to_voigt(mat: np.ndarray, engineering_strain: bool = False) -> np.ndarray:
    """
    Convert a 3x3 symmetric mat to a Voigt vector (6x1).

    Parameters:
        mat (np.ndarray): 3x3 symmetric matrix (stress/strain mat)
        engineering_strain (bool): If True, shear components are multiplied by 2 (engineering strain convention)

    Returns: np.ndarray: 6-element Voigt vector
    """
    v = np.zeros(6)
    v[0] = mat[0, 0]
    v[1] = mat[1, 1]
    v[2] = mat[2, 2]
    factor = 2 if engineering_strain else 1
    v[3] = factor * mat[1, 2]
    v[4] = factor * mat[0, 2]
    v[5] = factor * mat[0, 1]
    return v


def voigt_to_mat33(voigt: np.ndarray, engineering_strain: bool = True) -> np.ndarray:
    """
    Convert a Voigt vector (6x1) to a 3x3 symmetric tensor.

    Parameters:
        voigt (np.ndarray): 6-element vector
        engineering_strain (bool): If True, shear components are divided by 2 (engineering strain convention)

    Returns: np.ndarray: 3x3 symmetric matrix
    """
    mat = np.zeros((3, 3))
    mat[0, 0] = voigt[0]
    mat[1, 1] = voigt[1]
    mat[2, 2] = voigt[2]
    factor = 0.5 if engineering_strain else 1
    mat[1, 2] = mat[2, 1] = factor * voigt[3]
    mat[0, 2] = mat[2, 0] = factor * voigt[4]
    mat[0, 1] = mat[1, 0] = factor * voigt[5]
    return


def _dict_from_mat_npts(mat: np.ndarray, mat_comps: list[str], npts: int) -> dict:
    d = {"npts": npts}
    d.update({c: v for c, v in zip(mat_comps, mat.flatten(), strict=True)})
    if mat.shape[0] == mat.shape[1]:
        d["iso_avg"], d["det"] = np.trace(mat) / mat.shape[0] , np.linalg.det(mat)
    return d


class PertKind(StrEnum):
    DISPL = "displ"
    E = "E"
    H = "H"
    STRAIN = "strain"


@dataclasses.dataclass
class Perturbation:
    """
    Stores info on the perturbation.
    """
    kind: str
    #values: np.ndarray
    #ip0: int
    cart_dir: np.ndarray | None = None
    iatom: int | None = None
    strain: np.ndarray | None = None

    def __post_init__(self):
        """
        Validation logic.
        """
        if self.kind not in PertKind:
            raise ValueError(f"Invalid {self.kind=}")

        if self.kind == PertKind.DISPL:
            if self.iatom is None or self.cart_dir is None:
                raise ValueError("iatom and cart_dir must be specified for a `displ` perturbation.")

        if self.kind in (PertKind.E, PertKind.H):
            if self.cart_dir is None:
                raise ValueError("cart_dir must be specified for a `E` or `H` perturbations.")

        if self.kind == PertKind.STRAIN:
            if self.strain is None:
                raise ValueError("strain matrix must be specified for a `strain` perturbation.")

    #@lazy_property
    #def step(self) -> float:
    #    dx = np.zeros(len(self)-1)
    #    for (i, x) in enumerate(self.mesh[:-1]):
    #        dx[i] = self.mesh[i+1] - x
    #    return self.dx[0] if np.allclose(self.dx[0], self.dx) else None

    @lazy_property
    def label(self) -> str:
        return f"{self.name}: {dir2str(self.cart_dir)}"

    @lazy_property
    def dir_str(self) -> str:
        return f"{dir2str(self.cart_dir)}"

    @lazy_property
    def tex(self) -> str:
        return {
            PertKind.E: r"{\mathcal{E}}",
            PertKind.H: r"{\mathcal{H}}",
            PertKind.DISPL: r"{\Delta\tau}",
            PertKind.STRAIN: r"{\varepsilon}",
        }[self.kind]

    @lazy_property
    def name(self) -> str:
        return {
            PertKind.E: "Electric field",
            PertKind.H: "Magnetic field",
            PertKind.DISPL: "Atomic displacement",
            PertKind.STRAIN: "Strain",
        }[self.kind]


class _BaseFdWork(Work):

    def get_data_dict(self) -> dict:
        """
        Note: The suffix `_pv` stands for Field-direction and Field-Value.
        """
        natom = len(self[0].input.structure)
        npert, np_vals = len(self.perts), len(self.pert_values)
        has_pol = any(task.input.get("berryopt", 0) != 0 for task in self)
        # FIXME
        #has_mag = False
        #has_mag = any(task.input.get("zeeman", 0) != 0 for task in self)

        data = {
            "input_structure": self[0].input.structure,
            "relaxed_structure": self.relaxed_structure,
            "step_au": self.step_au,
            "perts": self.perts,
            "pert_values": self.pert_values,
            "pert_cart_dirs": self.pert_cart_dirs,
            "ipv0": self.ipv0,
            "params_p": [],
            "has_pol": has_pol,
            #"has_mag": has_mag,
        }

        data["etotals_pv"] = etotals_pv = np.empty((npert, np_vals))
        data["eterms_pv"] = eterms_pv = np.empty((npert, np_vals), dtype=object)
        data["cart_forces_pv"] = cart_forces_pv = np.empty((npert, np_vals, natom, 3))
        data["carts_stresses_pv"] = carts_stresses_pv = np.empty((npert, np_vals, 3, 3))

        if has_pol:
            data["cart_pol_pv"] = cart_pol_pv = np.empty((npert, np_vals, 3))
            data["cart_pole_pv"] = cart_pole_pv = np.empty((npert, np_vals, 3))
            data["cart_poli_pv"] = cart_poli_pv = np.empty((npert, np_vals, 3))

        # Read energy, forces and stress from the GSR files.
        for ip, pert in enumerate(self.perts):
            for ipv, p_val in enumerate(self.pert_values):
                task = self.tasks_pv[ip, ipv]
                with task.open_gsr() as gsr:
                    etotals_pv[ip, ipv] = gsr.r.read_value("etotal")
                    eterms_pv[ip, ipv] = gsr.r.read_energy_terms(unit="Ha")
                    cart_forces_pv[ip, ipv] = gsr.r.read_value("cartesian_forces") # Ha/Bohr units.
                    carts_stresses_pv[ip, ipv] = gsr.r.read_cart_stress_tensor(units="au")
                    # Add parameters that might be used for convergence studies.
                    data["params_p"].append(gsr.params)

                if has_pol:
                    with task.open_abo() as abo:
                        pol = abo.get_berry_phase_polarization()
                        cart_pol_pv[ip, ipv] = pol.total
                        cart_pole_pv[ip, ipv] = pol.electronic
                        cart_poli_pv[ip, ipv] = pol.ionic

        return data


class FiniteDisplWork(_BaseFdWork):
    """
    Work for the computation ... finite differences.
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       num_points: int,
                       step_au: float = 0.01,
                       pert_cart_dirs=None,
                       mask_iatom=None,
                       relax: bool = False,
                       relax_opts: dict | None = None,
                       manager=None):
        """
        Build a work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            num_points:
            step_au: Finite difference step for the displacement in Bohr (a.u.)
            pert_cart_dirs:
            mask_iatom:
            relax: False if the initial structural relaxation should not be performed.
            relax_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        work.scf_input = scf_input.deepcopy()
        structure = scf_input.structure
        natom = len(structure)

        work.step_au = float(step_au)
        work.pert_values, work.ipv0 = build_mesh(0.0, num_points, step_au, "=")
        work.num_deltas = len(work.pert_values)

        # Here we normalize the directions to 1. NB: pymatgen structures uses Ang and not Bohr.
        if pert_cart_dirs is not None:
            work.pert_cart_dirs = np.eye(3)

        work.pert_cart_dirs = np.reshape(work.pert_cart_dirs, (-1, 3))

        for idir, cart_dir in enumerate(work.pert_cart_dirs):
            norm = structure.lattice.norm(cart_dir, frac_coords=False)
            work.pert_cart_dirs[idir] = cart_dir / norm

        if mask_iatom is None:
            mask_iatom = np.ones(natom, dtype=bool)

        work.mask_iatom = np.array(mask_iatom).astype(bool)
        if len(work.mask_iatom) != natom:
            raise ValueError(f"{len(work.mask_iatom)=} != {natom=}")

        work.num_dirs = len(work.pert_cart_dirs)

        # Build list of perturbations.
        work.perts = []
        for iatom, mask in zip(range(natom), work.mask_iatom, strict=True):
            if not mask: continue
            for cart_dir in work.pert_cart_dirs:
                work.perts.append(Perturbation(kind=PertKind.DISPL, cart_dir=cart_dir, iatom=iatom))

        work.relax = relax
        if work.relax:
            relax_input = scf_input.make_relax_input(**relax_opts)
            work.initial_relax_task = work.register_relax_task(relax_input)
        else:
            work._add_tasks_with_displacements(scf_input.structure)
            work.relaxed_structure = scf_input.structure

        return work

    def _add_tasks_with_displacements(self, structure: Structure):
        nperts = len(self.perts)
        self.tasks_pv = np.empty((nperts, self.num_deltas), dtype=object)

        for ip, pert in enumerate(self.perts):
            iatom, cart_dir = pert.iatom, pert.cart_dir
            for iv, delta_au in enumerate(self.pert_values):
                new_structure = structure.copy()
                # Note Bohr --> Ang conversion.
                new_structure.translate_sites([iatom], delta_au * abu.Bohr_Ang * cart_dir,
                                               frac_coords=False, to_unit_cell=False)
                new_input = self.scf_input.new_with_structure(new_structure)
                task = self.register_scf_task(new_input)
                self.tasks_pv[ip, iv] = task

    def on_ok(self, sender):
        """This method is called when one task reaches status `S_OK`."""
        if self.relax and sender == self.initial_relax_task:
            # Get relaxed structure from GSR file.
            self.relaxed_structure = sender.get_final_structure()
            self._add_tasks_with_displacements(self.relaxed_structure)
            self.flow.allocate(build=True)

        return super().on_ok(sender)

    @lazy_property
    def natom(self) -> int:
        """Number of atoms in the unit cell."""
        return len(work.scf_input.structure)

    def get_data(self):
        d = self.get_data_dict()
        return DisplData(**d)

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        return super().on_all_ok()


class FiniteStrainWork(_BaseFdWork):
    """
    Work for the computation at finite Strain with finite differences.
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input,
                       num_points: int,
                       norm_step: float,
                       shear_step: float,
                       voigt_inds=None,
                       relax: bool = False,
                       relax_opts: dict | None = None,
                       manager=None):
        """
        Build a Work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            num_points: Number of points for finite difference.
            norm_step: Finite difference step for normal strain.
            shear_step: Finite difference step for shear strain.
            voigt_inds:
            relax: False if the initial structural relaxation should not be performed.
            relax_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default if None.
        """
        work = cls(manager=manager)
        work.scf_input = scf_input.deepcopy()

        if "ecutsm" not in scf_input:
            ecutsm = 0.5
            work.scf_input.set_vars(ecutsm=ecutsm)
            print("AbinitInput does not define ecutsm.\n A default value of %s will be added" % ecutsm)

        if voigt_inds is None:
            strains = []
            # Normal strain.
            for ind in [(0, 0), (1, 1), (2, 2)]:
                strains.append(Strain.from_index_amount(ind, amount=1))
            # Shear strain.
            for ind in [(0, 1), (0, 2), (1, 2)]:
                strains.append(Strain.from_index_amount(ind, amount=1))

        strains = np.reshape(strains, (-1, 3, 3))
        for strain in strains:
            if not np.array_equal(strain, strain.T):
                raise ValueError(f"The strain matrix should be symmetric but got: {strain}")

        work.pert_values, work.ipv0 = build_mesh(0.0, num_points, step, "=")
        work.num_points = len(work.pert_values)
        check_num_points_for_order(num_points=work.num_points, order=1, kind="=")

        # Build list of perturbations.
        work.perts = [Perturbation(kind=PertKind.STRAIN, strain=strain) for strain in work.strains]

        #from pymatgen.analysis.elasticity.strain import DeformedStructureSet
        #DeformedStructureSet(structure: Structure,
        #                     norm_strains: Sequence[float] = (-0.01, -0.005, 0.005, 0.01),
        #                     shear_strains: Sequence[float] = (-0.06, -0.03, 0.03, 0.06),
        #                     symmetry=False,

        # FIXME Different pert_values for normal and shear strain
        npert = len(work.perts)
        work.tasks_pv = np.empty((npert, len(work.pert_values)), dtype=object)
        for ip, pert in enumerate(work.perts):
            for ipv, pert_value in enumerate(work.pert_values):
                # Apply strain to the lattice.
                new_structure = scf_input.structure.apply_strain(pert_value * pert.strain, inplace=False)
                work.tasks_pv[ip, ipv] = work.register_scf_task(scf_input.new_with_structure(new_structure))

        return work

    #def on_all_ok(self):
    #    """This method is called when all tasks have reached S_OK."""
    #    data = self.get_data_dict()
    #    #data.pickle_dump(self.outdir.path)
    #    return super().on_all_ok()


class _FieldWork(_BaseFdWork):
    """Base class for finite field + finite difference Work."""

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       num_points: int,
                       step_au: float,
                       pert_cart_dirs: np.ndarray | None = None,
                       relax: bool = False,
                       relax_opts: dict | None = None,
                       manager=None):
        """
        Build the work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF calculation used as template to generate the other inputs.
            num_points: Number of points for finite difference.
            step_au: Finite difference step for the magnetic field in a.u.
            pert_cart_dirs:
            relax: False if the initial structural relaxation should not be performed.
            relax_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        work.step_au = step_au
        work.pert_values, work.ipv0 = build_mesh(0.0, num_points, step_au, "=")
        work.num_points = len(work.pert_values)
        check_num_points_for_order(num_points=work.num_points, order=1, kind="=")

        if pert_cart_dirs is None:
            work.pert_cart_dirs = np.eye(3)

        work.pert_cart_dirs = np.reshape(work.pert_cart_dirs, (-1, 3))

        work.scf_input_template = scf_input.deepcopy()

        if isinstance(work, FiniteEfieldWork):
            work.pert_kind = PertKind.E
        elif isinstance(work, FiniteHfieldWork):
            work.pert_kind = PertKind.H
        else:
            raise TypeError(f"Don't know how to handle {type(work)=}")

        # Build list of perturbations.
        work.perts = [Perturbation(kind=work.pert_kind, cart_dir=cart_dir) for cart_dir in work.pert_cart_dirs]

        work.relax = relax
        if work.relax:
            relax_input = scf_input.make_relax_input(**relax_opts)
            work.initial_relax_task = work.register_relax_task(relax_input)
        else:
            if work.pert_kind == PertKind.E:
                work._add_tasks_with_efield(scf_input.structure)
            elif work.pert_kind == PertKind.H:
                work._add_tasks_with_zeemanfield(scf_input.structure)
            else:
                raise TypeError(f"Don't know how to handle {work.pert_kind=}")

            work.relaxed_structure = scf_input.structure

        return work

    def on_ok(self, sender):
        """This method is called when one task reaches status `S_OK`."""
        if self.relax and sender == self.initial_relax_task:
            # Get relaxed structure from GSR file.
            self.relaxed_structure = sender.get_final_structure()

            if work.pert_kind == PertKind.E:
                self._add_tasks_with_efield(self.relaxed_structure)
            elif work.pert_kind == PertKind.H:
                self._add_tasks_with_zeemanfield(self.relaxed_structure)
            else:
                raise TypeError(f"Don't know how to handle {work.pert_kind}")

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        return super().on_all_ok()


class FiniteHfieldWork(_FieldWork):
    r"""
    Work for the computation of the dynamical magnetic charges with finite differences.

    The dynamical magnetic charges are defined as:

        Z_jv^m=Ω_0 (∂M_v)/(∂u_j ) = (∂F_j)/(∂H_v ) = Ω_0 (∂^2 E)/(∂H_β ∂u_i).

    Here we compute them as derivatives of forces wrt to the Zeeman magnetic field.
    """

    def _add_tasks_with_zeemanfield(self, structure: Structure) -> None:
        """Build new GS tasks with zeemanfield."""
        scf_input = self.scf_input_template.new_with_structure(structure)

        npert, np_vals = self.npert, len(self.pert_values)
        self.tasks_pv = np.empty((npert, np_vals), dtype=object)
        task_pv0 = None

        for ip, pert in enumerate(self.perts):
            for ipv, p_val in enumerate(self.pert_values):
                is_pv0 = abs(p_val) < 1e-16
                new_inp = scf_input.new_with_vars(zeemanfield=p_val * pert.cart_dir)
                if is_pv0:
                    # Avoid computing H=0 multiple times.
                    if task_pv0 is None:
                        task_pv0 = self.register_scf_task(new_inp)
                    self.tasks_pv[ip, ipv] = task_pv0
                else:
                    self.tasks_pv[ip, ipv] = self.register_scf_task(new_inp)

    def get_data(self) -> ZeemanData:
        """
        Read data from the GSR files, and compute Zm with finite differences.
        """
        d = super().get_data_dict()
        return ZeemanData(**d)


class FiniteEfieldWork(_FieldWork):
    r"""
    """
    def _add_tasks_with_efield(self, structure: Structure) -> None:
        """Build new GS tasks with finite electric field."""
        scf_input = self.scf_input_template.new_with_structure(structure)

        npert, np_vals = len(self.perts), len(self.pert_values)
        self.tasks_pv = np.empty((npert, np_vals), dtype=object)
        task_pv0 = None

        for ip, pert in enumerate(self.perts):
            for ipv, p_val in enumerate(self.pert_values):
                is_pv0 = abs(p_val) < 1e-16
                new_inp = scf_input.new_with_vars(efield=p_val * pert.cart_dir)
                if is_pv0:
                    # Avoid computing the zero-field case multiple times.
                    # Also the task at zero field uses berryopt -1 to get the polarization.
                    if task_pv0 is None:
                        new_inp.set_vars(berryopt=-1)
                        task_pv0 = self.register_berry_task(new_inp)
                    self.tasks_pv[ip, ipv] = task_pv0
                else:
                    # Finite E-field.
                    new_inp.set_vars(berryopt=4)
                    self.tasks_pv[ip, ipv] = self.register_berry_task(new_inp)

        # Now add dependencies to the tasks.
        for ip, pert in enumerate(self.perts):
            for ipv in range(0, self.ipv0):
                self.tasks_pv[ip, ipv].add_deps({self.tasks_pv[ip, ipv+1]: "WFK"})

            ntasks = len(self.tasks_pv[ip])
            for ipv in range(self.ipv0+1, ntasks):
                self.tasks_pv[ip, ipv].add_deps({self.tasks_pv[ip, ipv-1]: "WFK"})

    def get_data(self) -> ElectricFieldData:
        d = self.get_data_dict()
        return ElectricFieldData(**d)


@dataclasses.dataclass(kw_only=True)
class _FdData(HasPickleIO):
    """
    All values are in a.u.

    _pv stands for perturbation and perturbation Value.
    """
    input_structure: Structure
    relaxed_structure: Structure
    step_au: float
    ipv0: int
    has_pol: bool
    #has_mag: bool

    perts: list[Perturbation]
    pert_values: np.ndarray
    pert_cart_dirs: np.ndarray
    etotals_pv: np.ndarray
    eterms_pv: np.ndarray
    cart_forces_pv: np.ndarray
    carts_stresses_pv: np.ndarray
    params_p: list[dict]

    cart_pol_pv: Optional[np.ndarray] = None
    cart_pole_pv: Optional[np.ndarray] = None
    cart_poli_pv: Optional[np.ndarray] = None

    # npts -> dForce/dPert with shape [iat, 3, ip] in Cart. coords.
    dforces_dpert_npts: dict[int, np.array] = field(init=False)

    # npts -> dStress/dPert with shape [3, 3, ip] in Cart. coords.
    dstress_dpert_npts: dict[int, np.array] = field(init=False)

    dpol_dpert_npts: dict[int, np.array] = field(init=False)

    def __post_init__(self):
        """
        Compute quantities with finite differences.
        """
        natom, npert, np_vals = len(self.input_structure), self.npert, len(self.pert_values)

        # Use all stencils compatible with input num_points so that we can monitor the convergence afterwards.
        self.dforces_dpert_npts = {}
        self.dstress_dpert_npts = {}
        if self.has_pol:
            self.dpol_dpert_npts = {}

        for acc, weights in central_fdiff_weights[1].items():
            if np_vals < len(weights): continue
            nn = acc // 2
            fd_slice = slice(self.ipv0 - nn, self.ipv0 + nn + 1)
            npts = len(weights)

            # Finite differences for forces.
            dforce_dpert = np.empty((natom, 3, npert))
            for iat, iat_dir, ip in itertools.product(range(natom), range(3), range(npert)):
                fvals_f = self.cart_forces_pv[ip, :, iat, iat_dir]
                dforce_dpert[iat, iat_dir, ip] = np.sum(fvals_f[fd_slice] * weights) / self.step_au

            self.dforces_dpert_npts[npts] = dforce_dpert

            # Finite differences for stresses.
            dstress_dpert = np.empty((3, 3, npert))
            for ii, jj, ip in itertools.product(range(3), range(3), range(npert)):
                svals_f = self.carts_stresses_pv[ip, :, ii, jj] * self.relaxed_structure.volume # * abu.Ang_Bohr ** 3  TODO?
                dstress_dpert[ii, jj, ip] = np.sum(svals_f[fd_slice] * weights) / self.step_au

            self.dstress_dpert_npts[npts] = dstress_dpert

            # Finite differences for polarizations (if available).
            if self.has_pol:
                dpol_dpert = np.empty((3, npert))
                for ii, ip in itertools.product(range(3), range(npert)):
                    dpol_dpert[ii, ip] = np.sum(self.cart_pol_pv[ip, fd_slice, ii] * weights) / self.step_au
                self.dpol_dpert_npts[npts] = dpol_dpert
                #self.dpole_dpert_npts[npts] = dpole_dpert  TODO ?
                #self.dpoli_dpert_npts[npts] = dpoli_dpert  TODO ?

    @lazy_property
    def npert(self) -> int:
        """Number of field directions."""
        return len(self.perts)

    @lazy_property
    def pert_kind(self) -> str:
        """
        """
        pert_kind = self.perts[0].kind
        all_kinds = [p.kind for p in self.perts]
        if any(_ != pert_kind for _ in all_kinds):
            raise ValueError(f"Expecting perturbations of the same kind but got {all_kinds}")
        return pert_kind

    @lazy_property
    def natom(self) -> int:
        """Numbef of atoms in the unit cell."""
        return len(self.input_structure)

    def label_for_pdir(self, pdir: int) -> str:
        pdir_str = dir2str(self.pert_cart_dirs[pdir])
        return "$%s_{%s}$ (a.u.)" % (self.pert_tex, pdir_str)

    @lazy_property
    def pert_dir_comps(self) -> list[str]:
        """
        """
        if any(pert.cart_dir is None for pert in self.perts):
            raise TypeError("pert_dir_comps requires perturbation with directions!")

        return [dir2str(pert.cart_dir) for pert in self.perts]

    @add_fig_kwargs
    def plot_etotal(self, mode="diff", ax=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot energies as a function of the amplitude of the perturbation.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for ip, pert in enumerate(self.perts):
            e_values = self.etotals_pv[ip] * abu.Ha_meV / self.natom
            if mode == "diff":
                e_values -= e_values[self.ipv0]
            ax.plot(self.pert_values, e_values, marker="o", label=pert.label)

        set_grid_legend(ax, fontsize, xlabel=f"${pert.tex}$ (a.u.)",
                        ylabel=r"$\Delta$ Energy/atom (meV)" if mode == "diff" else "Energy/atom (meV)")
        return fig

    @add_fig_kwargs
    def plot_forces(self, elements=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot Cartesian forces as a function of the of the amplitude of the perturbation.
        """
        if elements is not None: elements = list_strings(elements)
        nrows, ncols = 3, self.npert
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        for iat_dir in range(3):
            for ip, pert in enumerate(self.perts):
                ax = ax_mat[iat_dir, ip]
                ax.set_title(f"{pert.label}, Atom_dir: {idir2s(iat_dir)}", fontsize=fontsize)
                for iat, site in enumerate(self.relaxed_structure):
                    if elements is not None and site.species_string not in elements: continue
                    ax.plot(self.pert_values, self.cart_forces_pv[ip, :, iat, iat_dir], marker="o",
                            label=site.species_string + r"$_{\text{%s}}$" % iat)

                ax.legend(loc="best", fontsize=fontsize, shadow=True)
                #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig

    @add_fig_kwargs
    def plot_stresses(self, fontsize=8, **kwargs) -> Figure:
        """
        Plot Cartesian stresses as a function of the finite external field for all the directions.
        """
        nrows, ncols = self.npert, 1
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        for ip, pert in enumerate(self.perts):
            ax = ax_mat[ip]
            ax.set_title(pert.label, fontsize=fontsize)
            for ii, jj in itertools.product(range(3), range(3)):
                ax.plot(self.pert_values, self.carts_stresses_pv[ip, :, ii, jj], marker="o",
                        label=r"$\sigma_{%s}$" % (f"{ii}{jj}"))

            ax.legend(loc="best", fontsize=fontsize, shadow=True)
            #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig

    def get_df_zeff_iatom(self, iatom: int) -> pd.Dataframe:
        """
        Return dataframe with effective charges for the given atom index and all the FD points.
        """
        field2zeff = {PertKind.E: "Ze", PertKind.H: "Zm"}

        if self.pert_kind in field2zeff:
            zeff_name, what_to_diff = [self.pert_kind], "forces"

        elif self.pert_kind == PertKind.DISPL:
            if self.has_pol:
                zeff_name, what_to_diff = "Ze", "polarization"
            #elif self.has_mag:
            #    zeff_name, what_to_diff = "Zm", "magnetization"
        else:
            raise ValueError(f"Don't know how to compute eff_charges with {self.pert_kind=}")

        xyz_comps = "x y z".split()
        rows = []

        if what_to_diff == "forces":
            zeff_comps = list(itertools.product(xyz_comps, self.pert_dir_comps))
            for npts, dforces_dpert in self.dforces_dpert_npts.items():
                zeff_atm = dforces_dpert[iatom]
                rows.append(_dict_from_mat_npts(zeff_atm, zeff_comps, npts))

        if what_to_diff == "polarization":
            for npts, dpol_dpert in self.dpol_dpert_npts.items():
                # dpol_dpert has shape (3, npert) where npert are atomic displacements.
                cnt = 0
                zeff_atm, atom_comps = np.empty((3, 3)), []
                for ip, pert in enumerate(self.perts):
                    if pert.iatom != iatom: continue
                    cnt += 1
                    iat_dir = ip % 3
                    zeff_atm[iat_dir,:] = dpol_dpert[:, ip]
                    atom_comps.append(pert.dir_str)

                if cnt != 3:
                    raise RuntimeError(f"Need all 3 directions for {iatom=} to compute Ze from polarization!")

                zeff_comps = list(itertools.product(atom_comps, xyz_comps))
                zeff_atm *= self.relaxed_structure.volume * abu.Ang_Bohr ** 3
                rows.append(_dict_from_mat_npts(zeff_atm, zeff_comps, npts))

        if what_to_diff == "magnetization":
            raise NotImplementedError("Zm from magnetization")

        # Build dataframe and add metadata.
        df = pd.DataFrame(rows)
        df.attrs["zeff_name"] = zeff_name

        return df

    def print_eff_charges(self, elements: None | list[str] = None, file=sys.stdout, verbose: int = 0) -> None:
        """
        Print effective charges to `file`. Show only elements in `elements` if not None.
        """
        if elements is not None: elements = list_strings(elements)

        def _p(*args, **kwargs):
            print(*args, file=file, **kwargs)

        if verbose:
            _p("Input structure:")
            _p(self.input_structure)
            _p("")
            _p("Relaxed structure:")
            _p(self.relaxed_structure)
            _p("")

        for iatom, site in enumerate(self.relaxed_structure):
            if elements is not None and site.species_string not in elements: continue
            df = self.get_df_zeff_iatom(iatom)
            zeff_name = df.attrs["zeff_name"]
            _p(f"{zeff_name}[atom_dir, {self.pert_kind}_dir] in Cart. coords for {iatom=}: element: {site.species_string}, frac_coords: {site.frac_coords}")
            _p(df)
            _p("")


@dataclasses.dataclass(kw_only=True)
class DisplData(_FdData):
    """
    """
    # TODO
    #def get_force_constant_df(self):


@dataclasses.dataclass(kw_only=True)
class StrainData(_FdData):
    """
    """
    # TODO
    #def get_elastic_df(self):


@dataclasses.dataclass(kw_only=True)
class ElectricFieldData(_FdData):
    """
    This object stores the dynamical magnetic charges Zm computed with finite differences.
    All values are in a.u. and tensors are in Cartesian coordinates.
    """
    def get_epsinf_df(self) -> pd.Dataframe:
        """
        Return dataframe with the eps_infinity tensor obtained with different FD points.
        """
        eps_comps = list(itertools.product(self.pert_dir_comps, self.pert_dir_comps))
        rows = []
        for npts, eps in self.dpol_dpert_npts.items():
            eps = 4.0 * np.pi * eps
            eps[np.diag_indices_from(eps)] += 1.0
            rows.append(_dict_from_mat_npts(eps, eps_comps, npts))

        return pd.DataFrame(rows)

    def get_piezo_df(self) -> pd.Dataframe:
        """
        Return dataframe with the piezo-electric tensor obtained with different FD points.
        """
        voigt_comps = [str(i) for i in range(1, 7)]
        piezo_comps = list(itertools.product(voigt_comps, self.pert_dir_comps))
        rows = []
        for npts, dstress_dpert in self.dstress_dpert_npts.items():
            # dStress/dPert has shape [3, 3, ip] in Cart. coords.
            piezo = np.empty((3, 6))
            for ip, pert in enumerate(self.perts):
                piezo[ip] = mat33_to_voigt(dstress_dpert[:,:,ip])
            rows.append(_dict_from_mat_npts(piezo, piezo_comps, npts))

        return pd.DataFrame(rows)

    @add_fig_kwargs
    def plot_polarization(self, what: str = "total", **kwargs) -> Figure:
        """
        Plot the polarization as a function of the Electric field.
        """
        nrows, ncols = 3, 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=True)

        # Select the quantity to plot depending on `what`.
        vals_pv = {
            "total": self.cart_pol_pv,
            "electronic": self.cart_pole_pv,
            "ionic": self.cart_poli_pv,
        }[what]

        for pol_dir in range(3):
            ax = ax_list[pol_dir]
            #ax.set_title(f"H_dir: {idir2s(pdir)}, Atom_dir: {idir2s(iat_dir)}", fontsize=fontsize)
            for ip, pert in enumerate(self.perts):
                ax.plot(self.pert_values, vals_pv[ip, :, pol_dir], marker="o")
                        #label=site.species_string + r"$_{\text{%s}}$" % iat)
            #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig


@dataclasses.dataclass(kw_only=True)
class ZeemanData(_FdData):
    """
    This object stores the dynamical magnetic charges Zm computed with finite differences.
    All values are in a.u. and tensors are in Cartesian coordinates.
    """

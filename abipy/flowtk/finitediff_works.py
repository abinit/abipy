# coding: utf-8
"""Work subclasses related to GS calculations."""
from __future__ import annotations

import sys
#import json
import itertools
#import pickle
import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from dataclasses import field
#from monty.json import MSONable
from monty.string import list_strings #, marquee
from monty.functools import lazy_property
from abipy.core.structure import Structure
from abipy.tools.numtools import build_mesh
from abipy.tools.derivatives import central_fdiff_weights, check_num_points_for_order # finite_diff
from abipy.tools.typing import Figure
from abipy.abio.inputs import AbinitInput
from abipy.abio.outputs import BerryPhasePolarization
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


@dataclasses.dataclass
class PertInfo:
    """
    Stores info on the perturbation.
    """
    field_type: str
    cart_dir: np.ndarray
    iatom: int | None = None
    voigt_idx: tuple | None = None

    def __post_init__(self):
        """Validation logic."""
        #if self.field_type == "Displ":
        #if not self.name:
        #    raise ValueError("Name cannot be empty")
        #if self.age < 0:
        #    raise ValueError("Age cannot be negative")

    @lazy_property
    def label(self) -> str:
        return f"{self.name}: {dir2str(self.cart_dir)}"

    @lazy_property
    def dir_str(self) -> str:
        return f"{dir2str(self.cart_dir)}"

    @lazy_property
    def tex(self) -> str:
        return {
            "E": r"\mathcal{E}",
            "H": r"\mathcal{H}",
        }[self.field_type]

    @lazy_property
    def name(self) -> str:
        return {
            "E": "Electric field",
            "H": "Magnetic field",
        }[self.field_type]


class FiniteDisplWork(Work):
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
                       manager=None):
        """
        Build a work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            num_points:
            pert_cart_dirs:
            step_au: Finite difference step for the displacement in Bohr (a.u.)
            pert_cart_dirs:
            mask_iatom:
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        work.scf_input = scf_input.deepcopy()
        structure = scf_input.structure
        natom = len(structure)

        work.step_au = float(step_au)
        work.displ_values, work.ix0 = build_mesh(0.0, num_points, step_au, "=")
        work.num_deltas = len(work.displ_values)

        # Here we normalize the directions to 1. NB: pymatgen structures uses Ang and not Bohr.
        if pert_cart_dirs is not None:
            work.pert_cart_dirs = np.reshape(pert_cart_dirs, (-1, 3))
        else:
            work.pert_cart_dirs = np.eye(3)

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
                work.perts.append(PertInfo(field_type="Displ", cart_dir=cart_dir, iatom=iatom))

        nperts = len(work.perts)
        work.scf_tasks_pv = np.empty((nperts, work.num_deltas), dtype=object)

        for ip, pert in enumerate(work.perts):
            iatom, cart_dir = pert.iatom, pert.cart_dir
            for iv, delta_au in enumerate(work.displ_values):
                new_structure = structure.copy()
                # Note Bohr --> Ang conversion.
                new_structure.translate_sites([iatom], delta_au * abu.Bohr_Ang * cart_dir,
                                               frac_coords=False, to_unit_cell=False)
                new_input = scf_input.new_with_structure(new_structure)
                task = work.register_scf_task(new_input)
                work.scf_tasks_pv[ip, iv] = task

        return work

    @lazy_property
    def natom(self) -> int:
        """Number of atoms in the unit cell."""
        return len(work.scf_input.structure)

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        #data = self.get_data_dict()
        #data.pickle_dump(self.outdir.path)
        return super().on_all_ok()


#class FiniteStrainWork(Work):
#    """
#    Work for the computation at finite Strain with finite differences.
#    """
#    @classmethod
#    def from_scf_input(cls, scf_input, delta=1e-4, voigt_inds=None, ecutsm=0.5, manager=None):
#        """
#        Build a Work from an AbinitInput representing a GS SCF calculation.
#
#        Args:
#           scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
#	        delta:
#           ecutsm: Value of ecutsm input variable. If `scf_input` does not provide ecutsm, this
#               value will be used else the vale in `scf_input`.
#           manager: TaskManager instance. Use default if None.
#        """
#        work = cls(manager=manager)
#	     work.scf_input = scf_input.deepcopy()
#
#        if "ecutsm" not in scf_input:
#            work.scf_input.set_vars(ecutsm=ecutsm)
#            print("Input does not define ecutsm.\n",
#                  "A default value of %s will be added to all the EOS inputs" % ecutsm)

#        work.unstrained_task = work.register_scf_task(work.scf_input)
#
#        if voigt_inds is None:
#           voigt_inds = [(0, 0)]

#        work.pert_values, work.ipv0 = build_mesh(0.0, num_points, step_au, "=")
#        work.tasks_voigt = {}
#
#        for voigt in voigt_inds:
#           work.tasks_voigt[voigt] = defaultdict(list)
#           for isign in (-1, +1):
#           for delta in work.pert_values:
#               # Apply strain to the lattice.
#               strain = np.zeros((3, 3))
#               strain[voigt] = float(isign) * delta
#               new_structure = scf_input.structure.deepcopy()
#               new_structure.apply_strain(strain)
#               new_input = scf_input.new_with_structure(new_structure)
#               # Perform GS calculations with strained cell.
#               task = work.register_scf_task(new_input)
#               work.tasks_voigt[voigt].append(task)
#
#        return work
#
#    def on_all_ok(self):
#        """This method is called when all tasks have reached S_OK."""
#        data = self.get_data()
#        data.pickle_dump(self.outdir.path)
#        return super().on_all_ok()


class _FieldWork(Work):
    """Base class for finite field + finite difference Work."""

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       num_points: int,
                       step_au: float,
                       relax: bool = False,
                       relax_opts: dict | None = None,
                       manager=None):
        """
        Build the work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF calculation used as template to generate the other inputs.
            num_points: Number of points for finite difference.
            step_au: Finite difference step for the magnetic field in a.u.
            relax: False if the initial structural relaxation should not be performed.
            relax_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        work.step_au = step_au
        work.pert_values, work.ipv0 = build_mesh(0.0, num_points, step_au, "=")
        work.num_points = len(work.pert_values)
        check_num_points_for_order(num_points=work.num_points, order=1, kind="=")

        work.pert_cart_dirs = np.array([
            #(1, 1, 1),  # This is the direction used in the tutorial.
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ])
        work.scf_input_template = scf_input.deepcopy()

        if isinstance(work, FiniteEfieldWork):
            work.field_type = "E"
        elif isinstance(work, FiniteHfieldWork):
            work.field_type = "H"
        else:
            raise TypeError(f"Don't know how to handle {type(work)=}")

        # Build list of perturbations.
        work.perts = [PertInfo(field_type=work.field_type, cart_dir=cart_dir) for cart_dir in work.pert_cart_dirs]

        work.relax = relax
        if work.relax:
            relax_input = scf_input.make_relax_input()
            work.initial_relax_task = work.register_relax_task(relax_input)
        else:
            if work.field_type == "E":
                work._add_tasks_with_efield(scf_input.structure)
            elif work.field_type == "H":
                work._add_tasks_with_zeemanfield(scf_input.structure)
            else:
                raise TypeError(f"Don't know how to handle {work.field_type}")

            work.relaxed_structure = scf_input.structure

        return work

    def on_ok(self, sender):
        """This method is called when one task reaches status `S_OK`."""
        if self.relax and sender == self.initial_relax_task:
            # Get relaxed structure from GSR file.
            self.relaxed_structure = sender.get_final_structure()

            if work.field_type == "E":
                self._add_tasks_with_efield(self.relaxed_structure)
            elif work.field_type == "H":
                self._add_tasks_with_zeemanfield(self.relaxed_structure)
            else:
                raise TypeError(f"Don't know how to handle {work.field_type}")

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        return super().on_all_ok()

    def get_data_dict(self) -> dict:
        """
        Note: The suffix `_pv` stands for Field-direction and field-Value.
        """
        natom = len(self[0].input.structure)
        npert, np_vals = len(self.perts), len(self.pert_values)
        has_pol = any(task.input.get("berryopt", 0) != 0 for task in self)

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

        self.tasks_pv = np.empty((len(self.pert_cart_dirs), len(self.pert_values)), dtype=object)
        task_pv0 = None
        for pdir, p_cart_dir in enumerate(self.pert_cart_dirs):
            for ipv, p_val in enumerate(self.pert_values):
                is_pv0 = abs(p_val) < 1e-16
                new_inp = scf_input.new_with_vars(zeemanfield=p_val * p_cart_dir)
                if is_pv0:
                    # Avoid computing H=0 multiple times.
                    if task_pv0 is None:
                        task_pv0 = self.register_scf_task(new_inp)
                    self.tasks_pv[pdir, ipv] = task_pv0
                else:
                    self.tasks_pv[pdir, ipv] = self.register_scf_task(new_inp)

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

        self.tasks_pv = np.empty((len(self.pert_cart_dirs), len(self.pert_values)), dtype=object)
        task_pv0 = None
        for pdir, p_cart_dir in enumerate(self.pert_cart_dirs):
            for ipv, p_val in enumerate(self.pert_values):
                is_pv0 = abs(p_val) < 1e-16
                new_inp = scf_input.new_with_vars(efield=p_val * p_cart_dir)
                if is_pv0:
                    # Avoid computing the zero-field case multiple times.
                    # Also the task at zero field uses berryopt -1.
                    if task_pv0 is None:
                        new_inp.set_vars(berryopt=-1)
                        task_pv0 = self.register_berry_task(new_inp)
                    self.tasks_pv[pdir, ipv] = task_pv0
                else:
                    new_inp.set_vars(berryopt=4)
                    self.tasks_pv[pdir, ipv] = self.register_berry_task(new_inp)

        # Now add dependencies to the tasks.
        for pdir, p_cart_dir in enumerate(self.pert_cart_dirs):
            for ipv in range(0, self.ipv0):
                self.tasks_pv[pdir, ipv].add_deps({self.tasks_pv[pdir, ipv+1]: "WFK"})
            ntasks = len(self.tasks_pv[pdir])
            for ipv in range(self.ipv0+1, ntasks):
                self.tasks_pv[pdir, ipv].add_deps({self.tasks_pv[pdir, ipv-1]: "WFK"})

    def get_data(self) -> ElectricFieldData:
        d = self.get_data_dict()
        return ElectricFieldData(**d)



def idir2s(idir: int):
    """Convert direction index to string."""
    return {0: "x", 1: "y", 2: "z"}[idir]


@dataclasses.dataclass(kw_only=True)
class _FiniteFieldDataMixin(HasPickleIO):

    input_structure: Structure
    relaxed_structure: Structure
    step_au: float
    ipv0: int
    has_pol: bool

    perts: list[PertInfo]
    pert_values: np.ndarray
    pert_cart_dirs: np.ndarray
    etotals_pv: np.ndarray
    eterms_pv: np.ndarray
    cart_forces_pv: np.ndarray
    carts_stresses_pv: np.ndarray
    params_p: list[dict]

    # Mapping npts -> dForce/dPert  [iat, 3, pdir] in Cart. coords.
    dforces_dpert_npts: dict[int, np.array] = field(init=False)
    # Mapping npts -> dStress/dPert [3, 3, pdir] in Cart. coords.
    dstress_dpert_npts: dict[int, np.array] = field(init=False)

    # Redefined in the subclasses.
    pert_type: str = "None"
    pert_name: str = "None"
    pert_tex: str = "None"

    def __post_init__(self):
        """
        Compute quantities with finite differences.
        """
        natom = len(self.input_structure)
        npert, np_vals = len(self.pert_cart_dirs), len(self.pert_values)

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
            for iat, iat_dir, pdir in itertools.product(range(natom), range(3), range(npert)):
                fvals_f = self.cart_forces_pv[pdir, :, iat, iat_dir]
                dforce_dpert[iat, iat_dir, pdir] = np.sum(fvals_f[fd_slice] * weights) / self.step_au
            self.dforces_dpert_npts[npts] = dforce_dpert

            # Finit differences for stresses.
            dstress_dpert = np.empty((3, 3, npert))
            for ii, jj, pdir in itertools.product(range(3), range(3), range(npert)):
                svals_f = self.carts_stresses_pv[pdir, :, ii, jj] * self.relaxed_structure.volume # * abu.Angs2Bohr ** 3  TODO?
                dstress_dpert[ii, jj, pdir] = np.sum(svals_f[fd_slice] * weights) / self.step_au
            self.dstress_dpert_npts[npts] = dstress_dpert

            # Finite differences for polarizations (if available).
            if self.has_pol:
                dpol_dpert = np.empty((3, npert))
                for pdir, ii in itertools.product(range(npert), range(3)):
                    dpol_dpert[ii, pdir] = np.sum(self.cart_pol_pv[pdir, fd_slice, ii] * weights) / self.step_au
                self.dpol_dpert_npts[npts] = dpol_dpert
                #self.dpole_dpert_npts[npts] = dpole_dpert  TODO ?
                #self.dpoli_dpert_npts[npts] = dpoli_dpert  TODO ?

    @lazy_property
    def npert(self) -> int:
        """Number of field directions."""
        return len(self.perts)

    @lazy_property
    def natom(self) -> int:
        """Numbef of atoms in the unit cell."""
        return len(self.input_structure)

    def label_for_pdir(self, pdir: int) -> str:
        pdir_str = dir2str(self.pert_cart_dirs[pdir])
        return "$%s_{%s}$ (a.u.)" % (self.pert_tex, pdir_str)

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
            ax = ax_mat[ip, 0]
            ax.set_title(pert.label, fontsize=fontsize)
            for ii, jj in itertools.product(range(3), range(3)):
                ax.plot(self.pert_values, self.carts_stresses_pv[ip, :, ii, jj], marker="o",
                        label=r"$\sigma_{%s}$" % (f"{ii}{jj}"))

            ax.legend(loc="best", fontsize=fontsize, shadow=True)
            #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig

    def get_df_iatom(self, iatom: int) -> pd.Dataframe:
        """
        Return dataframe with effective charges for the given atom index and all the FD points.
        """
        force_comps = "x y z".split()
        field_comps = [pert.dir_str for pert in self.perts]
        comps = list(itertools.product(force_comps, field_comps))
        rows = []
        for npts, dforces_dpert in self.dforces_dpert_npts.items():
            zeff = dforces_dpert[iatom]
            d = {"npts": npts}
            d.update({c: v for c, v in zip(comps, zeff.flatten(), strict=True)})
            if zeff.shape == (3, 3):
                d["isoavg"], d["det"] = np.trace(zeff) / 3, np.linalg.det(zeff)
            rows.append(d)

        return pd.DataFrame(rows)

    def print_eff_charges(self, elements: None | list[str] = None, file=sys.stdout) -> None:
        """
        Print effective charges to `file`. Show only elements in `elements` if not None.
        """
        if elements is not None: elements = list_strings(elements)
        zeff_str = {"E": "Ze", "H": "Zm"}[self.pert_type]

        def _p(*args, **kwargs):
            print(*args, file=file, **kwargs)

        _p("Input structure:")
        _p(self.input_structure)
        _p("")
        _p("Relaxed structure:")
        _p(self.relaxed_structure)
        _p("")

        for iatom, site in enumerate(self.relaxed_structure):
            if elements is not None and site.species_string not in elements: continue
            df = self.get_df_iatom(iatom)
            _p(f"{zeff_str}[atom_dir, {self.pert_type}_dir] in Cart. coords for {iatom=}: element: {site.species_string}, frac_coords: {site.frac_coords}")
            _p(df)
            _p("")


@dataclasses.dataclass(kw_only=True)
class ElectricFieldData(_FiniteFieldDataMixin):
    """
    This object stores the dynamical magnetic charges Zm computed with finite differences.
    All values are in a.u. and Ze are in Cartesian coordinates.

    _df stands for Direction, Field
    """
    cart_pol_pv: np.ndarray
    cart_pole_pv: np.ndarray
    cart_poli_pv: np.ndarray
    dpol_dpert_npts: dict[int, np.array] = field(init=False)

    pert_type: str = "E"
    pert_name: str = "Electric field"
    pert_tex: str = r"\mathcal{E}"

    def get_epsinf_df(self) -> pd.Dataframe:
        """
        Return dataframe with eps_infinity obtained with different FD points.
        """
        comps = "x y z".split()
        comps = list(itertools.product(comps, comps))
        rows = []
        for npts, eps in self.dpol_dpert_npts.items():
            eps = 4 * np.pi * eps
            eps[np.diag_indices_from(eps)] += 1.0
            d = {"npts": npts}
            d.update({c: v for c, v in zip(comps, eps.flatten(), strict=True)})
            if eps.shape == (3, 3):
                d["isoavg"], d["det"] = np.trace(eps) / 3 , np.linalg.det(eps)
            rows.append(d)

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
class ZeemanData(_FiniteFieldDataMixin):
    """
    This object stores the dynamical magnetic charges Zm computed with finite differences.
    All values are in a.u. and Ze are in Cartesian coordinates.
    """
    pert_type: str = "H"
    pert_name: str = "Magnetic field"
    pert_tex: str = r"\mathcal{H}"

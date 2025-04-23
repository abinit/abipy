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


@dataclasses.dataclass(kw_only=True)
class FiniteDisplData(HasPickleIO):

    #deltas: np.ndarray
    ix0: int
    #structures_fv: np.ndarray
    #etotals_fv: np.ndarray
    #cart_forces_fv: np.ndarray
    #cart_stresses_fv: np.ndarray


class FiniteDisplWork(Work):
    """
    Work for the computation of forces with finite differences.
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       num_points: int,
                       step_au: float = 0.01,
                       displ_cart_dirs=None,
                       mask_iatom=None,
                       manager=None):
        """
        Build a work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            num_points:
            displ_cart_dirs:
            step_au: Finite difference step for the displacement in Bohr (a.u.)
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

        # Here we normalize the directions to 1 (NB: pymatgen structures uses Ang and not Bohr.
        if displ_cart_dirs is not None:
            work.displ_cart_dirs = np.reshape(displ_cart_dirs, (-1, 3))
        else:
            work.displ_cart_dirs = np.eye(3)

        for idir, cart_dir in enumerate(work.displ_cart_dirs):
            norm = structure.lattice.norm(cart_dir, frac_coords=False)
            work.displ_cart_dirs[idir] = cart_dir / norm

        if mask_iatom is None:
            mask_iatom = np.ones(natom, dtype=bool)

        work.mask_iatom = np.array(mask_iatom).astype(bool)
        if len(work.mask_iatom) != natom:
            raise ValueError(f"{len(work.mask_iatom)=} != {natom=}")

        work.num_dirs = len(work.displ_cart_dirs)
        work.scf_tasks_dav = np.empty((work.num_dirs, natom, work.num_deltas), dtype=object)

        for iatom, mask in zip(range(natom), work.mask_iatom, strict=True):
            if not mask: continue
            for idir, cart_dir in enumerate(work.displ_cart_dirs):
                for iv, delta_au in enumerate(work.displ_values):
                    new_structure = structure.copy()
                    # Note Bohr --> Ang.
                    new_structure.translate_sites([iatom], delta_au * abu.Bohr_Ang * cart_dir,
                                                   frac_coords=False, to_unit_cell=False)
                    new_input = scf_input.new_with_structure(new_structure)
                    task = work.register_scf_task(new_input)
                    work.scf_tasks_dav[idir, iatom, iv] = task

        return work

    @lazy_property
    def natom(self) -> int:
        """Number of atoms in the unit cell."""
        return len(work.scf_input.structure)

    def get_data_dict(self):
        """
        Read data from the GSR files
        """
        natom, nf_dirs, nf_vals = self.natom, len(self.field_cart_dirs), len(self.field_values)

        data = {
            #"input_structure": self[0].input.structure,
            #"relaxed_structure": self.relaxed_structure,
            "step_au": self.step_au,
            "field_values": self.field_values,
            "field_cart_dirs": self.field_cart_dirs,
            "if0": self.if0,
        }

        data["params_f"] = []
        data["etotals_fv"] = etotals_fv = np.empty((nf_dirs, nf_vals))
        data["eterms_fv"] = eterms_fv = np.empty((nf_dirs, nf_vals), dtype=object)
        data["cart_forces_fv"] = cart_forces_fv = np.empty((nf_dirs, nf_vals, natom, 3))
        data["cart_stresses_fv"] = cart_stresses_fv = np.empty((nf_dirs, nf_vals, 3, 3))

        if has_pol := any(task.input.get("berryopt", 0) != 0 for task in self):
            data["cart_pol_fv"] = cart_pol_fv = np.empty((nf_dirs, nf_vals, 3), dtype=object)
            data["cart_pole_fv"] = cart_pole_fv = np.empty((nf_dirs, nf_vals, 3), dtype=object)
            data["cart_poli_fv"] = cart_poli_fv = np.empty((nf_dirs, nf_vals, 3), dtype=object)

        # Read energy, forces and stress from the GSR files.
        for iatom, mask in zip(range(self.natom), self.mask_iatom, strict=True):
            if not mask: continue
            for idir, cart_dir in enumerate(work.cart_dirs):
               for iv, delta in enumerate(work.deltas):
                   task = work.scf_tasks_dav[idir, iatom, iv]
                   with task.open_gsr() as gsr:
                       etotals_fv[fdir, ifv] = gsr.r.read_value("etotal")
                       eterms_fv[fdir, ifv] = gsr.r.read_energy_terms(unit="Ha")
                       cart_forces_fv[fdir, ifv] = gsr.r.read_value("cartesian_forces") # Ha/Bohr units.
                       cart_stresses_fv[fdir, ifv] = gsr.r.read_cart_stress_tensor(units="au")
                       # Add parameters that might be used for convergence studies.
                       data["params_f"].append(gsr.params)

                       if has_pol:
                           with task.open_abo() as abo:
                               pol = abo.get_berry_phase_polarization()
                               cart_pol_fv[fdir, ifv] = pol.total
                               cart_pole_fv[fdir, ifv] = pol.electronic
                               cart_poli_fv[fdir, ifv] = pol.ionic

        # Use all stencils compatible with input num_points so that we can monitor the convergence afterwards.
        data["dforces_dfield_npts"] = {}
        data["dstress_dfield_npts"] = {}
        if has_pol:
            data["dpol_dfield_npts"] = {}

        for acc, weights in central_fdiff_weights[1].items():
            if self.num_points < len(weights): continue
            nn = acc // 2
            fd_slice = slice(self.if0-nn, self.if0+nn+1)
            npts = len(weights)

            # FD for forces.
            dforce_dfield = np.empty((natom, 3, nf_dirs))
            for iat, iat_dir, fdir in itertools.product(range(natom), range(3), range(nf_dirs)):
                fvals_f = cart_forces_fv[fdir, :, iat, iat_dir]
                dforce_dfield[iat, iat_dir, fdir] = np.sum(fvals_f[fd_slice] * weights) / self.step_au
            data["dforces_dfield_npts"][npts] = dforce_dfield

            # FD for stresses.
            dstress_dfield = np.empty((3, 3, nf_dirs))
            for ii, jj, fdir in itertools.product(range(3), range(3), range(nf_dirs)):
                svals_f = cart_stresses_fv[fdir, :, ii, jj] * self.relaxed_structure.volume # * abu.Angs2Bohr ** 3
                dstress_dfield[ii, jj, fdir] = np.sum(svals_f[fd_slice] * weights) / self.step_au
            data["dstress_dfield_npts"][npts] = dstress_dfield

            # FD for polarizations (if available)
            if has_pol:
                dpol_dfield = np.empty((3, nf_dirs))
                for fdir, ii in itertools.product(range(nf_dirs), range(3)):
                    dpol_dfield[ii, fdir] = np.sum(cart_pol_fv[fdir, fd_slice, ii] * weights) / self.step_au
                data["dpol_dfield_npts"][npts] = dpol_dfield
                #data["dpole_dfield_npts"][npts] = dpole_dfield  TODO ?
                #data["dpoli_dfield_npts"][npts] = dpoli_dfield  TODO ?

        #return FiniteDisplData(**data)
        return data

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        return super().on_all_ok()


#class FiniteStrainWork(Work):
#    """
#    Work for the computation of the stress tensor with finite differences
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
#
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

#        work.delta = float(delta)
#        work.tasks_voigt = {}
#
#        for voigt in voigt_inds:
#           work.tasks_voigt[voigt] = defaultdict(list)
#           for isign in (-1, +1):
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
#    def get_data_dict(self):
#        """
#	     It reads the energies and the volumes from the GSR files
#        """
#        with self.unstrained_task.open_gsr() as gsr0:
#            e0, v0 = gsr0.energy, gsr0.structure.volume
#            cart_stress_tensor = gsr0.cart_stress_tensor
#
#        fd_tensor = np.empty((3, 3))
#        for voigt in voigt_inds:
#            tasks = self.tasks_voigt[voigt]
#            energies_ev, volumes = np.empty(len(tasks), np.empty(len(tasks))
#            for i, task in enumerate(task):
#                    with task.open_gsr() as gsr:
#                       energies_ev[i] = float(gsr.energy))
#
#            d = (e_plus - e_minus) / (2 * self.delta * v0)
#                fd_tensor[voigt] = d
#                fd_tensor[voigt[1], voigt[0]] = d
#
#        data = {
#                 "dfpt_cart_stress_tensor": cart_stress_tensor,
#                 "finite_diff_cart_stress_tensor": fd_cart_stress_tensor,
#            }
#
#            with open(self.outdir.path_in("stress.json"), "wt") as fh:
#                json.dump(data, fh, indent=4, sort_keys=True)
#
#            return data
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
        work.field_values, work.if0 = build_mesh(0.0, num_points, step_au, "=")
        work.num_points = len(work.field_values)
        check_num_points_for_order(num_points=work.num_points, order=1, kind="=")

        #nspinor = scf_input.get("nspinor", 1)
        #if nspinor != 2:
        #    raise ValueError("nspinor should be 2 to have non-zero dyn magnetic charges while it is {nspinor}")

        work.field_cart_dirs = np.array([
            #(1, 1, 1),  # This is the direction used in the tutorial.
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ])
        work.scf_input_template = scf_input.deepcopy()

        work.relax = relax
        if work.relax:
            relax_input = scf_input.make_relax_input()
            work.initial_relax_task = work.register_relax_task(relax_input)
        else:
            if isinstance(work, FiniteEfieldWork):
                work._add_tasks_with_efield(scf_input.structure)
            elif isinstance(work, FiniteHfieldWork):
                work._add_tasks_with_zeemanfield(scf_input.structure)
            else:
                raise TypeError(f"Don't know how to handle {type(work)}")

            work.relaxed_structure = scf_input.structure

        return work

    def on_ok(self, sender):
        """This method is called when one task reaches status `S_OK`."""
        if self.relax and sender == self.initial_relax_task:
            # Get relaxed structure from GSR file.
            self.relaxed_structure = sender.get_final_structure()

            if isinstance(work, FiniteEfieldWork):
                self._add_tasks_with_efield(self.relaxed_structure)
            elif isinstance(work, FiniteHfieldWork):
                self._add_tasks_with_zeemanfield(self.relaxed_structure)
            else:
                raise TypeError(f"Don't know how to handle {type(work)}")

            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        return super().on_all_ok()

    def get_data_dict(self) -> dict:
        """

        Note: The suffix `_fv` stands for Field-direction and field-Value.
        """
        natom = len(self[0].input.structure)
        nf_dirs, nf_vals = len(self.field_cart_dirs), len(self.field_values)

        data = {
            "input_structure": self[0].input.structure,
            "relaxed_structure": self.relaxed_structure,
            "step_au": self.step_au,
            "field_values": self.field_values,
            "field_cart_dirs": self.field_cart_dirs,
            "if0": self.if0,
        }

        data["params_f"] = []
        data["etotals_fv"] = etotals_fv = np.empty((nf_dirs, nf_vals))
        data["eterms_fv"] = eterms_fv = np.empty((nf_dirs, nf_vals), dtype=object)
        data["cart_forces_fv"] = cart_forces_fv = np.empty((nf_dirs, nf_vals, natom, 3))
        data["cart_stresses_fv"] = cart_stresses_fv = np.empty((nf_dirs, nf_vals, 3, 3))

        if has_pol := any(task.input.get("berryopt", 0) != 0 for task in self):
            data["cart_pol_fv"] = cart_pol_fv = np.empty((nf_dirs, nf_vals, 3))
            data["cart_pole_fv"] = cart_pole_fv = np.empty((nf_dirs, nf_vals, 3))
            data["cart_poli_fv"] = cart_poli_fv = np.empty((nf_dirs, nf_vals, 3))

        # Read energy, forces and stress from the GSR files.
        for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
            for ifv, f_val in enumerate(self.field_values):
                task = self.tasks_fdir_f[fdir, ifv]
                with task.open_gsr() as gsr:
                    etotals_fv[fdir, ifv] = gsr.r.read_value("etotal")
                    eterms_fv[fdir, ifv] = gsr.r.read_energy_terms(unit="Ha")
                    cart_forces_fv[fdir, ifv] = gsr.r.read_value("cartesian_forces") # Ha/Bohr units.
                    cart_stresses_fv[fdir, ifv] = gsr.r.read_cart_stress_tensor(units="au")
                    # Add parameters that might be used for convergence studies.
                    data["params_f"].append(gsr.params)

                if has_pol:
                    with task.open_abo() as abo:
                        pol = abo.get_berry_phase_polarization()
                        cart_pol_fv[fdir, ifv] = pol.total
                        cart_pole_fv[fdir, ifv] = pol.electronic
                        cart_poli_fv[fdir, ifv] = pol.ionic

        # Use all stencils compatible with input num_points so that we can monitor the convergence afterwards.
        data["dforces_dfield_npts"] = {}
        data["dstress_dfield_npts"] = {}
        if has_pol:
            data["dpol_dfield_npts"] = {}

        for acc, weights in central_fdiff_weights[1].items():
            if self.num_points < len(weights): continue
            nn = acc // 2
            fd_slice = slice(self.if0-nn, self.if0+nn+1)
            npts = len(weights)

            # FD for forces.
            dforce_dfield = np.empty((natom, 3, nf_dirs))
            for iat, iat_dir, fdir in itertools.product(range(natom), range(3), range(nf_dirs)):
                fvals_f = cart_forces_fv[fdir, :, iat, iat_dir]
                dforce_dfield[iat, iat_dir, fdir] = np.sum(fvals_f[fd_slice] * weights) / self.step_au
            data["dforces_dfield_npts"][npts] = dforce_dfield

            # FD for stresses.
            dstress_dfield = np.empty((3, 3, nf_dirs))
            for ii, jj, fdir in itertools.product(range(3), range(3), range(nf_dirs)):
                svals_f = cart_stresses_fv[fdir, :, ii, jj] * self.relaxed_structure.volume # * abu.Angs2Bohr ** 3  TODO?
                dstress_dfield[ii, jj, fdir] = np.sum(svals_f[fd_slice] * weights) / self.step_au
            data["dstress_dfield_npts"][npts] = dstress_dfield

            # FD for polarizations (if available)
            if has_pol:
                dpol_dfield = np.empty((3, nf_dirs))
                for fdir, ii in itertools.product(range(nf_dirs), range(3)):
                    dpol_dfield[ii, fdir] = np.sum(cart_pol_fv[fdir, fd_slice, ii] * weights) / self.step_au
                data["dpol_dfield_npts"][npts] = dpol_dfield
                #data["dpole_dfield_npts"][npts] = dpole_dfield  TODO ?
                #data["dpoli_dfield_npts"][npts] = dpoli_dfield  TODO ?

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

        self.tasks_fdir_f = np.empty((len(self.field_cart_dirs), len(self.field_values)), dtype=object)
        task_f0 = None
        for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
            for ifv, f_val in enumerate(self.field_values):
                is_f0 = abs(f_val) < 1e-16
                new_inp = scf_input.new_with_vars(zeemanfield=f_val * f_cart_dir)
                if is_f0:
                    # Avoid computing H=0 multiple times.
                    if task_f0 is None:
                        task_f0 = self.register_scf_task(new_inp)
                    self.tasks_fdir_f[fdir, ifv] = task_f0
                else:
                    self.tasks_fdir_f[fdir, ifv] = self.register_scf_task(new_inp)

    def get_data(self) -> ZeemanData:
        """
        Read data from the GSR files, and compute Zm with finite differences.
        """
        d = super().get_data_dict()
        return ZeemanData(**d)


def idir2s(idir: int):
    """Convert direction index to string."""
    return {0: "x", 1: "y", 2: "z"}[idir]


@dataclasses.dataclass(kw_only=True)
class _FiniteFieldDataMixin(HasPickleIO):

    input_structure: Structure
    relaxed_structure: Structure
    step_au: float
    if0: int
    field_values: np.ndarray
    field_cart_dirs: np.ndarray
    etotals_fv: np.ndarray
    eterms_fv: np.ndarray
    cart_forces_fv: np.ndarray
    cart_stresses_fv: np.ndarray
    dforces_dfield_npts: dict[int, np.array]  # Mapping npts -> dForce/dField  [iat, 3, fdir] in Cart. coords.
    dstress_dfield_npts: dict[int, np.array]  # Mapping npts -> dStress/dField [3, 3, fdir] in Cart. coords.
    params_f: list[dict]

    # Redefined in the subclasses.
    field_type: str = "None"
    field_name: str = "None"
    field_tex: str = "None"

    @lazy_property
    def nf_dirs(self) -> int:
        """Number of field directions."""
        return len(self.field_cart_dirs)

    @lazy_property
    def natom(self) -> int:
        """Numbef of atoms in the unit cell."""
        return len(self.input_structure)

    def label_for_fdir(self, fdir: int) -> str:
        fdir_str = dir2str(self.field_cart_dirs[fdir])
        return "$%s_{%s}$ (a.u.)" % (self.field_tex, fdir_str)

    @add_fig_kwargs
    def plot_etotal(self, mode="diff", ax=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot energies as a function of the finite external field for all the directions.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
            e_values = self.etotals_fv[fdir] * abu.Ha_meV / self.natom
            if mode == "diff":
                e_values -= e_values[self.if0]
            ax.plot(self.field_values, e_values, marker="o", label=f"H_dir: {idir2s(fdir)}")

        set_grid_legend(ax, fontsize, xlabel=f"${self.field_tex}$ (a.u.)",
                        ylabel=r"$\Delta$ Energy/atom (meV)" if mode == "diff" else "Energy/atom (meV)")
        return fig

    @add_fig_kwargs
    def plot_forces(self, elements=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot Cartesian forces as a function of the finite external field.
        """
        if elements is not None: elements = list_strings(elements)
        nrows, ncols = 3, self.nf_dirs
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        for iat_dir in range(3):
            for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
                ax = ax_mat[iat_dir, fdir]
                # TODO H_dir is wrong now
                ax.set_title(f"H_dir: {idir2s(fdir)}, Atom_dir: {idir2s(iat_dir)}", fontsize=fontsize)
                for iat, site in enumerate(self.relaxed_structure):
                    if elements is not None and site.species_string not in elements: continue
                    ax.plot(self.field_values, self.cart_forces_fv[fdir, :, iat, iat_dir], marker="o",
                            label=site.species_string + r"$_{\text{%s}}$" % iat)

                ax.legend(loc="best", fontsize=fontsize, shadow=True)
                #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig

    @add_fig_kwargs
    def plot_stresses(self, fontsize=8, **kwargs) -> Figure:
        """
        Plot Cartesian stresses as a function of the finite external field for all the directions.
        """
        nrows, ncols = self.nf_dirs, 1
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
            ax = ax_mat[fdir, 0]
            ax.set_title(f"H_dir: {idir2s(fdir)}", fontsize=fontsize)
            for ii, jj in itertools.product(range(3), range(3)):
                ax.plot(self.field_values, self.cart_stresses_fv[fdir, :, ii, jj], marker="o",
                        label=r"$\sigma_{%s}$" % (f"{ii}, {jj}"))

            ax.legend(loc="best", fontsize=fontsize, shadow=True)
            #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig

    def get_df_iatom(self, iatom: int) -> pd.Dataframe:
        """
        Return dataframe with effective charges for the given atom index and all the FD points.
        """
        force_comps = "x y z".split()
        field_comps = [dir2str(cart_dir) for cart_dir in self.field_cart_dirs]
        comps = list(itertools.product(force_comps, field_comps))
        rows = []
        for npts, dforces_dfield in self.dforces_dfield_npts.items():
            zeff = dforces_dfield[iatom]
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
        zeff_str = {"E": "Ze", "H": "Zm"}[self.field_type]

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
            _p(f"{zeff_str}[atom_dir, {self.field_type}_dir] in Cart. coords for {iatom=}: element: {site.species_string}, frac_coords: {site.frac_coords}")
            _p(df)
            _p("")


@dataclasses.dataclass(kw_only=True)
class ElectricFieldData(_FiniteFieldDataMixin):
    """
    This object stores the dynamical magnetic charges Zm computed with finite differences.
    All values are in a.u. and Ze are in Cartesian coordinates.

    _df stands for Direction, Field
    """
    cart_pol_fv: np.ndarray
    cart_pole_fv: np.ndarray
    cart_poli_fv: np.ndarray
    dpol_dfield_npts: dict[int, np.array]

    field_type: str = "E"
    field_name: str = "Electric field"
    field_tex: str = r"\mathcal{E}"

    def get_eps_df(self) -> pd.Dataframe:
        """
        Return dataframe with eps_infinity obtained with different FD points.
        """
        comps = "x y z".split()
        comps = list(itertools.product(comps, comps))
        rows = []
        for npts, eps in self.dpol_dfield_npts.items():
            eps = 4 * np.pi * eps
            eps[np.diag_indices_from(eps)] += 1.0
            d = {"npts": npts}
            d.update({c: v for c, v in zip(comps, eps.flatten(), strict=True)})
            if eps.shape == (3, 3):
                d["isoavg"], d["det"] = np.trace(eps) / 3 , np.linalg.det(eps)
            rows.append(d)

        return pd.DataFrame(rows)

    @add_fig_kwargs
    def plot_polarization(self, what="total", **kwargs) -> Figure:
        """
        Plot the polarization as a function of the Electric field.
        """
        nrows, ncols = 3, 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=True)

        # Select the quantity to plot depending on `what`.
        vals_fv = {
            "total": self.cart_pol_fv,
            "electronic": self.cart_pole_fv,
            "ionic": self.cart_poli_fv,
        }[what]

        for pol_dir in range(3):
            ax = ax_list[pol_dir]
            #ax.set_title(f"H_dir: {idir2s(fdir)}, Atom_dir: {idir2s(iat_dir)}", fontsize=fontsize)
            for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
                ax.plot(self.field_values, vals_fv[fdir, :, pol_dir], marker="o")
                        #label=site.species_string + r"$_{\text{%s}}$" % iat)
            #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig


@dataclasses.dataclass(kw_only=True)
class ZeemanData(_FiniteFieldDataMixin):
    """
    This object stores the dynamical magnetic charges Zm computed with finite differences.
    All values are in a.u. and Ze are in Cartesian coordinates.
    """
    field_type: str = "H"
    field_name: str = "Magnetic field"
    field_tex: str = r"\mathcal{H}"


class FiniteEfieldWork(_FieldWork):
    r"""
    """

    def _add_tasks_with_efield(self, structure: Structure) -> None:
        """Build new GS tasks with finite electric field."""
        scf_input = self.scf_input_template.new_with_structure(structure)

        self.tasks_fdir_f = np.empty((len(self.field_cart_dirs), len(self.field_values)), dtype=object)
        task_f0 = None
        for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
            for ifv, f_val in enumerate(self.field_values):
                is_f0 = abs(f_val) < 1e-16
                new_inp = scf_input.new_with_vars(efield=f_val * f_cart_dir)
                if is_f0:
                    # Avoid computing the zero-field case multiple times.
                    # Also the task at zero field uses berryopt -1.
                    if task_f0 is None:
                        new_inp.set_vars(berryopt=-1)
                        task_f0 = self.register_berry_task(new_inp)
                    self.tasks_fdir_f[fdir, ifv] = task_f0
                else:
                    new_inp.set_vars(berryopt=4)
                    self.tasks_fdir_f[fdir, ifv] = self.register_berry_task(new_inp)

        # Now add dependencies to the tasks.
        for fdir, f_cart_dir in enumerate(self.field_cart_dirs):
            ntasks = len(self.tasks_fdir_f[fdir])
            for ifv in range(0, self.if0):
                deps = {self.tasks_fdir_f[fdir, ifv+1]: "WFK"}
                self.tasks_fdir_f[fdir, ifv].add_deps(deps)
            for ifv in range(self.if0+1, ntasks):
                deps = {self.tasks_fdir_f[fdir, ifv-1]: "WFK"}
                self.tasks_fdir_f[fdir, ifv].add_deps(deps)

    def get_data(self) -> ElectricFieldData:
        d = self.get_data_dict()
        return ElectricFieldData(**d)

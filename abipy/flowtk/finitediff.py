# coding: utf-8
"""
This module provide Works for finite difference calculations and related post-processing tools.

IMPORTANT: In Abinit Stress is equal to dE/d_strain * (1/ucvol) See m_forstr.F90.
"""
from __future__ import annotations

# Handle
#=   KILLED BY SIGNAL: 9 (Killed)
# slurmstepd: error: Detected 1 oom_kill event in StepId=8214672.0. Some of the step tasks have been OOM Killed.
#srun: error: cns264: task 1: Out Of Memory

import sys
import pickle
import itertools
import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from dataclasses import field
from typing import Optional
from io import StringIO
from pathlib import Path
from monty.string import list_strings #, marquee
from monty.functools import lazy_property
from monty.termcolor import cprint
from pymatgen.analysis.elasticity.strain import Strain
from abipy.core.structure import Structure
#from abipy.core.mixins import NotebookWriter #
from abipy.tools.numtools import build_mesh
from abipy.tools.derivatives import central_fdiff_weights
from abipy.tools.tensors import DielectricDataList
from abipy.tools import duck
from abipy.tools.typing import Figure
from abipy.abio.inputs import AbinitInput
from abipy.abio.enums import StrEnum
from abipy.tools.serialization import HasPickleIO
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_grid_legend,
    rotate_ticklabels, set_visible, set_ax_xylabels, linear_fit_ax, quadratic_fit_ax)
#from abipy.tools.serialization import mjson_write #, pmg_serialize
from .works import Work
#from .flows import Flow

#def centered_indices(n):
#    half = n // 2
#    if n % 2 == 0:
#        return list(range(-half, half)), half
#    #else:
#    #return list(range(-half, half + 1)),

#VOIGT_TO_TUPLE = {
#    0: (0, 0),
#    1: (1, 1),
#    2: (2, 2),
#    3: (0, 1),
#    4: (0, 2),
#    5: (1, 2),
#}


NORMAL_STRAIN_INDS = [0, 1, 2]

SHEAR_STRAIN_INDS = [3, 4, 5]

ALL_STRAIN_INDS = NORMAL_STRAIN_INDS + SHEAR_STRAIN_INDS


def _mesh_for_fd_accuracy(acc, order, step) -> tuple:
    num_points = len(central_fdiff_weights[order][acc]) // 2
    values, ip0 = build_mesh(0.0, num_points, step, "=")
    return num_points, values, ip0


def vec2str(vec, variables: str = 'xyz') -> str:
    """
    >>> vec2str((1, 2, 3)))
    "x + 2y + 3z"
    >>> vec2str((0, -1, 4)))
    "-y + 4z"
    >>> vec2str((1, 0, 0)))
    "x"
    >>> vec2str((0, 0, 0)))   # Output: (empty string)
    ""
    """
    terms = []
    for coeff, var in zip(vec, variables, strict=True):
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
    Convert a 3x3 symmetric matrix to a Voigt vector (6x1).

    Parameters:
        mat: 3x3 symmetric matrix
        engineering_strain: If True, shear components are multiplied by 2 (engineering strain convention)

    Returns: 6-element Voigt vector
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
        voigt: 6-element vector
        engineering_strain: If True, shear components are divided by 2 (engineering strain convention)

    Returns: 3x3 symmetric matrix
    """
    mat = np.zeros((3, 3))
    mat[0, 0] = voigt[0]
    mat[1, 1] = voigt[1]
    mat[2, 2] = voigt[2]
    factor = 0.5 if engineering_strain else 1
    mat[1, 2] = mat[2, 1] = factor * voigt[3]
    mat[0, 2] = mat[2, 0] = factor * voigt[4]
    mat[0, 1] = mat[1, 0] = factor * voigt[5]
    return mat



@dataclasses.dataclass(kw_only=True)
class _FdData(HasPickleIO):
    """
    Base class storing energies, forces, stresses for the different perturbed configurations.
    Provides methods to visualize the results and compute tensors.
    Except for the Structure that uses Angstrom, all values are in a.u. and tensors
    are given in Cartesian coordinates.
    """
    ions_mode: str
    initial_structure: Structure

    has_pol: bool                     # True if polarization has been computed with Berry phase.
    has_mag: bool                     # True if magnetization has been computed.
    perts: list[Perturbation]         # List of perturbations.
    params_pv: np.ndarray             # Parameters used for each perturbation. Each entry is a dict.

    # The `_pv` suffix stands for perturbation and perturbation value.
    structures_pv: np.ndarray         # Array with structures: shape (npert, np_vals)
    etotals_pv: np.ndarray            # (npert, np_vals)
    eterms_pv: np.ndarray             # (npert, np_vals)
    cart_forces_pv: np.ndarray        # (npert, np_vals, natom, 3)
    carts_stresses_pv: np.ndarray     # (npert, np_vals, 6) Voigt form

    # Macro Polarization computed with Berry phase
    cart_pol_pv: Optional[np.ndarray] = None      # (npert, np_vals, 3)
    cart_pole_pv: Optional[np.ndarray] = None     # (npert, np_vals, 3)
    cart_poli_pv: Optional[np.ndarray] = None     # (npert, np_vals, 3)

    # Magnetization (spin part).
    cart_mag_pv: Optional[np.ndarray] = None      # (npert, np_vals, 3)

    # npts -> dForce/dPert with shape (natom, 3, npert) in Cart. coords.
    dforces_dpert_npts: dict[int, np.array] = field(init=False)

    # npts -> dStress/dPert with shape (6, npert) in Cart. coords (Voigt shape)
    dstress_dpert_npts: dict[int, np.array] = field(init=False)

    # npts -> dPol/dPert with shape (3, npert) in Cart. coords.
    dpol_dpert_npts: dict[int, np.array] = field(init=False)

    # npts -> dMag/dPert with shape (3, npert) in Cart. coords.
    dmag_dpert_npts: dict[int, np.array] = field(init=False)

    def __post_init__(self):
        """
        Compute derivatives wrt perturbations using finite differences.
        """
        npert = self.npert
        np_vals = max(len(pert.values) for pert in self.perts)

        # Use all stencils compatible with input num_points so that we can monitor the convergence afterwards.
        self.dforces_dpert_npts = {}
        self.dstress_dpert_npts = {}
        self.dpol_dpert_npts = {}
        self.dmag_dpert_npts = {}

        # FIXME Is this safe?
        # Do we want to allow for non-linear meshes.
        ipv0 = self.perts[0].ipv0

        for acc, weights in central_fdiff_weights[1].items():
            if np_vals < len(weights): continue
            nn = acc // 2
            npts = len(weights)
            # fd_slice is used to select the values for the Finite difference.
            fd_slice = slice(ipv0 - nn, ipv0 + nn + 1)

            # Finite differences for forces.
            dforce_dpert = np.empty((self.natom, 3, npert))
            for iat, iat_dir, ip in itertools.product(range(self.natom), range(3), range(npert)):
                pert = self.perts[ip]
                fvals_f = self.cart_forces_pv[ip, :, iat, iat_dir]
                dforce_dpert[iat, iat_dir, ip] = np.sum(fvals_f[fd_slice] * weights) / pert.step
            self.dforces_dpert_npts[npts] = dforce_dpert

            # Finite difference for stresses (Voigt form).
            dstress_dpert = np.empty((6, npert))
            for ivoigt, ip in itertools.product(range(6), range(npert)):
                pert = self.perts[ip]
                svals_f = self.carts_stresses_pv[ip, :, ivoigt]
                dstress_dpert[ivoigt, ip] = np.sum(svals_f[fd_slice] * weights) / pert.step
            self.dstress_dpert_npts[npts] = dstress_dpert

            # Finite differences for polarization (if available).
            if self.has_pol:
                dpol_dpert = np.empty((3, npert))
                for ii, ip in itertools.product(range(3), range(npert)):
                    pert = self.perts[ip]
                    # Shape is (npert, np_vals, 3)
                    dpol_dpert[ii, ip] = np.sum(self.cart_pol_pv[ip, fd_slice, ii] * weights) / pert.step
                self.dpol_dpert_npts[npts] = dpol_dpert
                #self.dpole_dpert_npts[npts] = dpole_dpert  TODO ?
                #self.dpoli_dpert_npts[npts] = dpoli_dpert  TODO ?

            # Finite differences for magnetization (if available).
            if self.has_mag:
                dmag_dpert = np.empty((3, npert))
                for ii, ip in itertools.product(range(3), range(npert)):
                    pert = self.perts[ip]
                    # Shape is (npert, np_vals, 3)
                    dmag_dpert[ii, ip] = np.sum(self.cart_mag_pv[ip, fd_slice, ii] * weights) / pert.step
                self.dmag_dpert_npts[npts] = dmag_dpert

    @lazy_property
    def natom(self) -> int:
        """Numbef of atoms in the unit cell."""
        return len(self.initial_structure)

    @lazy_property
    def npert(self) -> int:
        """Number of perturbations."""
        return len(self.perts)

    @lazy_property
    def pert_kind(self) -> str:
        """Kind of perturbation treated."""
        pert_kind = self.perts[0].kind
        all_kinds = [p.kind for p in self.perts]
        if any(k != pert_kind for k in all_kinds):
            raise ValueError(f"Expecting perturbations of the same kind but got {all_kinds}")
        return pert_kind

    @lazy_property
    def pert_dir_comps(self) -> list[str]:
        """
        List with the compononents associated to self.perts
        """
        if any(pert.cart_dir is None for pert in self.perts):
            raise TypeError("pert_dir_comps requires perturbation with directions!")

        return [vec2str(pert.cart_dir) for pert in self.perts]

    def __str__(self) -> str:
        return self.to_string()

    def get_elements_iatom_set(self, elements, iat_list) -> tuple:
        """
        Helper functions to convert elements and iat_list provided by users.
        """
        if elements is not None: elements = list_strings(elements)
        if iat_list is not None: iat_list = set(iat_list)
        if elements is not None and iat_list is not None:
            raise ValueError("elements and iat_list are mutually exclusive.")
        return elements, iat_list

    def add_geo_and_params(self, dct: dict, with_geo: bool, with_params: bool, with_spglib: bool = False, **kwargs):
        """
        Add info on structure and computational parameters to dct
        """
        # Add info on structure.
        if with_geo:
            geo_dict = self.structure.get_dict4pandas(with_geo=with_geo, with_spglib=with_spglib, **kwargs)
            dct.update(geo_dict)

        if with_params:
            dic.update(self.params)

    def print_relaxed_coords(self,
                             elements: list[str] | None = None,
                             iat_list: int | None = None) -> None:
        """
        Print relaxed atomic coordinates for each perturbation.

        Args:
            elements: String or list of strings with the chemical elements to select. None to select all.
            iat_list: List of atom indices to shown. None to select all.
        """
        elements, iat_set = self.get_elements_iatom_set(elements, iat_list)

        for iatom in range(self.natom):
            init_site = self.initial_structure[iatom]
            if elements is not None and init_site.species_string not in elements: continue
            if iat_set is not None and iatom not in iat_set: continue
            for ip, pert in enumerate(self.perts):
                print(init_site, "Initial site")
                for ipv, p_val in enumerate(pert.values):
                    site = self.structures_pv[ip, ipv][iatom]
                    print(site, f"{ip=}, {ipv=} value={pert.values[ipv]} {pert.name}")
                print("\n")

    def get_zeffs_list(self):
        """
        Dataframe with the effective charges for the given atom index and all the FD points.
        """
        # Zeff can be computed in different ways.
        # Here we find the kind of Zeff and the quantity to differentiate.
        field2zeff = {PertKind.E: "Ze", PertKind.H: "Zm"}

        if self.pert_kind in field2zeff:
            zeff_name, what_to_diff = field2zeff[self.pert_kind], "forces"

        elif self.pert_kind == PertKind.DISPL:
            if self.has_pol:
                zeff_name, what_to_diff = "Ze", "polarization"
            elif self.has_mag:
                zeff_name, what_to_diff = "Zm", "magnetization"
            else:
                raise ValueError(f"{self.pert_kind=} but neither polarization nor magnetization are available!")
        else:
            raise ValueError(f"Don't know how to compute eff_charges with {self.pert_kind=}")

        # ade stands from atom, direction and electric field.
        from abipy.dfpt.ddb import Zeffs, ZeffsList
        zeffs_list = ZeffsList()

        rows, xyz_comps = [], "x y z".split()
        natom = len(self.initial_structure)

        if what_to_diff == "forces":
            zeff_comps = list(itertools.product(xyz_comps, self.pert_dir_comps))
            for npts, dforces_dpert in self.dforces_dpert_npts.items():
                z_ade = dforces_dpert.copy()
                params = {"npts": npts}
                zeffs = Zeffs(zeff_name, z_ade, self.initial_structure, params=params)
                zeffs_list.append(zeffs)

        if what_to_diff in ("polarization", "magnetization"):
            if what_to_diff == "polarization":
                dvec_dpert_npts = self.dpol_dpert_npts
            if what_to_diff == "magnetization":
                dvec_dpert_npts = self.dmag_dpert_npts

            # dvec_dpert has shape (3, npert) where npert is 3*natom atomic displacements.
            for npts, dpol_dpert in dvec_dpert_npts.items():
                z_ade, atom_comps, cnt = np.empty((natom, 3, 3)), [], 0
                params = {"npts": npts}
                for ip, pert in enumerate(self.perts):
                    cnt += 1
                    iat_dir = ip % 3
                    z_ade[pert.iatom, iat_dir,:] = dpol_dpert[:, ip]
                    atom_comps.append(pert.dir_str)

                z_ade *= self.initial_structure.volume * abu.Ang_Bohr ** 3
                zeff_comps = list(itertools.product(atom_comps, xyz_comps))
                zeffs = Zeffs(zeff_name, z_ade, self.initial_structure, params=params)
                zeffs_list.append(zeffs)

        return zeffs_list

    def print_eff_charges(self,
                          elements: None | list[str] = None,
                          iat_list: None | list[int] = None,
                          file=sys.stdout,
                          verbose: int = 0) -> None:
        """
        Print effective charges to `file`.

        Args:
            elements: String or list of strings with the chemical elements to select. Default: All atoms are shown.
            iat_list: List of atom indices to shown. None to select all.
        """
        elements, iat_set = self.get_elements_iatom_set(elements, iat_list)

        def _p(*args, **kwargs):
            print(*args, file=file, **kwargs)

        if verbose:
            _p("Initial structure:")
            _p(self.initial_structure)
            _p("")

        zeffs_list = self.get_zeffs_list()
        zeff_name = zeffs_list[0].name
        _p(f"{zeff_name}[atom_dir, {self.pert_kind}_dir] in Cart. coords and a.u.")
        df = zeffs_list.concat(elements=elements)
        _p(df)

    @add_fig_kwargs
    def plot_etotal(self, mode="diff", fontsize=8, sharey=False, **kwargs) -> Figure:
        """
        Plot the total energy as a function of the amplitude of the perturbation.

        Args:
            mode: "diff" to plot the difference wrt to the unperturbed configuration.
        """
        nrows, ncols = self.npert, 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=True)

        for ip, pert in enumerate(self.perts):
            ax = ax_list[ip]
            ys = self.etotals_pv[ip] * abu.Ha_meV / self.natom
            if mode == "diff": ys -= ys[pert.ipv0]
            ax.plot(pert.values, ys, marker="o", label=pert.label)
            quadratic_fit_ax(ax, pert.values, ys, fontsize)

            ylabel = r"$\Delta$ Energy/atom (meV)" if mode == "diff" else "Energy/atom (meV)"
            set_grid_legend(ax, fontsize,
                            xlabel=f"${pert.tex}$ (a.u.)" if ip == len(self.perts) - 1 else None,
                            ylabel=ylabel if ip == 0 else None,
                            )
        return fig

    @add_fig_kwargs
    def plot_forces(self,
                    elements: None | list[str] = None,
                    iat_list: None | list[int] = None,
                    fontsize=8, sharey=False,
                    **kwargs) -> Figure:
        """
        Plot Cartesian forces as a function of the amplitude of the perturbation.

        Args:
            elements: String or list of strings with the chemical elements to select. None to select all.
            iat_list: List of atom indices to shown. None to select all.
        """
        elements, iat_set = self.get_elements_iat_set(elements, iat_list)

        nrows, ncols = 3, self.npert
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)

        for iat_dir in range(3):
            for ip, pert in enumerate(self.perts):
                ax = ax_mat[iat_dir, ip]
                ax.set_title(f"{pert.label}, Atom_dir: {pert.dir_str}", fontsize=fontsize)
                for iat, site in enumerate(self.initial_structure):
                    if elements is not None and site.species_string not in elements: continue
                    if iat_set is not None and iat not in iat_set: continue
                    ys = self.cart_forces_pv[ip, :, iat, iat_dir]
                    ax.plot(pert.values, ys, marker="o", label=site.species_string + r"$_{\text{%s}}$" % iat)

                ax.legend(loc="best", fontsize=fontsize, shadow=True)
                #set_grid_legend(ax, fontsize, xlabel=, ylabel=)

        return fig

    @add_fig_kwargs
    def plot_stresses(self, fontsize=8, sharey=False, **kwargs) -> Figure:
        """
        Plot Cartesian stresses as a function of the perturbation amplitude.
        """
        nrows, ncols = self.npert, 1
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)
        for ip, pert in enumerate(self.perts):
            ax = ax_mat[ip, 0]
            ax.set_title(pert.label, fontsize=fontsize)
            for ivoigt in range(6):
                ys = self.carts_stresses_pv[ip, :, ivoigt]
                ax.plot(pert.values, ys, marker="o", label=r"$\sigma_{%s}$" % (f"{ivoigt}"))

            ax.legend(loc="best", fontsize=fontsize, shadow=True)
            #set_grid_legend(ax, fontsize, xlabel=f"${pert.tex}$ (a.u.)", ylabel=)

        return fig

    @add_fig_kwargs
    def plot_polarization(self, what: str = "total", fontsize=8, sharey=False, **kwargs) -> Figure:
        """
        Plot the polarization as a function of the perturbation amplitude.
        """
        if self.cart_pol_pv is None:
            raise ValueError("The polarization has not been computed.")

        # Select the quantity to plot depending on `what`. Shape is (npert, np_vals, 3).
        vals_pv = {
            "total": self.cart_pol_pv,
            "electronic": self.cart_pole_pv,
            "ionic": self.cart_poli_pv,
        }[what]

        nrows, ncols = len(self.perts), 3
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)

        for ip, pert in enumerate(self.perts):
            for pol_dir in range(3):
                ax = ax_mat[ip, pol_dir]
                #ax.set_title(f"H_dir: {idir2s(pdir)}, Atom_dir: {idir2s(iat_dir)}", fontsize=fontsize)
                ys = vals_pv[ip, :, pol_dir]
                ax.plot(pert.values, ys, marker="o", label=pert.label)
                quadratic_fit_ax(ax, pert.values, ys, fontsize)

                set_grid_legend(ax, fontsize,
                                xlabel=f"${pert.tex}$ (a.u.)" if ip == len(self.perts) - 1 else None,
                                #ylabel=ylabel if ip == 0 else None,
                                )
        return fig

    @add_fig_kwargs
    def plot_magnetization(self, fontsize=8, sharey=False, **kwargs) -> Figure:
        """
        Plot the magnetization as a function of the perturbation amplitude.
        """
        if self.cart_mag_pv is None:
            raise ValueError("Polarization has not been computed.")

        nrows, ncols = len(self.perts), 3
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)

        for ip, pert in enumerate(self.perts):
            for mag_dir in range(3):
                # (npert, np_vals, 3)
                ys = self.cart_mag_pv[ip, :, mag_dir]
                ax = ax_mat[ip, mag_dir]
                #ax.set_title(f"H_dir: {idir2s(pdir)}, Atom_dir: {idir2s(iat_dir)}", fontsize=fontsize)
                ax.plot(pert.values, ys, marker="o", label=pert.label)
                quadratic_fit_ax(ax, pert.values, ys, fontsize)

                set_grid_legend(ax, fontsize,
                                xlabel=f"${pert.tex}$ (a.u.)" if ip == len(self.perts) - 1 else None,
                                #ylabel=ylabel if ip == 0 else None,
                                )
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """Generates figures common to the different subclasses."""
        yield self.plot_etotal(show=False)
        yield self.plot_forces(show=False)
        yield self.plot_stresses(show=False)
        if self.cart_pol_pv is not None:
            yield self.plot_polarization(show=False)
        if self.cart_mag_pv is not None:
            yield self.plot_magnetization(show=False)



@dataclasses.dataclass(kw_only=True)
class DisplData(_FdData):
    """
    Specialized class to handle finite differences wrt atomic displacement at q = 0.

    NB: Ze and Zm are implemented in get_zeffs_list of the super class.
    """

    def get_force_constants(self, npts: int) -> np.ndarray:
        # K_mn = d2E/{du_m du_n} = -dF_m/ du_n
        # dforces_dpert has shape (natom, 3, npert)
        # TODO: Singular value decomposition
        dforces_dpert = self.dforces_dpert_npts[npts]
        kmn = np.empty(self.natom, 3, self.natom, 3)
        for ip, pert in enumerate(self.perts):
            idir = ip % 3
            kmn[pert.iatom, idir] = - dforces_dpert[:,:,ip]
        return kmn

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose"""
        strio = StringIO()
        if self.has_pol or self.has_mag:
            self.print_eff_charges(file=strio)
        #print("piezoelectric tensor in Cartesian coords:\n", self.get_piezoel_df(), end=2*"\n", file=strio)

        strio.seek(0)
        return strio.read()


# WVH


@dataclasses.dataclass(kw_only=True)
class StrainData(_FdData):
    """
    Specialized class to handle finite diff. wrt strain.
    """

    def get_elastic(self, npts: int) -> np.ndarray:
        """
        Elastic tensor obtained with npts FD points. Eq (5) of WVH.
        """
        # dStress/dPert has shape (6, npert) in Cart. coords.
        dstress_dpert = self.dstress_dpert_npts[npts]
        cmat = np.empty((6, 6))
        for ip, pert in enumerate(self.perts):
            iv1 = pert.voigt_ind
            cmat[iv1] = dstress_dpert[:,ip]

        return cmat

    def get_elastic_df(self) -> pd.DataFrame:
        """
        Dataframe with the elastic tensor obtained with different FD points.
        """
        voigt_comps = [str(i) for i in range(1, 7)]
        cmat_comps = list(itertools.product(voigt_comps, voigt_comps))
        rows = []
        for npts in self.dstress_dpert_npts:
            cmat = self.get_elastic(npts)
            d = _dict_from_mat_npts(cmat, cmat_comps, npts)
            rows.append(d)

        return pd.DataFrame(rows)

    def get_internal_strain(self, npts: int) -> np.ndarray:
        """
        Internal-strain tensor obtained with npts FD points. Eq (7) of WVH.
        """
        # dForces/dPert has shape (natom, 3, npert) in Cart. coords.
        dforces_dpert = self.dforces_dpert_npts[npts]
        lmat = np.empty((self.natom, 3, 6))
        for ip, pert in enumerate(self.perts):
            iv1 = pert.voigt_ind
            print(iv1)
            lmat[:,:iv1] = dforces_dpert[:,ip]

        return np.reshape(lmat, (self.natom*3, 6))

    #def get_internal_strain_df(self) -> pd.DataFrame:

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level `verbose`"""
        strio = StringIO()
        #print("internal-strain tensor in Cartesian coords:\n", self.get_internal_strain_df(), end=2*"\n", file=strio)
        print(f"Elastic tensor in Cartesian coords and a.u. ({self.ions_mods}):\n",
               self.get_elastic_df(), end=2*"\n", file=strio)

        strio.seek(0)
        return strio.read()


class _HasExternalField:
    """
    Mixin class for calculations in which the perturbation is an external field.
    """

    def find_ip_pert_from_cart_dir(self, field_cart_dir) -> tuple[int, Perturbation]:
        """
        Find the perturbation from `field_cart_dir` that can be either a vector or an integer.
        Return perturbation index and perturbation.
        """
        if duck.is_intlike(field_cart_dir):
            ip = int(field_cart_dir)
            return ip, self.perts[ip]

        field_cart_dir = np.array(field_cart_dir)
        for ip, pert in enumerate(self.perts):
            if np.all(np.abs(pert.cart_dir - field_cart_dir) < 1e-6):
                return ip, pert

        raise ValueError(f"Cannot find perturbation with {field_cart_dir=}")

    @add_fig_kwargs
    def plot_forces_vs_field(self,
                             field_cart_dir,
                             elements: None | list[str] = None,
                             iat_list: None | list[int] = None,
                             fontsize=8, sharey=False,
                             **kwargs) -> Figure:
        """
        Plot Cartesian forces as a function of the amplitude of the perturbation.

        Args:
            elements: String or list of strings with the chemical elements to select. None to select all.
            iat_list: List of atom indices to shown. None to select all.
        """
        elements, iat_set = self.get_elements_iatom_set(elements, iat_list)
        nrows = self.natom
        if elements is not None:
            nrows = len(elements)
        if iat_set is not None:
            nrows = len(iat_set)

        ip, pert = self.find_ip_pert_from_cart_dir(field_cart_dir)

        nrows, ncols = nrows, 3
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)

        irow = -1
        for iat, site in enumerate(self.initial_structure):
            if elements is not None and site.species_string not in elements: continue
            if iat_set is not None and iat not in iat_set: continue
            irow += 1
            for iat_dir in range(3):
                ax = ax_mat[irow, iat_dir]
                #ax.set_title(f"{pert.label}, Atom_dir: {pert.dir_str}", fontsize=fontsize)
                ys = self.cart_forces_pv[ip, :, iat, iat_dir]
                ax.plot(pert.values, ys, marker="o", label=site.species_string + r"$_{\text{%s}}$" % iat)
                quadratic_fit_ax(ax, pert.values, ys, fontsize)

                ax.legend(loc="best", fontsize=fontsize, shadow=True)
                #set_grid_legend(ax, fontsize,
                #                xlabel=f"${pert.tex}$ (a.u.)" if ip == len(self.perts) - 1 else None,
                #                ylabel=ylabel if ip == 0 else None,
                #                )

        return fig


@dataclasses.dataclass(kw_only=True)
class ElectricFieldData(_FdData, _HasExternalField):
    """
    Specialized class to handle finite diff. wrt the electric field.
    """

    def get_epsinf(self, npts: int) -> np.ndarray:
        """eps_infinity obtained with npts FD points."""
        # FIXME: FD and DFPT do not agree here!
        eps = self.dpol_dpert_npts[npts] * 4.0 * np.pi
        eps[np.diag_indices_from(eps)] += 1.0
        return eps

    def get_epsinf_data(self, with_geo=False, with_params=False) -> DielectricDataList:
        """
        Dataframe with the components of eps_infinity obtained with different FD points.

        Args:
            with_geo: True to add info on structure.
            with_params: True to add calculations parameters.
        """
        #comps2inds = {"xx": (0,0), "yy": (1,1), "zz": (2,2),
        #              "xy": (0, 1), "xz": (0, 2), "yx": (1, 0), "yz": (1, 2), "zx": (2, 0), "zy": (2, 1)}

        #eps_inf_comps = list(itertools.product(self.pert_dir_comps, self.pert_dir_comps))

        diel_data = DielectricDataList()
        for npts in self.dpol_dpert_npts:
            eps_inf = self.get_epsinf(npts)
            meta = {"npts": npts}
            if with_params:
                meta.update(self.params)
            diel_data.append((eps_inf, self.initial_structure, meta))

        return diel_data

    def get_piezoel(self, npts: int, proper: bool) -> np.ndarray:
        """Improper piezo-electric tensor obtained with npts FD points. Eq (8) of WVH."""
        # dstress_dpert has shape (6, npert) in Cart. coords.
        dstress_dpert = self.dstress_dpert_npts[npts]
        piezoel = np.empty((3, 6))
        for ip, pert in enumerate(self.perts):
            piezoel[ip] = -dstress_dpert[:,ip]

        if proper:
            # Compute proper tensor. Eq (A9) of WVH.
            if not self.has_pol:
                raise RuntimeError("Polarization is needed to compute the proper piezoelectric tensor.")
            raise NotImplementedError()
            # Go from Voigt to (3,3)
            # Add polarization terms
            # Shape: (npert, np_vals, 3)
            cart_pol = self.cart_pol_pv[ip, ipv0]
            # From (3, 3) to Voigt.

        return piezoel

    def get_piezoel_df(self, proper, with_geo=False, with_params=False, **kwargs) -> pd.Dataframe:
        """
        Dataframe with the components of the piezo-electric tensor obtained with different FD points.

        Args:
            proper
            with_geo: True to add info on structure.
            with_params: True to add calculations parameters.
            kwargs: Optional kwargs passed to add_geo_and_params.
        """
        voigt_comps = [str(i) for i in range(1, 7)]
        piezoel_comps = list(itertools.product(self.pert_dir_comps, voigt_comps))
        rows = []
        for npts in self.dpol_dpert_npts:
            piezoel = self.get_piezoel(npts, proper)
            d = _dict_from_mat_npts(piezoel, piezoel_comps, npts)
            self.add_geo_and_params(d, with_geo, with_params, **kwargs)
            d["proper"] = proper
            rows.append(d)

        return pd.DataFrame(rows)

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose"""
        strio = StringIO()
        print(f"Epsilon_inf tensor in Cartesian coords and a.u. ({self.ions_mode}):\n",
              self.get_epsinf_data(), end=2*"\n", file=strio)
        proper = False
        print(f"Piezoelectric tensor in Cartesian coords and a.u. ({self.ions_mode}):\n",
              self.get_piezoel_df(proper), end=2*"\n", file=strio)
        self.print_eff_charges(file=strio)

        strio.seek(0)
        return strio.read()

    #def yield_figs(self, **kwargs):  # pragma: no cover
    #    """This function *generates* a predefined list of matplotlib figures with minimal input from the user."""
    #    # First, yield everything from the superclass
    #    yield from super().yield_figs()
    #    #self.plot_forces_vs_field([1, 0, 0], elements="Al")
    #    yield self.plot_polarization(show=False)



@dataclasses.dataclass(kw_only=True)
class ZeemanData(_FdData, _HasExternalField):
    """
    Specialized class to handle finite diff. wrt the Zeeman magnetic field.
    """

    def get_piezomag(self, npts: int) -> np.ndarray:
        """Piezo-magnetic tensor obtained with npts FD points."""
        # dstress_dpert has shape (6, npert) in Cart. coords.
        dstress_dpert = self.dstress_dpert_npts[npts]
        piezomag = np.empty((3, 6))
        for ip, pert in enumerate(self.perts):
            piezomag[ip] = -dstress_dpert[:,ip]
        return piezomag

    def get_piezomag_df(self, with_geo=False, with_params=False, **kwargs) -> pd.Dataframe:
        """
        Dataframe with the components of the piezo-magnetic tensor obtained with different FD points.

        Args:
            with_geo: True to add info on structure.
            with_params: True to add calculations parameters.
            kwargs: Optional kwargs passed to add_geo_and_params.
        """
        voigt_comps = [str(i) for i in range(1, 7)]
        piezomag_comps = list(itertools.product(self.pert_dir_comps, voigt_comps))
        rows = []
        for npts in self.dstress_dpert_npts:
            piezomag = self.get_piezomag(npts)
            d = _dict_from_mat_npts(piezomag, piezomag_comps, npts)
            self.add_geo_and_params(d, with_geo, with_params, **kwargs)
            rows.append(d)

        return pd.DataFrame(rows)

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose"""
        strio = StringIO()
        print(f"Piezomagnetic tensor in Cartesian coords and a.u. ({self.ions_mode}):\n",
              self.get_piezomag_df(), end=2*"\n", file=strio)
        self.print_eff_charges(file=strio)

        strio.seek(0)
        return strio.read()



def _dict_from_mat_npts(mat: np.ndarray, mat_comps: list[str], npts: int, with_info: bool = True,
                        comps2inds: dict | None = None) -> dict:
    """
    Convert a numpy array to a dict that can be used to construct a pandas DataFrame.
    """
    d = {"npts": npts}
    if comps2inds is None:
        d.update({c: v for c, v in zip(mat_comps, mat.flatten(), strict=True)})
    else:
        for k, ind in comps2inds.items():
            d[k] = mat[ind]

    if with_info and mat.shape[0] == mat.shape[1]:
        d["iso_avg"]= np.trace(mat) / mat.shape[0]
        #d["det"] = np.linalg.det(mat)
        tmp_mat = (mat + mat.T) / 2.0
        eigvals = np.linalg.eigvalsh(mat)  # Assuming Hermitian matrix
        d["det"] = np.prod(eigvals)
        d["posdef"] = np.all(eigvals > 0)
        #d["min_eig"] = np.min(eigvals)

    return d



class PertKind(StrEnum):
    """String enumerator for the different kind of perturbations."""
    DISPL = "displ"
    E = "E"
    H = "H"
    STRAIN = "strain"


class IonsMode(StrEnum):
    CLAMPED = "clamped_ions"
    RELAXED = "relaxed_ions"


@dataclasses.dataclass
class Perturbation:
    """
    This object stores info on the perturbation and its amplitude.
    """
    kind: str
    values: np.ndarray

    # Optional arguments.
    cart_dir: np.ndarray | None = None
    iatom: int | None = None
    voigt_ind: int | None = None
    strain: np.ndarray | None = None

    def __post_init__(self):
        """Implement validation logic."""
        if self.kind not in PertKind:
            raise ValueError(f"Invalid {self.kind=}")

        if self.kind == PertKind.DISPL:
            if self.iatom is None or self.cart_dir is None:
                raise ValueError("iatom and cart_dir must be specified for `displ` perturbations.")

        if self.kind in (PertKind.E, PertKind.H):
            if self.cart_dir is None:
                raise ValueError("cart_dir must be specified for `E` or `H` perturbations.")

        if self.kind == PertKind.STRAIN:
            if self.strain is None or self.voigt_ind is None:
                raise ValueError("strain matrix and voigt_ind must be specified for `strain` perturbations.")
            if self.voigt_ind not in ALL_STRAIN_INDS:
                raise ValueError(f"Invalid {self.voigt_ind=}")

    @lazy_property
    def step(self) -> float:
        """Step of the linear mesh. Raises ValueError if mesh is not linear."""
        dx = np.zeros(len(self.values) - 1)
        for i, x in enumerate(self.values[:-1]):
            dx[i] = self.values[i+1] - x

        if np.allclose(dx[0], dx):
            return float(dx[0])

        raise ValueError(f"Mesh is not homogenous: {dx=}")

    # TODO: Is this safe to use?
    @lazy_property
    def ipv0(self) -> int:
        """Index of the """
        return np.argmin(np.abs(self.values))

    @lazy_property
    def label(self) -> str:
        """Label string used in the plots."""
        if self.kind == PertKind.STRAIN:
            return "${%s}_{%s}$" % (self.tex, self.voigt_ind)

        return "${%s}_{%s}$" % (self.tex, vec2str(self.cart_dir))

    @lazy_property
    def dir_str(self) -> str:
        """String with the direction of the perturbation."""
        return "" if self.kind == PertKind.STRAIN else f"{vec2str(self.cart_dir)}"

    @lazy_property
    def tex(self) -> str:
        """Latex symbol"""
        return {
            PertKind.E: r"{\mathcal{E}}",
            PertKind.H: r"{\mathcal{H}}",
            PertKind.DISPL: r"{u}",
            PertKind.STRAIN: r"{\varepsilon}",
        }[self.kind]

    @lazy_property
    def name(self) -> str:
        """Name of the perturbation."""
        return {
            PertKind.E: "Electric field",
            PertKind.H: "Magnetic field",
            PertKind.DISPL: "Atomic displacement",
            PertKind.STRAIN: "Strain",
        }[self.kind]


class _BaseFdWork(Work):
    """
    Base class for finite difference Works.
    """

    @lazy_property
    def natom(self) -> int:
        """Number of atoms in the unit cell."""
        return len(self[0].input.structure)

    @lazy_property
    def npert(self) -> int:
        """Number of perturbations."""
        return len(self.perts)

    @lazy_property
    def all_ions_modes(self) -> list[str]:
        return [IonsMode.CLAMPED, IonsMode.RELAXED] if self.relax_ions else [IonsMode.CLAMPED]

    @lazy_property
    def gs_tasks_ids(self) -> set:
        """Set with the ids of the gs tasks."""
        return {task.node_id for task in self.gs_tasks_pv.flat}

    def allocate_tasks_pv(self, relax_ions: bool, relax_ions_opts) -> None:
        """
        Allocate arrays with GS tasks and relax tasks. Also, set input variables for ionic relaxations.
        """
        # Provide default values for relaxation algorithm, and allow user to override settings.
        self.relax_ions = relax_ions
        self.relax_ions_opts = {"ionmov": 2, "tolmxf": 1e-6}
        if relax_ions_opts:
            self.relax_ions_opts.update(**relax_ion_opts)

        np_vals = max(len(pert.values) for pert in self.perts)
        self.gs_tasks_pv = np.empty((self.npert, np_vals), dtype=object)
        self.relax_tasks_pv = np.empty((self.npert, np_vals), dtype=object)

    def get_data_dict(self, ions_mode: str) -> dict:
        """
        This method is shared by all the Finite difference works.
        It reads energies, forces, stresses, polarizations and magnetizations from the GSR.nc files
        and the main output files of the tasks and builds a dictionary that can be used
        to instanciate the appropriate subclass of _FdData.
        """
        npert, np_vals = len(self.perts), max(len(pert.values) for pert in self.perts)

        has_pol = all(task.input.get("berryopt", 0) != 0 for task in self)

        # Detect if a magnetic calculation is being performed by looking at nsppol and then nspinor.
        has_mag = False
        if all(task.input.get("nsppol", 1) == 2 for task in self):
            has_mag = True

        if all(task.input.get("nspinor", 1) == 2 for task in self):
            if all(task.input.get("nspden", 4) == 4 for task in self):
                has_mag = True

        data = {
            "initial_structure": self[0].input.structure,
            "ions_mode": ions_mode,
            "perts": self.perts,
            "has_pol": has_pol,
            "has_mag": has_mag,
        }

        # The suffix `_pv` stands for perturbation index and perturbation value.
        data["etotals_pv"] = etotals_pv = np.empty((npert, np_vals))
        data["eterms_pv"] = eterms_pv = np.empty((npert, np_vals), dtype=object)
        data["cart_forces_pv"] = cart_forces_pv = np.empty((npert, np_vals, self.natom, 3))
        data["carts_stresses_pv"] = carts_stresses_pv = np.empty((npert, np_vals, 6))
        data["params_pv"] = params_pv = np.empty((npert, np_vals), dtype=object)

        if has_pol:
            data["cart_pol_pv"] = cart_pol_pv = np.empty((npert, np_vals, 3))
            data["cart_pole_pv"] = cart_pole_pv = np.empty((npert, np_vals, 3))
            data["cart_poli_pv"] = cart_poli_pv = np.empty((npert, np_vals, 3))

        if has_mag:
            data["cart_mag_pv"] = cart_mag_pv = np.empty((npert, np_vals, 3))

        data["structures_pv"] = structures_pv = np.empty((npert, np_vals), dtype=object)

        # Read energy, forces and stress from the GSR files.
        for ip, pert in enumerate(self.perts):
            for ipv, p_val in enumerate(pert.values):

                # Select the appropriate task.
                if ions_mode == IonsMode.CLAMPED:
                    task = self.gs_tasks_pv[ip, ipv]
                elif ions_mode == IonsMode.RELAXED:
                    task = self.relax_tasks_pv[ip, ipv]
                else:
                    raise ValueError(f"Invalid {ions_mode=}")

                if task is None:
                    raise RuntimeError(f"Got None instead of task for {ip=}, {ipv=}, {ions_mode=} and work type: {type(self)}")
                #print(f"{task=}")

                with task.open_gsr() as gsr:
                    structures_pv[ip, ipv] = gsr.structure.copy()
                    etotals_pv[ip, ipv] = gsr.r.read_value("etotal")
                    eterms_pv[ip, ipv] = gsr.r.read_energy_terms(unit="Ha")
                    cart_forces_pv[ip, ipv] = gsr.r.read_value("cartesian_forces") # Ha/Bohr units.
                    # Abinit stores 6 unique components of this symmetric 3x3 tensor:
                    # Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
                    carts_stresses_pv[ip, ipv] = gsr.r.read_value("cartesian_stress_tensor")
                    # Add parameters that might be used for convergence studies afterwards.
                    params_pv[ip, ipv] = gsr.params.copy()
                    if has_mag:
                        # Get magnetization from the GSR file.
                        cart_mag_pv[ip, ipv] = gsr.get_magnetization()

                if has_pol:
                    # Read polarization from the abo file.
                    with task.open_abo() as abo:
                        pol = abo.get_berry_phase_polarization()
                        cart_pol_pv[ip, ipv] = pol.total
                        cart_pole_pv[ip, ipv] = pol.electronic
                        cart_poli_pv[ip, ipv] = pol.ionic

        return data

    def get_data(self) -> dict:
        return {ions_mode: self.DataCls(**self.get_data_dict(ions_mode)) for ions_mode in self.all_ions_modes}

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        data = self.get_data()
        for ions_mode, obj in data.items():
            with open(Path(self.outdir.path) / f"{ions_mode}_{obj.__class__.__name__}.pickle", "wb") as fh:
                pickle.dump(obj, fh)
            #mjson_write(obj, Path(self.outdir.path ) / f"{ions_mode}_{obj.__class__.__name__}.json")

        return super().on_all_ok()


class FiniteDisplWork(_BaseFdWork):
    """
    This work displaces the atoms in unit cell by a finite amount and performs
    GS calculations to get forces and stresses.
    """
    DataCls = DisplData

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       fd_accuracy: int,
                       step_au: float = 0.01,
                       pert_cart_dirs=None,
                       mask_iatom=None,
                       extra_abivars: dict | None = None,
                       manager=None):
        """
        Build a work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            fd_accuracy:
            step_au: Finite difference step for the displacement in Bohr (a.u.)
            pert_cart_dirs:
            mask_iatom:
            extra_abivars: dictionary with extra Abinit variables to be added to scf_input.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        work.scf_input = scf_input.deepcopy()
        if extra_abivars is not None:
            work.scf_input.set_vars(**extra_abivars)

        structure = scf_input.structure
        natom = len(structure)

        num_points, work.pert_values, _ipv0 = _mesh_for_fd_accuracy(fd_accuracy, 1, step_au)

        # Here we normalize the directions to 1. NB: pymatgen structures uses Ang and not Bohr.
        if pert_cart_dirs is None:
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

        # Build list of perturbations.
        work.perts = []
        for iatom, mask in zip(range(natom), work.mask_iatom, strict=True):
            if not mask: continue
            for cart_dir in work.pert_cart_dirs:
                work.perts.append(Perturbation(kind=PertKind.DISPL, values=work.pert_values, cart_dir=cart_dir, iatom=iatom))

        work.relax_ions = False
        work.allocate_tasks_pv(work.relax_ions, {})
        work._add_tasks_with_displacements(work.scf_input.structure)

        return work

    def _add_tasks_with_displacements(self, structure: Structure):
        """
        """
        for ip, pert in enumerate(self.perts):
            iatom, cart_dir = pert.iatom, pert.cart_dir
            for iv, delta_au in enumerate(pert.values):
                # Note Bohr --> Ang conversion.
                new_structure = structure.copy()
                new_structure.translate_sites([iatom], delta_au * abu.Bohr_Ang * cart_dir,
                                               frac_coords=False, to_unit_cell=False)
                new_input = self.scf_input.new_with_structure(new_structure)
                task = self.register_scf_task(new_input)
                self.gs_tasks_pv[ip, iv] = task


class FiniteStrainWork(_BaseFdWork):
    """
    This work deforms the initial unit cell by a finite amout and performs
    GS calculations to get forces and stresses.
    """
    DataCls = StrainData

    @classmethod
    def from_scf_input(cls,
                       scf_input,
                       fd_accuracy: int,
                       norm_step: float,
                       shear_step: float,
                       voigt_inds=None,
                       extra_abivars: dict | None = None,
                       relax_ions: bool = False,
                       relax_ions_opts: dict | None = None,
                       manager=None):
        """
        Build a Work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            fd_accuracy:
            norm_step: Finite difference step for normal strain.
            shear_step: Finite difference step for shear strain.
            voigt_inds:
            extra_abivars: dictionary with extra Abinit variables to be added to scf_input.
            relax_ions: False if the initial structural relaxation should not be performed.
            relax_ions_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default if None.
        """
        work = cls(manager=manager)
        work.scf_input = scf_input.deepcopy()
        if extra_abivars is not None:
            work.scf_input.set_vars(**extra_abivars)

        if "ecutsm" not in work.scf_input:
            ecutsm = 0.5
            work.scf_input.set_vars(ecutsm=ecutsm)
            cprint("AbinitInput does not define ecutsm.\nA default value of %s will be added" % ecutsm, color="yellow")

        # Build list of strain matrices.
        if voigt_inds is None:
            voigt_inds = list(range(6))

        strains = []
        for vind in voigt_inds:
            strains.append(Strain.from_index_amount(vind, amount=1.0))
        strains = np.reshape(strains, (-1, 3, 3))

        for strain in strains:
            if not np.array_equal(strain, strain.T):
                raise ValueError(f"The strain matrix should be symmetric but got: {strain}")

        # Use different list of pert_values for normal and shear strain.
        num_points, norm_values, _ipv0 = _mesh_for_fd_accuracy(fd_accuracy, 1, norm_step)
        num_points, shear_values, _ipv0 = _mesh_for_fd_accuracy(fd_accuracy, 1, shear_step)

        # Build list of perturbations.
        work.perts = []
        for vind, strain in zip(voigt_inds, strains, strict=True):
            work.perts.append(Perturbation(kind=PertKind.STRAIN, voigt_ind=vind, strain=strain,
                              values=norm_values if vind in NORMAL_STRAIN_INDS else shear_values))

        work.allocate_tasks_pv(relax_ions, relax_ions_opts)
        for ip, pert in enumerate(work.perts):
            work._add_tasks_with_strains_ipv(ip, None, work.scf_input.structure, IonsMode.CLAMPED)

        return work

        #from pymatgen.analysis.elasticity.strain import DeformedStructureSet
        #DeformedStructureSet(structure: Structure,
        #                     norm_strains: Sequence[float] = (-0.01, -0.005, 0.005, 0.01),
        #                     shear_strains: Sequence[float] = (-0.06, -0.03, 0.03, 0.06),
        #                     symmetry=False,

    def _add_tasks_with_strains_ipv(self,
                                    ip: int, ipv_select: int | None,
                                    structure: Structure,
                                    ions_mode: str) -> None:
        """Build new GS tasks with strained cells."""
        scf_input = self.scf_input
        pert = self.perts[ip]
        task_pv0 = None # TODO

        if ions_mode == IonsMode.CLAMPED:
            if ipv_select is not None:
                raise ValueError(f"ipv_select should be None if {ions_mode=} but got {ipv_select=}")

        for ipv, p_val in enumerate(pert.values):
            if ipv_select is not None and ipv != ipv_select: continue
            is_pv0 = abs(p_val) < 1e-16

            if ions_mode == IonsMode.CLAMPED:
                # Apply strain to the lattice and build new SCF tasks.
                strained_structure = scf_input.structure.apply_strain(p_val * pert.strain, inplace=False)
                task = self.register_scf_task(scf_input.new_with_structure(strained_structure))
                # Add the (ip, iv) indices as attribute of the task.
                task.attrs["ip_ipv"] = (ip, ipv)
                self.gs_tasks_pv[ip, ipv] = task

            elif ions_mode == IonsMode.RELAXED:
                # Get strained structure from gs_tasks_pv and relax atoms.
                relax_inp = self.gs_tasks_pv[ip,ipv].input.new_with_vars(**self.relax_ions_opts)
                task = self.register_relax_task(relax_inp)
                task.add_deps({self.gs_tasks_pv[ip,ipv]: "WFK"})
                # Add the (ip, iv) indices as attribute of the task.
                task.attrs["ip_ipv"] = (ip, ipv)
                self.relax_tasks_pv[ip, ipv] = task
            else:
                raise ValueError(f"Invalid {ions_mode=}")

    def on_ok(self, sender):
        """This method is called when one task reaches status `S_OK`."""
        # NB: Only gs tasks should trigger relaxed ions calculations.

        if self.relax_ions and sender.node_id in self.gs_tasks_ids:
            # Get structure from the GSR file.
            relaxed_structure = sender.get_final_structure()
            ip, iv = sender.attrs["ip_ipv"]
            self._add_tasks_with_strains_ipv(ip, iv, relaxed_structure, IonsMode.RELAXED)
            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        """This method is called when all tasks have reached S_OK."""
        data = self.get_data()
        for ions_mode, obj in data.items():
            with open(Path(self.outdir.path) / f"{ions_mode}_{obj.__class__.__name__}.pickle", "wb") as fh:
                pickle.dump(obj, fh)
            #mjson_write(obj, Path(self.outdir.path ) / f"{ions_mode}_{obj.__class__.__name__}.json")

        # TODO: Use ElasticData from dfpt.elastic module
        #from abipy.dfpt.elastic import ElasticData
        #el_data = ElasticData(
        #             self.initial_structure,
        #             params=self.params_pv[0,0]
        #             elastic_clamped=None,
        #             elastic_relaxed=None,
        #             elastic_stress_corr=None,
        #             elastic_relaxed_fixed_D=None,
        #             piezo_clamped=None,
        #             piezo_relaxed=None,
        #             d_piezo_relaxed=None,
        #             g_piezo_relaxed=None,
        #             h_piezo_relaxed=None)

        return super().on_all_ok()



class _FieldWork(_BaseFdWork):
    """Base class for finite field + finite difference Work."""

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       fd_accuracy: int,
                       step_au: float,
                       pert_cart_dirs: np.ndarray | None = None,
                       extra_abivars: dict | None = None,
                       relax_ions: bool = False,
                       relax_ions_opts: dict | None = None,
                       manager=None):
        """
        Build the work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF calculation used as template to generate the other inputs.
            fd_accuracy:
            step_au: Finite difference step for the magnetic field in a.u.
            pert_cart_dirs:
            extra_abivars: dictionary with extra Abinit variables to be added to scf_input.
            relax_ions: False if the initial structural relaxation should not be performed.
            relax_ions_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        num_points, work.pert_values, _ipv0 = _mesh_for_fd_accuracy(fd_accuracy, 1, step_au)

        if pert_cart_dirs is None:
            work.pert_cart_dirs = np.eye(3)

        work.pert_cart_dirs = np.reshape(work.pert_cart_dirs, (-1, 3))
        work.scf_input = scf_input.deepcopy()
        if extra_abivars is not None:
            work.scf_input.set_vars(**extra_abivars)

        if isinstance(work, FiniteEfieldWork):
            work.pert_kind = PertKind.E
        elif isinstance(work, FiniteHfieldWork):
            work.pert_kind = PertKind.H
        else:
            raise TypeError(f"Don't know how to handle {type(work)=}")

        # Build list of perturbations.
        work.perts = [Perturbation(kind=work.pert_kind, values=work.pert_values, cart_dir=cart_dir)
                      for cart_dir in work.pert_cart_dirs]

        work.allocate_tasks_pv(relax_ions, relax_ions_opts)

        for ip, pert in enumerate(work.perts):
            if work.pert_kind == PertKind.E:
                work._add_tasks_with_efield_ipv(ip, None, work.scf_input.structure, IonsMode.CLAMPED)
            elif work.pert_kind == PertKind.H:
                work._add_tasks_with_zeeman_field_ipv(ip, None, work.scf_input.structure, IonsMode.CLAMPED)
            else:
                raise TypeError(f"Don't know how to handle {work.pert_kind=}")

        return work

    def on_ok(self, sender):
        """This method is called when one task reaches status `S_OK`."""
        # NB: Only gs tasks should trigger relaxed ions calculations.

        if self.relax_ions and sender.node_id in self.gs_tasks_ids:
            # Get structure from the GSR file.
            # FIXME: Do we really need get_final_structure eventhough this is GS?
            relaxed_structure = sender.get_final_structure()
            ip, iv = sender.attrs["ip_ipv"]

            if self.pert_kind == PertKind.E:
                self._add_tasks_with_efield_ipv(ip, iv, relaxed_structure, IonsMode.RELAXED)
            elif self.pert_kind == PertKind.H:
                self._add_tasks_with_zeeman_field_ipv(ip, iv, relaxed_structure, IonsMode.RELAXED)
            else:
                raise TypeError(f"Don't know how to handle {self.pert_kind}")

            self.flow.allocate(build=True)

        return super().on_ok(sender)


class FiniteHfieldWork(_FieldWork):
    r"""
    This work performs GS calculations with a finite Zeeman field.
    """
    DataCls = ZeemanData

    def _add_tasks_with_zeeman_field_ipv(self, ip: int, ipv_select: int | None,
                                         structure: Structure, ions_mode: str) -> None:
        """Build new GS tasks with zeemanfield."""
        scf_input = self.scf_input.new_with_structure(structure)
        pert = self.perts[ip]
        task_pv0 = None

        if ions_mode == IonsMode.CLAMPED:
            register_func = self.register_scf_task
            tasks_pv = self.gs_tasks_pv
            relax_ions_opts = {}
            if ipv_select is not None:
                raise ValueError(f"ipv_select should be None if {ions_mode=} but got {ipv_select=}")

        elif ions_mode == IonsMode.RELAXED:
            register_func = self.register_relax_task
            tasks_pv = self.relax_tasks_pv
            relax_ions_opts = self.relax_ions_opts
        else:
            raise ValueError(f"Invalid {ions_mode=}")

        for ipv, p_val in enumerate(pert.values):
            if ipv_select is not None and ipv != ipv_select: continue
            is_pv0 = abs(p_val) < 1e-16

            new_inp = scf_input.new_with_vars(zeemanfield=p_val * pert.cart_dir, **relax_ions_opts)
            if is_pv0:
                # Avoid computing the zero-field case multiple times.
                # FIXME Not sure this works as expected.
                if task_pv0 is None:
                    task_pv0 = register_func(new_inp)
                tasks_pv[ip, ipv] = task_pv0
            else:
                tasks_pv[ip, ipv] = register_func(new_inp)


class FiniteEfieldWork(_FieldWork):
    r"""
    This work performs GS calculations with a finite Electric field.
    """
    DataCls = ElectricFieldData

    def _add_tasks_with_efield_ipv(self,
                                   ip: int,
                                   ipv_select: int | None,
                                   structure: Structure, ions_mode: str) -> None:
        """Build new GS tasks with finite electric field."""
        scf_input = self.scf_input.new_with_structure(structure)
        pert = self.perts[ip]
        task_pv0 = None

        if ions_mode == IonsMode.CLAMPED:
            register_func = self.register_berry_task
            tasks_pv = self.gs_tasks_pv
            relax_ions_opts = {}
            if ipv_select is not None:
                raise ValueError(f"ipv_select should be None if {ions_mode=} but got {ipv_select=}")

        elif ions_mode == IonsMode.RELAXED:
            register_func = self.register_relax_task
            tasks_pv = self.relax_tasks_pv
            relax_ions_opts = self.relax_ions_opts
        else:
            raise ValueError(f"Invalid {ions_mode=}")

        for ipv, p_val in enumerate(pert.values):
            if ipv_select is not None and ipv != ipv_select: continue
            is_pv0 = abs(p_val) < 1e-16
            new_inp = scf_input.new_with_vars(efield=p_val * pert.cart_dir, **relax_ions_opts)

            if tasks_pv[ip, ipv] is not None:
                raise RuntimeError(f"Expecting None for {ip=}, {ipv=} but got {str(tasks_pv[ip, ipv])}")

            if is_pv0:
                # Avoid computing the zero-field case multiple times.
                # Also, the task at zero field uses berryopt -1 to get the polarization.
                if task_pv0 is None:
                    new_inp.set_vars(berryopt=-1)
                    task_pv0 = register_func(new_inp)
                tasks_pv[ip, ipv] = task_pv0
            else:
                # Finite electric field E computation.
                new_inp.set_vars(berryopt=4)
                tasks_pv[ip, ipv] = register_func(new_inp)

            if ions_mode == IonsMode.RELAXED:
                # Will relaxed ions at finite E using the WFK(E) of the gs_task
                tasks_pv[ip, ipv].add_deps({self.gs_tasks_pv[ip, ipv]: "WFK"})

            # Add the (ip, iv) indices as attributes of the task
            tasks_pv[ip, ipv].attrs["ip_ipv"] = (ip, ipv)

        # Now add dependencies for GS tasks: connect tasks with +E and -E starting from E = 0.
        if ions_mode == IonsMode.CLAMPED:
            for ipv in range(0, pert.ipv0):
                tasks_pv[ip, ipv].add_deps({tasks_pv[ip, ipv+1]: "WFK"})

            for ipv in range(pert.ipv0+1, len(tasks_pv[ip])):
                tasks_pv[ip, ipv].add_deps({tasks_pv[ip, ipv-1]: "WFK"})

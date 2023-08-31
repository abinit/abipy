"""
Objects to perform ASE calculations with machine-learned potentials.
"""
from __future__ import annotations

import sys
import os
import io
import time
import contextlib
import json
import pickle
import warnings
import dataclasses
import stat
import numpy as np
import pandas as pd
try:
    import ase
except ImportError as exc:
    raise ImportError("ase not installed. Try `pip install ase`.") from exc

from pathlib import Path
from inspect import isclass
from multiprocessing import Pool
from typing import Type, Any, Optional, Union
from monty.string import marquee, list_strings # is_string,
from monty.functools import lazy_property
from monty.collections import AttrDict #, dict2namedtuple
from pymatgen.core import Structure as PmgStructure
from pymatgen.io.ase import AseAtomsAdaptor
from ase import units
from ase.atoms import Atoms
from ase.io.trajectory import write_traj, Trajectory
from ase.optimize.optimize import Optimizer
from ase.calculators.calculator import Calculator
from ase.io.vasp import write_vasp_xdatcar, write_vasp
from ase.neb import NEB
from ase.md.nptberendsen import NPTBerendsen, Inhomogeneous_NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from abipy.core import Structure
import abipy.core.abinit_units as abu
from abipy.tools.iotools import workdir_with_prefix, PythonScript
from abipy.tools.printing import print_dataframe
from abipy.abio.enums import StrEnum, EnumMixin
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_grid_legend,
    set_visible, set_ax_xylabels)


###################
# Helper functions
###################

def nprocs_for_ntasks(nprocs, ntasks, title=None) -> int:
    """
    Return the number of procs to be used in a multiprocessing Pool.
    If negative or None, use hlaf the procs in the system.
    """
    if nprocs is None or nprocs <= 0:
        nprocs = max(1, os.cpu_count() // 2)
    else:
        nprocs = int(nprocs)

    nprocs = min(nprocs, ntasks)
    if title is not None:
        print(title)
        print(f"Using multiprocessing pool with {nprocs=} for {ntasks=} ...")
    return nprocs


_CELLPAR_KEYS = ["a", "b", "c", "angle(b,c)", "angle(a,c)", "angle(a,b)"]


ASENEB_METHODS = ['aseneb', 'eb', 'improvedtangent', 'spline', 'string']


class RX_MODE(EnumMixin, StrEnum):  # StrEnum added in 3.11
    """
    Relaxation mode string flags.
    """
    no   = "no"
    ions = "ions"
    cell = "cell"

    #@classmethod
    #def from_ionmov_optcell(cls, ionmov: int, optcell: int) -> RX_MODE:
    #    cls.no
    #    cls.ions
    #    cls.cell


def to_ase_atoms(structure: PmgStructure, calc=None) -> Atoms:
    """Convert pymatgen structure to ASE atoms. Optionally, attach a calculator."""
    structure = Structure.as_structure(structure)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    if calc:
        atoms.calc = calc
    return atoms


def get_atoms(obj: Any) -> Atoms:
    """Return ASE Atoms from object."""
    if isinstance(obj, str):
        return to_ase_atoms(Structure.from_file(obj))
    if isinstance(obj, PmgStructure):
        return to_ase_atoms(obj)
    if isinstance(obj, Atoms):
        return obj
    raise TypeError(f"Don't know how to construct Atoms object from {type(obj)}")


def abisanitize_atoms(atoms: Atoms, **kwargs) -> Atoms:
    """
    Call abisanitize, return new Atoms instance.
    """
    structure = Structure.as_structure(atoms)
    new_structure = structure.abi_sanitize(**kwargs)
    return to_ase_atoms(get_atoms(new_structure), calc=atoms.calc)


def fix_atoms(atoms: Atoms,
              fix_inds: list[int] | None = None,
              fix_symbols: list[str] | None = None) -> None:
    """
    Fix atoms by indices and by symbols.

    Args:
        atoms: ASE atoms
        fix_inds: List of site indices to fix. None to ignore constraint.
        fix_symbols: List of chemical elements to fix. None to ignore the constraint.
    """
    from ase.constraints import FixAtoms
    cs = []; app = cs.append
    if fix_inds is not None:
        app(FixAtoms(indices=fix_inds))
    if fix_symbols is not None:
        fix_symbols = set(fix_symbols)
        app(FixAtoms(mask=[atom.symbol in fix_symbols for atom in atoms]))
    if cs:
        atoms.set_constraint(constraint=cs)


_FMT2FNAME = {
    "poscar": "POSCAR",
    "abinit": "run.abi",
    #"qe": "qe.in",
}

def write_atoms(atoms: Atoms, workdir, verbose: int,
                formats=None, prefix=None, postfix=None) -> list[tuple[Path, str]]:
    """
    Write atoms to file(s), return list with (Path, fmt) tuples.

    Args:
        atoms: ASE atoms
        workdir: Working directory.
        verbose: Verbosity level.
        formats: List of strings with file formats. If None all known formats are used.
        prefix: String to be prepended to filenames.
        prefix: String to be appended to filenames.
    """
    workdir = Path(workdir)
    structure = Structure.as_structure(atoms)
    fmt2fname = _FMT2FNAME
    if formats is not None:
        fmt2fname = {k: _FMT2FNAME[k] for k in list_strings(formats)}

    outpath_fmt = []
    for fmt, fname in fmt2fname.items():
        if prefix: fname = prefix + fname
        if postfix: fname = fname + postfix
        outpath = workdir / fname
        if verbose > 1: print(f"Writing atoms to: {outpath:} with {fmt=}")
        with open(outpath, "wt") as fh:
            fh.write(structure.convert(fmt=fmt))
        outpath_fmt.append((outpath, fmt))
    return outpath_fmt


def print_atoms(atoms: Atoms, title=None, cart_forces=None, stream=sys.stdout) -> None:
    """
    Print atoms object to stream.

    Args:
        atoms: ASE atoms.
        title: Optional string with the title.
        cart_forces: np.array with cart_forces to print.
        stream: Output stream
    """
    def pf(*args):
        print(*args, file=stream)

    scaled_positions = atoms.get_scaled_positions()
    if title is not None:
        pf(title)
    if cart_forces is None:
        pf("Frac coords:")
    else:
        pf("Frac coords and cart forces:")

    for ia, (atom, frac_coords) in enumerate(zip(atoms, scaled_positions)):
        if cart_forces is None:
            pf("\t", frac_coords)
        else:
            pf("\t", frac_coords, cart_forces[ia])


def diff_two_structures(label1, structure1, label2, structure2, fmt, file=sys.stdout):
    """
    Diff two structures using format `fmt`and print results to `file`.
    """
    lines1 = Structure.as_structure(structure1).convert(fmt=fmt).splitlines()
    lines2 = Structure.as_structure(structure2).convert(fmt=fmt).splitlines()
    pad = max(max(len(l) for l in lines1), len(label1), len(label2))
    print(label1.ljust(pad), " | ", label2, file=file)
    for l1, l2 in zip(lines1, lines2):
        print(l1.ljust(pad), " | ", l2, file=file)


@dataclasses.dataclass
class AseResults:
    """
    Container with the results produced by the ASE calculator.
    """
    atoms: Atoms
    ene: float
    stress: np.ndarray  # 3x3 matrix with stress
    forces: np.ndarray

    @classmethod
    def from_traj_inds(cls, trajectory, *inds) -> AseResults:
        """Build list of AseResults from a trajectory and list of indices."""
        return [cls.from_atoms(trajectory[i]) for i in inds]

    @classmethod
    def from_atoms(cls, atoms: Atoms, calc=None) -> AseResults:
        """Build the object from an atoms instance with a calculator."""
        if calc is not None:
            atoms.calc = calc
        results = cls(atoms=atoms,
                      ene=float(atoms.get_potential_energy()),
                      stress=atoms.get_stress(voigt=False),
                      forces=atoms.get_forces())
        if calc is not None:
            atoms.calc = None
        return results

    @property
    def pressure(self) -> float:
        return -self.stress.trace() / 3

    def get_voigt_stress(self):
        """xx, yy, zz, yz, xz, xy"""
        from ase.stress import full_3x3_to_voigt_6_stress
        return full_3x3_to_voigt_6_stress(self.stress)

    @property
    def volume(self) -> float:
        """Volume of the unit cell in Ang^3."""
        return self.atoms.get_volume()

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level `verbose`."""
        lines = []; app = lines.append

        app(f"Energy: {self.ene} (eV)")
        app(f"Pressure: {self.pressure} ")
        fstats = self.get_fstats()
        for k, v in fstats.items():
            app(f"{k} = {v}")
        #app('Stress tensor:', r.stress)
        if verbose:
            app('Forces (eV/Ang):')
            positions = self.atoms.get_positions()
            df = pd.DataFrame(dict(
                x=positions[:,0],
                y=positions[:,1],
                z=positions[:,2],
                fx=self.forces[:,0],
                fy=self.forces[:,1],
                fz=self.forces[:,2],
            ))
            app(str(df))

        return "\n".join(lines)

    def get_fstats(self) -> dict:
        """
        Return dictionary with statistics on forces.
        """
        fmods = np.array([np.linalg.norm(force) for force in self.forces])
        #fmods = np.sqrt(np.einsum('ij, ij->i', forces, forces))
        #return AttrDict(
        return dict(
            fmin=fmods.min(),
            fmax=fmods.max(),
            fmean=fmods.mean(),
            #fstd=fmods.std(),
            drift=np.linalg.norm(self.forces.sum(axis=0)),
        )

    def get_dict4pandas(self, with_geo=True, with_fstats=True) -> dict:
        """
        Dictionary with results used to build pandas dataframe.
        """
        d = {k: getattr(self, k) for k in ["ene", "volume", "pressure"]}
        if with_geo:
            d.update(dict(zip(_CELLPAR_KEYS, self.atoms.cell.cellpar())))
        if with_fstats:
            d.update(self.get_fstats())

        return d


def zip_sort(xs, ys):
    isort = xs.argsort()
    return xs[isort].copy(), ys[isort].copy()


def diff_stats(xs, ys):
    abs_diff = np.abs(ys - xs)
    return AttrDict(
       MAE=abs_diff.mean(),
       ADIFF_MIN=abs_diff.min(),
       ADIFF_MAX=abs_diff.max(),
       ADIFF_STD=abs_diff.std(),
    )


def linear_fit_ax(ax, xs, ys, fontsize, with_label=True, with_ideal_line=False) -> tuple[float]:
    """
    """
    from scipy.stats import linregress
    result = linregress(xs, ys)
    label = r"Linear fit $\alpha={:.2f}$, $r^2$={:.2f}".format(result.slope, result.rvalue**2)
    ax.plot(xs, result.slope*xs + result.intercept, 'r', label=label if with_label else None)
    if with_ideal_line:
        # Plot y = x line
        ax.plot([xs[0], xs[-1]], [ys[0], ys[-1]], color='k', linestyle='-',
                linewidth=1, label='Ideal' if with_label else None)
    return result


def make_square_axes(ax_mat):
    """
    Make an axes square in screen units.
    Should be called after plotting.
    """
    return
    for ax in ax_mat.flat:
        #ax.set_aspect(1 / ax.get_data_ratio())
        #ax.set(adjustable='box', aspect='equal')
        ax.set(adjustable='datalim', aspect='equal')
    #ax.set_aspect(1 / ax.get_data_ratio())


class AseResultsComparator:
    """
    This object allows one to compare energies, forces and stressee computed
    for the same structure but with different methods e.g. results obtained
    with different ML potentials.
    """

    ALL_VOIGT_COMPS = "xx yy zz yz xz xy".split()

    @classmethod
    def pickle_load(cls, workdir):
        """
        Reconstruct the object from a pickle file located in workdir.
        """
        with open(Path(workdir) / f"{cls.__name__}.pickle", "rb") as fh:
            return pickle.load(fh)

    @classmethod
    def from_ase_results(cls, keys: list[str], results_list: list[list[AseResults]]):
        """
        Build object from list of keys and list of AseResults.
        """
        if len(keys) != len(results_list):
            raise ValueError(f"{len(keys)=} != {len(results_list)=}")

        # Extract energies.
        ene_list = []
        for i in range(len(keys)):
            ene_list.append(np.array([r.ene for r in results_list[i]]))

        # Extract forces.
        forces_list = []
        for i in range(len(keys)):
            forces_list.append(np.array([r.forces for r in results_list[i]]))

        # Extract stress.
        stress_list = []
        for i in range(len(keys)):
            stress_list.append(np.array([r.get_voigt_stress() for r in results_list[i]]))

        structure = Structure.as_structure(results_list[0][0].atoms)

        return cls(structure, keys, np.array(ene_list), np.array(forces_list), np.array(stress_list))

    def __init__(self, structure, keys, ene_list, forces_list, stress_list):
        """
        Args:
            structure: Structure object.
            keys: List of strings, each string is associated to a different set of energies/forces/stresses.
            ene_list: array of shape (nkeys, nsteps)
            forces_list: array of shape (nkeys,nsteps,natom,3) with Cartesian forces.
            stress_list: array of shape (nkeys, nsteps, 6) with stress in Voigt notation.
        """
        self.structure = structure
        self.keys = keys
        self.ene_list = ene_list         # [nkeys, nsteps]
        self.forces_list = forces_list   # [nkeys, nsteps, natom, 3]
        self.stress_list = stress_list   # [nkeys, nsteps, 6]

        # Consistency check
        nkeys = len(self)
        if self.ene_list.shape != (nkeys, self.nsteps):
            raise ValueError(f"{self.ene_list.shape=} != ({nkeys=}, {self.nsteps=})")
        if self.forces_list.shape != (nkeys, self.nsteps, self.natom, 3):
            raise ValueError(f"{self.forces_list.shape=} != ({nkeys=}, {self.nsteps=}, {self.natom=}, 3)")
        if self.stress_list.shape != (nkeys, self.nsteps, 6):
            raise ValueError(f"{self.stress_list.shape=} != ({nkeys=}, {self.nsteps=}, 6)")

        # Index of the reference key.
        self.iref = 0

    def __len__(self):
        return len(self.keys)

    @lazy_property
    def nsteps(self) -> int:
        """Number of steps in the trajectory."""
        return self.forces_list.shape[1]

    @lazy_property
    def natom(self) -> int:
        """Number of atoms."""
        return len(self.structure)

    def get_key_pairs(self) -> list[tuple]:
        """
        Return list with (key_ref, key) tuple.
        """
        return [(self.keys[self.iref], key) for ik, key in enumerate(self.keys) if ik != self.iref]

    def inds_of_keys(self, key1: str, key2: str) -> tuple[int,int]:
        """Tuple with the indices of key1, key2"""
        ik1 = self.keys.index(key1)
        ik2 = self.keys.index(key2)
        return ik1, ik2

    def idir_from_direction(self, direction: str) -> int:
        """Index from direction string."""
        idir = {"x": 0, "y": 1, "z": 2}[direction]
        return idir

    def ivoigt_from_comp(self, voigt_comp: str) -> int:
        iv = "xx yy zz yz xz xy".split().index(voigt_comp)
        return iv

    def get_aligned_energies_traj(self, istep=-1) -> np.ndarray:
        """
        Return energies in eV aligned with respect to self.iref key.
        Use the energy at the `istep` step index.
        """
        out_ene_list = self.ene_list.copy()
        for i in range(len(self)):
            if i == self.iref: continue
            shift = self.ene_list[i,istep] - self.ene_list[self.iref,istep]
            out_ene_list[i] -= shift

        return out_ene_list

    def xy_energies_for_keys(self, key1: str, key2: str, sort=True) -> tuple:
        """
        Return (xs, ys) sorted arrays with aligned energies for (key1, key2).
        """
        aligned_ene_list = self.get_aligned_energies_traj()
        ik1, ik2 = self.inds_of_keys(key1, key2)
        xs = aligned_ene_list[ik1]
        ys = aligned_ene_list[ik2]

        return zip_sort(xs, ys) if sort else (xs, ys)

    def xy_forces_for_keys(self, key1, key2, direction) -> tuple:
        """
        Return (xs, ys), sorted arrays with forces along the cart direction for (key1, key2).
        """
        idir = self.idir_from_direction(direction)
        ik1, ik2 = self.inds_of_keys(key1, key2)
        xs = self.forces_list[ik1,:,:,idir].flatten()
        ys = self.forces_list[ik2,:,:,idir].flatten()

        return zip_sort(xs, ys)

    def traj_forces_for_keys(self, key1, key2) -> tuple:
        """
        Return arrays with the cart direction of forces along the trajectory for (key1, key2).
        """
        ik1, ik2 = self.inds_of_keys(key1, key2)
        xs, ys = self.forces_list[ik1], self.forces_list[ik2]

        return xs, ys

    def xy_stress_for_keys(self, key1, key2, voigt_comp, sort=True) -> tuple:
        """
        Return xs, ys sorted arrays with the stress along the voigt component for (key1, key2).
        """
        # (nkeys, self.nsteps, 6)
        iv = self.ivoigt_from_comp(voigt_comp)
        ik1, ik2 = self.inds_of_keys(key1, key2)
        xs = self.stress_list[ik1,:,iv].flatten()
        ys = self.stress_list[ik2,:,iv].flatten()

        return zip_sort(xs, ys) if sort else (xs, ys)

    def get_forces_dataframe(self) -> pd.DataFrame:
        """
        Return dataFrame with columns (fx, fy, fz, isite, istep, key)
        """
        # [nkeys, nsteps, natom, 3]
        d_list = []
        for ik, key in enumerate(self.keys):
            f_traj = self.forces_list[ik]
            for istep in range(self.nsteps):
                for iatom in range(self.natom):
                    fx, fy, fz = f_traj[istep, iatom,:]
                    d = dict(fx=fx, fy=fy, fz=fz, iatom=iatom, istep=istep, key=key)
                    d_list.append(d)

        df = pd.DataFrame(d_list).sort_values(by=["istep", "iatom"], ignore_index=True)
        return df

    def get_stress_dataframe(self) -> pd.DataFrame:
        """
        Return DataFrame with columns [sxx,syy,szz, ... ,istep,key]
        """
        # [nkeys, nsteps, 6]
        d_list = []
        for ik, key in enumerate(self.keys):
            stress_traj = self.stress_list[ik]
            for istep in range(self.nsteps):
                d = {comp: stress_traj[istep, ic] for ic, comp in enumerate(self.ALL_VOIGT_COMPS)}
                d.update(istep=istep, key=key)
                d_list.append(d)

        df = pd.DataFrame(d_list).sort_values(by=["istep"], ignore_index=True)
        return df

    @add_fig_kwargs
    def plot_energies(self, fontsize=8, **kwargs):
        """
        Compare energies aligned wrt to self.iref entry
        """
        key_pairs = self.get_key_pairs()
        nrows, ncols = 1, len(key_pairs)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False,
                                               #subplot_kw=dict(box_aspect=1) #, layout="constrained",
                                               #aspect='equal', adjustable='box',
                                               )
        irow = 0
        for icol, (key1, key2) in enumerate(key_pairs):
           xs, ys = self.xy_energies_for_keys(key1, key2)
           stats = diff_stats(xs, ys)
           ax = ax_mat[irow, icol]
           ax.scatter(xs, ys, marker="o")
           ax.grid(True)
           ax.set_xlabel(f"{key1} energy", fontsize=fontsize)
           ax.set_ylabel(f"{key2} energy", fontsize=fontsize)
           linear_fit_ax(ax, xs, ys, fontsize=fontsize, with_label=True)
           ax.legend(loc="best", shadow=True, fontsize=fontsize)
           if irow == 0:
               ax.set_title(f"{key1}/{key2} MAE: {stats.MAE:.6f}", fontsize=fontsize)

        if "title" not in kwargs: fig.suptitle(f"Energies in eV for {self.structure.latex_formula}")
        #make_square_axes(ax_mat)
        return fig

    @add_fig_kwargs
    def plot_forces(self, fontsize=8, **kwargs):
        """
        Compare forces.
        """
        key_pairs = self.get_key_pairs()
        nrows, ncols = 3, len(key_pairs)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False,
                                               #subplot_kw=dict(box_aspect=1), # layout="constrained",
                                               #aspect='equal', adjustable='box',
                                               )

        for icol, (key1, key2) in enumerate(key_pairs):
            for irow, direction in enumerate(("x", "y", "z")):
                xs, ys = self.xy_forces_for_keys(key1, key2, direction)
                stats = diff_stats(xs, ys)
                ax = ax_mat[irow, icol]
                ax.scatter(xs, ys, marker="o")
                ax.grid(True)
                linear_fit_ax(ax, xs, ys, fontsize=fontsize, with_label=True)
                ax.legend(loc="best", shadow=True, fontsize=fontsize)
                f_tex = f"$F_{direction}$"
                if icol == 0:
                    ax.set_ylabel(f"{key2} {f_tex}", fontsize=fontsize)
                if irow == 2:
                    ax.set_xlabel(f"{key1} {f_tex}", fontsize=fontsize)
                ax.set_title(f"{key1}/{key2} MAE: {stats.MAE:.6f}", fontsize=fontsize)

        if "title" not in kwargs: fig.suptitle(f"Cartesian forces in ev/Ang for {self.structure.latex_formula}")
        #make_square_axes(ax_mat)
        return fig

    @add_fig_kwargs
    def plot_stresses(self, fontsize=6, **kwargs):
        """
        Compare stress components.
        """
        key_pairs = self.get_key_pairs()
        nrows, ncols = 6, len(key_pairs)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False,
                                               #subplot_kw=dict(box_aspect=1), # layout="constrained",
                                               #aspect='equal', adjustable='box',
                                               )

        for icol, (key1, key2) in enumerate(key_pairs):
            for irow, voigt_comp in enumerate(self.ALL_VOIGT_COMPS):
                xs, ys = self.xy_stress_for_keys(key1, key2, voigt_comp)
                stats = diff_stats(xs, ys)
                ax = ax_mat[irow, icol]
                ax.scatter(xs, ys, marker="o")
                linear_fit_ax(ax, xs, ys, fontsize=fontsize, with_label=True)
                ax.legend(loc="best", shadow=True, fontsize=fontsize)
                ax.grid(True)
                s_tex = "$\sigma_{%s}$" % voigt_comp
                if icol == 0:
                    ax.set_ylabel(f"{key2} {s_tex}", fontsize=fontsize)
                if irow == (len(self.ALL_VOIGT_COMPS) - 1):
                    ax.set_xlabel(f"{key1} {s_tex}", fontsize=fontsize)
                ax.set_title(f"{key1}/{key2} MAE: {stats.MAE:.6f}", fontsize=fontsize)

        if "title" not in kwargs: fig.suptitle(f"Stresses in (eV/Ang^2) for {self.structure.latex_formula}")
        #make_square_axes(ax_mat)
        return fig

    @add_fig_kwargs
    def plot_energies_traj(self, delta_mode=True, fontsize=6, markersize=2, **kwargs):
        """
        Plot energies along the trajectory.

        Args:
            delta_mode: True to plot differences instead of absolute values.
        """
        key_pairs = self.get_key_pairs()
        nrows, ncols = 1, len(key_pairs)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False,
                                               #subplot_kw=dict(box_aspect=1), layout="constrained",
                                               )

        for icol, (key1, key2) in enumerate(key_pairs):
            e1, e2 = self.xy_energies_for_keys(key1, key2, sort=False)
            stats = diff_stats(e1, e2)

            ax = ax_mat[0, icol]
            if delta_mode:
                # Plot delta energy along the trajectory.
                ax.plot(np.abs(e1 - e2), marker="o", markersize=markersize)
                ax.set_yscale("log")
            else:
                ax.plot(e1, marker="o", color="red",  label=key1, markersize=markersize)
                ax.plot(e2, marker="o", color="blue", label=key2, markersize=markersize)

            set_grid_legend(ax, fontsize, xlabel='trajectory',
                            ylabel="$|\Delta_E|$" if delta_mode else "$E$",
                            grid=True, legend_loc="upper left",
                            title=f"{key1}/{key2} MAE: {stats.MAE:.6f} eV")

        head = "$\Delta$-Energy in eV" if delta_mode else "Energy in eV"
        if "title" not in kwargs: fig.suptitle(f"{head} for {self.structure.latex_formula}")

        return fig

    @add_fig_kwargs
    def plot_forces_traj(self, delta_mode=True, fontsize=6, markersize=2, **kwargs):
        """
        Plot forces along the trajectory.

        Args:
            delta_mode: True to plot differences instead of absolute values.
        """
        # Fx,Fy,Fx along rows, pairs along columns.
        key_pairs = self.get_key_pairs()
        nrows, ncols = 3, len(key_pairs)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False,
                                               #subplot_kw=dict(box_aspect=1), layout="constrained",
                                               )

        atom1_cmap = plt.get_cmap("viridis")
        atom2_cmap = plt.get_cmap("jet")
        marker_idir = {0: ">", 1: "<", 2: "^"}

        for icol, (key1, key2) in enumerate(key_pairs):
            # Arrays of shape: [nsteps, natom, 3]
            f1_tad, f2_tad = self.traj_forces_for_keys(key1, key2)
            for idir, direction in enumerate(("x", "y", "z")):
                last_row = idir == 2
                fp_tex = f"F_{direction}"
                xs, ys = self.xy_forces_for_keys(key1, key2, direction)
                stats = diff_stats(xs, ys)
                ax = ax_mat[idir, icol]
                ax.set_title(f"{key1}/{key2} MAE: {stats.MAE:.6f}", fontsize=fontsize)

                zero_values = False
                for iatom in range(self.natom):
                    if delta_mode:
                        # Plot delta of forces along the trajectory.
                        style = dict(marker=marker_idir[idir], markersize=markersize,
                                     color=atom1_cmap(float(iatom) / self.natom))
                        abs_delta = np.abs(f1_tad[:,iatom,idir] - f2_tad[:,iatom,idir])
                        zero_values = zero_values or np.any(abs_delta == 0.0)
                        ax.plot(abs_delta, **style, label=f"$\Delta {fp_tex}$" if iatom == 0 else None)
                    else:
                        f1_style = dict(marker=marker_idir[idir], markersize=markersize,
                                        color=atom1_cmap(float(iatom) / self.natom))
                        f2_style = dict(marker=marker_idir[idir], markersize=markersize,
                                        color=atom2_cmap(float(iatom) / self.natom))
                        with_label = (iatom, idir) == (0,0)
                        ax.plot(f1_tad[:,iatom,idir], **f1_style,
                                label=f"${key1}\, {fp_tex}$" if with_label else None)
                        ax.plot(f2_tad[:,iatom,idir], **f2_style,
                                label=f"${key2}\, {fp_tex}$" if with_label else None)

                if delta_mode:
                    ax.set_yscale("log" if not zero_values else "symlog")

                set_grid_legend(ax, fontsize, xlabel='trajectory' if last_row else None,
                                grid=True, legend=not delta_mode, legend_loc="upper left",
                                ylabel=f"$|\Delta {fp_tex}|$" if delta_mode else f"${fp_tex}$")

        head = "$\Delta$-forces in eV/Ang" if delta_mode else "Forces in eV/Ang"
        if "title" not in kwargs: fig.suptitle(f"{head} for {self.structure.latex_formula}")

        return fig

    @add_fig_kwargs
    def plot_stress_traj(self, delta_mode=True, markersize=2, fontsize=6, **kwargs):
        """
        Plot stresses along the trajectory.

        Args:
            delta_mode: True to plot differences instead of absolute values.
        """
        # Sxx,Syy,Szz,... along rows, pairs along columns.
        key_pairs = self.get_key_pairs()
        nrows, ncols = 6, len(key_pairs)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False,
                                               #subplot_kw=dict(box_aspect=1), layout="constrained",
                                               )

        marker_voigt = {"xx": ">", "yy": "<", "zz": "^", "yz": 1, "xz": 2, "xy":3}

        for icol, (key1, key2) in enumerate(key_pairs):
            # Plot (delta of) stresses along the trajectory.
            for iv, voigt_comp in enumerate(self.ALL_VOIGT_COMPS):
                xs, ys = self.xy_stress_for_keys(key1, key2, voigt_comp)
                stats = diff_stats(xs, ys)
                last_row = iv == len(self.ALL_VOIGT_COMPS) - 1
                ax = ax_mat[iv, icol]
                ax.set_title(f"{key1}/{key2} MAE: {stats.MAE:.6f}", fontsize=fontsize)
                voigt_comp_tex = "{" + voigt_comp + "}"
                s1, s2 = self.xy_stress_for_keys(key1, key2, voigt_comp, sort=False)
                s_style = dict(marker=marker_voigt[voigt_comp], markersize=markersize)

                if delta_mode:
                    abs_delta_stress = np.abs(s1 - s2)
                    ax.plot(abs_delta_stress, **s_style, label=f"$|\Delta \sigma_{voigt_comp_tex}|$")
                    ax.set_yscale("log")
                else:
                    ax.plot(s1, **s_style, label=f"${key1}\,\sigma_{voigt_comp_tex}$" if iv == 0 else None)
                    ax.plot(s2, **s_style, label=f"${key2}\,\sigma_{voigt_comp_tex}$" if iv == 0 else None)
                    #ax.set_ylim(-1, +1)

                set_grid_legend(ax, fontsize, xlabel='trajectory' if last_row else None,
                                grid=True, legend=not delta_mode, legend_loc="upper left",
                                ylabel=f"$|\Delta \sigma_{voigt_comp_tex}|$ " if delta_mode else "$\sigma$ ")

        head = r"$\Delta \sigma$ (eV/Ang$^2$)" if delta_mode else "Stress tensor (eV/Ang$^2$)"
        if "title" not in kwargs: fig.suptitle(f"{head} for {self.structure.latex_formula}")

        return fig


class AseRelaxation:
    """
    Container with the results produced by the ASE calculator.
    """
    def __init__(self, dyn, traj_path):
        self.dyn = dyn
        self.traj_path = str(traj_path)

    @lazy_property
    def traj(self):
        """ASE trajectory."""
        if self.traj_path is None:
            raise RuntimeError("Cannot read ASE traj as traj_path is None")
        from ase.io import read
        return read(self.traj_path, index=":")

    #def __str__(self):
    #def to_string(self, verbose=0)

    def summarize(self, tags=None, mode="smart", stream=sys.stdout):
        """"""
        if self.traj_path is None: return
        r0, r1 = AseResults.from_traj_inds(self.traj, 0, -1)
        if tags is None: tags = ["unrelaxed", "relaxed"],
        df = dataframe_from_results_list(tags, [r0, r1], mode=mode)
        print_dataframe(df, end="\n", file=stream)

    #def plot(self, **kwargs):


def dataframe_from_results_list(index: list, results_list: list[AseResults],
                                mode="smart") -> pd.DataFrame:
    assert len(index) == len(results_list)
    df = pd.DataFrame([r.get_dict4pandas() for r in results_list], index=index)

    if mode == "smart":
        # Remove columns with the same values e.g. geometry params.
        def is_unique(s):
            a = s.to_numpy()
            return (a[0] == a).all()

        for k in (["volume",] + _CELLPAR_KEYS):
            if k in df and is_unique(df[k]):
                df.drop(columns=k, inplace=True)

    return df


def ase_optimizer_cls(s: str | Optimizer) -> Type | list[str]:
    """
    Return an ASE Optimizer subclass from string `s`.
    If s == "__all__", return list with all Optimizer subclasses supported by ASE.
    """
    from ase import optimize
    def is_ase_optimizer(key: str) -> bool:
        return isclass(obj := getattr(optimize, key)) and issubclass(obj, Optimizer)

    valid_keys = [key for key in dir(optimize) if is_ase_optimizer(key)]

    if s == "__all__":
        return valid_keys

    if isinstance(s, Optimizer):
        return s

    if s not in valid_keys:
        raise ValueError(f"Unknown optimizer {s}, must be one of {valid_keys}")

    return getattr(optimize, s)


def relax_atoms(atoms: Atoms, relax_mode: str, optimizer: str, fmax: float, pressure: float,
                verbose: int, steps: int = 500,
                opt_kwargs=None, traj_path=None, calculator=None) -> AseRelaxation:
    """
    Relax atoms using an ASE calculator and ASE algorithms.

    Args:
        atoms: ASE atoms.
        relax_mode: "ions" to relax ions only, "cell" for ions + cell, "no" for no relaxation.
        optimizer: name of the ASE optimizer to use.
        fmax: tolerance for relaxation convergence. Here fmax is a sum of force and stress forces.
        pressure: Target pressure.
        verbose: whether to print stdout.
        steps: max number of steps for relaxation.
        opt_kwargs (dict): kwargs for the ASE optimizer class.
        traj_path:
        calculator:
    """
    from ase.constraints import ExpCellFilter
    from ase.io import read

    RX_MODE.validate(relax_mode)
    if relax_mode == RX_MODE.no:
        raise ValueError(f"Invalid {relax_mode:}")

    opt_kwargs = opt_kwargs or {}
    if traj_path is not None:
        opt_kwargs["trajectory"] = str(traj_path)

    if calculator is not None:
        atoms.calc = calculator

    # Run relaxation
    opt_class = ase_optimizer_cls(optimizer)
    stream = sys.stdout if verbose else io.StringIO()
    def pf(*args, **kwargs):
        print(*args, file=stream, **kwargs)

    with contextlib.redirect_stdout(stream):
        pf(f"Relaxation parameters: fmax: {fmax}, relax_mode: {relax_mode}, steps: {steps}, optimizer: {optimizer}")
        if atoms.constraints and verbose > 1:
            # Print constraints.
            pf(f"Number of constraints: {len(atoms.constraints)}")
            for c in atoms.constraints:
                pf("\t", c)
            pf("")

        dyn = opt_class(ExpCellFilter(atoms, scalar_pressure=pressure), **opt_kwargs) if relax_mode == RX_MODE.cell else \
              opt_class(atoms, **opt_kwargs)

        t_start = time.time()
        converged = dyn.run(fmax=fmax, steps=steps)
        t_end = time.time()
        pf("Converged:", converged)
        pf('Relaxation completed in %2.4f sec\n' % (t_end - t_start))

    return AseRelaxation(dyn, traj_path)


def silence_tensorflow() -> None:
    """
    Silence every unnecessary warning from tensorflow.
    """
    # https://stackoverflow.com/questions/35911252/disable-tensorflow-debugging-information
    import logging
    logging.getLogger('tensorflow').setLevel(logging.ERROR)
    os.environ["KMP_AFFINITY"] = "noverbose"
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    try:
        import tensorflow as tf
        tf.get_logger().setLevel('ERROR')
        tf.autograph.set_verbosity(3)
    except (ModuleNotFoundError, ImportError):
        pass


class _MyCalculatorMixin:
    """
    Add _delta_forces_list and _delta_stress_list internal attributes to an ASE calculator.
    Extend `calculate` method so that forces and stresses are corrected accordingly.
    """

    def set_delta_forces(self, delta_forces):
        """F_Abinitio - F_ML"""
        if getattr(self, "_delta_forces_list", None) is None:
            self._delta_forces_list = []
        self._delta_forces_list.append(delta_forces)

    def get_delta_forces(self):
        delta_forces_list = getattr(self, "_delta_forces_list", None)
        if delta_forces_list is not None: return delta_forces_list[-1]
        return None

    def set_delta_stress(self, delta_stress):
        """S_Abinitio - S_ML"""
        if getattr(self, "_delta_stress_list", None) is None:
            self._delta_stress_list = []
        self._delta_stress_list.append(delta_stress)

    def get_delta_stress(self):
        delta_stress_list = getattr(self, "_delta_stress_list", None)
        if delta_stress_list is not None: return delta_stress_list[-1]
        return None

    def calculate(
         self,
         atoms: Atoms | None = None,
         properties: list | None = None,
         system_changes: list | None = None,
     ):
        """
        Perform calculation for an input Atoms.

        Args:
            atoms (ase.Atoms): ase Atoms object
            properties (list): list of properties to calculate
            system_changes (list): monitor which properties of atoms were
                changed for new calculation. If not, the previous calculation
                results will be loaded.
        """
        super().calculate(atoms=atoms, properties=properties, system_changes=system_changes)

        # Apply delta correction to forces.
        forces = self.results["forces"]
        delta_forces = self.get_delta_forces()
        if delta_forces is not None:
            #print("Updating forces with delta_forces:\n", forces)
            forces += delta_forces
            self.results.update(
                forces=forces,
            )

        # Apply delta correction to stress.
        stress = self.results["stress"]
        delta_stress = self.get_delta_stress()
        if delta_stress is not None:
            #print("Updating stress with delta_stress:\n", stress)
            stress += delta_stress
            self.results.update(
                stress=stress,
            )


def as_calculator(obj) -> Calculator:
    """Build a calculator."""
    if isinstance(obj, Calculator):
        return obj

    # Assume string
    return CalcBuilder(obj).get_calculator()


class CalcBuilder:
    """
    Factory class to build an ASE calculator with ML potential
    Supports different backends defined by `name` string.
    """

    ALL_NN_TYPES = [
        "m3gnet",
        "matgl",
        "chgnet",
        "alignn",
        #"quip",
    ]

    def __init__(self, name: str, **kwargs):
        self.name = name

        # Extract nn_type and model_name from name
        self.nn_type, self.model_name = name, None
        if ":" in name:
            self.nn_type, self.model_name = name.split(":")

        if  self.nn_type not in self.ALL_NN_TYPES:
            raise ValueError(f"Invalid {name=}, it should be in {self.ALL_NN_TYPES=}")

        self._model = None

    def __str__(self):
        if self.model_name is not None:
            return f"{self.__class__.__name__} nn_type: {self.nn_type}, model_name: {self.model_name}"
        else:
            return f"{self.__class__.__name__} nn_type: {self.nn_type}"

    # pickle support.
    def __getstate__(self):
        return dict(name=self.name)

    def __setstate__(self, d):
        self.name = d["name"]
        self._model = None

    def get_calculator(self) -> Calculator:
        """
        Return ASE calculator with ML potential.
        """

        if self.nn_type == "m3gnet":
            # m3gnet legacy version.
            if self._model is None:
                silence_tensorflow()
            try:
                from m3gnet.models import Potential, M3GNet, M3GNetCalculator
            except ImportError as exc:
                raise ImportError("m3gnet not installed. Try `pip install m3gnet`.") from exc

            if self._model is None:
                assert self.model_name is None
                self._model = Potential(M3GNet.load())

            class MyM3GNetCalculator(M3GNetCalculator, _MyCalculatorMixin):
                """Add delta_forces and delta_stress"""

            return MyM3GNetCalculator(potential=self._model)

        if self.nn_type == "matgl":
            # See https://github.com/materialsvirtuallab/matgl
            try:
                import matgl
                from matgl.ext.ase import M3GNetCalculator
            except ImportError as exc:
                raise ImportError("matgl not installed. Try `pip install matgl`.") from exc

            if self._model is None:
                model_name = "M3GNet-MP-2021.2.8-PES" if self.model_name is None else self.model_name
                self._model = matgl.load_model(model_name)

            class MyM3GNetCalculator(M3GNetCalculator, _MyCalculatorMixin):
                """Add delta_forces and delta_stress"""

            return MyM3GNetCalculator(potential=self._model)

        if self.nn_type == "chgnet":
            try:
                from chgnet.model.dynamics import CHGNetCalculator
                from chgnet.model.model import CHGNet
            except ImportError as exc:
                raise ImportError("chgnet not installed. Try `pip install chgnet`.") from exc

            if self._model is None:
                assert self.model_name is None
                self._model = CHGNet.load()

            class MyCHGNetCalculator(CHGNetCalculator, _MyCalculatorMixin):
                """Add delta_forces and delta_stress"""

            return MyCHGNetCalculator(model=self._model)

        if self.nn_type == "alignn":
            try:
                from alignn.ff.ff import AlignnAtomwiseCalculator, default_path
            except ImportError as exc:
                raise ImportError("alignn not installed. See https://github.com/usnistgov/alignn") from exc

            class MyAlignnCalculator(AlignnAtomwiseCalculator, _MyCalculatorMixin):
                """Add delta_forces and delta_stress"""

            model_name = default_path() if self.model_name is None else self.model_name
            return AlignnAtomwiseCalculator(path=model_name)

        #if self.nn_type == "quip":
        #    try:
        #        from quippy.potential import Potential
        #    except ImportError as exc:
        #        raise ImportError("quippy not installed. Try `pip install quippy-ase`.\n" +
        #                          "See https://github.com/libAtoms/QUIP") from exc

        #    class MyQuipPotential(Potential, _MyCalculatorMixin):
        #        """Add delta_forces and delta_stress"""

        #    assert self.model_name is None
        #    args_str = ""
        #    return MyQuipPotential(args_str="")

        raise ValueError(f"Invalid {self.name=}")


class _MlBase:
    """
    Base class for all Ml subclasses providing helper methods to
    perform typical tasks such as writing files in the workdir
    and object persistence via pickle.
    """
    @classmethod
    def pickle_load(cls, workdir):
        """
        Reconstruct the object from a pickle file located in workdir.
        """
        with open(Path(workdir) / f"{cls.__name__}.pickle", "rb") as fh:
            return pickle.load(fh)

    def __init__(self, workdir, prefix=None):
        """
        Build directory with `prefix` if `workdir` is None else create it.
        Raise RuntimeError if workdir already exists.
        """
        self.workdir = workdir_with_prefix(workdir, prefix)
        self.basename_info = []
        self.delta_forces = None
        self.delta_stress = None

    def pickle_dump(self):
        """Write pickle file for object persistence."""
        with open(self.workdir / f"{self.__class__.__name__}.pickle", "wb") as fh:
            pickle.dump(self, fh)

    def set_delta_forces_stress(self, delta_forces, delta_stress) -> None:
        """Set the value of the delta corrections."""
        self.delta_forces = delta_forces
        self.delta_stress = delta_stress

    def __str__(self):
        # Delegated to the subclass.
        return self.to_string()

    @lazy_property
    def calc_builder(self):
        return CalcBuilder(self.nn_name)

    def get_calculator_name(self) -> tuple[Calculator, str]:
        return self.get_calculator(), self.calc_builder.name

    def get_calculator(self) -> Calculator:
        """Return ASE calculator."""
        calc = self.calc_builder.get_calculator()

        if self.delta_forces is not None:
            #print("Setting delta_forces:\n", self.delta_forces)
            calc.set_delta_forces(self.delta_forces)

        if self.delta_stress is not None:
            #print("Setting delta_stress:\n", self.delta_stress)
            calc.set_delta_stress(self.delta_stress)

        return calc

    def add_basename_info(self, basename: str, info: str) -> None:
        """
        Register basename with info in the internal buffer used to generate
        the README.md file in _finalize. Print WARNING if basename is already registered.
        """
        if any(basename == t[0] for t in self.basename_info):
            print(f"WARNING: {basename:} already in basename_info:")
        self.basename_info.append((basename, info))

    def mkdir(self, basename: str, info: str) -> Path:
        """Create directory in workdir, return Path object."""
        self.add_basename_info(basename, info)
        dirpath = self.workdir / basename
        dirpath.mkdir()
        return dirpath

    def get_path(self, basename: str, info: str) -> Path:
        """Return Path in workdir."""
        self.add_basename_info(basename, info)
        return self.workdir / str(basename)

    def savefig(self, basename: str, fig, info: str) -> None:
        """Save matplotlib figure in workdir."""
        self.add_basename_info(basename, info)
        fig.savefig(self.workdir / basename)

    def write_traj(self, basename: str, traj, info: str) -> None:
        """Write ASE trajectory in workdir."""
        self.add_basename_info(basename, info)
        with open(self.workdir / basename, "wb") as fd:
            write_traj(fd, traj)

    def write_json(self, basename: str, data, info: str,
                   indent=4, stream=None, **kwargs) -> None:
        """Write data in JSON format and mirror output to `stream`."""
        self.add_basename_info(basename, info)
        with open(self.workdir / basename, "wt") as fh:
            json.dump(data, fh, indent=indent, **kwargs)

        if stream is not None:
            # Print JSON to stream as well.
            print("", file=stream)
            print(marquee(info, mark="="), file=stream)
            print(json.dumps(data, indent=4), file=stream, end="\n")

    def write_df(self, df, basename: str, info: str, fmt="csv") -> None:
        """Write dataframe to file."""
        self.add_basename_info(basename, info)
        filepath = self.workdir / basename
        if fmt == "csv":
            df.to_csv(filepath)
        else:
            raise ValueError(f"Invalid format {fmt=}")

    def write_script(self, basename: str, text: str, info: str) -> Path:
        """
        Write text script to basename file.
        """
        self.add_basename_info(basename, info)
        _, ext = os.path.splitext(basename)
        shebang = {
            ".py": "#!/usr/bin/env python",
            ".sh": "#!/bin/bash",
        }[ext]

        header = ""
        if "python" in shebang:
            header = """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
"""

        path = self.workdir / basename
        with path.open("wt") as fh:
            fh.write(f"""\
{shebang}

# {info}

{header}

{text}
""")
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        return path

    def _finalize(self) -> None:
        """Called at the end of the `run` method to write the README.md file in the workdir."""
        if self.basename_info:
            # Generate README.md file.
            md_lines = ["## Directory content\n",]
            for path, info in self.basename_info:
                path = os.path.basename(str(path))
                md_lines.append(f"- `{path}`: {info}")

            md_str = "\n".join(md_lines)
            with open(self.workdir / "README.md", "wt") as fh:
                fh.write(md_str)
            print("\n", md_str, end=2*"\n")

            # Print WARNINGs if files do not exist.
            for basename, _ in self.basename_info:
                p = self.workdir / basename
                if not p.exists():
                    print(f"WARNING: Cannot find `{basename}` in {self.workdir}")

        print("\nResults available in directory:", self.workdir)


class MlRelaxer(_MlBase):
    """
    Relax structure with ASE and ML-potential.
    """

    @classmethod
    def from_abinit_yaml_file(cls, filepath: str, workdir=None, prefix=None) -> MlRelaxer:
        """
        Build object from a YAML file produced by ABINIT in hybrid relaxation mode.
        """
        # Read yaml file produced by Abinit:
        #
        #  iteration_state: {dtset: 1, itime: 1, icycle: 1, }
        #  comment   : Summary of ground state results
        #  lattice_vectors:
        #  - [ -5.1690735,  -5.1690735,   0.0000000, ]
        #  - [ -5.1690735,   0.0000000,  -5.1690735, ]
        #  - [  0.0000000,  -5.1690735,  -5.1690735, ]
        #  lattice_lengths: [   7.31017,    7.31017,    7.31017, ]
        #  lattice_angles: [ 60.000,  60.000,  60.000, ] # degrees, (23, 13, 12)
        #  lattice_volume:   2.7622826E+02
        #  convergence: {deltae: -1.926E-11, res2:  3.761E-10, residm:  1.588E-05, diffor: null, }
        #  etotal    :  -8.46248947E+00
        #  entropy   :   0.00000000E+00
        #  fermie    :   1.42500714E-01
        #  cartesian_stress_tensor: # hartree/bohr^3
        #  - [  3.06384355E-06,   0.00000000E+00,   0.00000000E+00, ]
        #  - [  0.00000000E+00,   3.06384355E-06,   0.00000000E+00, ]
        #  - [  0.00000000E+00,   0.00000000E+00,   3.06384355E-06, ]
        #  pressure_GPa:  -9.0141E-02
        #  xred      :
        #  - [  5.0000E-01,   5.0000E-01,   5.0000E-01, Si]
        #  - [  2.5000E-01,   2.5000E-01,   2.5000E-01, Si]
        #  cartesian_forces: # hartree/bohr
        #  - [  5.08705549E-32,   5.08705549E-32,  -1.52611665E-31, ]
        #  - [ -5.08705549E-32,  -5.08705549E-32,   1.52611665E-31, ]
        #  force_length_stats: {min:   1.68718543E-31, max:   1.68718543E-31, mean:   1.68718543E-31, }
        #  format_version: 1
        #  natom: 2
        #  ionmov: 1
        #  optcell: 2
        #  nn_name: matgl
        #  prtvol: 1

        from ruamel.yaml import YAML
        with open(os.path.expanduser(filepath), "rt") as fh:
            doc = YAML().load(fh) #; print(doc)

        format_version = doc.pop("format_version")
        natom = doc.pop("natom")
        ntypat = doc.pop("ntypat")
        typat = np.array(doc.pop("typat"), dtype=int)
        znucl = np.array(doc.pop("znucl"), dtype=float)
        rprim = np.array(doc.pop("lattice_vectors"))
        xred = np.array(doc.pop("xred"))
        xred = np.array(xred[:,:3], dtype=float)
        abi_cart_forces = np.array(doc.pop("cartesian_forces")) * abu.Ha_eV / abu.Bohr_Ang
        abi_cart_stresses = np.array(doc.pop("cartesian_stress_tensor")) * abu.Ha_eV / (abu.Bohr_Ang**3)

        ionmov = doc.pop("ionmov")
        optcell = doc.pop("optcell")
        iatfix = doc.pop("iatfix") # [3,natom] array or None if unconstrained.
        strtarget = np.array(doc.pop("strtarget"), dtype=float)
        nn_name = doc.pop("nn_name", default="chgnet")
        verbose = doc.pop("prtvol")

        structure = Structure.from_abivars(
            acell=3*[1.0],
            rprim=rprim,
            typat=typat,
            xred=xred,
            ntypat=ntypat,
            znucl=znucl,
        ) #; print(structure)

        atoms = structure.to_ase_atoms()
        if iatfix is not None:
            raise NotImplementedError()
            #aseml.fix_atoms(atoms, fix_inds=fix_inds, fix_symbols=fix_symbols)

        ######################################################################
        # Consistency check as not all the Abinit options are supported by ASE
        ######################################################################
        relax_mode = RX_MODE.cell if optcell != 0 else RX_MODE.ions

        allowed_optcells = (0, 2)
        if optcell not in allowed_optcells:
            raise ValueError(f"{optcell=} not in {allowed_optcells=}")

        # Target pressure is taken from strtarget.
        # The components of the stress tensor are stored in a.u. according to:
        # (1,1) → 1; (2,2) → 2; (3,3) → 3; (2,3) → 4; (3,1) → 5; (1,2) → 6.
        pressure = -strtarget[0] * abu.HaBohr3_GPa
        if np.any(strtarget[:3] != strtarget[0]):
            raise ValueError(f"Only hydrostatic stress in strtarget is supported. {strtarget=}")
        if np.any(strtarget[3:] != 0.0):
            raise ValueError(f"Off diagonal components in strtarget are not supported. {strtarget=}")

        # Set internal parameters according to yaml file and build object.
        fmax, steps, optimizer = 0.01, 500, "BFGS"

        new = cls(atoms, relax_mode, fmax, pressure, steps, optimizer, nn_name, verbose,
                  workdir=workdir, prefix=prefix)

        # Set delta forces and delta stress if script is called by ABINIT.

        #new.set_delta_forces_stress_from_abiml_nc(filepath)
        ## Compute ML forces for this structure.
        #atoms = to_ase_atoms(structure, calc=new.get_calculator())
        #ml_cart_forces = atoms.get_forces()
        #ml_cart_stress = atoms.get_stress()

        ## Set delta forces/stresses so that the next invokation to get_calculator include the deltas
        #new.set_delta_forces_stress(abi_cart_forces - ml_cart_forces, abi_stress - ml_cart_stress)

        return new

    def write_output_file_for_abinit(self) -> Path:
        """
        Write output file with results in a format that can be parsed by ABINIT.
        Return path to output file.
        """
        filepath = self.get_path("ABI_MLRELAXER.out", "Output file for hybrid relaxation with ABINIT.")
        format_version = 1

        def fmt_vec3(vec) -> str:
            return "{:.12e} {:.12e} {:.12e}".format(*vec)

        with open(filepath, "wt") as fh:
            fh.write("%i # format_version\n" % format_version)
            fh.write("%i # natom\n" % len(self.atoms))
            # Write lattice vectors.
            rprimd = self.atoms.cell.array * abu.Ang_Bohr
            for i in range(3):
                fh.write("%s # lattice vector %i\n" % (fmt_vec3(rprimd[i]), i+1))
            # Write relaxed fractional coordinates.
            fh.write("xred\n")
            for atom in self.atoms:
                fh.write(fmt_vec3(atom.scaled_position) + "\n")

        return filepath

    def __init__(self, atoms: Atoms, relax_mode, fmax, pressure, steps, optimizer, nn_name, verbose,
                 workdir, prefix=None):
        """
        Args:
            atoms: ASE atoms to relax.
            relax_mode:
            fmax: tolerance for relaxation convergence. Here fmax is a sum of force and stress forces.
            pressure: Target pressure.
            steps: max number of steps for relaxation.
            optimizer: name of the ASE optimizer to use.
            nn_name:
            verbose: Verbosity level.
        """
        super().__init__(workdir, prefix)
        self.atoms = atoms
        self.relax_mode = relax_mode
        RX_MODE.validate(relax_mode)
        self.fmax = fmax
        self.steps = steps
        self.optimizer = optimizer
        self.pressure = pressure
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     relax_mode  = {self.relax_mode}
     fmax        = {self.fmax}
     steps       = {self.steps}
     optimizer   = {self.optimizer}
     pressure    = {self.pressure}
     nn_name     = {self.nn_name}
     workdir     = {self.workdir}
     verbose     = {self.verbose}

=== ATOMS ===

{self.atoms}

"""

    def run(self):
        """Run structural relaxation."""
        #self.pickle_dump()
        workdir = self.workdir
        self.atoms.calc = self.get_calculator()

        print(f"Relaxing structure with relax mode: {self.relax_mode} ...")
        relax_kws = dict(calculator=self.atoms.calc,
                         optimizer=self.optimizer,
                         relax_mode=self.relax_mode,
                         fmax=self.fmax,
                         pressure=self.pressure,
                         steps=self.steps,
                         traj_path=self.get_path("relax.traj", "ASE relaxation trajectory"),
                         verbose=1,
                        )

        relax = relax_atoms(self.atoms, **relax_kws)
        relax.summarize(tags=["unrelaxed", "relaxed"])

        # Write files with final structure and dynamics.
        formats = ["poscar",]
        outpath_fmt = write_atoms(self.atoms, workdir, self.verbose, formats=formats)
        for outp, fmt in outpath_fmt:
            self.add_basename_info(outp.name, f"Final structure in {fmt} format.")

        label = "xdatcar with structural relaxation"
        write_vasp_xdatcar(self.get_path("XDATCAR", label), relax.traj, label=label)

        # Write output file for Abinit
        self.write_output_file_for_abinit()

        self._finalize()
        return relax


class MlMd(_MlBase):
    """Perform MD calculations with ASE and ML potential."""

    def __init__(self, atoms: Atoms, temperature, timestep, steps, loginterval,
                 ensemble, nn_name, verbose, workdir, prefix=None):
        """
        Args:
            atoms:
            temperature:
            timestep:
            steps:
            loginterval:
            ensemble:
            nn_name:
            verbose: Verbosity level.
            workdir:
            prefix:
        """
        super().__init__(workdir, prefix)
        self.atoms = atoms
        self.temperature = temperature
        self.timestep = timestep
        self.steps = steps
        self.loginterval = loginterval
        self.ensemble = ensemble
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

    temperature = {self.temperature} K
    timestep    = {self.timestep} fs
    steps       = {self.steps}
    loginterval = {self.loginterval}
    ensemble    = {self.ensemble}
    calculator  = {self.calc_builder}
    workdir     = {self.workdir}
    verbose     = {self.verbose}

=== ATOMS ===

{self.atoms}

"""

    def run(self) -> None:
        """Run MD"""
        #self.pickle_dump()
        workdir = self.workdir
        self.atoms.calc = self.get_calculator()

        traj_file = self.get_path("md.traj", "ASE MD trajectory")
        logfile = self.get_path("md.log", "ASE MD log file")

        md = MolecularDynamics(
            atoms=self.atoms,
            ensemble=self.ensemble,
            temperature=self.temperature,   # K
            timestep=self.timestep,         # fs,
            #pressure,
            trajectory=str(traj_file),      # save trajectory to md.traj
            logfile=str(logfile),           # log file for MD
            loginterval=self.loginterval,   # interval for record the log
            #append_trajectory,
        )

        self.write_script("diffusion_coeff.py", text=f"""\
from ase.md.analysis import DiffusionCoefficient
from ase.io import read

# For an MD simulation with timestep of N, and images written every M iterations, our timestep here is N * M.
timestep = {self.timestep} * {self.loginterval}
traj = read("{str(traj_file)}", index=":")
dc = DiffusionCoefficient(traj, timestep, atom_indices=None, molecule=False)
dc.calculate(ignore_n_images=0, number_of_segments=1)
dc.print_data()
dc.plot(ax=None, show=True)
""", info="Python script to compute and visualize diffusion coefficients.")

        self.write_script("plot_energies.py", text=f"""\
df = pd.read_csv("{str(logfile)}", sep="\s+")
print(df)
xname = "Time[ps]"
ynames = [k for k in df.keys() if k != xname]
print("=== Summary statistics ===")
print(df[ynames].describe())

axes = df.plot.line(x=xname, y=ynames, subplots=True)
fig = axes[0].get_figure()
plt.show()
""", info="Python script to visualize energies vs Time.")

        md.run(steps=self.steps)

        #trajectory = read(traj_file, index=":")
        #write_vasp_xdatcar(workdir / "XDATCAR", trajectory,
        #                   label=f"xdatcar with relaxation generated by {self.__class__.__name__}")


class _MlNebBase(_MlBase):
    """
    Base class for Neb calculations
    """

    def postprocess_images(self, images):
        """
        post-process ASE NEB calculation.
        See <https://wiki.fysik.dtu.dk/ase/tutorials/neb/diffusion.html>
        """
        from ase.neb import NEBTools
        nebtools = NEBTools(images)

        # get the actual maximum force at this point in the simulation.
        max_force = nebtools.get_fmax()
        # get the calculated barrier and the energy change of the reaction.
        ef, de = nebtools.get_barrier()

        neb_data = dict(max_force=float(max_force),
                        energies_images=[float(image.get_potential_energy()) for image in images],
                        barrier_with_fit=float(ef),
                        energy_change_with_fit=float(de),
                        )

        # get the barrier without any interpolation between highest images.
        ef, de = nebtools.get_barrier(fit=False)
        neb_data.update(barrier_without_fit=float(ef),
                        energy_change_without_fit=float(de),
        )

        self.write_json("neb_data.json", neb_data, info="JSON document with NEB results",
                        stream=sys.stdout if self.verbose else None)

        # create a figure like that coming from ase-gui.
        self.savefig("neb_barrier.png", nebtools.plot_band(), info="Figure with NEB barrier")
        return neb_data

    def read_neb_data(self) -> dict:
        """
        Read results from the JSON file produced by postprocess_images
        """
        with open(self.workdir / 'neb_data.json', "rt") as fh:
            return json.load(fh)


class MlGsList(_MlNebBase):
    """
    Perform ground-state calculations for a list of atoms with ASE and ML-potential.
    Inherits from _MlNebBase so that we can reuse postprocess_images and read_neb_data.
    """

    def __init__(self, atoms_list: list[Atoms], nn_name, verbose,
                 workdir, prefix=None):
        """
        Args:
            atoms_list: List of ASE atoms
            nn_name:
            verbose: Verbosity level.
        """
        super().__init__(workdir, prefix)
        self.atoms_list = atoms_list
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     nn_name  = {self.nn_name}
     workdir  = {self.workdir}
     verbose  = {self.verbose}

"""

    def run(self) -> None:
        """Run list of GS calculations."""
        #self.pickle_dump()
        workdir = self.workdir

        results = []
        for ind, atoms in enumerate(self.atoms_list):
            write_vasp(self.workdir / f"IND_{ind}_POSCAR", atoms, label=None)
            atoms.calc = self.get_calculator()
            results.append(AseResults.from_atoms(atoms))

        write_vasp_xdatcar(self.workdir / "XDATCAR", self.atoms_list,
                           label=f"XDATCAR with list of atoms.")

        self.postprocess_images(self.atoms_list)
        self._finalize()


class MlNeb(_MlNebBase):
    """
    Perform NEB calculation with ASE and ML potential.
    """

    def __init__(self, initial_atoms: Atoms, final_atoms: Atoms,
                 nimages, neb_method, climb, optimizer, relax_mode, fmax, pressure,
                 nn_name, verbose, workdir, prefix=None):
        """
        Args:
            initial_atoms
            final_atoms:
            nimages:
            neb_method:
            climb:
            optimizer:
            relax_mode:
            fmax:
            pressure:
            nn_name:
            verbose:
            workdir:
            prefix:
        """
        super().__init__(workdir, prefix)
        self.initial_atoms = get_atoms(initial_atoms)
        self.final_atoms = get_atoms(final_atoms)
        self.nimages = nimages
        self.neb_method = neb_method
        if self.neb_method not in ASENEB_METHODS:
            raise ValueError(f"{self.neb_method} not in {ASENEB_METHODS}")
        self.climb = climb
        self.optimizer = optimizer
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.pressure = pressure
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        s = f"""\

{self.__class__.__name__} parameters:

     nimages     = {self.nimages}
     neb_method  = {self.neb_method}
     climb       = {self.climb}
     optimizer   = {self.optimizer}
     pressure    = {self.pressure}
     relax_mode  = {self.relax_mode}
     fmax        = {self.fmax}
     nn_name     = {self.nn_name}
     workdir     = {self.workdir}
     verbose     = {self.verbose}

=== INITIAL ATOMS ===

{self.initial_atoms}

=== FINAL ATOMS ===

{self.final_atoms}

"""
        if verbose:
            #s += scompare_two_atoms("initial image", self.initial_atoms, "final image", self.final_atoms)
            file = io.StringIO()
            fmt = "poscar"
            diff_two_structures("initial image", self.initial_atoms,
                                "final image", self.final_atoms, fmt, file=file)
            s += "\n" + file.getvalue()
        return s

    def run(self) -> None:
        """Run NEB"""
        #self.pickle_dump()
        workdir = self.workdir
        initial_atoms, final_atoms = self.initial_atoms, self.final_atoms

        if self.relax_mode != RX_MODE.no:
            relax_kws = dict(calculator=self.get_calculator(),
                             optimizer=self.optimizer,
                             relax_mode=self.relax_mode,
                             fmax=self.fmax,
                             pressure=self.pressure,
                             verbose=self.verbose,
                             )

            print(f"Relaxing initial image with relax mode: {self.relax_mode} ...")
            relax = relax_atoms(initial_atoms,
                                traj_path=self.get_path("initial_relax.traj", "ASE Relaxation of the first image"),
                                **relax_kws)

            relax.summarize(tags=["initial_unrelaxed", "initial_relaxed"])

            print(f"Relaxing final image with relax mode: {self.relax_mode} ...")
            relax = relax_atoms(final_atoms,
                                traj_path=self.get_path("final_relax.traj", "ASE Relaxation of the last image"),
                                **relax_kws)

            relax.summarize(tags=["final_unrelaxed", "final_relaxed"])

        # Generate several instances of the calculator. It is probably fine to have just one, but just in case...
        calculators = [self.get_calculator() for i in range(self.nimages)]
        neb = make_ase_neb(initial_atoms, final_atoms, self.nimages, calculators, self.neb_method, self.climb,
                           method='linear', mic=False)

        write_vasp_xdatcar(workdir / "INITIAL_NEB_XDATCAR", neb.images,
                           label=f"XDATCAR with initial NEB images.")

        # Optimize
        opt_class = ase_optimizer_cls(self.optimizer)
        nebtraj_file = str(workdir / "neb.traj")
        logfile = self.get_path("neb.log", "Log file of NEB calculation.")
        optimizer = opt_class(neb, trajectory=nebtraj_file, logfile=logfile)
        print("Starting NEB algorithm with optimizer:", opt_class, "...")
        optimizer.run(fmax=self.fmax)

        # To read the last nimages atoms e.g. 5: read('neb.traj@-5:')
        images = ase.io.read(f"{str(nebtraj_file)}@-{self.nimages}:")
        write_vasp_xdatcar(workdir / "FINAL_NEB_XDATCAR", images,
                           label=f"XDATCAR with final NEB images.")

        # write vasp poscar files for each image in vasp_neb
        dirpath = self.mkdir("VASP_NEB", info="Directory with POSCAR files for each NEB image.")
        for im, image in enumerate(images):
            subdir = dirpath / str(im).zfill(2)
            subdir.mkdir()
            ase.io.write(subdir / "POSCAR", image, format="vasp")

        neb_data = self.postprocess_images(images)

        self.write_script("ase_gui.sh", text=f"""\
# To visualize the results, use:

ase gui {nebtraj_file}@-{self.nimages}

# then select `tools->neb` in the gui.
""", info="Shell script to visualize NEB results with ase gui")

        self.write_script("ase_nebplot.sh", text=f"""\
# This command create a series of plots showing the progression of the neb relaxation

ase nebplot --share-x --share-y --nimages {self.nimages} {nebtraj_file}
""", info="Shell script to create a series of plots showing the progression of the neb relaxation")

        self._finalize()


class MultiMlNeb(_MlNebBase):
    """
    Perform a multi-NEB calculation with ASE and ML potential.
    """

    def __init__(self, atoms_list: list[Atoms], nimages, neb_method, climb, optimizer, relax_mode, fmax, pressure,
                 nn_name, verbose, workdir, prefix=None):
        """
        Args:
            atoms_list:
            nimages:
            neb_method:
            climb:
            optimizer:
            relax_mode:
            fmax:
            pressure:
            nn_name:
            verbose:
            workdir:
            prefix:
        """
        super().__init__(workdir, prefix)
        self.atoms_list = atoms_list
        self.nimages = nimages
        self.neb_method = neb_method
        self.climb = climb
        self.optimizer = optimizer
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.pressure = pressure
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        s = f"""\

{self.__class__.__name__} parameters:

     nimages     = {self.nimages}
     neb_method  = {self.neb_method}
     climb       = {self.climb}
     optimizer   = {self.optimizer}
     relax_mode  = {self.relax_mode}
     fmax        = {self.fmax}
     pressure    = {self.pressure}
     nn_name     = {self.nn_name}
     workdir     = {self.workdir}
     verbose     = {self.verbose}
"""
        return s

    def run(self) -> None:
        """
        Run multi NEB calculations.
        """
        #self.pickle_dump()
        workdir = self.workdir
        atoms_list = self.atoms_list
        camp_dirs = [workdir / f"CAMP_{i}" for i in range(len(atoms_list) - 1)]

        energies = []
        for i in range(len(atoms_list) - 1):
            ml_neb = MlNeb(atoms_list[i], atoms_list[i+1],
                           self.nimages, self.neb_method, self.climb, self.optimizer,
                           self.relax_mode, self.fmax, self.pressure, self.nn_name, self.verbose, camp_dirs[i])
            ml_neb.run()

            # Read energies from json files and remove first/last point depending on CAMP index..
            data = ml_neb.read_neb_data()
            enes = data['energies_images']
            if i == 0: enes = enes[:-1]
            if i == len(camp_dirs) - 1: enes = enes[1:]
            energies.extend(enes)

        #print("energies", energies)
        ax, fig, plt = get_ax_fig_plt()
        ax.plot(energies, marker="o")
        ax.set_xlabel('Path index')
        ax.set_ylabel('Energy [eV]')
        ef = max(energies) - energies[0]
        er = max(energies) - energies[-1]
        de = energies[-1] - energies[0]
        ax.set_title(r'$E_\mathrm{{f}} \approx$ {:.3f} eV; '
                     r'$E_\mathrm{{r}} \approx$ {:.3f} eV; '
                     r'$\Delta E$ = {:.3f} eV'.format(ef, er, de))
        self.savefig("neb_barrier.png", fig, info="Figure with NEB barrier")

        self._finalize()


def make_ase_neb(initial: Atoms, final: Atoms, nimages: int,
                 calculators: list, neb_method: str, climb: bool,
                 method='linear', mic=False) -> NEB:
    """
    Make a NEB band consisting of nimages. See https://databases.fysik.dtu.dk/ase/ase/neb.html

    Args:
        initial: First point.
        final: Last point.
        nimages: Number of images.
        calculators: List of ASE calculators.
        neb_method: String defining NEB algorithm.
        climb: True to use a climbing image.
        method: str
            Method by which to interpolate: 'linear' or 'idpp'.
            linear provides a standard straight-line interpolation, while
            idpp uses an image-dependent pair potential.
        mic: Map movement into the unit cell by using the minimum image convention.
    """
    images = [initial]
    images += [initial.copy() for i in range(nimages - 2)]
    images += [final]

    apply_constraint = None
    if initial.constraints:
        if not final.constraints:
            raise RuntimeError("Both initial and final points should have constraints!")
        if len(initial.constraints) != len(final.constraints):
            raise RuntimeError("different number of constraints in initial and final")
        for ci, cf in zip(initial.constraints, final.constraints):
            if ci.__class__ != cf.__class__:
                raise RuntimeError(f"Constraints in initial and final points should belong to the same class: {ci}, {cf}")
        apply_constraint = True

    # Set calculators
    for image, calculator in zip(images, calculators): #, strict=True):
        image.calc = calculator

    # Compute energy/forces for the extrema in order to have them in the trajectory.
    _ = AseResults.from_traj_inds(images, 0, -1)

    neb = NEB(images, method=neb_method, climb=climb)
    # Interpolate linearly the positions of the middle images
    neb.interpolate(method=method, mic=mic, apply_constraint=apply_constraint)

    return neb


class MlOrderer(_MlBase):
    """
    Order a disordered structure using pymatgen and ML potential.
    """
    def __init__(self, structure, max_ns, optimizer, relax_mode, fmax, pressure, steps, nn_name, verbose,
                 workdir, prefix=None):
        """
        Args:
            structure:
            max_ns:
            optimizer:
            relax_mode:
            fmax:
            pressure:
            steps:
            nn_name:
            verbose:
            workdir:
            prefix:
        """
        super().__init__(workdir, prefix)
        self.structure = Structure.as_structure(structure)
        self.max_ns = max_ns
        self.optimizer = optimizer
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.pressure = pressure
        self.steps = steps
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        s = f"""\

{self.__class__.__name__} parameters:

     max_ns      = {self.max_ns}
     optimizer   = {self.optimizer}
     relax_mode  = {self.relax_mode}
     fmax        = {self.fmax}
     pressure    = {self.pressure}
     steps       = {self.steps}
     nn_name     = {self.nn_name}
     workdir     = {self.workdir}
     verbose     = {self.verbose}


=== STRUCTURE ===

{self.structure}

"""
        return s

    def run(self) -> None:
        """
        Run MlOrderer.
        """
        #self.pickle_dump()
        workdir = self.workdir
        from pymatgen.core import Lattice
        specie = {"Cu0+": 0.5, "Au0+": 0.5}
        structure = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.677), [specie], [[0, 0, 0]])
        #structure = self.structure
        #print(structure)

        # Each dict in d_list contains the following entries:
        #{
        #    "energy": output[0],
        #    "energy_above_minimum": (output[0] - lowest_energy) / num_atoms,
        #    "structure": s_copy.get_sorted_structure(),
        #}
        from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation
        trans = OrderDisorderedStructureTransformation()
        d_list = trans.apply_transformation(structure, return_ranked_list=max(self.max_ns, 2))
        print("Number of structures after OrderedDisordered:", len(d_list))
        if self.verbose > 2:
            for d in d_list:
                print(d)

        # Note that the OrderDisorderedTransformation (with a sufficiently large return_ranked_list parameter)
        # returns all orderings, including duplicates without accounting for symmetry.
        # A computed ewald energy is returned together with each structure.
        # To eliminate duplicates, the best way is to use StructureMatcher's group_structures method
        from pymatgen.analysis.structure_matcher import StructureMatcher
        matcher = StructureMatcher()

        # Add ew_pos index to structures to faciliatate reindexing after sorting.
        ew_structures = [d["structure"] for d in d_list]
        for ew_pos, s in enumerate(ew_structures):
            s.ew_pos = ew_pos
        ew_energies = [d["energy"] for d in d_list]
        ew_energies_above_minimum = [d["energy_above_minimum"] for d in d_list]

        groups = matcher.group_structures(ew_structures)
        print("Number of structures after StructureMatcher:", len(groups))
        if self.verbose > 2:
            for group in groups:
                print(group[0])

        if self.relax_mode != RX_MODE.no:
            print(f"Relaxing structures with relax mode: {self.relax_mode}")
            relax_kws = dict(calculator=self.get_calculator(),
                             optimizer=self.optimizer,
                             relax_mode=self.relax_mode,
                             fmax=self.fmax,
                             pressure=self.pressure,
                             return_trajectory=True,
                             verbose=self.verbose,
                            )

            rows = []
            for group in groups:
                s = group[0]
                #print("s.ew_pos:", s.ew_pos)
                #relax = relax_atoms(self.atoms, **relax_kws)
                rel_s, trajectory = s.relax(**relax_kws)
                r0, r1 = AseResults.from_traj_inds(trajectory, 0, -1)
                df = dataframe_from_results_list(["unrelaxed", "relaxed"], [r0, r1])
                print(df, end=2*"\n")

                rows.append(dict(
                    ew_unrelaxed_energy=ew_energies[s.ew_pos],
                    unrelaxed_energy=r0.ene,
                    relaxed_energy=r1.ene,
                    relaxed_structure=rel_s,
                    #unrelaxed_pressure=r0.pressure
                    #relaxed_pressure=r1.pressure
                ))

            df = pd.DataFrame(rows).sort_values("relaxed_energy")
            print(df.drop("relaxed_structure", axis=1))

        # TODO: Post-process
        self._finalize()


class MlPhonons(_MlBase):
    """Compute phonons with ASE and ML potential."""

    def __init__(self, atoms: Atoms, supercell, kpts, asr, nqpath,
                 relax_mode, fmax, pressure, steps, optimizer, nn_name,
                 verbose, workdir, prefix=None):
        """
        Args:
            atoms: ASE atoms.
            supercell: tuple with supercell dimension.
            kpts:
            asr: Enforce acoustic sum-rule.
            nqpath: Number of q-point along the q-path.
            relax_mode:
            fmax:
            steps:
            optimizer:
            nn_name,
            verbose:
            workdir:
            prefix:
        """
        super().__init__(workdir, prefix)
        self.atoms = get_atoms(atoms)
        self.supercell = supercell
        self.kpts = kpts
        self.asr = asr
        self.nqpath = nqpath
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.pressure = pressure
        self.steps = steps
        self.optimizer = optimizer
        self.nn_name = nn_name
        self.verbose = verbose

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        s = f"""\

{self.__class__.__name__} parameters:

     supercell  = {self.supercell}
     kpts       = {self.kpts}
     asr        = {self.asr}
     nqpath     = {self.nqpath}
     relax_mode = {self.relax_mode}
     fmax       = {self.fmax}
     steps      = {self.steps}
     optimizer  = {self.optimizer}
     pressure   = {self.pressure}
     nn_name    = {self.nn_name}
     workdir    = {self.workdir}
     verbose    = {self.verbose}

=== ATOMS ===

{self.atoms}
"""
        return s

    def run(self) -> None:
        """Run MlPhonons."""
        #self.pickle_dump()
        workdir = self.workdir
        calculator = self.get_calculator()
        atoms = self.atoms

        if self.relax != RX_MODE.no:
            print(f"Relaxing atoms with relax mode: {self.relax_mode}.")
            relax_kws = dict(calculator=calculator,
                             optimizer=self.optimizer,
                             relax_mode=self.relax_mode,
                             fmax=self.fmax,
                             pressure=self.pressure,
                             steps=self.steps,
                             traj_path=self.get_path("relax.traj", "ASE relax trajectory"),
                             verbose=self.verbose,
                            )

            relax = relax_atoms(atoms, **relax_kws)
            #self.write_traj("relax.traj", traj, info="")

            r0, r1 = AseResults.from_traj_inds(relax.traj, 0, -1)
            df = dataframe_from_results_list(["initial_unrelaxed", "initial_relaxed"], [r0, r1])
            print(df, end=2*"\n")

        # Phonon calculator
        from ase.phonons import Phonons
        ph = Phonons(atoms, calculator, supercell=self.supercell, delta=0.05)
        #ph.read_born_charges(name=, neutrality=True)
        ph.run()
        #print("Phonons Done")

        # Read forces and assemble the dynamical matrix
        ph.read(acoustic=self.asr, born=False)
        ph.clean()

        # Calculate phonon dispersion along a path in the Brillouin zone.
        path = atoms.cell.bandpath(npoints=self.nqpath)
        bs = ph.get_band_structure(path, born=False)
        dos = ph.get_dos(kpts=self.kpts).sample_grid(npts=100, width=1e-3)

        # Plot the band structure and DOS
        import matplotlib.pyplot as plt
        plt.rc("figure", dpi=150)
        fig = plt.figure(1, figsize=(7, 4))
        bs_ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
        emax = 0.035
        bs.plot(ax=bs_ax, emin=0.0, emax=emax)
        dos_ax = fig.add_axes([0.8, 0.07, 0.17, 0.85])
        dos_ax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color="grey", edgecolor="black", lw=1)
        dos_ax.set_ylim(0, emax)
        dos_ax.set_yticks([])
        dos_ax.set_xticks([])
        dos_ax.set_xlabel("PH DOS", fontsize=14)

        #title = f"Phonon band structure and DOS of {atoms.symbols} with supercell: {self.supercell}",
        title = f"Phonon band structure and DOS with supercell: {self.supercell}",
        fig.suptitle(title, fontsize=8, y=1.02)
        self.savefig("phonons.png", fig, info=title)

        self._finalize()



class MlCompareWithAbinitio(_MlNebBase):
    """
    Compare ab-initio energies, forces and stresses with ML results.
    """

    def __init__(self, filepaths, nn_names, traj_range, verbose, workdir, prefix=None):
        """
        Args:
            filepaths: List of file produced by the ab-initio code with energies, forces and stresses.
            nn_names: String or list of strings defining the NN potential.
            traj_range: Trajectory range. None to include all steps.
            verbose: Verbosity level.
            workdir: Working directory.
        """
        super().__init__(workdir, prefix)
        self.filepaths = list_strings(filepaths)
        self.traj_range = traj_range
        self.nn_names = list_strings(nn_names)
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     filepaths   = {self.filepaths}
     traj_range  = {self.traj_range}
     nn_names    = {self.nn_names}
     workdir     = {self.workdir}
     verbose     = {self.verbose}

"""

    def get_abinitio_results(self) -> list[AseResults]:
        results = []
        for filepath in self.filepaths:
            results.extend(self._get_results_filepath(filepath))
        return results

    def _get_results_filepath(self, filepath) -> list[AseResults]:
        """
        Extract ab-initio results from self.filepath according to the file extension.
        """
        basename = os.path.basename(filepath)
        abi_results = []
        from fnmatch import fnmatch

        if basename.endswith("_HIST.nc"):
            # Abinit HIST file produced by a structural relaxation.
            from abipy.dynamics.hist import HistFile
            with HistFile(filepath) as hist:
                # etotals in eV units.
                etotals = hist.etotals
                if self.traj_range is None: self.traj_range = range(0, len(hist.etotals), 1)
                forces_hist = hist.r.read_cart_forces(unit="eV ang^-1")
                # GPa units.
                stress_cart_tensors, pressures = hist.reader.read_cart_stress_tensors()
                for istep, (structure, ene, stress, forces) in enumerate(zip(hist.structures, etotals, stress_cart_tensors, forces_hist)):
                    if not istep in self.traj_range: continue
                    r = AseResults(atoms=get_atoms(structure), ene=float(ene), forces=forces, stress=stress)
                    abi_results.append(r)
                return abi_results

        elif fnmatch(basename, "vasprun*.xml*"):
            # Assume Vasprun file with structural relaxation or MD results.
            def get_energy_step(step: dict) -> float:
                """Copied from final_energy property in vasp.outputs."""
                final_istep = step
                total_energy = final_istep["e_0_energy"]
                # Addresses a bug in vasprun.xml. See https://www.vasp.at/forum/viewtopic.php?f=3&t=16942
                final_estep = final_istep["electronic_steps"][-1]
                electronic_energy_diff = final_estep["e_0_energy"] - final_estep["e_fr_energy"]
                total_energy_bugfix = np.round(electronic_energy_diff + final_istep["e_fr_energy"], 8)
                if np.abs(total_energy - total_energy_bugfix) > 1e-7:
                    return total_energy_bugfix
                return total_energy

            from pymatgen.io.vasp.outputs import Vasprun
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                vasprun = Vasprun(filepath)

            num_steps = len(vasprun.ionic_steps)
            if self.traj_range is None: self.traj_range = range(0, num_steps, 1)
            for istep, step in enumerate(vasprun.ionic_steps):
                #print(step.keys())
                if not istep in self.traj_range: continue
                structure, forces, stress = step["structure"], step["forces"], step["stress"]
                ene = get_energy_step(step)
                r = AseResults(atoms=get_atoms(structure), ene=float(ene), forces=forces, stress=stress)
                abi_results.append(r)
            return abi_results

        raise ValueError(f"Don't know how to extract data from: {filepath=}")

    def run(self, nprocs, print_dataframes=True) -> AseResultsComparator:
        """
        Run calculation with nprocs processes.
        """
        #self.pickle_dump()
        workdir = self.workdir

        labels = ["abinitio"]
        abi_results = self.get_abinitio_results()
        results_list = [abi_results]

        ntasks = len(abi_results)
        #nprocs = nprocs_for_ntasks(nprocs, ntasks, title="Begin ML computation")
        nprocs = 1

        for nn_name in self.nn_names:
            labels.append(nn_name)
            # Use ML to compute quantities with the same ab-initio trajectory.
            if nprocs == 1:
                calc = as_calculator(nn_name)
                items = [AseResults.from_atoms(res.atoms, calc=calc) for res in abi_results]
            else:
                raise NotImplementedError("run with multiprocessing!")
                args_list = [(nn_name, res) for res in abi_results]
                with Pool(processes=nprocs) as pool:
                    items = pool.map(_map_run_compare, args_list)

            results_list.append(items)

        comp = AseResultsComparator.from_ase_results(labels, results_list)

        # Write pickle file for object persistence.
        with open(self.workdir / f"{comp.__class__.__name__}.pickle", "wb") as fh:
            pickle.dump(comp, fh)

        py_path = self.workdir / "analyze.py"
        print("Writing python script to analyze the results in:", py_path.name)
        with PythonScript(py_path) as script:
            script.add_text("""
def main():
    from abipy.ml.aseml import AseResultsComparator
    c = AseResultsComparator.pickle_load(".")
    with_stress = True
    from abipy.tools.plotting import Exposer
    exposer = "mpl" # or "panel"
    with Exposer.as_exposer(exposer) as e:
        e(c.plot_energies(show=False))
        e(c.plot_forces(delta_mode=True, show=False))
        e(c.plot_energies_traj(delta_mode=True, show=False))
        e(c.plot_energies_traj(delta_mode=False, show=False))
        if with_stress:
            e(c.plot_stresses(delta_mode=True, show=False))
        e(c.plot_forces_traj(delta_mode=True, show=False))
        e(c.plot_stress_traj(delta_mode=True, show=False))
""").add_main()

        if print_dataframes:
            # Write dataframes to disk in CSV format.
            forces_df = comp.get_forces_dataframe()
            self.write_df(forces_df, "cart_forces.csv", info="CSV file with cartesian forces.")
            stress_df = comp.get_stress_dataframe()
            self.write_df(stress_df, "voigt_stress.csv", info="CSV file with cartesian stresses in Voigt notation")

        self._finalize()
        return comp


#__COMPARE_CALC = None
#
#def _map_run_compare(args: tuple) -> AseResults:
#    """Function passed to pool.map."""
#    nn_name, res = args
#    global __COMPARE_CALC
#    if __COMPARE_CALC is None:
#        #__COMPARE_CALC = as_calculator(nn_name)
#    return AseResults.from_atoms(res.atoms, calc=__COMPARE_CALC)


class MolecularDynamics:
    """
    Molecular dynamics class

    Based on https://github.com/materialsvirtuallab/m3gnet/blob/main/m3gnet/models/_dynamics.py
    """

    def __init__(
        self,
        atoms: Atoms,
        ensemble: str = "nvt",
        temperature: int = 300,
        timestep: float = 1.0,
        pressure: float = 1.01325 * units.bar,
        taut: Optional[float] = None,
        taup: Optional[float] = None,
        compressibility_au: Optional[float] = None,
        trajectory: Optional[Union[str, Trajectory]] = None,
        logfile: Optional[str] = None,
        loginterval: int = 1,
        append_trajectory: bool = False,
    ):
        """
        Args:
            atoms (Atoms): atoms to run the MD
            ensemble (str): choose from 'nvt' or 'npt'. NPT is not tested,
                use with extra caution
            temperature (float): temperature for MD simulation, in K
            timestep (float): time step in fs
            pressure (float): pressure in eV/A^3
            taut (float): time constant for Berendsen temperature coupling
            taup (float): time constant for pressure coupling
            compressibility_au (float): compressibility of the material in A^3/eV
            trajectory (str or Trajectory): Attach trajectory object
            logfile (str): open this file for recording MD outputs
            loginterval (int): write to log file every interval steps
            append_trajectory (bool): Whether to append to prev trajectory
        """
        self.atoms = atoms

        if taut is None:
            taut = 100 * timestep * units.fs
        if taup is None:
            taup = 1000 * timestep * units.fs

        ensemble = ensemble.lower()
        if ensemble == "nvt":
            self.dyn = NVTBerendsen(
                self.atoms,
                timestep * units.fs,
                temperature_K=temperature,
                taut=taut,
                trajectory=trajectory,
                logfile=logfile,
                loginterval=loginterval,
                append_trajectory=append_trajectory,
            )

        elif ensemble == "npt":
            """
            NPT ensemble default to Inhomogeneous_NPTBerendsen thermo/barostat
            This is a more flexible scheme that fixes three angles of the unit
            cell but allows three lattice parameter to change independently.
            """
            self.dyn = Inhomogeneous_NPTBerendsen(
                self.atoms,
                timestep * units.fs,
                temperature_K=temperature,
                pressure_au=pressure,
                taut=taut,
                taup=taup,
                compressibility_au=compressibility_au,
                trajectory=trajectory,
                logfile=logfile,
                loginterval=loginterval,
                # append_trajectory=append_trajectory,
                # this option is not supported in ASE at this point (I have sent merge request there)
            )

        elif ensemble == "npt_berendsen":
            """
            This is a similar scheme to the Inhomogeneous_NPTBerendsen.
            This is a less flexible scheme that fixes the shape of the
            cell - three angles are fixed and the ratios between the three
            lattice constants.
            """
            self.dyn = NPTBerendsen(
                self.atoms,
                timestep * units.fs,
                temperature_K=temperature,
                pressure_au=pressure,
                taut=taut,
                taup=taup,
                compressibility_au=compressibility_au,
                trajectory=trajectory,
                logfile=logfile,
                loginterval=loginterval,
                append_trajectory=append_trajectory,
            )

        else:
            raise ValueError(f"{ensemble=} not supported")

        self.trajectory = trajectory
        self.logfile = logfile
        self.loginterval = loginterval
        self.timestep = timestep

    def run(self, steps: int):
        """
        Thin wrapper of ase MD run

        Args:
            steps (int): number of MD steps
        """
        from ase.md import MDLogger
        self.dyn.attach(MDLogger(self.dyn, self.atoms, '-', header=True, stress=False,
                        peratom=True, mode="a"), interval=self.loginterval)
        self.dyn.run(steps)


def traj_to_qepos(traj_filepath: str, pos_filepath: str) -> None:
    """
    Convert ASE trajectory file to QE POS file.

    Args:
        traj_filepath: Name of ASE trajectory file
        pos_filepath: Name of output POS file.
    """
    traj = Trajectory(traj_filepath)
    nStepsTraj = len(traj)
    nAtoms = len(traj[0])
    #print(nStepsTraj)
    posArray = np.zeros((nStepsTraj, nAtoms, 3), dtype=float)
    # positionsArray = np.zeros((nStepsTraj), dtype=float)

    count = -1
    for atoms in traj:
        count = count + 1
        #print(atoms.positions)
        posArray[count,:,:] = atoms.positions
    #print(posArray.shape)

    with open(pos_filepath, 'w+') as posFile:
        for i in range(nStepsTraj):
            posFile.write(str(i)+ '\n')
            for j in range(nAtoms):
                posFile.write(str(posArray[i,j,0]) + ' ' + str(posArray[i,j,1]) + ' ' + str(posArray[i,j,2]) + '\n')


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
import abipy.core.abinit_units as abu
try:
    import ase
except ImportError as exc:
    raise ImportError("ase not installed. Try `pip install ase`.") from exc
from pathlib import Path
from inspect import isclass
from multiprocessing import Pool
from typing import Type, Any, Optional, Union
from enum import IntEnum
from tabulate import tabulate
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from monty.json import MontyEncoder
from monty.collections import AttrDict
from pymatgen.core import Structure as PmgStructure
from pymatgen.io.ase import AseAtomsAdaptor
from ase import units
from ase.atoms import Atoms
from ase.io.trajectory import write_traj, Trajectory
from ase.io import read
from ase.optimize.optimize import Optimizer
from ase.calculators.calculator import Calculator
from ase.io.vasp import write_vasp_xdatcar, write_vasp
from ase.neb import NEB
from ase.md.nptberendsen import NPTBerendsen, Inhomogeneous_NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution, Stationary, ZeroRotation)
from abipy.core import Structure
from abipy.tools.iotools import workdir_with_prefix, PythonScript, yaml_safe_load_path
from abipy.tools.typing import Figure, PathLike
from abipy.tools.printing import print_dataframe
from abipy.tools.serialization import HasPickleIO
from abipy.tools.context_managers import Timer
from abipy.abio.enums import StrEnum, EnumMixin
from abipy.core.mixins import TextFile, NotebookWriter
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_grid_legend,
    set_visible, set_ax_xylabels, linear_fit_ax)


_CELLPAR_KEYS = ["a", "b", "c", "angle(b,c)", "angle(a,c)", "angle(a,b)"]


ASENEB_METHODS = ['aseneb', 'eb', 'improvedtangent', 'spline', 'string']


class RX_MODE(EnumMixin, StrEnum):  # StrEnum added in 3.11
    """
    Relaxation mode string flags.
    """
    no   = "no"
    ions = "ions"
    cell = "cell"


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
    Fix atoms by indices and/or by symbols.

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


class AseTrajectoryPlotter:
    """
    Plot an ASE trajectory with matplotlib.
    """
    def __init__(self, traj: Trajectory):
        self.traj = traj
        self.natom = len(traj[0])
        self.traj_size = len(traj)

    @classmethod
    def from_file(cls, filepath: PathLike) -> AseTrajectoryPlotter:
        """Initialize an instance from file filepath"""
        return cls(read(filepath, index=":"))

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level verbose."""
        lines = [f"ASE trajectory with {len(self.traj)} configuration(s)."]
        app = lines.append
        if len(self.traj) == 1:
            first = AseResults.from_atoms(self.traj[0])
            app(first.to_string(verbose=verbose))
        else:
            first, last = AseResults.from_atoms(self.traj[0]), AseResults.from_atoms(self.traj[-1])
            app("First configuration:")
            app(first.to_string(verbose=verbose))
            app("Last configuration:")
            app(last.to_string(verbose=verbose))

        return "\n".join(lines)

    @add_fig_kwargs
    def plot(self, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot energies, force stats, and pressure as a function of the trajectory index.
        """
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=3, ncols=1,
                                                sharex=True, sharey=False, squeeze=True)

        # Plot total energy in eV.
        energies = [float(atoms.get_potential_energy()) for atoms in self.traj]
        ax = ax_list[0]
        marker = "o"
        ax.plot(energies, marker=marker)
        ax.set_ylabel('Energy (eV)')

        # Plot Force stats.
        forces_traj = np.reshape([atoms.get_forces() for atoms in self.traj], (self.traj_size, self.natom, 3))
        fmin_steps, fmax_steps, fmean_steps, fstd_steps = [], [], [], []
        for forces in forces_traj:
            fmods = np.sqrt([np.dot(force, force) for force in forces])
            fmean_steps.append(fmods.mean())
            fstd_steps.append(fmods.std())
            fmin_steps.append(fmods.min())
            fmax_steps.append(fmods.max())

        markers = ["o", "^", "v", "X"]
        ax = ax_list[1]
        ax.plot(fmin_steps, label="min |F|", marker=markers[0])
        ax.plot(fmax_steps, label="max |F|", marker=markers[1])
        ax.plot(fmean_steps, label="mean |F|", marker=markers[2])
        #ax.plot(fstd_steps, label="std |F|", marker=markers[3])
        ax.set_ylabel('F stats (eV/A)')
        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        # Plot pressure.
        voigt_stresses_traj = np.reshape([atoms.get_stress() for atoms in self.traj], (self.traj_size, 6))
        pressures = [-sum(vs[0:3])/3 for vs in voigt_stresses_traj]
        ax = ax_list[2]
        ax.plot(pressures, marker=marker)
        ax.set_ylabel('Pressure (GPa)')

        for ix, ax in enumerate(ax_list):
            set_axlims(ax, xlims, "x")
            ax.grid(True)
            if ix == len(ax_list) - 1:
                ax.set_xlabel('Trajectory index', fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_lattice(self, ax_list=None,
                     fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot lattice lengths/angles/volume as a function the of the trajectory index.

        Args:
            ax_list: List of axis or None if a new figure should be created.
            fontsize: fontsize for legends and titles
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
        """
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=3, ncols=1,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        def cell_dict(atoms):
            return dict(zip(_CELLPAR_KEYS, atoms.cell.cellpar()))

        cellpar_list = [cell_dict(atoms) for atoms in self.traj]
        df = pd.DataFrame(cellpar_list)
        #print(df)

        # plot lattice parameters.
        ax = ax_list[0]
        markers = ["o", "^", "v"]
        for i, label in enumerate(["a", "b", "c"]):
            ax.plot(df[label].values, label=label, marker=markers[i])
        ax.set_ylabel("abc (A)")

        # plot lattice angles.
        ax = ax_list[1]
        for i, label in enumerate(["angle(b,c)", "angle(a,c)", "angle(a,b)"]):
            ax.plot(df[label].values, label=label, marker=markers[i])
        ax.set_ylabel(r"$\alpha\beta\gamma$ (degree)")

        # plot lattice volume.
        ax = ax_list[2]
        volumes = [atoms.get_volume() for atoms in self.traj]
        marker = "o"
        ax.plot(volumes, label="Volume", marker=marker)
        ax.set_ylabel(r'$V\, (A^3)$')

        for ix, ax in enumerate(ax_list):
            set_axlims(ax, xlims, "x")
            if ix == len(ax_list) - 1:
                ax.set_xlabel('Trajectory index', fontsize=fontsize)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig


def get_fstats(cart_forces: np.ndarray) -> dict:
    """
    Return dictionary with statistics on cart_forces.
    """
    fmods = np.array([np.linalg.norm(f) for f in cart_forces])
    #fmods = np.sqrt(np.einsum('ij, ij->i', cart_forces, cart_forces))
    #return AttrDict(
    return dict(
        fmin=fmods.min(),
        fmax=fmods.max(),
        fmean=fmods.mean(),
        fstd=fmods.std(),
        drift=np.linalg.norm(cart_forces.sum(axis=0)),
    )


@dataclasses.dataclass
class AseResults(HasPickleIO):
    """
    Container with the results produced by an ASE calculator.
    """
    atoms: Atoms
    ene: float
    stress: np.ndarray   # 3x3 matrix with stress tensor.
    forces: np.ndarray
    magmoms: np.ndarray  # None if calculator does not provide magmoms.

    @classmethod
    def from_traj_inds(cls, trajectory: Trajectory, *inds) -> AseResults:
        """Build list of AseResults from a trajectory and list of indices."""
        return [cls.from_atoms(trajectory[i]) for i in inds]

    @classmethod
    def from_atoms(cls, atoms: Atoms, calc=None) -> AseResults:
        """Build the object from an atoms instance with a calculator."""
        if calc is not None:
            atoms.calc = calc

        from ase.stress import voigt_6_to_full_3x3_strain
        stress_voigt = atoms.get_stress()
        stress = voigt_6_to_full_3x3_strain(stress_voigt)

        from ase.calculators.calculator import PropertyNotImplementedError
        try:
            magmoms = atoms.get_magnetic_moments()
        except PropertyNotImplementedError:
            magmoms = None

        results = cls(atoms=atoms.copy(),
                      ene=float(atoms.get_potential_energy()),
                      stress=stress,
                      forces=atoms.get_forces(),
                      magmoms=magmoms)

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
        app(f"Pressure: {self.pressure} (Gpa)")

        fstats = self.get_fstats()
        for k, v in fstats.items():
            app(f"{k} = {v} (eV/Ang)")

        if True:
        #if verbose:
            app('Forces (eV/Ang):')
            positions = self.atoms.get_positions()
            data = dict(
                x=positions[:,0],
                y=positions[:,1],
                z=positions[:,2],
                fx=self.forces[:,0],
                fy=self.forces[:,1],
                fz=self.forces[:,2],
            )
            # Add magmoms if available.
            if self.magmoms is not None:
                data["magmoms"] = self.magmoms

            df = pd.DataFrame(data)
            app(df.to_string())

        app('Stress tensor:')
        for row in self.stress:
            app(str(row))

        return "\n".join(lines)

    def get_fstats(self) -> dict:
        """
        Return dictionary with statistics on forces.
        """
        return get_fstats(self.forces)

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


class AseResultsComparator(HasPickleIO):
    """
    This object allows one to compare energies, forces and stresses computed
    for the same structure but with different methods e.g. results obtained
    with different ML potentials.
    """

    ALL_VOIGT_COMPS = "xx yy zz yz xz xy".split()

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

    def pickle_dump_and_write_script(self, workdir: PathLike) -> None:
        """
        Write pickle file for object persistence and python script.
        """
        workdir = Path(str(workdir))
        self.pickle_dump(workdir)

        py_path = workdir / "analyze.py"
        print("Writing python script to analyze the results in:", str(py_path))
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

    #def get_rdf(self):
    #    try:
    #        from pymatgen.analysis.diffusion.aimd import RadialDistributionFunction, RadialDistributionFunctionFast
    #    except ImportError as exc:
    #        raise ImportError("pymatgen-analysis-diffusion. Try `pip install pymatgen-analysis-diffusion`.") from exc

    #    rdf = RadialDistributionFcuntion.from_species(
    #        structures: list,
    #        ngrid: int = 101,
    #        rmax: float = 10.0,
    #        cell_range: int = 1,
    #        sigma: float = 0.1,
    #        species: tuple | list = ("Li", "Na"),
    #        reference_species: tuple | list | None = None,
    #    )

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
    def __init__(self, dyn, r0, r1, traj_path):
        self.dyn = dyn
        self.r0, self.r1 = r0, r1
        self.traj_path = str(traj_path)

    @lazy_property
    def traj(self):
        """ASE trajectory."""
        if self.traj_path is None:
            raise RuntimeError("Cannot read ASE traj as traj_path is None")
        return read(self.traj_path, index=":")

    #def __str__(self):
    #def to_string(self, verbose=0)

    def summarize(self, tags=None, mode="smart", stream=sys.stdout):
        """"""
        if self.traj_path is None: return
        r0, r1 = self.r0, self.r1
        #r0, r1 = AseResults.from_traj_inds(self.traj, 0, -1)
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

        r0 = AseResults.from_atoms(atoms)

        dyn = opt_class(ExpCellFilter(atoms, scalar_pressure=pressure), **opt_kwargs) if relax_mode == RX_MODE.cell else \
              opt_class(atoms, **opt_kwargs)

        t_start = time.time()
        converged = dyn.run(fmax=fmax, steps=steps)
        t_end = time.time()
        pf('Relaxation completed in %2.4f sec\n' % (t_end - t_start))
        if not converged:
            raise RuntimeError("ASE relaxation didn't converge")

        r1 = AseResults.from_atoms(atoms)

    return AseRelaxation(dyn, r0, r1, traj_path)


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



class CORRALGO(IntEnum):
    """
    Enumerate the different algorithms used to correct the ML forces/stresses.
    """
    none = 0
    delta = 1
    one_point = 2
    #two_points = 3

    @classmethod
    def from_string(cls, string: str):
        """Build instance from string."""
        try:
           enum = getattr(cls, string)
           return enum
        except AttributeError as exc:
           raise ValueError(f'Error: {string} is not a valid value')


class _MyCalculator:
    """
    Add __abi_forces_list and __abi_stress_list internal attributes to an ASE calculator.
    Extend `calculate` method so that ML forces and stresses can be corrected.
    """

    def __init__(self, *args, **kwargs):
        #print("In _MyMlCalculatorInit with args:", args, ", kwargs:", kwargs)
        super().__init__(*args, **kwargs)

        from collections import deque
        maxlen = 2
        self.__correct_forces_algo = CORRALGO.none
        self.__correct_stress_algo = CORRALGO.none
        self.__abi_forces_list = deque(maxlen=maxlen)
        self.__abi_stress_list = deque(maxlen=maxlen)
        self.__abi_atoms_list = deque(maxlen=maxlen)
        self.__ml_forces_list = deque(maxlen=maxlen)
        self.__ml_stress_list = deque(maxlen=maxlen)
        self.__verbose = 0

    def set_correct_forces_algo(self, new_algo: int) -> int:
        """Set the correction algorithm for forces."""
        assert new_algo in CORRALGO
        old_algo = self.__correct_forces_algo
        self.__correct_forces_algo = new_algo
        return old_algo

    @property
    def correct_forces_algo(self) -> int:
        """Correction algorithm for forces."""
        return self.__correct_forces_algo

    def set_correct_stress_algo(self, new_algo: int) -> int:
        """Set the correction algorithm for the stress."""
        assert new_algo in CORRALGO
        old_algo = self.__correct_stress_algo
        self.__correct_stress_algo = new_algo
        return old_algo

    @property
    def correct_stress_algo(self) -> int:
        """Correction algorithm for the stress."""
        return self.__correct_stress_algo

    def store_abi_forstr_atoms(self, abi_forces, abi_stress, atoms):
        """
        Stores a copy of the ab-initio forces, stress tensor and atoms
        in the internal buffers. Also compute and store the corresponding ML values.
        """
        # Store copies in internal buffers.
        abi_forces = np.asarray(abi_forces).copy()
        self.__abi_forces_list.append(abi_forces)
        abi_stress = np.asarray(abi_stress).copy()
        self.__abi_stress_list.append(abi_stress)
        self.__abi_atoms_list.append(atoms.copy())

        # Compute ML forces and stresses for the input atoms.
        old_forces_algo = self.set_correct_forces_algo(CORRALGO.none)
        old_stress_algo = self.set_correct_stress_algo(CORRALGO.none)
        self.reset()
        ml_forces = self.get_forces(atoms=atoms)
        ml_stress = self.get_stress(atoms=atoms)
        self.reset()
        self.__ml_forces_list.append(ml_forces)
        self.__ml_stress_list.append(ml_stress)
        self.reset()
        self.set_correct_forces_algo(old_forces_algo)
        self.set_correct_stress_algo(old_stress_algo)

        if self.__verbose:
            def fmt_vec3(vec3) -> str:
                return "{:.6e} {:.6e} {:.6e}".format(*vec3)
            def fmt_vec6(vec6) -> str:
                return "{:.6e} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e}".format(*vec6)

            from ase.stress import full_3x3_to_voigt_6_stress
            print("abi_stress6:", fmt_vec6(full_3x3_to_voigt_6_stress(abi_stress)))
            print("ml_stress6: ", fmt_vec6(full_3x3_to_voigt_6_stress(ml_stress)))
            for iat in range(len(atoms)):
                print(f"abi_fcart_{iat=}:", fmt_vec3(abi_forces[iat]))
                print(f"ml_fcart_{iat=} :", fmt_vec3(ml_forces[iat]))

    def get_abi_ml_forces(self, pos=-1):
        """Return the ab-initio and the ML forces or (None, None) if not available."""
        if not self.__abi_forces_list: return (None, None)
        return self.__abi_forces_list[pos], self.__ml_forces_list[pos]

    def get_abi_ml_stress(self, pos=-1):
        """Return the ab-initio and the ML stress or (None, None) if not available."""
        if not self.__abi_stress_list: return (None, None)
        return self.__abi_stress_list[pos], self.__ml_stress_list[pos]

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
        #print("In super.calculate")
        super().calculate(atoms=atoms, properties=properties, system_changes=system_changes)

        if self.correct_forces_algo != CORRALGO.none:
            if self.__verbose: print(f"Applying ab-initio correction to the ml_forces with {self.correct_forces_algo=}")
            forces = self.results["forces"]
            abi_forces, ml_forces = self.get_abi_ml_forces()
            if abi_forces is not None:
                # Change forces only if have invoked store_abi_forstr_atoms
                if self.correct_forces_algo == CORRALGO.delta:
                    # Apply delta correction to forces.
                    alpha = 1.0
                    delta_forces = (abi_forces - ml_forces)/alpha
                    if self.__verbose > 1: print("Updating forces with delta_forces:\n", abi_forces)
                    forces += delta_forces
                    print(f"{delta_forces=}")
                    #AA: TODO: save the delta in list and call method...
                    dict ={'delta_forces': delta_forces,}
                    with open('delta_forces.json', 'a') as outfile:
                        json.dump(dict, outfile,indent=1,cls=MontyEncoder)

                elif self.correct_forces_algo == CORRALGO.one_point:
                    forces += abi_forces
                else:
                    raise ValueError(f"Invalid {self.correct_forces_algo=}")

                self.results.update(forces=forces)

        if self.correct_stress_algo != CORRALGO.none:
            if self.__verbose: print(f"Applying ab-initio correction to the ml_stress with {self.correct_stress_algo=}")
            stress = self.results["stress"]
            abi_stress, ml_stress = self.get_abi_ml_stress()
            if abi_stress is not None:
                # Change stresses only if have invoked store_abi_forstr_atoms
                if self.correct_stress_algo == CORRALGO.delta:
                    # Apply delta correction to stress.
                    delta_stress = abi_stress - ml_stress
                    if self.__verbose > 1: print("Updating stress with delta_stress:\n", delta_stress)
                    stress += delta_stress
                elif self.correct_stress_algo == CORRALGO.one_point:
                    stress += abi_stress
                else:
                    raise ValueError(f"Invalid {self.correct_stress_algo=}")
                self.results.update(stress=stress)


def as_calculator(obj) -> Calculator:
    """Build an ASE calculator."""
    if isinstance(obj, Calculator):
        return obj

    # Assume string
    return CalcBuilder(obj).get_calculator()


def find_spec_nn() -> list[str]:
    # https://stackoverflow.com/questions/14050281/how-to-check-if-a-python-module-exists-without-importing-it
    import importlib.util
    installed = []
    for nn_name in CalcBuilder.ALL_NN_TYPES:
        nn_spec = importlib.util.find_spec(nn_name)
        if nn_spec is not None:
            installed.append(nn_name)

    return installed


def get_installed_nn_names(verbose=0, printout=True) -> tuple[list[str], list[str]]:
    """
    Return list of strings with the names of the nn potentials installed in the environment
    """
    installed, versions = [], []
    from importlib import import_module

    for name in CalcBuilder.ALL_NN_TYPES:
        try:
            mod = import_module(name)
            installed.append(name)
            try:
                v = mod.__version__
            except AttributeError:
                v = "0.0.0"
            versions.append(v)

        except ImportError as exc:
            if verbose: print(exc)

    if printout:
        print("The following NN potentials are installed in the environment:", sys.executable, end=2*"\n")
        table = [t for t in zip(installed, versions)]
        print(tabulate(table, headers=["Package", "Version"]))
        print("")

    return installed, versions


def install_nn_names(nn_names="all", update=False, verbose=0) -> None:
    """
    Install NN potentials in the environment using pip.

    Args:
        nn_names: List of NN potentisl to install.
        update: True if packages should be updated.
        verbose: Verbosity level.
    """
    def pip_install(nn_name) -> int:
        from subprocess import call
        cmd = f"python -m pip install {nn_name}"
        #print("About to execute", cmd)
        try:
            retcode = call(cmd, shell=True)
            if retcode < 0:
                print(f"`{cmd=}` was terminated by signal", -retcode, file=sys.stderr)
            else:
                print(f"`{cmd=}` returned", retcode, file=sys.stderr)
        except OSError as exc:
            print(f"`{cmd=}` failed:", exc, file=sys.stderr)

        return retcode

    black_list = [
        "alignn",
    ]

    nn_names == list_strings(nn_names)
    if "all" in nn_names: nn_names = CalcBuilder.ALL_NN_TYPES
    installed, versions = get_installed_nn_names(verbose=verbose, printout=False)

    for name in nn_names:
        print(f"About to install nn_name={name} ...")
        if name in black_list:
            print("Cannot install {name} with pip!")
            continue
        if name in installed:
            if not update:
                print("{name} is already installed!. Use update to update the package")
                continue
            print("Upgrading: {name}")
        if name not in CalcBuilder.ALL_NN_TYPES:
            print(f"Ignoring unknown {name=}")
            continue

        pip_install(name)


class CalcBuilder:
    """
    Factory class to build an ASE calculator with a ML potential as backend.
    Supports different backends defined by `name` string.
    Possible formats are:

        1) nn_type e.g. m3gnet. See ALL_NN_TYPES for available keys.
        2) nn_type:model_name
        3) nn_type@filepath
    """

    ALL_NN_TYPES = [
        "m3gnet",
        "matgl",
        "chgnet",
        "alignn",
        "mace",
        "pyace",
        "nequip",
        "metatensor",
        "deepmd",
    ]


    def __init__(self, name: str, dftd3_args=None, **kwargs):
        self.name = name

        # Extract nn_type and model_name from name
        self.nn_type, self.model_name, self.model_path = name, None, None

        if ":" in name:
            self.nn_type, self.model_name = name.split(":")
        elif "@" in name:
            self.nn_type, self.model_path = name.split("@")

        if self.nn_type not in self.ALL_NN_TYPES:
            raise ValueError(f"Invalid {name=}, it should be in {self.ALL_NN_TYPES=}")

        # Handle DFTD3.
        self.dftd3_args = dftd3_args
        if self.dftd3_args and not isinstance(dftd3_args, dict):
            # Load parameters from Yaml file.
            self.dftd3_args = yaml_safe_load_path(self.dftd3_args)

        if self.dftd3_args:
            print("Activating dftd3 with arguments:", self.dftd3_args)

        self._model = None

    def __str__(self):
        if self.model_name is not None:
            return f"{self.__class__.__name__} nn_type: {self.nn_type}, model_name: {self.model_name}"

        return f"{self.__class__.__name__} nn_type: {self.nn_type}"

    # pickle support.
    def __getstate__(self):
        return dict(name=self.name)

    def __setstate__(self, d):
        self.name = d["name"]
        self._model = None

    def reset(self) -> None:
        self._model = None

    def get_calculator(self, with_delta: bool = True, reset: bool = False) -> Calculator:
        """
        Return an ASE calculator with ML potential.

        Args:
            with_delta: False if the calculator should not include delta corrections.
            reset: True if the internal cache for the model should be reset.
        """
        if reset: self.reset()

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

            class MyM3GNetCalculator(_MyCalculator, M3GNetCalculator):
                """Add abi_forces and abi_stress"""

            cls = MyM3GNetCalculator if with_delta else M3GNetCalculator
            calc = cls(potential=self._model)

        elif self.nn_type == "matgl":
            # See https://github.com/materialsvirtuallab/matgl
            try:
                import matgl
                from matgl.ext.ase import M3GNetCalculator
            except ImportError as exc:
                raise ImportError("matgl not installed. Try `pip install matgl`.") from exc

            if self._model is None:
                if self.model_path is not None:
                    self._model = matgl.load_model(self.model_path)
                else:
                    model_name = "M3GNet-MP-2021.2.8-PES" if self.model_name is None else self.model_name
                    self._model = matgl.load_model(model_name)

            class MyM3GNetCalculator(_MyCalculator, M3GNetCalculator):
                """Add abi_forces and abi_stress"""

            cls = MyM3GNetCalculator if with_delta else M3GNetCalculator
            calc = cls(potential=self._model)

        elif self.nn_type == "chgnet":
            try:
                from chgnet.model.dynamics import CHGNetCalculator
                from chgnet.model.model import CHGNet
            except ImportError as exc:
                raise ImportError("chgnet not installed. Try `pip install chgnet`.") from exc

            if self._model is None:
                assert self.model_name is None
                if self.model_path is not None:
                    self._model = CHGNet.from_file(self.model_path)
                elif self.model_name is not None:
                    self._model = CHGNet.load(model_name=self.model_name)
                else:
                    # Default model with model_name="MPtrj-efsm"
                    self._model = CHGNet.load()

            class MyCHGNetCalculator(_MyCalculator, CHGNetCalculator):
                """Add abi_forces and abi_stress"""

            cls = MyCHGNetCalculator if with_delta else CHGNetCalculator
            calc = cls(model=self._model)

        elif self.nn_type == "alignn":
            try:
                from alignn.ff.ff import AlignnAtomwiseCalculator, default_path, get_figshare_model_ff
            except ImportError as exc:
                raise ImportError("alignn not installed. See https://github.com/usnistgov/alignn") from exc

            class MyAlignnCalculator(_MyCalculator, AlignnAtomwiseCalculator):
                """Add abi_forces and abi_stress"""

            #if self.model_path is not None:
            #    get_figshare_model_ff(model_name=self.model_path)

            model_name = default_path() if self.model_name is None else self.model_name
            cls = MyAlignnCalculator if with_delta else AlignnAtomwiseCalculator
            calc = cls(path=model_name)

        elif self.nn_type == "pyace":
            try:
                from pyace import PyACECalculator
            except ImportError as exc:
                raise ImportError("pyace not installed. See https://pacemaker.readthedocs.io/en/latest/pacemaker/install/") from exc

            class MyPyACECalculator(_MyCalculator, PyACECalculator):
                """Add abi_forces and abi_stress"""

            if self.model_path is None:
                raise RuntimeError("PyACECalculator requires model_path e.g. nn_name='pyace@FILEPATH'")

            cls = MyPyACECalculator if with_delta else PyACECalculator
            calc = cls(basis_set=self.model_path)

        elif self.nn_type == "mace":
            try:
                from mace.calculators import MACECalculator
            except ImportError as exc:
                raise ImportError("mace not installed. See https://github.com/ACEsuit/mace") from exc

            class MyMACECalculator(_MyCalculator, MACECalculator):
                """Add abi_forces and abi_stress"""

            self.model_path = os.path.expanduser("~/NN_MODELS/2023-08-14-mace-universal.model")
            print("Using MACE model_path:", self.model_path)

            if self.model_path is None:
                raise RuntimeError("MACECalculator requires model_path e.g. nn_name='mace@FILEPATH'")

            cls = MyMACECalculator if with_delta else MACECalculator
            calc = cls(model_paths=self.model_path, device="cpu") #, default_dtype='float32')

        elif self.nn_type == "nequip":
            try:
                from nequip.ase.nequip_calculator import NequIPCalculator
            except ImportError as exc:
                raise ImportError("nequip not installed. See https://github.com/mir-group/nequip") from exc

            class MyNequIPCalculator(_MyCalculator, NequIPCalculator):
                """Add abi_forces and abi_stress"""

            if self.model_path is None:
                raise RuntimeError("NequIPCalculator requires model_path e.g. nn_name='nequip:FILEPATH'")

            cls = MyNequIPCalculator if with_delta else NequIPCalculator
            calc = cls.from_deployed_model(modle_path=self.model_path, species_to_type_name=None)

        elif self.nn_type == "metatensor":
            try:
                from metatensor.torch.atomistic.ase_calculator import MetatensorCalculator
            except ImportError as exc:
                raise ImportError("metatensor not installed. See https://github.com/lab-cosmo/metatensor") from exc

            class MyMetatensorCalculator(_MyCalculator, MetatensorCalculator):
                """Add abi_forces and abi_stress"""

            if self.model_path is None:
                raise RuntimeError("MetaTensorCalculator requires model_path e.g. nn_name='metatensor:FILEPATH'")

            cls = MyMetaTensorCalculator if with_delta else MetatensorCalculator
            calc = cls(self.model_path)

        elif self.nn_type == "deepmd":
            try:
                from deepmd.calculator import DP
            except ImportError as exc:
                raise ImportError("deepmd not installed. See https://tutorials.deepmodeling.com/") from exc

            class MyDpCalculator(_MyCalculator, DP):
                """Add abi_forces and abi_stress"""

            if self.model_path is None:
                raise RuntimeError("DeepMD calculator requires model_path e.g. nn_name='deepmd:FILEPATH'")

            cls = MyDp if with_delta else Dp
            calc = cls(self.model_path)

        else:
            raise ValueError(f"Invalid {self.nn_type=}")

        # Include DFTD3 vDW corrections on top of ML potential.
        if self.dftd3_args is not None:
            from ase.calculators.dftd3 import DFTD3
            calc = DFTD3(dft=calc, **self.dftd3_args)

        return calc


class MlBase(HasPickleIO):
    """
    Base class for all Ml subclasses providing helper methods to
    perform typical tasks such as writing files in the workdir
    and object persistence via pickle.
    """

    def __init__(self, workdir, prefix=None, exist_ok=False):
        """
        Build directory with `prefix` if `workdir` is None else create it.
        If exist_ok is False (the default), a FileExistsError is raised if the target directory already exists.
        """
        self.workdir = workdir_with_prefix(workdir, prefix, exist_ok=exist_ok)
        self.basename_info = []

    def __str__(self):
        # Delegated to the subclass.
        return self.to_string()

    def add_basename_info(self, basename: str, info: str) -> None:
        """
        Register basename with info in the internal buffer used to generate
        the README.md file in _finalize. Print WARNING if basename is already registered.
        """
        if any(basename == t[0] for t in self.basename_info):
            print(f"WARNING: {basename=} already in basename_info!")
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
            json.dump(data, fh, indent=indent, cls=MontyEncoder, **kwargs)

        if stream is not None:
            # Print JSON to stream as well.
            print("", file=stream)
            print(marquee(info, mark="="), file=stream)
            print(json.dumps(data, cls=MontyEncoder, indent=4), file=stream, end="\n")

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
        """
        Called at the end of the `run` method to write the README.md file in the workdir.
        """
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

        print("\nResults available in:", self.workdir)


class MlRelaxer(MlBase):
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

        doc = yaml_safe_load_path(filepath) #; print(doc)

        format_version = doc.pop("format_version")
        natom = doc.pop("natom")
        ntypat = doc.pop("ntypat")
        typat = np.array(doc.pop("typat"), dtype=int)
        znucl = np.array(doc.pop("znucl"), dtype=float)
        rprim = np.array(doc.pop("lattice_vectors"))
        xred = np.array(doc.pop("xred"))
        xred = np.array(xred[:,:3], dtype=float)
        # Read forces and stress in a.u. and convert.
        abi_cart_forces = np.array(doc.pop("cartesian_forces")) * abu.Ha_eV / abu.Bohr_Ang
        abi_cart_stresses = np.array(doc.pop("cartesian_stress_tensor")) * abu.Ha_eV / (abu.Bohr_Ang**3)

        ionmov = doc.pop("ionmov")
        optcell = doc.pop("optcell")
        iatfix = doc.pop("iatfix") # [3,natom] array or None if unconstrained.
        strtarget = np.array(doc.pop("strtarget"), dtype=float)
        nn_name = doc.get("nn_name", "chgnet")
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

        # Set internal parameters according to YAML file and build object.
        fmax, steps, optimizer = 0.01, 500, "BFGS"

        new = cls(atoms, relax_mode, fmax, pressure, steps, optimizer, nn_name, verbose,
                  workdir=workdir, prefix=prefix)

        # Set delta forces and delta stress if the script is called by ABINIT.

        return new

    def write_output_file_for_abinit(self) -> Path:
        """
        Write output file with results in a format that can be parsed by ABINIT.
        Return path to the output file.
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
            optimizer: name of the ASE optimizer to use for relaxation.
            nn_name: String defining the NN potential. See also CalcBuilder.
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
        workdir = self.workdir
        self.atoms.calc = CalcBuilder(self.nn_name).get_calculator()
        # TODO: Here I should add the ab-initio forces/stress to the calculator to correct the ML ones

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


def restart_md(traj_filepath, atoms, verbose) -> tuple[bool, int]:
    """
    Try to restart a MD run from an existent trajectory file.
    Return: (restart_bool, len_traj)
    """
    traj_filepath = str(traj_filepath)
    if not os.path.exists(traj_filepath):
        if verbose: print(f"Starting MD run from scratch.")
        return False, 0

    print(f"Restarting MD run from the last image of the trajectory file: {traj_filepath}")
    traj = read(traj_filepath, index=":")
    last = traj[-1]
    if verbose:
        print("input_cell:\n", atoms.cell)
        print("last_cell:\n", last.cell)
        print("input_positions:\n", atoms.positions)
        print("last_positions:\n", last.positions)
        print("input_velocites:\n", atoms.get_velocities())
        print("last_velocites:\n", last.get_velocities())

    atoms.set_cell(last.cell, scale_atoms=False, apply_constraint=True)
    atoms.set_positions(last.positions)
    atoms.set_velocities(last.get_velocities())
    #atoms.set_initial_magnetic_moments(magmoms=None)
    return True, len(traj)


class AseMdLog(TextFile):
    """
    Postprocessing tool for the log file produced by ASE MD.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AseMdLog
    """

    time_key = "Time[ps]"

    @lazy_property
    def df(self) -> pd.DataFrame:
        """
        DataFrame with the results.
        """
        # Here we implement the logic required to take into account a possible restart.
        begin_restart = False
        d = {}
        add_time = 0.0
        with open(self.filepath, mode="rt") as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    # Extract column names from the header and init dict.
                    columns = line.split()
                    d = {c: [] for c in columns}
                    continue

                if line.startswith(self.time_key):
                    # Found new header denoting restart.
                    begin_restart = True
                    continue

                if begin_restart:
                    # First iteration after restart. Set add_time and kkip it.
                    begin_restart = False
                    add_time = d[self.time_key][-1]
                    continue

                tokens = [float(tok) for tok in line.split()]
                tokens[0] += add_time
                for c, v in zip(columns, tokens):
                    d[c].append(v)

        return pd.DataFrame(d)

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose."""
        return "=== Summary statistics ===\n" + self.df.describe().to_string()

    @add_fig_kwargs
    def plot(self, **kwargs) -> Figure:
        """
        """
        ynames = [k for k in self.df.keys() if k != self.time_key]
        axes = self.df.plot.line(x=self.time_key, y=ynames, subplots=True)
        return axes[0].get_figure()

    @add_fig_kwargs
    def histplot(self, **kwargs) -> Figure:
        """
        """
        ynames = [k for k in self.df.keys() if k != self.time_key]
        axes = self.df.plot.hist(column=ynames, subplots=True)
        return axes[0].get_figure()

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        Generate a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)
        yield self.histplot(show=False)



class MlMd(MlBase):
    """
    Perform MD calculations with ASE and ML potential.
    """

    def __init__(self, atoms: Atoms, temperature, timestep, steps, loginterval,
                 ensemble, nn_name, verbose, workdir, prefix=None):
        """
        Args:
            atoms: ASE atoms.
            temperature: Temperature in K
            timestep:
            steps: Number of steps.
            loginterval:
            ensemble:
            nn_name: String defining the NN potential. See also CalcBuilder.
            verbose: Verbosity level.
            workdir: Working directory.
            prefix:
        """
        super().__init__(workdir, prefix, exist_ok=True)
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
    nn_name     = {self.nn_name}
    workdir     = {self.workdir}
    verbose     = {self.verbose}

=== ATOMS ===

{self.atoms}

"""

    def run(self) -> None:
        """Run MD"""
        workdir = self.workdir
        traj_file = self.get_path("md.traj", "ASE MD trajectory")
        logfile = self.get_path("md.aselog", "ASE MD log file")

        # Write JSON files with parameters.
        md_dict = dict(
            temperature=self.temperature,
            timestep   =self.timestep,
            steps      =self.steps,
            loginterval=self.loginterval,
            ensemble   =self.ensemble,
            nn_name    =self.nn_name,
            workdir    =str(self.workdir),
            verbose    =self.verbose,
        )
        self.write_json("md.json", md_dict, info="JSON file with ASE MD parameters")

        append_trajectory, len_traj = restart_md(traj_file, self.atoms, self.verbose)
        self.atoms.calc = CalcBuilder(self.nn_name).get_calculator()
        #forces = self.atoms.get_forces()

        if not append_trajectory:
            print("Setting momenta corresponding to the input temperature using MaxwellBoltzmannDistribution.")
            MaxwellBoltzmannDistribution(self.atoms, temperature_K=self.temperature)
        else:
            prev_steps = self.steps
            self.steps = self.steps - (len_traj * self.loginterval)
            if self.steps <= 0:
                raise RuntimeError(f"Have already performed {prev_steps} iterations!")

        md = MolecularDynamics(
            atoms=self.atoms,
            ensemble=self.ensemble,
            temperature=self.temperature,   # K
            timestep=self.timestep,         # fs,
            #pressure,
            trajectory=str(traj_file),      # save trajectory to md.traj
            logfile=str(logfile),           # log file for MD
            loginterval=self.loginterval,   # interval for record the log
            append_trajectory=append_trajectory, # If True, the new structures are appended to the trajectory
        )

        self.write_script("ase_diffusion_coeff.py", text=f"""\
from ase.md.analysis import DiffusionCoefficient
from ase.io import read

# For an MD simulation with timestep of N, and images written every M iterations, our timestep here is N * M.
timestep = {self.timestep} * {self.loginterval}
traj = read("{str(traj_file)}", index=":")
dc = DiffusionCoefficient(traj, timestep, atom_indices=None, molecule=False)
dc.calculate(ignore_n_images=0, number_of_segments=1)
dc.print_data()
dc.plot(ax=None, show=True)
""", info="Python script to compute and visualize diffusion coefficients with ASE.")

        self.write_script("plot_energies.py", text=f"""\
from abipy.ml.aseml import AseMdLog
log = AseMdLog("{str(logfile)}")
log.plot(savefig=None)
""", info="Python script to visualize energies vs time.")

        md.run(steps=self.steps)


class _MlNebBase(MlBase):
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
                        nn_name=str(self.nn_name),
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
            nn_name: String defining the NN potential. See also CalcBuilder.
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
        workdir = self.workdir

        results = []
        calc = CalcBuilder(self.nn_name).get_calculator()
        for ind, atoms in enumerate(self.atoms_list):
            write_vasp(self.workdir / f"IND_{ind}_POSCAR", atoms, label=None)
            atoms.calc = calc
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
            initial_atoms: initial ASE atoms.
            final_atoms: final ASE atoms.
            nimages: Number of images.
            neb_method:
            climb:
            optimizer: name of the ASE optimizer to use for relaxation.
            relax_mode: "ions" to relax ions only, "cell" for ions + cell, "no" for no relaxation.
            fmax: tolerance for relaxation convergence. Here fmax is a sum of force and stress forces.
            pressure: Target pressure
            nn_name: String defining the NN potential. See also CalcBuilder.
            verbose: Verbosity level.
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
        workdir = self.workdir
        initial_atoms, final_atoms = self.initial_atoms, self.final_atoms
        calculator = CalcBuilder(self.nn_name).get_calculator()

        if self.relax_mode != RX_MODE.no:
            relax_kws = dict(calculator=calculator,
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
        calculators = [CalcBuilder(self.nn_name).get_calculator() for i in range(self.nimages)]
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

ase gui "{nebtraj_file}@-{self.nimages}"

# then select `tools->neb` in the gui.
""", info="Shell script to visualize NEB results with ase gui")

        self.write_script("ase_nebplot.sh", text=f"""\
# This command create a series of plots showing the progression of the neb relaxation

ase nebplot --share-x --share-y --nimages {self.nimages} "{nebtraj_file}"
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
            atoms_list: List of ASE atoms.
            nimages: Number of NEB images.
            neb_method:
            climb:
            optimizer:
            relax_mode: "ions" to relax ions only, "cell" for ions + cell, "no" for no relaxation.
            fmax: tolerance for relaxation convergence. Here fmax is a sum of force and stress forces.
            pressure:
            nn_name: String defining the NN potential. See also CalcBuilder.
            verbose: Verbosity level.
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
            linear provides a standard straight-line interpolation, while idpp uses an image-dependent pair potential.
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


class MlOrderer(MlBase):
    """
    Order a disordered structure using pymatgen and ML potential.
    """
    def __init__(self, structure, max_ns, optimizer, relax_mode, fmax, pressure, steps, nn_name, verbose,
                 workdir, prefix=None):
        """
        Args:
            structure: Abipy Structure object or any object that can be converted to structure.
            max_ns:
            optimizer:
            relax_mode:
            fmax:
            pressure:
            steps:
            nn_name: String defining the NN potential. See also CalcBuilder.
            verbose: Verbosity level.
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

        calculator = CalcBuilder(self.nn_name).get_calculator()

        if self.relax_mode != RX_MODE.no:
            print(f"Relaxing structures with relax mode: {self.relax_mode}")
            relax_kws = dict(calculator=calculator,
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


class MlValidateWithAbinitio(_MlNebBase):
    """
    Compare ab-initio energies, forces and stresses with ML results.
    """

    def __init__(self, filepaths, nn_names, traj_range, verbose, workdir, prefix=None):
        """
        Args:
            filepaths: List of file produced by the ab-initio code with energies, forces and stresses.
            nn_names: String or list of strings defining the NN potential. See also CalcBuilder.
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
            from abipy.ml.tools import get_energy_step
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
        workdir = self.workdir

        labels = ["abinitio"]
        abi_results = self.get_abinitio_results()
        results_list = [abi_results]

        ntasks = len(abi_results)
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
        comp.pickle_dump_and_write_script(self.workdir)

        if print_dataframes:
            # Write dataframes to disk in CSV format.
            forces_df = comp.get_forces_dataframe()
            self.write_df(forces_df, "cart_forces.csv", info="CSV file with cartesian forces.")
            stress_df = comp.get_stress_dataframe()
            self.write_df(stress_df, "voigt_stress.csv", info="CSV file with cartesian stresses in Voigt notation")

        self._finalize()
        return comp


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
            if append_trajectory:
                raise NotImplementedError("append_trajectory with Inhomogeneous_NPTBerendsen")

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
        Thin wrapper of ase MD run.

        Args:
            steps (int): number of MD steps
        """
        from ase.md import MDLogger
        self.dyn.attach(MDLogger(self.dyn, self.atoms, '-', header=True, stress=False,
                        peratom=True, mode="a"), interval=self.loginterval)
        self.dyn.run(steps)


class GsMl(MlBase):
    """
    Single point calculation of energy, forces and stress with ML potential.
    """
    def __init__(self, atoms, nn_name, verbose, workdir, prefix=None):
        super().__init__(workdir, prefix)
        self.atoms = atoms
        self.nn_name = nn_name
        self.verbose = verbose

    def run(self):
        calc = CalcBuilder(self.nn_name).get_calculator()
        self.atoms.calc = calc
        res = AseResults.from_atoms(self.atoms)
        print(res.to_string(verbose=self.verbose))

        # Write ASE trajectory file with results.
        with open(self.workdir / "gs.traj", "wb") as fd:
            write_traj(fd, [self.atoms])

        return 0


class MlCompareNNs(MlBase):
    """
    Compare energies, forces and stresses obtained with different ML potentials.
    Also profile the time required.
    """
    def __init__(self, atoms, nn_names, num_tests, rattle, stdev_rvol, verbose, workdir, prefix=None):
        """
        Args:
            atoms: ASE atoms
            nn_names: String or list of strings defining the NN potential. See also CalcBuilder.
            num_tests: Number of tests.
            rattle: Displace atoms randomly with this stdev.
            stdev_rvol: Scale volumes randomly around input v0 with stdev: v0*value
            verbose: Verbosity level.
            workdir: Working directory.
        """
        super().__init__(workdir, prefix)
        self.initial_atoms = atoms
        self.nn_names = list_strings(nn_names)
        self.num_tests = num_tests
        self.rattle = rattle
        self.stdev_rvol = stdev_rvol
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     nn_names    = {self.nn_names}
     num_tests   = {self.num_tests}
     rattle      = {self.rattle}
     stdev_rvol  = {self.stdev_rvol}
     workdir     = {self.workdir}
     verbose     = {self.verbose}

=== ATOMS ===

{self.initial_atoms}

"""

    def run(self, print_dataframes=True) -> AseResultsComparator:
        """
        Run calculations.
        """
        workdir = self.workdir

        # Build list of atoms
        v0 = self.initial_atoms.cell.volume
        vols = np.random.normal(loc=v0, scale=self.stdev_rvol * v0, size=self.num_tests)
        atoms_list = []
        for it in range(self.num_tests):
            atoms = self.initial_atoms.copy()
            atoms.rattle(stdev=abs(self.rattle))
            atoms_list.append(Structure.as_structure(atoms).scale_lattice(vols[it]).to_ase_atoms())

        labels, results_list = [], []
        for nn_name in self.nn_names:
            labels.append(nn_name)
            calc = as_calculator(nn_name)
            with Timer(f"Computing GS properties for {len(atoms_list)} configurations with {nn_name=}...") as timer:
                items = [AseResults.from_atoms(atoms, calc=calc) for atoms in atoms_list]

            results_list.append(items)

        comp = AseResultsComparator.from_ase_results(labels, results_list)
        comp.pickle_dump_and_write_script(self.workdir)

        if print_dataframes:
            # Write dataframes to disk in CSV format.
            forces_df = comp.get_forces_dataframe()
            self.write_df(forces_df, "cart_forces.csv", info="CSV file with cartesian forces.")
            stress_df = comp.get_stress_dataframe()
            self.write_df(stress_df, "voigt_stress.csv", info="CSV file with cartesian stresses in Voigt notation")

        self._finalize()
        return comp



class MlCwfEos(MlBase):
    """
    Compute EOS with different ML potentials.
    """
    def __init__(self, elements, nn_names, verbose, workdir, prefix=None):
        """
        Args:
            atoms: ASE atoms
            elements: String or List of strings with the chemical symbols to consider.
            nn_names: String or list of strings defining the NN potential. See also CalcBuilder.
            verbose: Verbosity level.
            workdir: Working directory.
        """
        super().__init__(workdir, prefix)
        self.elements = list_strings(elements)
        self.nn_names = list_strings(nn_names)
        self.verbose = verbose

        self.configurations_set_name = {
            "unaries-verification-PBE-v1": ["X/SC", "X/FCC", "X/BCC", "X/Diamond"],
            "oxides-verification-PBE-v1": ["XO", "XO2", "XO3", "X2O", "X2O3", "X2O5"],
        }

        root = Path("/Users/giantomassi/git_repos/acwf-verification-scripts/0-preliminary-do-not-run")
        self.dirpath_set_name = {
            "unaries-verification-PBE-v1": root / "unaries" / "xsfs-unaries-verification-PBE-v1",
            "oxides-verification-PBE-v1": root / "oxides" / "xsfs-oxides-verification-PBE-v1",
        }

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     elements    = {self.elements}
     nn_names    = {self.nn_names}
     workdir     = {self.workdir}
     verbose     = {self.verbose}
"""

    def run(self):
        for nn_name in self.nn_names:
            for set_name in self.configurations_set_name:
                self.run_nn_name_set_name(nn_name, set_name)

    def run_nn_name_set_name(self, nn_name: str, set_name: str) -> dict:
        print(f"Computing ML EOS with {nn_name=}, {set_name=} ...")
        # This piece is taken from get_results.py
        try:
            from acwf_paper_plots.eosfit_31_adapted import BM, echarge
        except ImportError as exc:
            raise ImportError("ase not installed. Try `pip install ase`.") from exc

        calc = as_calculator(nn_name)
        warning_lines = []

        uuid_mapping = {}
        all_missing_outputs = {}
        completely_off = []
        failed_wfs = []
        all_eos_data = {}
        all_stress_data = {}
        all_BM_fit_data = {}
        num_atoms_in_sim_cell = {}

        import uuid
        for element in self.elements:
            for configuration in self.configurations_set_name[set_name]:
                my_uuid = uuid.uuid4().hex
                uuid_mapping[f'{element}-{configuration}'] = {
                    'structure': my_uuid,
                    'eos_workflow': my_uuid
                }

                #count = 7
                #increment = 0.02
                #return tuple(float(1 + i * increment - (count - 1) * increment / 2) for i in range(count))

                conf = configuration.replace("X/", "")
                filepath = self.dirpath_set_name[set_name] / f"{element}-{conf}.xsf"
                v0_atoms = read(filepath, format='xsf')
                v0 = v0_atoms.get_volume()
                volumes = (v0 * np.array([0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06])).tolist()
                energies, stresses = [], []

                # Initialize to None if the outputs are not there
                eos_data = None
                stress_data = None
                BM_fit_data = None
                num_atoms = len(v0_atoms)

                try:
                    for volume in volumes:
                        #ase = structure.get_ase().copy()
                        #ase.set_cell(ase.get_cell() * float(scale_factor)**(1 / 3), scale_atoms=True)
                        atoms = to_ase_atoms(Structure.as_structure(v0_atoms).scale_lattice(volume))
                        r = AseResults.from_atoms(atoms, calc=calc)
                        energies.append(r.ene)
                        stresses.append(r.stress.tolist())

                    eos_data = list(zip(volumes, energies))
                    stress_data = list(zip(volumes, stresses))
                    # This line disables the visualization of stress
                    stress_data = None
                    #for v, e in eos_data: print(v, e)
                    #except IndexError:

                    # Check if the central point was completely off (i.e. the minimum of the energies is
                    # on the very left or very right of the volume range)
                    min_loc = np.array(energies).argmin()
                    if min_loc == 0:
                        # Side is whether the minimum occurs on the left side (small volumes) or right side (large volumes)
                        completely_off.append({'element': element, 'configuration': configuration, 'side': 'left'})
                    elif min_loc == len(energies) - 1:
                        completely_off.append({'element': element, 'configuration': configuration, 'side': 'right'})

                    try:
                        min_volume, E0, bulk_modulus_internal, bulk_deriv, residuals = BM(np.array(eos_data))
                        bulk_modulus_GPa = bulk_modulus_internal * echarge * 1.0e21
                        #1 eV/Angstrom3 = 160.21766208 GPa
                        bulk_modulus_ev_ang3 = bulk_modulus_GPa / 160.21766208
                        BM_fit_data = {
                            'min_volume': min_volume,
                            'E0': E0,
                            'bulk_modulus_ev_ang3': bulk_modulus_ev_ang3,
                            'bulk_deriv': bulk_deriv,
                            'residuals': residuals[0]
                        }
                        if residuals[0] > 1.e-3:
                            warning_lines.append(f"WARNING! High fit residuals: {residuals[0]} for {element} {configuration}")
                    except ValueError as exc:
                        # If we cannot find a minimum
                        # Note that BM_fit_data was already set to None at the top
                        warning_lines.append(f"WARNING! Unable to fit for {element=} {configuration=}")
                        #print(str(exc))

                except Exception as exc:
                    warning_lines.append(f"WARNING! Unable to compute E(V) for {element=} {configuration=}")
                    #print(str(exc))

                all_eos_data[f'{element}-{configuration}'] = eos_data
                num_atoms_in_sim_cell[f'{element}-{configuration}'] = num_atoms
                if stress_data is None:
                    stress_data = list(zip(volumes, [None for _ in range(len(volumes))]))
                all_stress_data[f'{element}-{configuration}'] = stress_data
                all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data

        data = {
            'script_version': "0.0.4",
            'set_name': set_name,
            # Mapping from strings like "He-X2O" to a dictionary with the UUIDs of the structure and the EOS workflow
            'uuid_mapping': uuid_mapping,
            # A list of dictionaries with information on the workchains that did not finish with a 0 exit code
            'failed_wfs': failed_wfs,
            # A dictionary that indicate for which elements and configurations there are missing outputs,
            # (only for the workchains that still had enough volumes to be considered for a fit)
            'missing_outputs': all_missing_outputs,
            # A list of dictionaries that indicate which elements and configurations have been computed completely
            # off-centre (meaning that the minimum of all computed energies is on either of the two edges, i.e. for
            # the smallest or largest volume)
            'completely_off': completely_off,
            # Dictionary with the EOS data (volumes and energies datapoints). The keys are the same as the `uuid_mapping`.
            # Values can be None.
            'eos_data': all_eos_data,
            'stress_data': all_stress_data,
            # Birch-Murnaghan fit data. See above for the keys. Can be None.
            'BM_fit_data': all_BM_fit_data,
            'num_atoms_in_sim_cell': num_atoms_in_sim_cell
        }

        # Print some statistics on the results
        warning_lines.append("")
        #warning_lines.append("Counter of states: " + str(Counter(states)))
        good_cnt = len([eos_data for eos_data in data['eos_data'] if eos_data is not None])
        warning_lines.append("")
        warning_lines.append(f"Minimum completely off for {len(completely_off)}/{good_cnt}")
        warning_lines.append("Completely off systems (symbol indicates if the minimum is on the very left or right):")
        for system in data['completely_off']:
            warning_lines.append(
                f"- {system['element']} {system['configuration']} "
                f"({'<' if system['side'] == 'left' else '>'})"
            )

        SET_NAME = set_name
        PLUGIN_NAME = nn_name

        fname = self.workdir / f"warnings-{SET_NAME}-{PLUGIN_NAME}.txt"
        with open(fname, 'w') as fhandle:
            for line in warning_lines:
                fhandle.write(f"{line}\n")
                print(line)
        print(f"Warning log written to: '{fname}'.")

        # Output results to file
        fname = self.workdir / f"results-{SET_NAME}-{PLUGIN_NAME}.json"
        with open(fname, 'w') as fhandle:
            json.dump(data, fhandle, indent=2, sort_keys=True)
        print(f"Output results written to: '{fname}'.")

        return data

"""
Objects to perform ASE calculations with machine-learned potentials.
"""
from __future__ import annotations

import sys
import os
import pickle
import json
import itertools
import numpy as np
import pandas as pd
import abipy.tools.duck as duck

from pathlib import Path
from dataclasses import dataclass
#from typing import Type, Any, Optional, Union
from monty.json import MontyEncoder
from monty.functools import lazy_property
from monty.collections import dict2namedtuple # AttrDict,
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord import pbc_shortest_vectors
from abipy.core import Structure
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt #, get_axarray_fig_plt,
from abipy.tools.iotools import workdir_with_prefix, PythonScript, ShellScript
#from abipy.tools.printing import print_dataframe
from abipy.ml.aseml import relax_atoms, get_atoms, as_calculator, ase_optimizer_cls, RX_MODE, fix_atoms


@dataclass
class Entry:
    """
    An Entry stores the relaxed structure with the associated energy
    and the Cartesian forces.
    """
    isite: int             # Index of the site being relaxed.
    structure: Structure   # pymatgen Structure
    energy: float          # Energy in eV
    forces: np.ndarray     # Forces in eV/Ang

    @classmethod
    def from_atoms_and_calculator(cls, isite, atoms, calculator) -> Entry:
        """
        Build an Entry from an Atoms object by attaching a calculator
        to compute energy and forces at fixed geometry.
        """
        structure = Structure.as_structure(atoms)
        # Keep sites within the unit cell.
        structure.translate_sites(range(len(structure)), np.zeros(3), to_unit_cell=True)
        entry = cls(isite=isite,
                    structure=structure,
                    energy=float(calculator.get_potential_energy(atoms=atoms)),
                    forces=calculator.get_forces(atoms=atoms),
                   )

        return entry

    @classmethod
    def from_traj(cls, isite, traj) -> Entry:
        """
        Build an Entry by taking the last Atoms object in an ASE trajectory.
        """
        atoms = traj[-1]
        structure = Structure.as_structure(atoms)
        # Keep sites within the unit cell.
        structure.translate_sites(range(len(structure)), np.zeros(3), to_unit_cell=True)
        return cls(isite=isite,
                   structure=structure,
                   energy=float(atoms.get_potential_energy()),
                   forces=atoms.get_forces(),
                   )

    def __eq__(self, other: Entry) -> bool:
        """
        Invoked by python to evaluate `self == other`.
        """
        return abs(self.energy - other.energy) / len(self.structure) < 1e-4 and \
               np.abs(self.structure.lattice.matrix - other.structure.lattice.matrix).max() < 1e-3 and \
               all(np.abs(site1.frac_coords - site2.frac_coords).max() < 1e-3
                   for site1, site2 in zip(self.structure, other.structure))


class Entries(list):
    """
    A list of Entry objects.
    """

    @classmethod
    def merge_from_topdir(cls, topdir: Path) -> Entries:
        """
        Merge all entries starting from directory `topdir`.
        """
        topdir = Path(topdir)
        top = str(topdir)
        dirpaths = [name for name in os.listdir(top) if os.path.isdir(os.path.join(top, name))
                    and name.startswith("start:")]

        entries = cls()
        for dpath in dirpaths:
            pickle_path = topdir / Path(dpath) / "data.pickle"
            completed_path = topdir / Path(dpath) / "COMPLETED"
            if not (pickle_path.exists() and completed_path.exists()): continue
            with open(pickle_path, "rb") as fh:
                for e in pickle.load(fh):
                    if e in entries: continue
                    entries.append(e)

        if not entries:
            raise ValueError(f"Entries is empty. Cannot find `data.pickle` files in {topdir=}")

        return entries

    def get_data(self) -> ScannerData:
        """
        Build and return a ScannerData instance containing a pandas dataframe
        with the total energy in eV and the relaxed coords of the atomic site
        that has been relaxed.
        """
        dict_list = []
        for entry in self:
            site = entry.structure[entry.isite]
            x, y, z, fx, fy, fz = (*site.coords, *site.frac_coords)
            fmods = np.array([np.linalg.norm(force) for force in entry.forces])
            dict_list.append(dict(
                energy=entry.energy,
                xcart0=x, xcart1=y, xcart2=z,
                xred0=fx, xred1=fy, xred2=fz,
                fmin=fmods.min(), fmax=fmods.max(), fmean=fmods.mean(),
            ))

        df = pd.DataFrame(dict_list).sort_values(by="energy", ignore_index=True)

        # Add medata
        df.attrs['lattice_matrix'] = self[0].structure.lattice.matrix
        return ScannerData(df)


class ScannerData:

    @classmethod
    def from_json(cls, filepath: str):
        """Reconstruct object from json file."""
        with open(filepath, 'rt') as fh:
            data = json.load(fh)

        df = pd.DataFrame().as_dict(data["df"])
        df.attrs['lattice_matrix'] = np.reshape(data["lattice_matrix"], (3,3))

        return cls(df)

    def __init__(self, df):
        self.df = df

    def __str__(self):
        return self.to_string()

    #def to_string(self, verbose=0) -> str:

    @lazy_property
    def lattice(self):
        return Lattice(self.df.attrs["lattice_matrix"])

    def get_emin_emax(self, fact=0.001):
        """
        Return min and max energy. The range is extended by .1% on each side.
        """
        emin, emax = self.df['energy'].min(), self.df['energy'].max()
        emin -= fact * abs(emin)
        emax += fact * abs(emax)
        return emin, emax

    def search_enediff_dist(self, ene_diff, dist_tol,
                            emin=None, emax=None, verbose=1):
        """
        Find configurations that differ in energy by less than ene_diff in eV
        an relaxed sites that are less than dist_tol Angstrom apart.

        Args:
            ene_diff:
            dist_tol:
            verbose: Verbosity level.

        Return named tuple three arrays: pair_indices with the index of the pair in the df
            dist_list with distance between the sites and ediff_list with the energy difference.
        """
        def difference_matrix(a):
            x = np.reshape(a, (len(a), 1))
            return x - x.transpose()

        energy = self.df["energy"].to_numpy(dtype=float)
        ediff_mat = difference_matrix(energy)
        xreds = self.df[["xred0", "xred1", "xred2"]].to_numpy(dtype=float)

        pair_indices, dist_list, ediff_list = [], [], []
        dist_tol2 = dist_tol ** 2
        inds = np.triu_indices_from(ediff_mat)
        for i, j in zip(inds[0], inds[1]):
            if ediff_mat[i,j] > ene_diff or i == j: continue
            _, d2 = pbc_shortest_vectors(self.lattice, xreds[i], xreds[j], return_d2=True)
            dist = np.sqrt(float(d2))
            if dist > dist_tol: continue
            print(f"{i=}, {j=}, ediff:", ediff_mat[i, j], f"{dist=}", "x_i:", xreds[i], "x_j:", xreds[j])
            pair_indices.append((i, j))
            dist_list.append(dist)
            ediff_list.append(energy[j] - energy[i])

        print(f"Found {len(pair_indices)} configurations that differ less than {ene_diff=} eV and {dist_tol=} Ang")
        if verbose > 2:
            print(f"{pair_indices=}")

        return dict2namedtuple(pair_indices=np.array(pair_indices),
                               dist_list=np.array(dist_list),
                               ediff_list=np.array(ediff_list))

    #def binsearch_enediff_dist(self, ene_diff, dist_tol,
    #                        emin=None, emax=None, estep=1e-3, bins=None, verbose=0):
    #    """
    #    Find configurations that differ in energy by less than ene_diff in eV
    #    an relaxed sites that are less than dist_tol Angstrom apart.

    #    Args:
    #        ene_diff:
    #        dist_tol:
    #        bins:
    #        verbose: Verbosity level.

    #    Return named tuple three arrays: pair_indices with the index of the pair in the df
    #        dist_list with distance between the sites and ediff_list with the energy difference.
    #    """
    #    df = self.df #.copy()

    #    # Bin the 'energy' column with equal-width intervals.
    #    if bins is None:
    #        _emin, _emax = self.get_emin_emax()
    #        emin = emin if emin is not None else _emin
    #        emax = emax if emax is not None else _emax
    #        bins = np.arange(emin, emax, estep)

    #    df['energy_bin'] = pd.cut(df['energy'], bins=bins, retbins=False)
    #    #df['energy_bin'] = pd.cut(df['energy'], bins=10, retbins=False)

    #    print(df["energy_bin"].value_counts())
    #    if verbose:
    #        print(df[["energy", "energy_bin"]])

    #    pair_indices, dist_list, ediff_list = [], [], []
    #    dist_tol2 = dist_tol ** 2
    #    for energy_bin, group in df.groupby("energy_bin"):
    #        if group.empty: continue      # I don't understand why group can be empty
    #        if len(group) == 1: continue  # Need more than one configuration in the group.
    #        #print(f"{energy_bin=}", "group:\n", group)
    #        for index1, row1 in group.iterrows():
    #            fcoords1 = row1[["xred0", "xred1", "xred2"]].to_numpy(dtype=float)
    #            ene1 = float(row1["energy"])
    #            #cart1 = row1[["xcart0", "xcart1", "xcart2"]].to_numpy(dtype=float)
    #            #print(f"{index1=}", row1, cart1)
    #            for index2, row2 in group.iterrows():
    #                if index2 >= index1: continue
    #                fcoords2 = row2[["xred0", "xred1", "xred2"]].to_numpy(dtype=float)
    #                ene2 = float(row2["energy"])
    #                # Compute distance with minimum-image convention
    #                _, d2 = pbc_shortest_vectors(self.lattice, fcoords1, fcoords2, return_d2=True)
    #                if d2 > dist_tol2: continue
    #                pair_indices.append((index1, index2))
    #                dist_list.append(np.sqrt(d2))
    #                ediff_list.append(e2-e1)

    #    print(f"Found {len(pair_indices)} configurations that differ less than {ene_diff=} eV and {dist_tol=} Ang")
    #    if verbose:
    #        print(f"{pair_indices=}")

    #    return dict2namedtuple(pair_indices=np.array(pair_indices),
    #                           dist_list=np.array(dist_list),
    #                           ediff_list=np.array(ediff_list))

    @add_fig_kwargs
    def histplot(self, **kwargs):
        """
        Plot histogram to show distribution of energies.
        """
        ax, fig, plt = get_ax_fig_plt(ax=None)
        import seaborn as sns
        sns.histplot(self.df, x="energy", ax=ax) #, binwidth=1e-3)
        return fig


class RelaxScanner:
    """
    This object employs an ASE calculator to perform many different
    structural relaxations by placing ...
    """

    def __init__(self, structure, isite, nx, ny, nz, nn_name,
                 relax_mode: str = "ions", fmax: float = 1e-3, steps=500,
                 verbose: int = 0, optimizer_name="BFGS", pressure=0.0,
                 workdir=None, prefix=None):
        """
        Args:
            structure: Structure object or any file supported by pymatgen providing a structure.
            isite: Index of the site to be displaced or string with chemical element to be added.
            nx, ny, nz: Size of the mesh
            nn_name: String with the name of the NN potential e.g. "chgnet" or "m3gnet".
            relax_mode: "ions" to relax ions only, "cell" for ions + cell, "no" for no relaxation.
            fmax: total force tolerance for relaxation convergence.
                Here fmax is a sum of force and stress forces. Defaults to 0.1.
            steps: max number of steps for relaxation.
            verbose: whether to print stdout.
            optimizer_name: name of the ASE optimizer class to use
            pressure: Target pressure.
            workdir: Working directory.
        """
        self.initial_structure = Structure.as_structure(structure).copy()

        if str(isite).isdigit():
            self.isite = int(isite)
        else:
            # Assume isite is a chemical element to be added at the origin.
            self.initial_structure.insert(0, isite, [0, 0, 0])
            self.isite = 0

        assert len(self.initial_structure) > self.isite >= 0

        self.nx, self.ny, self.nz = nx, ny, nz
        self.nn_name = nn_name
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.steps = steps
        self.verbose = verbose
        self.optimizer_name = optimizer_name
        self.pressure = pressure

        self.workdir = workdir_with_prefix(workdir, prefix)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     site index  = {self.isite}
     nx          = {self.nx}
     ny          = {self.ny}
     nz          = {self.nz}
     relax_mode  = {self.relax_mode}
     fmax        = {self.fmax}
     steps       = {self.steps}
     optimizer   = {self.optimizer_name}
     pressure    = {self.pressure}
     calculator  = {self.nn_name}
     verbose     = {self.verbose}
     workdir     = {str(self.workdir)}

"""

#=== INITIAL STRUCTURE ===
#
#{self.initial_structure}


    def _make_directory(self, start, stop) -> Path:
        """
        Create sub-directory inside workdir to perform relaxations for a subset of structures
        defined by `start` and `stop`. Return path to directory.
        """
        # Make sure we are not already running similar ranges.
        top = str(self.workdir)
        dirpaths = [name for name in os.listdir(top) if os.path.isdir(os.path.join(top, name))
                    and name.startswith("start:")]
        for dpath in dirpaths:
            # Directory name has pattern: f"start:{start}-stop:{stop}"
            tokens = dpath.split("-")
            this_start = int(tokens[0].split(":")[1])
            this_stop = int(tokens[1].split(":")[1])
            if start >= this_start and start < this_stop:
                raise RuntimeError(f"{start=}, {stop=} range overlaps with {dpath=}")

        # Make sure we are not already running the same range.
        newdir = self.workdir / f"start:{start}-stop:{stop}"
        if newdir.exists():
            raise RuntimeError(f"{newdir=} already exists!")

        newdir.mkdir()
        return newdir

    def run_start_count(self, start: int, count: int) -> Path:
        """
        Invoke ASE to relax `count` structures starting from index `start`.
        Stores unique results in an internal list and save results to file in pickle format.
        Return Path to the directory with the results.

        Args:
            start: Index of initial structure.
            count: Number of structures to relax. If < 0, all structures are considered.
        """
        if count < 0:
            count = self.nx * self.ny * self.nz

        stop = start + count
        directory = self._make_directory(start, stop)
        calculator = as_calculator(self.nn_name)

        entries, errors = Entries(), []
        for cnt, (ix, iy, iz) in enumerate(itertools.product(range(self.nx), range(self.ny), range(self.nz))):
            if not (stop > cnt >= start): continue

            structure = self.initial_structure.copy()
            structure._sites[self.isite].frac_coords = np.array((ix/self.nx, iy/self.ny, iz/self.nz))
            atoms = get_atoms(structure)

            if self.relax_mode == RX_MODE.no:
                # Just GS energy, no relaxation.
                entry = Entry.from_atoms_and_calculator(self.isite, atoms, calculator)

            else:
                try:
                    # Relax atoms with constraints.
                    fix_atoms(atoms, fix_inds=[i for i in range(len(atoms)) if i != self.isite])
                    relax = relax_atoms(atoms,
                                        relax_mode=self.relax_mode,
                                        optimizer=self.optimizer_name,
                                        fmax=self.fmax,
                                        pressure=self.pressure,
                                        verbose=self.verbose,
                                        steps=self.steps,
                                        traj_path=directory / "relax.traj",
                                        calculator=calculator,
                    )

                    entry = Entry.from_traj(self.isite, relax.traj)

                except Exception:
                    entry = None
                    errors.append((ix, iy, iz))

            if entry is not None and entry not in entries:
                entries.append(entry)

        if errors:
            print(f"ASE relaxation raised {len(errors)} exception(s)!")

        # Save results in pickle format.
        with open(directory / "data.pickle", "wb") as fh:
            pickle.dump(entries, fh)

        with open(directory / "COMPLETED", "wt") as fh:
            fh.write("completed")

        return directory

    def run(self, nprocs) -> None:
        """
        Run relaxations using a multiprocessing Pool with nprocs processors.
        If nprocs is None or negative, half the number of CPUs in the system is used.
        """
        # Write python scrip to analyze the results.
        with PythonScript(self.workdir / "analyze.py") as script:
            script.add_text("""
from abipy.ml.relax_scanner import RelaxScanner, Entries
from abipy.tools.printing import print_dataframe

data = Entries.merge_from_topdir(".").get_data()
print_dataframe(data.df)
#data.histplot()
""")

        if nprocs is None or nprocs < 0:
            nprocs = max(1, os.cpu_count() // 2)

        if nprocs == 1:
            # Sequential execution.
            self.run_start_count(0, -1)
        else:
            # Split ntasks across nprocs.
            args_list = []
            ntasks = self.nx * self.ny * self.nz
            div, rest = divmod(ntasks, nprocs)
            for iproc in range(nprocs):
                start, stop = iproc * div, (iproc + 1) * div
                count = div
                if iproc == nprocs - 1:
                    stop += rest
                    count = div + rest
                args_list.append((self, start, count))

            from multiprocessing import Pool
            print(f"Using multiprocessing pool with {nprocs=} ...")
            with Pool(processes=nprocs) as pool:
                pool.map(_func, args_list)

        entries = Entries.merge_from_topdir(self.workdir)
        data = entries.get_data()
        print(data.df)
        data.df.to_json()

        data = dict(df=df.to_dict(),
                    lattice_matrix=df.attrs["lattice_matrix"].to_list())
        with open(self.workdir / "ScannerData.json", 'wt') as fh:
            json.dump(data, fh)


def _func(args):
    """Function passed to pool.map."""
    self, start, count = args
    self.run_start_count(start, count)


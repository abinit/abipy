"""
Objects to perform ASE calculations with machine-learned potentials.
"""
from __future__ import annotations

import sys
import os
import pickle
import json
import itertools
import warnings
import dataclasses
import numpy as np
import pandas as pd
#import abipy.tools.duck as duck

from pathlib import Path
from multiprocessing import Pool
#from typing import Type, Any, Optional, Union
from monty.json import MontyEncoder
from monty.functools import lazy_property
from monty.collections import dict2namedtuple
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord import pbc_shortest_vectors

from abipy.core import Structure
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt #, get_axarray_fig_plt,
from abipy.tools.iotools import workdir_with_prefix, PythonScript, ShellScript
from abipy.tools.printing import print_dataframe
from abipy.ml.aseml import (relax_atoms, get_atoms, as_calculator, ase_optimizer_cls, RX_MODE, fix_atoms,
                            MlNeb, MlGsList, CalcBuilder, make_ase_neb)


@dataclasses.dataclass
class Entry:
    """
    Stores the relaxed structure with the associated energy and Cartesian forces.
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
        # NB: Keep sites within the unit cell so that we can compare entries in __eq__
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
        Build an Entry by taking the last Atoms object from an ASE trajectory.
        """
        atoms = traj[-1]
        structure = Structure.as_structure(atoms)
        # NB: Keep sites within the unit cell so that we can compare entries in __eq__
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


@dataclasses.dataclass
class Pair:
    """
    """
    index1: int
    index2: int
    ediff: float  # Energy difference in eV
    dist:  float
    frac_coords1: np.ndarray
    frac_coords2: np.ndarray

    def get_dict4pandas(self) -> dict:
        """Useful to construct pandas DataFrames"""
        return dataclasses.asdict(self)


class RelaxScannerAnalyzer:
    """
    """

    @classmethod
    def from_topdir(cls, topdir: Path) -> RelaxScannerAnalyzer:
        """
        Merge all entries starting from directory `topdir`.
        """
        topdir = Path(topdir)
        top = str(topdir)
        dirpaths = [name for name in os.listdir(top) if os.path.isdir(os.path.join(top, name))
                    and name.startswith("start:")]

        entries = []
        for dpath in dirpaths:
            pickle_path = topdir / Path(dpath) / "entries.pickle"
            completed_path = topdir / Path(dpath) / "COMPLETED"
            if not (pickle_path.exists() and completed_path.exists()): continue
            with open(pickle_path, "rb") as fh:
                for e in pickle.load(fh):
                    if e in entries: continue
                    entries.append(e)

        if not entries:
            raise ValueError(f"Entries is empty. Cannot find `entries.pickle` files in {topdir=}")

        # Reconstruct RelaxScanner instance from pickle file.
        scanner = RelaxScanner.pickle_load(topdir)
        return cls(entries, scanner)

    def __init__(self, entries: list[Entry], scanner: RelaxScanner, verbose: int = 0):
        self.entries = entries
        self.scanner = scanner
        self.verbose = verbose

    @property
    def workdir(self):
        return self.scanner.workdir

    #@property
    #def initial_structure(self):
    #    return self.scanner.initial_structure

    @lazy_property
    def df(self) -> pd.DataFrame:
        """
        pandas dataframe with the total energy in eV and the relaxed coords
        of the atomic site that has been relaxed.
        """
        dict_list = []
        entries = self.entries
        for entry in entries:
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
        df.attrs['lattice_matrix'] = entries[0].structure.lattice.matrix
        df.attrs['structure'] = entries[0].structure
        return df

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        s = self.scanner.to_string(verbose=verbose)
        return s

    @lazy_property
    def lattice(self):
        return Lattice(self.df.attrs["lattice_matrix"])

    #def get_emin_emax(self, fact=0.001) -> tuple[float, float]:
    #    """
    #    Return min and max energy. The range is extended by .1% on each side.
    #    """
    #    emin, emax = self.df['energy'].min(), self.df['energy'].max()
    #    emin -= fact * abs(emin)
    #    emax += fact * abs(emax)
    #    return emin, emax

    def pairs_enediff_dist(self, ene_diff=1e-3, dist_tol=3.0) -> list[Pair]:
        """
        Find pairs that differ in energy less than ene_diff
        and with relaxed sites that are less than dist_tol Angstrom apart
        (minimum-image convention is applied).

        Args:
            ene_diff: Energy difference in eV.
            dist_tol: Tolerance on site distance in Ang

        Return list of Pair objects.
        """
        def adiff_matrix(vec):
            """Return matrix A_ij with all the differences |vec_i - vec_j|."""
            x = np.reshape(vec, (len(vec), 1))
            return np.abs(x - x.transpose())

        energy = self.df["energy"].to_numpy(dtype=float)
        aediff_mat = adiff_matrix(energy)
        xreds = self.df[["xred0", "xred1", "xred2"]].to_numpy(dtype=float)

        pairs = []
        inds = np.triu_indices_from(aediff_mat)
        for i, j in zip(inds[0], inds[1]):
            if aediff_mat[i,j] > ene_diff or i == j: continue
            _, d2 = pbc_shortest_vectors(self.lattice, xreds[i], xreds[j], return_d2=True)
            dist = np.sqrt(float(d2))
            if dist > dist_tol: continue
            pair = Pair(index1=i, index2=j, ediff=energy[j] - energy[i], dist=dist,
                        frac_coords1=xreds[i], frac_coords2=xreds[j])
            pairs.append(pair)
            if self.verbose > 1: print(pair)

        print(f"Found {len(pairs)} pairs with {ene_diff=} eV and {dist_tol=} Ang.")

        # Here we compute the transition energy for each pair either
        # by performing NEB or single point calculations along a linear path connecting the two sites.
        nprocs = -1
        nprocs = _get_nprocs(nprocs)

        #if nprocs == 1:
        d_list = [self.run_pair(pair, neb_method="None") for pair in pairs]
        #else:
        #    print(f"Using multiprocessing pool with {nprocs=} ...")
        #    kws_list = [dict(self=self, pair=pair, neb_method="None") for pair in pairs]
        #    with Pool(processes=nprocs) as pool:
        #        d_list = pool.map(_map_run_pair, kws_list)

        df = pd.DataFrame(d_list)
        df["ene_diff"], df["dist_tol"] = ene_diff, dist_tol
        print_dataframe(df)

        return pairs

    def run_pair(self, pair: Pair, neb_method="None", nimages=14, climb=False) -> dict:
        """
        Return dictionary with results.

        Args:
            pair: Info on Pair.
            frac_coords1: Fractional coordinates of the first site.
            frac_coords2: Fractional coordinates of the second site.
            neb_method:
            nimages:
            climb:
        """
        # Get atoms with constraints.
        scanner = self.scanner
        initial_atoms = scanner.get_atoms_with_frac_coords(pair.frac_coords1)
        final_atoms = scanner.get_atoms_with_frac_coords(pair.frac_coords2)

        relax_mode = "no"
        calc_builder = CalcBuilder(scanner.nn_name)
        pair_str = f"{pair.index1}-{pair.index2}"

        # TODO: Check mic
        if neb_method == "None":
            calculators = [calc_builder.get_calculator() for i in range(nimages)]
            neb = make_ase_neb(initial_atoms, final_atoms, nimages,
                               calculators, "aseneb", climb,
                               method='linear', mic=False)

            # NB: It's not a NEB ML object although it provides a similar API.
            ml_neb = MlGsList(neb.images, calc_builder, self.verbose,
                              workdir=self.workdir / f"GSLIST/pair:{pair_str}")

        else:
            ml_neb = MlNeb(initial_atoms, final_atoms,
                           nimages, neb_method, climb, scanner.optimizer_name,
                           relax_mode, scanner.fmax, scanner.pressure,
                           calc_builder, self.verbose,
                           workdir=self.workdir / f"NEB/pair:{pair_str}")

        if self.verbose: print(ml_neb.to_string(verbose=self.verbose))

        ml_neb.run()
        neb_data = ml_neb.read_neb_data()

        # Build out_data dict.
        out_data = pair.get_dict4pandas()

        # Add keys from neb_data.
        keys = ["barrier_with_fit", "energy_change_with_fit",
                "barrier_without_fit", "energy_change_without_fit"
               ]
        out_data.update({k: neb_data[k] for k in keys})
        out_data.update(dict(neb_method=neb_method, nimages=nimages, climb=climb))

        return out_data

    @add_fig_kwargs
    def histplot(self, ax=None, **kwargs):
        """
        Plot histogram to show energy distribution.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        import seaborn as sns
        sns.histplot(self.df, x="energy", ax=ax)
        return fig


class RelaxScanner:
    """
    This object employs an ASE calculator to perform many different
    structural relaxations by placing ...
    """

    @classmethod
    def pickle_load(cls, workdir) -> RelaxScanner:
        """
        Reconstruct the object from a pickle file located in workdir.
        """
        with open(Path(workdir) / f"{cls.__name__}.pickle", "rb") as fh:
            return pickle.load(fh)

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

        # Write pickle file for object persistence.
        with open(self.workdir / f"{self.__class__.__name__}.pickle", "wb") as fh:
            pickle.dump(self, fh)

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

=== INITIAL STRUCTURE ===

{self.initial_structure}
"""

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

    def get_atoms_with_frac_coords(self, frac_coords, with_fixex_atoms=True) -> Atoms:
        """
        Return Atoms instance with frac_coords at site index `isite`.
        """
        structure = self.initial_structure.copy()
        structure._sites[self.isite].frac_coords = np.array(frac_coords)
        atoms = get_atoms(structure)
        if with_fixex_atoms:
            fix_atoms(atoms, fix_inds=[i for i in range(len(atoms)) if i != self.isite])
        return atoms

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

        #with warnings.catch_warnings():
        #    warnings.simplefilter("ignore")
        calculator = as_calculator(self.nn_name)

        entries, errors = [], []
        for cnt, (ix, iy, iz) in enumerate(itertools.product(range(self.nx), range(self.ny), range(self.nz))):
            if not (stop > cnt >= start): continue

            atoms = self.get_atoms_with_frac_coords((ix/self.nx, iy/self.ny, iz/self.nz))

            if self.relax_mode == RX_MODE.no:
                # Just GS energy, no relaxation.
                entry = Entry.from_atoms_and_calculator(self.isite, atoms, calculator)

            else:
                try:
                    # Relax atoms with constraints.
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
        with open(directory / "entries.pickle", "wb") as fh:
            pickle.dump(entries, fh)

        with open(directory / "COMPLETED", "wt") as fh:
            fh.write("completed")

        return directory

    def run(self, nprocs) -> None:
        """
        Main entry point for client code.
        Run relaxations using a multiprocessing Pool with nprocs processors.
        If nprocs is None or negative, the total number of CPUs in the system is used.
        """
        # Write python scrip to analyze the results.
        with PythonScript(self.workdir / "analyze.py") as script:
            script.add_text("""
from abipy.ml.relax_scanner import RelaxScannerAnalyzer
rsa = RelaxScannerAnalyzer.from_topdir(".")
print_dataframe(rsa.df)
#rsa.histplot()

rsa.pairs_enediff_dist(ene_diff=1e-3, dist_tol=3.0)
""")

        nprocs = _get_nprocs(nprocs)
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

            print(f"Using multiprocessing pool with {nprocs=} ...")
            with Pool(processes=nprocs) as pool:
                pool.map(_map_run_start_count, args_list)

        entries = Entries.merge_from_topdir(self.workdir)
        data = entries.get_data()
        print(data.df)

        data = dict(df=data.df.to_dict(),
                    lattice_matrix=data.df.attrs["lattice_matrix"].tolist())
        with open(self.workdir / "RelaxScannerAnalyzer.json", 'wt') as fh:
            json.dump(data, fh)


def _get_nprocs(nprocs) -> int:
    """
    Return the number of procs to be used for multiprocessing Pool.
    If negative or None, use all procs in the system
    """
    if nprocs is None or nprocs <= 0:
        return max(1, os.cpu_count())
    return int(nprocs)


def _map_run_start_count(args):
    """Function passed to pool.map."""
    scanner, start, count = args
    return scanner.run_start_count(start, count)


def _map_run_pair(kwargs):
    """Function passed to pool.map."""
    #self = kwargs.pop("self")
    return {}
    #return self.run_pair(**kwargs)

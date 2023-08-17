"""
Objects to perform ASE calculations with machine-learned potentials.
"""
from __future__ import annotations

import os
import pickle
import itertools
import numpy as np
import pandas as pd
import abipy.tools.duck as duck

from pathlib import Path
from dataclasses import dataclass
#from typing import Type, Any, Optional, Union
from monty.json import MontyEncoder
from abipy.core import Structure
#from abipy.tools.plotting import get_ax_fig_plt #, get_axarray_fig_plt,
from abipy.tools.iotools import workdir_with_prefix, PythonScript, ShellScript
#from abipy.tools.printing import print_dataframe
from abipy.ml.aseml import relax_atoms, get_atoms, as_calculator, ase_optimizer_cls, RX_MODE, fix_atoms


@dataclass
class Entry:
    """
    An Entry stores the relaxed structure with the associated energy
    and the Cartesian forces.
    """
    isite: int
    structure: Structure   # pymatgen Structure
    energy: float          # Energy in eV
    forces: np.ndarray     # Forces in eV/Ang

    @classmethod
    def from_atoms_and_calculator(cls, isite, atoms, calculator) -> Entry:
        """
        Build an Entry from an Atoms object by attaching a calculator
        to compute energy and forces with fixed geometry.
        """
        entry = cls(isite=isite,
                    structure=Structure.as_structure(atoms),
                    energy=float(calculator.get_potential_energy(atoms=atoms)),
                    forces=calculator.get_forces(atoms=atoms),
                   )

        return entry

    @classmethod
    def from_traj(cls, isite, traj) -> Entry:
        """
        Build an Entry by reading the last Atoms object stored an ASE trajectory.
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
        """Invoked by python to evaluate `self == other`."""
        return abs(self.energy - other.energy) / len(self.structure) < 1e-4 and \
               np.abs(self.structure.lattice.matrix - other.structure.lattice.matrix).max() < 1e-3 and \
               all(np.abs(site1.frac_coords - site2.frac_coords).max() < 1e-3
                   for site1, site2 in zip(self.structure, other.structure))
               #self.structure == other.structure


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

        return entries

    def get_dataframe(self, with_cart_coords=True, with_frac_coords=True) -> pd.DataFrame:
        """
        Build and return a pandas dataframe with the total energy in eV
        and the relaxed coords of the atomic site that has been displaced.
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
        if not with_cart_coords: df = df.drop(labels=["xcart0", "xcart1", "xcart2"], axis=1)
        if not with_frac_coords: df = df.drop(labels=["xred0", "xred1", "xred2"], axis=1)
        return df


class RelaxScanner:
    """
    This object employs an ASE calculator to perform a structural relaxation.
    for ...
    """

    def __init__(self, structure, isite, nx, ny, nz, calculator,
                 relax_mode: str = "ions", fmax: float = 1e-3, steps=500,
                 verbose: int = 0, optimizer_name="BFGS", pressure=0.0,
                 workdir=None, prefix=None):
        """
        Args:
            structure: Structure object or any file supported by pymatgen providing a structure.
            isite: Index of the site to be displaced or string with chemical element to be added.
            nx, ny, nz: Size of the mesh giving the initial position of isite.
            calculator: ASE calculator or string with the name of the nn potential e.g. "chgnet" or "m3gnet".
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
        self.calculator_name = calculator
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.steps = steps
        self.verbose = verbose
        self.optimizer_name = optimizer_name
        self.pressure = pressure

        #if workdir is not None:
        #    self.workdir = Path(workdir)
        #else:
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
     calculator  = {self.calculator_name}
     verbose     = {self.verbose}
     workdir     = {str(self.workdir)}

"""

#=== INITIAL STRUCTURE ===
#
#{self.initial_structure}



    def _make_directory(self, start, stop) -> Path:
        """
        Create sub-directory in workdir to perform relaxations for a subset of structures
        defined by `start` and `stop`. Return path to directory.
        """
        # Make sure we are not already running similar ranges
        top = str(self.workdir)
        dirpaths = [name for name in os.listdir(top) if os.path.isdir(os.path.join(top, name))
                    and name.startswith("start:")]
        for dpath in dirpaths:
            tokens = dpath.split("-")
            # Directory name has pattern: f"start:{start}-stop:{stop}"
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
        Invoke ASE to relax `ncount` structures starting from index `start`.
        Stores unique results in an internal list, finally save results to file in pickle format.

        Args:
            start: Index of initial structure.
            count: Number of structures to relax. If < 0, all structures are considered.

        Return: Path to directory.
        """
        if count < 0:
            count = self.nx * self.ny * self.nz

        stop = start + count
        directory = self._make_directory(start, stop)

        calculator = as_calculator(self.calculator_name)

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

    def run(self, nprocs):
        """

        Args:
            nprocs: Number of processors to user.
        """
        with PythonScript(self.workdir / "analyze.py") as script:
            script.add_text("""

from abipy.ml.relax_scanner import RelaxScanner, Entries
from abipy.tools.printing import print_dataframe

entries = Entries.merge_from_topdir(".")
assert entries
df = entries.get_dataframe()
print_dataframe(df)
""")

        if nprocs == 1:
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
            with Pool(processes=nprocs) as pool:
                pool.map(func, args_list)

        entries = Entries.merge_from_topdir(self.workdir)
        df = entries.get_dataframe()
        print(df)
        df.to_csv(self.workdir / "configurations.csv")


def func(args):
    self, start, count = args
    self.run_start_count(start, count)



if __name__ == "__main__":
    #from abipy.ml.aseml import silence_tensorflow
    #silence_tensorflow()
    import sys
    filepath = sys.argv[1]
    isite = 0
    #isite = 1
    nx, ny, nz = 4, 4, 4
    nx, ny, nz = 2, 2, 2
    nx, ny, nz = 3, 3, 3
    optimizer_name = "chgnet"
    optimizer_name = "m3gnet"
    structure = Structure.from_file(filepath)
    #structure.perturb(distance=0.01)
    relax_mode = "ions"
    #relax_mode = "no"
    scanner = RelaxScanner(structure, isite, nx, ny, nz, optimizer_name, relax_mode=relax_mode, verbose=0)
    print(scanner)
    scanner.run(nprocs=2)

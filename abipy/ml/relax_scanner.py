"""
Objects to perform ASE calculations with machine-learned potentials.
"""
from __future__ import annotations

import os
import pickle
import itertools
import numpy as np
import pandas as pd

from pathlib import Path
from dataclasses import dataclass
#from typing import Type, Any, Optional, Union
from monty.json import MontyEncoder
from abipy.core import Structure
#from abipy.tools.plotting import get_ax_fig_plt #, get_axarray_fig_plt,
from abipy.tools.iotools import workdir_with_prefix
#from abipy.tools.printing import print_dataframe
from abipy.ml.aseml import relax_atoms, get_atoms, as_calculator, ase_optimizer_cls, RX_MODE, fix_atoms



@dataclass
class Entry:
    """
    Stores the relaxed structure with the associated energy and Cartesian forces.
    """
    isite: int
    structure: Structure   # pymatgen Structure
    energy: float          # Energy in eV
    forces: np.ndarray     # Forces in eV/Ang

    @classmethod
    def from_atoms_and_calculator(cls, isite, atoms, calculator) -> Entry:
        """
        Build an Entry from an Atoms object and a calculator.
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
        Build an Entry by reading the last Atoms object from an ASE trajectory.
        """
        atoms = traj[-1]
        return cls(isite=isite,
                   structure=Structure.as_structure(atoms),
                   energy=float(atoms.get_potential_energy()),
                   forces=atoms.get_forces(),
                   )

    def __eq__(self, other: Entry) -> bool:
        """Invoked by python to evaluate `self == other`."""
        return abs(self.energy - other.energy) / len(self.structure) < 1e-4 and \
               self.structure == other.structure


class Entries(list):
    """
    A list of Entry objects.
    """

    @classmethod
    def merge_from_topdir(cls, topdir=".") -> Entries:
        """
        Merge all entries starting from directory `topdir`.
        """
        top = str(topdir)
        dirpaths = [name for name in os.listdir(top) if os.path.isdir(os.path.join(top, name))
                    and name.startswith("start:")]

        entries = cls()
        for dpath in dirpaths:
            #print(f"{dpath=}")
            #tokens = dpath.split("-")
            #this_start = int(tokens[0].split(":"))
            #this_stop = int(tokens[1].split(":"))
            with open(Path(dpath) / "data.pickle", "rb") as fh:
                for e in pickle.load(fh):
                    if e in entries: continue
                    entries.append(e)

        return entries

    def get_dataframe(self) -> pd.DataFrame:
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
                #xcart0=x, xcart1=y, xcart2=z,
                xred0=fx, xred1=fy, xred2=fz,
                fmin=fmods.min(),
                fmax=fmods.max(),
                fmean=fmods.mean(),
                #fstd=fmods.std(),
                #drift=np.linalg.norm(self.forces.sum(axis=0)),
            ))

        return pd.DataFrame(dict_list).sort_values(by="energy", ignore_index=True)


class RelaxScanner:
    """
    This object employs an ASE calculator to perform a structural relaxation.
    for ...
    """

    def __init__(self, structure, isite, nx, ny, nz, calculator,
                 relax_mode: str = "ions", fmax: float = 1e-3, steps=500,
                 verbose: int = 0, optimizer_name="BFGS", pressure=0.0):
        """
        Args:
            structure: Structure object or any file supported by pymatgen providing a structure.
            isite: Index of the site to be displaced.
            nx, ny, nz: Size of the mesh giving the initial position of isite.
            calculator: ASE calculator or string with the name of the nn potential e.g. "chgnet" or "m3gnet".
            relax_mode: "ions" to relax ions only, "cell" for ions + cell, "no" for no relaxation.
            fmax: total force tolerance for relaxation convergence.
                Here fmax is a sum of force and stress forces. Defaults to 0.1.
            steps: max number of steps for relaxation.
            verbose: whether to print stdout.
            optimizer_name: name of the ASE optimizer class to use
            pressure: Target pressure.
        """
        self.initial_structure = Structure.as_structure(structure).copy()
        self.isite = isite
        assert len(self.initial_structure) > self.isite >= 0
        self.nx, self.ny, self.nz = nx, ny, nz
        self.calculator = as_calculator(calculator)
        self.relax_mode = relax_mode
        RX_MODE.validate(self.relax_mode)
        self.fmax = fmax
        self.steps = steps
        self.verbose = verbose
        self.optimizer_name = optimizer_name
        self.pressure = pressure

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
     calculator  = {self.calculator}
     verbose     = {self.verbose}

=== INITIAL STRUCTURE ===

{self.initial_structure}

"""

    def _mkdir(self, start, stop) -> Path:
        """
        """
        # Make sure we are not already running similar ranges
        top = "."
        dirpaths = [name for name in os.listdir(top) if os.path.isdir(os.path.join(top, name))
                    and name.startswith("start:")]
        for dpath in dirpaths:
            tokens = dpath.split("-")
            this_start = int(tokens[0].split(":")[1])
            this_stop = int(tokens[1].split(":")[1])
            if start >= this_start and start < this_stop:
                raise RuntimeError(f"{start=}, {stop=} range overlaps with {dpath=}")

        # Make sure we are not already running the same range.
        workdir = Path(f"start:{start}-stop:{stop}-nx:{self.nx}-ny:{self.ny}-nz:{self.nz}")
        if workdir.exists():
            raise RuntimeError(f"{workdir=} already exists!")

        workdir.mkdir()
        return workdir

    def build(self, nprocs: int, template_filepath: str, workdir=None, prefix=None):
        """
        Build directory with `prefix` if `workdir` is None else create it.
        Raise RuntimeError if workdir already exists.
        """
        workdir = workdir_with_prefix(workdir, prefix)

        # Read the Slurm template from the external file.
        template_lines = open(template_filepath ,"rt").readlines()

        # Split ntasks across nprocs.
        ntask = self.nx * self.ny * self.nz
        div, rest = divmod(ntask, nprocs)
        for iproc in range(nprocs):
            start, stop = iproc * div, (iproc + 1) * div
            count = div
            if iproc == nprocs - 1:
                stop += rest
                count = div + rest

            new_dir = self._mkdir(start, stop)
            with open(new_dir / "job.sh'", "wt")  as fh:
                text = "foobar"
                template = template_lines.append(text)
                fh.write(template)

        return workdir

    def run_start_count(self, start: int, count: int) -> None:
        """
        Invoke ASE to relax `ncount` structures starting from index `start`.
        Stores unique results in an internal list, finally save results to file in pickle format.

        Args:
            start: Index of initial structure.
            count: Number of structures to relax. If < 0, all structures are considered.
        """
        if count < 0:
            count = self.nx * self.ny * self.nz

        stop = start + count
        workdir = self._mkdir(start, stop)

        entries, errors = Entries(), []
        for cnt, (ix, iy, iz) in enumerate(itertools.product(range(self.nx), range(self.ny), range(self.nz))):
            if not (stop > cnt >= start): continue

            structure = self.initial_structure.copy()
            structure._sites[self.isite].frac_coords = np.array((ix/self.nx, iy/self.ny, iz/self.nz))
            atoms = get_atoms(structure)

            if self.relax_mode == RX_MODE.no:
                entry = Entry.from_atoms_and_calculator(self.isite, atoms, self.calculator)

            else:
                try:
                    fix_atoms(atoms, fix_inds=[i for i in range(len(atoms)) if i != self.isite])
                    relax = relax_atoms(atoms,
                                        relax_mode=self.relax_mode,
                                        optimizer=self.optimizer_name,
                                        fmax=self.fmax,
                                        pressure=self.pressure,
                                        #verbose=self.verbose,
                                        verbose=1,
                                        steps=self.steps,
                                        #opt_kwargs=None,
                                        traj_path=workdir / "relax.traj",
                                        calculator=self.calculator,
                    )

                    entry = Entry.from_traj(self.isite, relax.traj)

                except Exception:
                    entry = None
                    errors.append((ix, iy, iz))

            if errors:
                print(f"ASE relaxation raised {len(errors)} exceptions")

            if entry is not None and entry not in entries:
                entries.append(entry)

        # Save results in pickle format.
        with open(workdir / "data.pickle", "wb") as fh:
            pickle.dump(entries, fh)


if __name__ == "__main__":
    import sys
    filepath = sys.argv[1]
    isite = 0
    isite = 1
    nx, ny, nz = 4, 4, 4
    nx, ny, nz = 2, 2, 2
    nx, ny, nz = 3, 3, 3
    optimizer_name = "chgnet"
    optimizer_name = "m3gnet"
    structure = Structure.from_file(filepath)
    #structure.perturb(distance=0.01)
    relax_mode = "ions"
    #relax_mode = "no"
    scanner = RelaxScanner(structure, isite, nx, ny, nz, optimizer_name, relax_mode=relax_mode)
    print(scanner)
    scanner.run_start_count(0, -1)
    #scanner.build(ntasks: int, template_file: str, workdir=None, prefix=None):
    entries = Entries.merge_from_topdir(".")
    df = entries.get_dataframe()
    print(df)

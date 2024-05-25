"""
"""
from __future__ import annotations

import os
#import numpy as np
import abipy.core.abinit_units as abu
#try:
#    import ase
#except ImportError as exc:
#    raise ImportError("ase not installed. Try `pip install ase`.") from exc
from pathlib import Path
#from multiprocessing import Pool
#from typing import Type, Any, Optional, Union
#from enum import IntEnum
#from tabulate import tabulate
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from monty.string import list_strings  # marquee,
#from monty.functools import lazy_property
#from monty.json import MontyEncoder
#from monty.collections import AttrDict
#from ase import units
#from ase.atoms import Atoms
#from ase.io.trajectory import write_traj, Trajectory
from ase.io import read
#from ase.calculators.calculator import Calculator
#from ase.io.vasp import write_vasp_xdatcar, write_vasp
#from ase.stress import voigt_6_to_full_3x3_strain
#from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from abipy.core import Structure
from abipy.electrons.gsr import GsrFile
#from abipy.tools.iotools import workdir_with_prefix, PythonScript, yaml_safe_load_path
from abipy.tools.typing import PathLike
import abipy.flowtk.qutils as qu
#from abipy.tools.serialization import HasPickleIO
#from abipy.tools.context_managers import Timer
#from abipy.tools.parallel import get_max_nprocs, pool_nprocs_pmode
#from abipy.abio.enums import StrEnum, EnumMixin
#from abipy.core.mixins import TextFile, NotebookWriter
#from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_grid_legend,
#    set_visible, set_ax_xylabels, linear_fit_ax)


class ExtxyzIOWriter:
    """
    Example:

        # To find all the vasprun.xml files starting from the a top-level directory, use:

            xyz_writer = ExtxyzIOWriter.from_top(".", "vasprun.xml")
            print(xyz_writer)
            xyz_writer.write("foo.xyz")

        # A similar syntax can be used for Abinit GSR files:

            ExtxyzIOWriter.from_top(".", "GSR.nc").write("foo.xyz")

        # To specify an explicit list of files, use:

            xyz_writer = ExtxyzIOWriter(["dir1/vasprun.xml", "dir2/vasprun.xml"])
            xyz_writer.write("foo.xyz")
    """

    SUPPORTED_EXTS = [
        "vasprun.xml",
        "GSR.nc",
    ]

    @classmethod
    def from_top(cls, top: PathLike, ext: str):
         from monty.os.paths import find_exts
         filepaths = find_exts(str(top), ext)
         return cls(filepaths)

    def __init__(self, filepaths: list[PathLike]):
        self.filepaths = list_strings(filepaths)
        if not self.filepaths:
            raise RuntimeError("Empty list of filepaths!")

        for ext in self.SUPPORTED_EXTS:
            if all(f.endswith(ext) for f in self.filepaths):
                self.ext = ext
                break
        else:
            raise ValueError(f"Cannot detect extension from filepaths, should be in: {self.SUPPORTED_EXTS}")

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []
        for i, path in enumerate(self.filepaths):
            lines.append(f"[{i}]: {path}")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.to_string()

    def write(self, xyz_filepath: PathLike, overwrite: bool=False):
        """
        """
        if not overwrite and os.path.isfile(xyz_filepath):
            raise RuntimeError(f"Cannot overwrite pre-existent file: {xyz_filepath=}, use overwrite=True to allow overwriting.")

        with open(xzy_filepath, "wt") as fh:
            for atoms in self.yield_atoms():
                write(fh, atoms, format='extxyz', append=True)

    def yield_atoms(self):
        """
        """
        for filepath in self.filepaths:
            if self.ext == "vasprun.xml":
                vasprun = Vasprun(filepath)
                last_step = vasprun.ionic_steps[-1]
                structure, forces, stress = last_step["structure"], last_step["forces"], last_step["stress"]
                energy = get_energy_step(last_step)

            elif self.ext == "GSR.nc":
                with GsrFile(filepath) as gsr:
                      if not gsr.is_scf_run:
                          raise RuntimeError("GSR file was not produced by a SCF run!")
                      structure, forces, stress_gpa = gsr.structure, gsr.cart_forces, gsr.cart_stress_tensor
                      stress = stress_gpa / abu.eVA3_GPa
                      energy = float(gsr.energy)

            else:
                raise ValueError(f"Format {self.ext=} is not supported!")

            atoms = structure.to_ase_atoms()

            # Attach calculator with results.
            atoms.calc = SinglePointCalculator(atoms,
                                               energy=energy,
                                               free_energy=energy,
                                               forces=forces,
                                               stress=stress,
                                               )
            yield atoms


def check_vasp_success(vasprun, outcar, verbose: int = 1) -> bool:
    """
    Check if a VASP calculation completed successfully.

    Returns: True if the calculation completed successfully, False otherwise.
    """
    def my_print(*args, **kwargs):
        if verbose: print(*args, **kwargs)

    try:
        if not vasprun.converged:
            my_print("Calculation did not converge.")
            return False

        #outcar = Outcar(f"{directory}/OUTCAR")
        if outcar is not None:
            if outcar.run_stats.get("Elapsed time (sec)"):
                my_print("Calculation completed in {} seconds.".format(outcar.run_stats["Elapsed time (sec)"]))
            else:
                my_print("Elapsed time not found in OUTCAR.")
                return False

        my_print("Calculation completed successfully.")
        return True

    except Exception as e:
        my_print(f"Error checking calculation status: {e}")
        return False


class SinglePointRunner:
    """

    Usage example:

    .. code-block:: python

        runner = SinglePointRunner("out.traj", "outdir", traj_range=(0,-1,100), "vasp")
        runner.sbatch()
        runner.collect_xyz("foo.xyz")
    """

    def __init__(self, traj_path: PathLike, topdir: PathLike, traj_range: range, code: str = "vasp",
                 slurm_script: PathLike = "run.sh",
                 custodian_script: PathLike = "run_custodian.py",
                 verbose=0, **kwargs):
        """
        """
        self.traj_path = traj_path
        self.topdir = Path(str(topdir)).absolute()
        self.traj_range = traj_range
        if not isinstance(traj_range, range):
            raise TypeError(f"Got type{traj_range} instead of range")
        self.code = code

        err_msgs = []
        if not os.path.exists(slurm_script):
            open(slurm_script, "wt").write(qu.get_slurm_script()
            err_msgs.append("""\
No template for slurm submission script has been found. A default template that requires customization has been generated for you!""")
        else:
            self.slurm_script = open(slurm_script, "rt").read()

        if code == "vasp":
            if not os.path.exists(custodian_script):
                open(custodian_script, "wt").write(qu.get_custodian_script()
                err_msgs.append("""\
No template for custodian script has been found. A default template that requires customization has been generated for you!""")
            else:
                self.custodian_script = open(slurm_script, "rt").read()

        if err_msgs:
            raise RuntimeError("\n".join(err_msgs))

        self.verbose = int(verbose)
        self.kwargs = kwargs

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []
        app = lines.append

        return "\n".join(lines)

    #def get_custom_incar_settings():
    #    # Define custom INCAR settings
    #    custom_incar_settings = {
    #        'ENCUT': 600,         # Override plane-wave cutoff energy
    #        'ISMEAR': -5,         # Use tetrahedron method with Blöchl corrections
    #        'SIGMA': 0.01,        # Smearing width
    #        'EDIFF': 1E-8,        # Electronic energy convergence criterion
    #        'NSW': 0,             # Number of ionic steps (static calculation)
    #        'ISIF': 2,            # Stress tensor calculation
    #        'LREAL': 'Auto',      # Projection operators (automatic)
    #        'LWAVE': False,       # Do not write WAVECAR
    #        'LCHARG': True        # Write CHGCAR
    #    }
    #    return custom_incar_settings

    def sbatch(self, max_jobs: int=100) -> int:
        """
        """
        if not self.topdir.exists(): self.topdir.mkdir()

        num_jobs = 0
        for index in self.traj_range:
            workdir = self.topdir / f"SINGLEPOINT_{index}"
            if workdir.exists():
                print("{workdir=} already exists. Ignoring it")
                continue

            try:
                atoms = read(self.traj_path, index=index)
            except StopIteration as exc:
                print(f"ASE trajectory does not have more that {index=} configurations. Exiting sbatch loop!")
                break

            structure = Structure.as_structure(atoms)
            workdir.mkdir()
            script_filepath = workdir / "run.sh"

            if self.code == "vasp":
                # Generate VASP input files using the Materials Project settings for a single-point calculation
                from pymatgen.io.vasp.sets import MPStaticSet
                vasp_input_set = MPStaticSet(structure, **self.kwargs)
                vasp_input_set.write_input(workdir)
                with open(workdir / "run_custodian.py", wt) as fh:
                    fh.write(self.custodian_script)

            else:
                raise ValueError(f"Unsupported {self.code=}")

            with open(script_filepath, "wt") as fh:
                fh.write(self.slurm_template)

            queue_id = qu.slurm_sbatch(script_filepath)
            num_jobs += 1
            if num_jobs == max_jobs:
                print(f"Reached {max_jobs=}, will stop firing new jobs!")

        return num_jobs

    def write_xyz(self, xyz_filepath: PathLike, dryrun=False) -> None:
        """
        """
        ext = {
            "vasp": "vasprun.xml",
            "abinit": "GSR.nc",
        }[self.code]

        writer = ExtxyzIOWriter.from_top(self.workdir, ext)
        writer.write(xyz_filepath)


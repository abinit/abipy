"""
"""
from __future__ import annotations

import os
import shutil
import abipy.core.abinit_units as abu
from pathlib import Path
from ase.calculators.singlepoint import SinglePointCalculator
from monty.string import list_strings  # marquee,
from monty.termcolor import cprint
from ase.io import read
from ase.stress import full_3x3_to_voigt_6_stress
from ase.io import write
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.io.vasp.sets import MatPESStaticSet # , MPStaticSet
from abipy.core import Structure
from abipy.electrons.gsr import GsrFile
from abipy.tools.typing import PathLike
import abipy.flowtk.qutils as qu
from abipy.ml.tools import get_energy_step


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
        """
        Scan for files with extension ext starting from the top directory top.
        """
        from monty.os.path import find_exts
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

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []
        for i, path in enumerate(self.filepaths):
            lines.append(f"[{i}]: {path}")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.to_string()

    def write(self, xyz_filepath: PathLike, overwrite: bool = False):
        """
        """
        if not overwrite and os.path.isfile(xyz_filepath):
            raise RuntimeError(f"Cannot overwrite pre-existent file: {xyz_filepath=}, use overwrite=True to allow overwriting.")

        with open(xyz_filepath, "wt") as fh:
            for atoms in self.yield_atoms():
                write(fh, atoms, format='extxyz', append=True)

    def yield_atoms(self):
        """
        """
        for filepath in self.filepaths:
            if self.ext == "vasprun.xml":
                vasprun = Vasprun(filepath)
                dirname = os.path.dirname(filepath)
                outcar_path = os.path.join(dirname, "OUTCAR")
                outcar = Outcar(outcar_path) if os.path.exists(outcar_path) else None

                ok = check_vasp_success(vasprun, outcar, verbose=1)

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
            stress = full_3x3_to_voigt_6_stress(stress)

            # Attach calculator with results.
            atoms.calc = SinglePointCalculator(atoms,
                                               energy=energy,
                                               free_energy=energy,
                                               forces=forces,
                                               stress=stress,
                                               )
            yield atoms


def check_vasp_success(vasprun: Vasprun, outcar: Outcar | None, verbose: int = 1) -> bool:
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

        traj_range = range(0, -1, 100)
        runner = SinglePointRunner("out.traj", "outdir", traj_range)
        runner.sbatch()
        runner.collect_xyz("foo.xyz")
    """
    slurm_script_name = "run.sh"

    custodian_script_name = "run_custodian.py"

    def __init__(self, traj_path: PathLike, traj_range: range,
                 topdir: PathLike = ".", code: str = "vasp",
                 vasp_set_cls=MatPESStaticSet,
                 verbose: int = 0):
        """
        Args:
            traj_path: Path to ASE trajectory file.
            traj_range:
            topdir:
            code:
            verbose:
        """
        self.traj_path = traj_path
        self.topdir = Path(str(topdir)).absolute()
        self.traj_range = traj_range
        if not isinstance(traj_range, range):
            raise TypeError(f"Got {type(traj_range)} instead of range")
        self.code = code

        err_lines = []
        slurm_body = ""

        if code == "vasp":
            self.vasp_set_cls = vasp_set_cls
            slurm_body = f"python {self.custodian_script_name}"
            if not os.path.exists(self.custodian_script_name):
                open(self.custodian_script_name, "wt").write(qu.get_custodian_template())
                err_lines.append(f"""\
No custodian script: {self.custodian_script_name} has been found in {str(self.topdir)}.
A template that requires customization has been generated for you!""")
            else:
                self.custodian_script_str = open(self.custodian_script_name, "rt").read()

        if not os.path.exists(self.slurm_script_name):
            open(self.slurm_script_name, "wt").write(qu.get_slurm_template(slurm_body))
            err_lines.append(f"""\
No slurm submission script: {self.slurm_script_name} has been found in {str(self.topdir)}.
A template that requires customization has been generated for you!""")
        else:
            self.slurm_script_str = open(self.slurm_script_name, "rt").read()

        if err_lines:
            raise RuntimeError("\n".join(err_lines))

        self.verbose = int(verbose)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []
        app = lines.append

        return "\n".join(lines)

    def sbatch(self, max_jobs: int = 100) -> list[int]:
        """
        Submit max_jobs SinglePoint calculations with structures taken from the ASE trajectory file.
        """
        if not self.topdir.exists(): self.topdir.mkdir()

        job_ids = []
        for index in self.traj_range:
            workdir = self.topdir / f"SINGLEPOINT_{index}"
            if workdir.exists():
                print(f"{str(workdir)} already exists. Ignoring it")
                continue

            try:
                atoms = read(self.traj_path, index=index)
            except StopIteration as exc:
                print(f"ASE trajectory does not have more that {index+1} configurations. Exiting sbatch loop!")
                break

            structure = Structure.as_structure(atoms)
            workdir.mkdir()

            if self.code == "vasp":
                # Generate VASP input files using the Materials Project settings for a single-point calculation

                user_incar_settings = {
                    "NCORE": 2,
                    'LWAVE': False,       # Do not write WAVECAR
                    'LCHARG': False,      # Do not Write CHGCAR
                }
                vasp_input_set = self.vasp_set_cls(structure, user_incar_settings=user_incar_settings)
                vasp_input_set.write_input(workdir)
                with open(workdir / self.custodian_script_name, "wt") as fh:
                    fh.write(self.custodian_script_str)

            else:
                raise ValueError(f"Unsupported {self.code=}")

            try:
                job_id = qu.slurm_write_and_sbatch(workdir / "run.sh", self.slurm_script_str)

            except Exception as exc:
                cprint(exc, "red")
                cprint("Job sumbission failed. Will remove directory and exit sbatch loop.", color="red")
                shutil.rmtree(workdir)
                break

            job_ids.append(job_id)
            if len(job_ids) == max_jobs:
                print(f"Reached {max_jobs=}, will stop firing new jobs!")

        return job_ids

    def write_xyz(self, xyz_filepath: PathLike, dry_run=False) -> None:
        """
        """
        ext = {
            "vasp": "vasprun.xml",
            "abinit": "GSR.nc",
        }[self.code]

        writer = ExtxyzIOWriter.from_top(self.topdir, ext)
        writer.write(xyz_filepath)

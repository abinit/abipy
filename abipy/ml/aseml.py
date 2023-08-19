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
import warnings
import stat
import numpy as np
import pandas as pd
try:
    import ase
except ImportError as exc:
    raise ImportError("ase not installed. Try `pip install ase`.") from exc

from pathlib import Path
from dataclasses import dataclass
from inspect import isclass
from typing import Type, Any, Optional, Union
from monty.string import marquee, list_strings # is_string,
from monty.functools import lazy_property
from pymatgen.core import Structure as PmgStructure
from pymatgen.io.ase import AseAtomsAdaptor
from ase import units
from ase.atoms import Atoms
from ase.io.trajectory import write_traj, Trajectory
from ase.optimize.optimize import Optimizer
from ase.calculators.calculator import Calculator
from ase.io.vasp import write_vasp_xdatcar
from ase.neb import NEB
from ase.md.nptberendsen import NPTBerendsen, Inhomogeneous_NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from abipy.core import Structure
from abipy.tools.plotting import get_ax_fig_plt #, get_axarray_fig_plt,
from abipy.tools.iotools import workdir_with_prefix
from abipy.tools.printing import print_dataframe
from abipy.abio.enums import StrEnum, EnumMixin

###################
# Helper functions
###################

_CELLPAR_KEYS = ["a", "b", "c", "angle(b,c)", "angle(a,c)", "angle(a,b)"]


ASENEB_METHODS = ['aseneb', 'eb', 'improvedtangent', 'spline', 'string']


class RX_MODE(EnumMixin, StrEnum):  # StrEnum added in 3.11
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
    new_structure = structure.abisanitize(**kwargs)
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


@dataclass
class AseResults:
    """
    Container with the results produced by the ASE calculator.
    """
    atoms: Atoms
    ene: float
    stress: np.ndarray
    forces: np.ndarray

    @classmethod
    def from_inds(cls, trajectory, *inds) -> AseResults:
        return [cls.from_atoms(trajectory[i]) for i in inds]

    @classmethod
    def from_atoms(cls, atoms: Atoms) -> AseResults:
        """Build the object from atoms with a calculator."""
        return cls(atoms=atoms,
                   ene=float(atoms.get_potential_energy()),
                   stress=atoms.get_stress(voigt=False),
                   forces=atoms.get_forces())

    @property
    def pressure(self) -> float:
        return -self.stress.trace() / 3

    @property
    def volume(self) -> float:
        """Volume of unit cell."""
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
        """Dictionary with statistics on forces."""
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
        r0, r1 = AseResults.from_inds(self.traj, 0, -1)
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
        optimizer: name of the ASE optimizer class to use
        fmax: total force tolerance for relaxation convergence.
            Here fmax is a sum of force and stress forces. Defaults to 0.1.
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
    """Silence every unnecessary warning from tensorflow."""
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
    Add _delta_forces and _delta_stress attributes to an ASE calculator.
    Extend `calculate` so that forces and stresses are corrected accordingly.
    """
    def set_delta_forces(self, delta_forces):
        """F_Abinitio - F_ML"""
        self._delta_forces = delta_forces

    def get_delta_forces(self):
        return getattr(self, "_delta_forces", None)

    def set_delta_stress(self, delta_stress):
        """S_Abinitio - S_ML"""
        self._delta_stress = delta_stress

    def get_delta_stress(self):
        return getattr(self, "_delta_stress", None)

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
            forces += delta_forces
            #print("Updating forces with delta_forces:\n", forces)
            self.results.update(
                forces=forces,
            )

        # Apply delta correction to stress.
        stress = self.results["stress"]
        delta_stress = self.get_delta_stress()
        if delta_stress is not None:
            stress += delta_stress
            #print("Updating stress with delta_stress:\n", stress)
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
    Supports different backends defined by `name`.
    """
    def __init__(self, name: str, **kwargs):
        self.name = name
        self._model = None

    def __str__(self):
        return f"{self.__class__.__name__}: name: {self.name}"

    def get_calculator(self) -> Calculator:
        """Return ASE Calculator with ML potential."""
        if self.name == "m3gnet_old":
            # Legacy version.
            if self._model is None:
                silence_tensorflow()
            try:
                from m3gnet.models import Potential, M3GNet, M3GNetCalculator
            except ImportError as exc:
                raise ImportError("m3gnet not installed. Try `pip install m3gnet`.") from exc

            if self._model is None:
                self._model = Potential(M3GNet.load())

            class MyM3GNetCalculator(M3GNetCalculator, _MyCalculatorMixin):
                """Add delta_forces"""

            return MyM3GNetCalculator(potential=self._model)

        if self.name == "m3gnet":
            # See https://github.com/materialsvirtuallab/matgl
            try:
                import matgl
                from matgl.ext.ase import M3GNetCalculator
            except ImportError as exc:
                raise ImportError("matgl not installed. Try `pip install matgl`.") from exc

            if self._model is None:
                self._model = matgl.load_model("M3GNet-MP-2021.2.8-PES")

            class MyM3GNetCalculator(M3GNetCalculator, _MyCalculatorMixin):
                """Add delta_forces"""

            return MyM3GNetCalculator(potential=self._model)

        if self.name == "chgnet":
            try:
                from chgnet.model.dynamics import CHGNetCalculator
                from chgnet.model.model import CHGNet
            except ImportError as exc:
                raise ImportError("chgnet not installed. Try `pip install chgnet`.") from exc

            if self._model is None:
                self._model = CHGNet.load()

            class MyCHGNetCalculator(CHGNetCalculator, _MyCalculatorMixin):
                """Add delta_forces"""

            return MyCHGNetCalculator(model=self._model)

        raise ValueError(f"Invalid {self.name=}")


class _MlBase:
    """
    Base class for all Ml subclasses. Provides helper functions to
    perform typical tasks such as writing files in the workdir.
    """
    def __init__(self, workdir, prefix=None):
        """
        Build directory with `prefix` if `workdir` is None else create it.
        Raise RuntimeError if workdir already exists.
        """
        self.workdir = workdir_with_prefix(workdir, prefix)
        self.basename_info = []
        self.delta_forces = None
        self.delta_stress = None

    def set_delta_forces_stress_from_abiml_nc(self, filepath: str) -> None:
        """
        Read ab-initio forces from a netcdf file produced by ABINIT and use
        these values to set the delta corrections in the calculator.
        """
        if os.path.basename(filepath) != "ABIML_RELAX_IN.nc": return

        # Read structure, forces and stresses from the nc file produced by ABINIT.
        from abipy.iotools import ETSF_Reader
        with ETSF_Reader(filepath) as reader:
            abi_forces = reader.read_value("cart_forces")
            abi_stress = reader.read_value("cart_stress")
            structure = reader.read_structure()

        # Compute ML forces for this structure.
        atoms = to_ase_atoms(structure, calc=self.get_calculator())
        ml_forces = atoms.get_forces()
        ml_stress = atoms.get_stress()

        # Set delta forces/stresses so that the next invokation to get_calculator include the deltas
        self.set_delta_forces_stress(abi_forces - ml_forces, abi_stress - ml_stress)

    def set_delta_forces_stress(self, delta_forces, delta_stress) -> None:
        """Set the value of the delta corrections."""
        self.delta_forces = delta_forces
        self.delta_stress = delta_stress

    def __str__(self):
        # Delegated to the subclass.
        return self.to_string()

    def get_calculator(self) -> Calculator:
        """Return ASE calculator."""
        calc = self.calc_builder.get_calculator()

        if self.delta_forces is not None:
            print("Setting delta_forces:\n", self.delta_forces)
            calc.set_delta_forces(self.delta_forces)
        if self.delta_stress is not None:
            print("Setting delta_stress:\n", self.delta_stress)
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

    def write_script(self, basename: str, text: str, info: str) -> Path:
        """Write text script to basename."""
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
    Relax structure with ASE and ML-potential
    """

    def __init__(self, atoms: Atoms, relax_mode, fmax, pressure, steps, optimizer, calc_builder, verbose,
                 workdir, prefix=None):
        """
        Args:
            atoms: ASE atoms to relax.
            relax_mode:
            fmax:
            pressure:
            steps:
            optimizer:
            calc_builder:
            verbose:
        """
        super().__init__(workdir, prefix)
        self.atoms = atoms
        self.relax_mode = relax_mode
        RX_MODE.validate(relax_mode)
        self.fmax = fmax
        self.steps = steps
        self.optimizer = optimizer
        self.pressure = pressure
        self.calc_builder = calc_builder
        self.verbose = verbose

        self.atoms.calc = self.get_calculator()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     relax_mode  = {self.relax_mode}
     fmax        = {self.fmax}
     steps       = {self.steps}
     optimizer   = {self.optimizer}
     pressure    = {self.pressure}
     calculator  = {self.calc_builder}
     workdir     = {self.workdir}
     verbose     = {self.verbose}

=== ATOMS ===

{self.atoms}

"""

    def run(self) -> None:
        """Run structural relaxation."""
        workdir = self.workdir

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

        #from abipy.core.structure import StructDiff
        #diff = StructDiff(["INPUT", "ML RELAXED"], [initial_atoms, self.atoms])

        self._finalize()


class MlMd(_MlBase):
    """Perform MD calculations with ASE and ML potential."""

    def __init__(self, atoms: Atoms, temperature, timestep, steps, loginterval,
                 ensemble, calc_builder, verbose, workdir, prefix=None):
        """
        Args:
            atoms:
            temperature:
            timestep:
            steps:
            loginterval:
            ensemble:
            calc_builder:
            verbose:
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
        self.calc_builder = calc_builder
        self.verbose = verbose

        self.atoms.calc = self.get_calculator()

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
        workdir = self.workdir

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
""", info="Python script to compute visualize energies vs Time.")

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

    def __init__(self, atoms_list: list[Atoms], calc_builder, verbose,
                 workdir, prefix=None):
        """
        Args:
            atoms_list: List of ASE atoms
            calc_builder:
            verbose:
        """
        super().__init__(workdir, prefix)
        self.atoms_list = atoms_list
        self.calc_builder = calc_builder
        self.verbose = verbose

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level `verbose`."""
        return f"""\

{self.__class__.__name__} parameters:

     calculator  = {self.calc_builder}
     workdir     = {self.workdir}
     verbose     = {self.verbose}

"""

    def run(self) -> None:
        """Run list of GS calculations."""
        workdir = self.workdir

        results = []
        for atoms in self.atoms_list:
            atoms.calc = self.get_calculator()
            results.append(AseResults.from_atoms(atoms))

        #dict_list = [r.get_dict4pandas() for r in results]
        #df = pd.DataFrame(dict_list)

        neb_data = self.postprocess_images(self.atoms_list)

        self._finalize()


class MlNeb(_MlNebBase):
    """
    Perform NEB calculation with ASE and ML potential.
    """

    def __init__(self, initial_atoms: Atoms, final_atoms: Atoms,
                 nimages, neb_method, climb, optimizer, relax_mode, fmax, pressure,
                 calc_builder, verbose, workdir, prefix=None):
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
            calc_builder:
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
        self.calc_builder = calc_builder
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
     calculator  = {self.calc_builder}
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

        if self.relax_mode != RX_MODE.no:
            relax_kws = dict(calculator=self.get_calculator(),
                             optimizer=self.optimizer,
                             relax_mode=self.relax_mode,
                             fmax=self.fmax,
                             pressure=self.pressure,
                             verbose=self.verbose,
                             )

            print(f"Relaxing initial image with relax mode: {self.relax} ...")
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

        # Optimize
        opt_class = ase_optimizer_cls(self.optimizer)
        nebtraj_file = str(workdir / "neb.traj")
        logfile = self.get_path("neb.log", "Log file of NEB calculation.")
        optimizer = opt_class(neb, trajectory=nebtraj_file, logfile=logfile)

        #print("Starting NEB algorithm with optimizer:", optimizer, "...")
        optimizer.run(fmax=self.fmax)

        # To read the last nimages atoms e.g. 5: read('neb.traj@-5:')
        images = ase.io.read(f"{str(nebtraj_file)}@-{self.nimages}:")
        write_vasp_xdatcar(workdir / "XDATCAR", images,
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
                 calc_builder, verbose, workdir, prefix=None):
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
            calc_builder:
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
        self.calc_builder = calc_builder
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
     calculator  = {self.calc_builder}
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
                           self.relax_mode, self.fmax, self.pressure, self.calc_builder, self.verbose, camp_dirs[i])
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
        neb_method:
        climb: True to use a climbing image.
        method:
    """
    images = [initial]
    images += [initial.copy() for i in range(nimages - 2)]
    images += [final]

    apply_constraint = None
    if initial.constraints:
        if not final.constraints:
            raise RuntimeError("Both initial and final points should have constraints!")
        for ci, cf in zip(initial.constraints, final.constraints): #, strict=True):
            if ci.__class__ != cf.__class__:
                raise RuntimeError(f"Constraints in initial and final points should belong to the same class: {ci}, {cf}")
        apply_constraint = True

    # Set calculators
    for image, calculator in zip(images, calculators): #, strict=True):
        image.calc = calculator

    # Compute energy/forces for the extrema in order to have them in the trajectory.
    _ = AseResults.from_inds(images, 0, -1)

    neb = NEB(images, method=neb_method, climb=climb)
    # Interpolate linearly the positions of the middle images
    neb.interpolate(method=method, mic=mic, apply_constraint=apply_constraint)

    return neb


class MlOrderer(_MlBase):
    """
    Order a disordered structure using pymatgen and ML potential.
    """
    def __init__(self, structure, max_ns, optimizer, relax_mode, fmax, pressure, steps, calc_builder, verbose,
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
            calc_builder:
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
        self.calc_builder = calc_builder
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
     calculator  = {self.calc_builder}
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
                r0, r1 = AseResults.from_inds(trajectory, 0, -1)
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
                 relax_mode, fmax, pressure, steps, optimizer, calc_builder,
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
            calc_builder,
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
        self.calc_builder = calc_builder
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
     calculator = {self.calc_builder}
     workdir    = {self.workdir}
     verbose    = {self.verbose}

=== ATOMS ===

{self.atoms}
"""
        return s

    def run(self) -> None:
        """Run MlPhonons."""
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

            r0, r1 = AseResults.from_inds(relax.traj, 0, -1)
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

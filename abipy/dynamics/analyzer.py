"""
Tools to analyze MD trajectories and compute diffusion coefficients.
"""
from __future__ import annotations

import dataclasses
import warnings
import json
import os
import numpy as np
import pandas as pd
import pymatgen.core.units as units

from pathlib import Path
from scipy.stats import linregress
from scipy import optimize
from matplotlib.offsetbox import AnchoredText
from monty.functools import lazy_property
from monty.bisect import find_le
from monty.string import list_strings, marquee
from monty.collections import AttrDict #, dict2namedtuple
from monty.termcolor import cprint
from pymatgen.util.string import latexify
from pymatgen.core.lattice import Lattice
from pymatgen.core.composition import Composition
#from abipy.core.mixins import TextFile # , NotebookWriter
from abipy.core.structure import Structure
from abipy.tools.typing import Figure, PathLike
from abipy.tools.serialization import HasPickleIO
from abipy.tools.iotools import try_files # , file_with_ext_indir
from abipy.tools.context_managers import Timer
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, get_color_symbol,
                                  set_ticks_fontsize, set_logscale, linear_fit_ax)
from abipy.tools.parallel import pool_nprocs_pmode #, get_max_nprocs
from abipy.dynamics.cpx import parse_file_with_blocks, EvpFile


__author__ = "Giuliana Materzanini, Tommaso Chiarotti, Matteo Giantomassi"


Ang2PsTocm2S = 0.0001
e2s = 1.602188**2 # electron charge in Coulomb scaled by 10.d-19**2
kbs = 1.38066     # Boltzmann constant in Joule/K scaled by 10.d-23

#bohr2a = 0.529177
#kBoltzAng10m9 = 1.38066e-02
#kBoltz = 1.38066e-23
kBoltzEv = 8.617333e-05
#nCar = 56  # FIXME: Harcoded


def common_oxidation_states() -> dict:
    oxi2symbols = {
        +1: ["Li", "Na", "K", "Rb", "Cs", "Fr"],
        +2: ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"],
    }

    # map: symbol to oxidation state.
    symb2oxi = {}
    for oxi, slist in oxi2symbols.items():
        for s in slist:
            symb2oxi[s] = oxi
    return symb2oxi


def read_structure_postac_ucmats(traj_filepath: PathLike, step_skip: int) -> tuple[Structure, np.ndarray, np.ndarray, int]:
    """
    Read all configurations from an ASE trajectory file.

    Args:
        traj_filepath: File path.
        step_skip: Sampling frequency.
            time_step should be multiplied by this number to get the real time between measurements.

    Returns:
        tuple with: initial Structure, (nsteps, natom, 3) array with the Cartesian coords,
        (nsteps,3,3) array with cell vectors.
    """
    from ase.io import read
    traj = read(traj_filepath, index=":")
    traj_len = len(traj)
    structure = Structure.as_structure(traj[0])

    pos_tac, ucmats = [], []
    for it in range(0, traj_len, step_skip):
        atoms = traj[it]
        pos_tac.append(atoms.positions)
        ucmats.append(atoms.cell.array)

    del traj
    pos_tac = np.array(pos_tac, dtype=float)
    ucmats = np.array(ucmats, dtype=float)

    return structure, pos_tac, ucmats, traj_len


def parse_lammps_input(filepath: PathLike, verbose: int = 0):
    """
    Extract parameters, atoms and structure from a LAMMPS input file.
    Returns namedtuple with results.
    """
    dirpath = Path(os.path.dirname(filepath)).absolute()
    from pymatgen.io.lammps.inputs import LammpsInputFile
    inp = LammpsInputFile.from_file(filepath, ignore_comments=False)

    atom_style = inp.get_args("atom_style")
    units = inp.get_args("units")
    timestep = float(inp.get_args("timestep"))
    # dump ID group-ID style N file attribute1 attribute2 ...
    dump = inp.get_args("dump")
    if not dump:
        raise ValueError(f"Cannot find dump command in LAMMPS input file: {filepath}")
    tokens = dump.split()
    loginterval = int(tokens[3])
    traj_filepath = dirpath / tokens[4]

    read_data = dirpath / inp.get_args("read_data")
    log = inp.get_args("log")
    if not log:
        log = "log.lammps"
    log = dirpath / log

    if verbose:
        print(f"{atom_style=}, {units=}, {timestep=}, {read_data=}")

    # Build atoms from read_data.
    from ase.io.lammpsdata import read_lammps_data
    atoms = read_lammps_data(
        read_data,
        Z_of_type=None,
        sort_by_id=False,
        units=units,
        style=atom_style,
    )

    # Extract temperature from the fix command.
    temperature = None
    fix_list = inp.get_args("fix")
    for fix in fix_list:
        tokens = fix.split()
        try:
            i = tokens.index("temp")
        except ValueError:
            continue
        temp = float(tokens[i+1])
        if temperature is not None and temp != temperature:
            raise ValueError(f"Found two different temperatures {temp=} and {temperature=} in LAMMPS input file: {filepath}")
        temperature = temp

    if temperature is None:
        raise ValueError(f"Cannot detect temperature in LAMMPS input file: {filepath}")

    return AttrDict(atoms=atoms, structure=Structure.as_structure(atoms),
                    units=units, temperature=temperature, timestep=timestep, loginterval=loginterval,
                    traj_filepath=traj_filepath, log=log, lammps_input=inp)


class MdAnalyzer(HasPickleIO):
    """
    High-level interface to read MD trajectories and metadata from external files,
    compute the MSQD and plot the results.
    """

    @classmethod
    def from_abiml_dir(cls, directory: PathLike, step_skip: int = 1) -> MdAnalyzer:
        """
        Build an instance from a directory containing an ASE trajectory file and
        a JSON file with the MD parameters as produced by the `abiml.py md` script.
        """
        directory = Path(str(directory))
        structure, pos_tac, ucmats, traj_len = read_structure_postac_ucmats(directory / "md.traj", step_skip)

        # Read metadata from the JSON file.
        with open(directory / "md.json", "rt") as fh:
            meta = json.load(fh)

        temperature = meta["temperature"]
        # Convert from ASE fs to ps
        timestep = meta["timestep"] * 1e-3
        loginterval = meta["loginterval"]
        engine = meta["nn_name"]
        times = (np.arange(0, traj_len) * timestep * loginterval)[::step_skip].copy()

        log_path = try_files([directory / "md.aselog", directory / "md.log"])
        from abipy.ml.aseml import AseMdLog
        with AseMdLog(log_path) as log:
            evp_df = log.df.copy()

        return cls(structure, temperature, times, pos_tac, ucmats, engine, evp_df=evp_df)

    @classmethod
    def from_hist_file(cls, hist_filepath: PathLike, step_skip: int=1) -> MdAnalyzer:
        """
        Build an instance from an ABINIT HIST.nc file.
        """
        from abipy.dynamics.hist import HistFile
        with HistFile(hist_filepath) as hist:
            structure = hist.structure.copy()
            #hist.r.read_dimvalue("time")
            pos_tac = hist.r.read_value("xcart") * units.bohr_to_ang
            ucmats = hist.r.read_value("rprimd") * units.bohr_to_ang
            #temperature = None
            #evp_df = None

        #times = (np.arange(0, traj_len) * r.timestep * r.loginterval)[::step_skip].copy()
        if step_skip != 1:
            ucmats = ucmats[::step_skip].copy()
            pos_tac = pos_tac[::step_skip].copy()

        raise NotImplementedError()
        #return cls(structure, temperature, times, pos_tac, ucmats, "abinit", evp_df=evp_df)

    @classmethod
    def from_vaspruns(cls, filepaths: list) -> MdAnalyzer:
        """
        Build an instance from a list of Vasprun files (must be ordered in sequence of MD simulation).
        """
        def get_structures(vaspruns):
            # This piece of code is shamelessy taken from
            # https://github.com/materialsvirtuallab/pymatgen-analysis-diffusion/blob/master/pymatgen/analysis/diffusion/analyzer.py
            for i, vr in enumerate(vaspruns):
                if i == 0:
                    step_skip = vr.ionic_step_skip or 1
                    final_structure = vr.initial_structure
                    temperature = vr.parameters["TEEND"]
                    timestep = vr.parameters["POTIM"] # fs
                    yield step_skip, temperature, timestep

                # check that the runs are continuous
                from pymatgen.util.coord import pbc_diff
                fdist = pbc_diff(
                    vr.initial_structure.frac_coords, final_structure.frac_coords
                )
                if np.any(fdist > 0.001):
                    raise ValueError("initial and final structures do not match.")
                final_structure = vr.final_structure

                assert (vr.ionic_step_skip or 1) == step_skip
                for s in vr.ionic_steps:
                    yield s["structure"]

        from pymatgen.io.vasp.outputs import Vasprun
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vaspruns = [Vasprun(path) for path in list_strings(filepaths)]

        s = get_structures(vaspruns)
        step_skip, temperature, timestep = next(s)

        # Extract Cartesian positions.
        pos_tac = []
        for i, strc in enumerate(s):
            if i == 0: structure = strc
            pos_tac.append(strc.coords)

        nsteps, natom = i + 1, len(structure)
        pos_tac = np.reshape(pos_tac, (nsteps, natom, 3))
        times = np.arange(0, nsteps) * timestep * step_skip
        evp_df = None

        return cls(structure, temperature, times, pos_tac, "vasp", evp_df=evp_df)

    #@classmethod
    #def from_qe_dir(cls, directory: PathLike, step_skip: int=1):
    #    traj_filepath = file_with_ext_indir(".lammpstrj", directory)
    #    return cls.from_qe_input(filepath: PathLike, step_skip=step_skip):

    @classmethod
    def from_qe_input(cls, filepath: PathLike, step_skip: int=1):
        """
        Build an instance from a CP/PW input file.

        Conventions for Quantum ESPRESSO input files:
           ".pwi" -> pw.x
           ".cpi" -> cp.x
        """
        filepath = Path(str(filepath)).absolute()
        engine = None
        if filepath.name.endswith(".pwi"):
            engine = "qepw"
            default_prefix = "pwscf"
        if filepath.name.endswith(".cpi"):
            engine = "qecp"
            default_prefix = "cp"

        if engine is None:
            raise ValueError("QE input file should end with .pwi for pw.x or .cpi for cp.x")

        # Use ASE to parse QE input file.
        from ase.io.espresso import read_espresso_in, read_fortran_namelist
        with open(filepath, "rt") as fh:
            atoms = read_espresso_in(fh)
            fh.seek(0)
            sections, card_lines = read_fortran_namelist(fh)

        structure = Structure.as_structure(atoms)
        natoms = len(structure)

        # https://www.quantum-espresso.org/Doc/INPUT_CP.html
        control = sections["control"]
        calculation = control["calculation"]

        # time step for molecular dynamics, in Hartree atomic units
        # (1 a.u.=2.4189 * 10^-17 s : beware, PW code use Rydberg atomic units, twice that much!!!)
        timestep = control.get("dt", 1.0)
        timestep = timestep * 2.4189 * 1e-5
        if engine == "qepw": timestep *= 2
        isave = control.get("isave", 100)
        prefix = control.get("prefix", default_prefix)
        outdir = Path(control.get("outdir", filepath.cwd()))

        # ELECTRONS section
        #electrons = sections["electrons"]
        # IONS section
        ions = sections["ions"]
        temperature = ions.get("tempw", 300)
        # CELL section
        #cell = sections["cell"]

        # Get atomic positions from the pos file.
        pos_filepath = outdir / (prefix + ".pos")
        pos_nsteps, pos_tac, pos_headers = parse_file_with_blocks(pos_filepath, natoms)
        pos_tac = np.reshape(pos_tac, (pos_nsteps, natoms, 3))

        # Get lattices from the cell file.
        cel_filepath = outdir / (prefix + ".cel")
        cel_nsteps, ucmats, pos_headers = parse_file_with_blocks(cel_filepath, 3)
        ucmats = np.reshape(ucmats, (cel_nsteps, 3, 3))

        if pos_nsteps != cel_nsteps:
            raise RuntimeError(f"Found differrent no. iterations in pos and cel file: {pos_nsteps=} != {cel_nsteps=}")

        times = np.arange(0, pos_nsteps) * timestep * step_skip
        if step_skip != 1:
            ucmats = ucmats[::step_skip].copy()
            pos_tac = pos_tac[::step_skip].copy()

        # Extract energies from CP EVP file
        with EvpFile(outdir / (prefix + ".evp")) as evp:
            evp_df = evp.df.copy()

        return cls(structure, temperature, times, pos_tac, ucmats, engine, evp_df=evp_df)

    @classmethod
    def from_lammps_dir(cls, directory: PathLike, step_skip: int=1, basename="in.lammps") -> MdAnalyzer:
        """
        Build an instance from a directory containing a LAMMPS input file.
        """
        return cls.from_lammp_input(Path(str(directory)) / basename, step_skip=step_skip)

    @classmethod
    def from_lammps_input(cls, input_filepath: PathLike, step_skip: int=1) -> MdAnalyzer:
        """
        Build an instance from a LAMMPS input file.

        Args:
            input_filepath: LAMMPS input file.
        """
        r = parse_lammps_input(input_filepath, verbose=0)
        #structure = r.structure

        structure, pos_tac, ucmats, traj_len = read_structure_postac_ucmats(r.traj_filepath, step_skip)

        if r.units != "metal":
            from ase.calculators.lammps import convert
            pos_tac = convert(pos_tac, "distance", r.units, "metal")
            ucmats = convert(ucmats, "distance", r.units, "metal")

        times = (np.arange(0, traj_len) * r.timestep * r.loginterval)[::step_skip].copy()

        # Extract dataframe from the log file.
        evp_df = None
        if r.log.exists():
            from pymatgen.io.lammps.outputs import parse_lammps_log
            df_list = parse_lammps_log(r.log)
            evp_df = df_list[-1]

        return cls(structure, r.temperature, times, pos_tac, ucmats, "lammps", evp_df=evp_df)

    @classmethod
    def from_lammpstrj(cls, traj_filepath: PathLike, input_filepath: PathLike, step_skip: int = 1) -> MdAnalyzer:
        """
        Build an instance from a LAMMPS trajectory file and a log file.

        Args:
            traj_filepath:
            input_filepath:
        """
        structure, pos_tac, ucmats, traj_len = read_structure_postac_ucmats(traj_filepath, step_skip)

        if input_filepath.endswith(".evp"):
            # Extract times from CP EVP file (Giuliana's way)
            temperature = None
            with EvpFile(input_filepath) as evp:
                evp_df = evp.df.copy()
                timestep = evp.times[1] - evp.times[0]
                loginterval = 1
        else:
            raise NotImplementedError()

        times = (np.arange(0, traj_len) * timestep * loginterval)[::step_skip].copy()

        return cls(structure, temperature, times, pos_tac, ucmats, "lammps", evp_df=evp_df)

    def __init__(self,
                structure: Structure,
                temperature: float,
                times: np.ndarray,
                cart_positions: np.ndarray,
                ucmats: np.ndarray,
                engine: str,
                pos_order: str = "tac",
                evp_df=None | pd.DataFrame,
                ):
        """
        Args:
            structure: Structure object (first geometry of the MD run).
            temperature: Temperature in Kelvin.
            times: Array with times in ps units.
            cart_positions: Cartesian positions in Ang. Default shape: (nt, natom, 3).
            ucmats: Array of lattice matrix of every step. Used for NPT.
                For NVT-AIMD, the lattice at each time step is set to the lattice in the "structure" argument.
            pos_order: "tac" if cart_positions has shape (nt, natom, 3).
                       "atc" if cart_positions has shape (natom, nt, 3).
            evp_df:
        """
        self.structure = structure
        self.times = times
        self.engine = engine

        if pos_order == "tac":
            self.pos_atc = cart_positions.transpose(1, 0, 2).copy()
        elif pos_order == "atc":
            self.pos_atc = cart_positions
        else:
            raise ValueError(f"Invalid {pos_order=}")

        self.evp_df = evp_df

        if np.all(ucmats[it] == ucmats[0] for it in range(len(ucmats))):
            self.lattices = None
        else:
            self.lattices = np.array([Lattice(mat) for mat in ucmats])

        self.latex_formula = self.structure.latex_formula
        self.temperature = temperature
        self.set_color_symbol("VESTA")
        self.verbose = 0

        self.consistency_check()

    def consistency_check(self) -> None:
        """
        Perform internal consistency check.
        """
        if self.pos_atc.shape != (self.natom, self.nt, 3):
            raise ValueError(f"Invalid shape {self.pos_atc.shape=}, expecting: {(self.natom, self.nt, 3)}")
        if len(self.times) != self.nt:
            raise ValueError(f"{len(self.times)=} != {self.nt=}")

        # Check times mesh.
        ierr = 0
        for it in range(self.nt-1):
            dt = self.times[it+1] - self.times[it]
            if abs(dt - self.timestep) > 1e-3:
                ierr += 1
                if ierr < 10: print(f"{dt=} != {self.timestep=}")
        if ierr:
            raise ValueError(f"Time-mesh is not linear. There are {ierr} points with wrong timestep")

        if self.lattices is not None and len(self.lattices) != self.nt:
            raise ValueError(f"{len(self.lattices)=} != {self.nt=}")

    def get_params_dict(self) -> dict:
        """Dictionary with the most important parameters."""
        attr_names = [
            "latex_formula", "temperature", "timestep", "nt", "max_time", "natom", "avg_volume", "engine",
        ]
        d = {aname: getattr(self, aname) for aname in attr_names}
        return d

    def deepcopy(self) -> MdAnalyzer:
        """Deep copy of the object."""
        import copy
        return copy.deepcopy(self)

    def iter_structures(self):
        """Generate pymatgen structures."""
        species = self.structure.species
        pos_tac = self.pos_atc.transpose(1, 0, 2).copy()

        if self.lattices is None:
            # Same lattice.
            const_lattice = self.structure.lattice
            for coords in pos_tac:
                yield Structure(const_lattice, species, coords, coords_are_cartesian=True)
        else:
            for coords, lattice in zip(pos_tac, self.lattices):
                yield Structure(lattice.copy(), species, coords, coords_are_cartesian=True)

    #def iter_atoms(self):
    #    """Generate ASE atoms"""
    #    Atoms(symbols=None,
    #          positions=None, numbers=None,
    #          tags=None, momenta=None, masses=None,
    #          magmoms=None, charges=None,
    #          scaled_positions=None,
    #          cell=None, pbc=None, celldisp=None,
    #          constraint=None,
    #          calculator=None,
    #          info=None,
    #          velocities=None)

    def resample_step(self, start_at_step: int, take_every: int) -> MdAnalyzer:
        """
        Resample the trajectory. Start at iteration start_at_step and increase
        the timestep by taking every `take_every` iteration.
        """
        return self.resample_time(start_time=start_at_step * self.timestep, new_timestep=take_every * self.timestep)

    def resample_time(self, start_time: float, new_timestep: float) -> MdAnalyzer:
        """
        Resample the trajectory. Start at time `start_time` and use new timestep `new_timestep`.
        """
        # TODO: This is not true anymore!
        # NB: Cannot change the object in place as SigmaBerend and DiffusionData keep a reference to self.
        new = self.deepcopy()

        old_timestep = new.times[1] - new.times[0]
        if not (new.times[-1] > start_time > new.times[0]):
            raise ValueError(f"Invalid start_time should be between {new.times[0]} and {new.times[-1]})")

        it0 = int(start_time / old_timestep)
        if it0 != 0:
            new.pos_atc = new.pos_atc[:,it0:,:]
            new.times = new.times[it0:] - new.times[it0]
            if new.lattices is not None:
                new.lattices = new.lattices[it0:]
            if self.evp_df is not None:
                self.evp_df = self.evp_df.iloc[it0:]

        if new_timestep < old_timestep:
            raise ValueError(f"Invalid {new_timestep=} should be >= {old_timestep}")

        istep = int(new_timestep / old_timestep)
        if istep != 1:
            new.pos_atc = new.pos_atc[:,::istep,:].copy()
            new.times = new.times[::istep] - new.times[0]
            if new.lattices is not None:
                new.lattices = new.lattices[::istep].copy()
            if self.evp_df is not None:
                self.evp_df = self.evp_df.iloc[::istep]

        new.consistency_check()
        return new

    @property
    def timestep(self) -> float:
        """Timestep in ps."""
        return self.times[1] - self.times[0]

    @property
    def max_time(self) -> float:
        """Maximum simulation time in ps."""
        return self.times[-1]

    @property
    def nt(self) -> int:
        """Number of points in the MD trajectory."""
        return self.pos_atc.shape[1]

    @lazy_property
    def natom(self) -> int:
        """Number of atoms."""
        return len(self.structure)

    @property
    def temperature(self) -> float:
        """Temperature in Kelvin."""
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        """Set temperature in Kelvin."""
        self._temperature = value

    @property
    def verbose(self) -> int:
        """Verbosity level."""
        return self._verbose

    @verbose.setter
    def verbose(self, value: int):
        """Set temperature in Kelvin."""
        self._verbose = value

    @property
    def engine(self) -> str:
        return self._engine

    @engine.setter
    def engine(self, value):
        """Set engine string."""
        self._engine = value

    @property
    def formula(self) -> str:
        """Returns the formula as a string."""
        return self.structure.formula

    @property
    def latex_formula(self) -> str:
        return self._latex_formula

    @latex_formula.setter
    def latex_formula(self, value):
        """LaTeX formatted formula. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$."""
        self._latex_formula = latexify(value)

    @property
    def latex_formula_n_temp(self) -> str:
        return f"{self.latex_formula}\nT = {self.temperature} K"

    @property
    def latex_avg_volume(self) -> str:
        return "V$_{\mathrm{ave}}$ = " + f"{self.avg_volume:.2f}" + '$\mathrm{{\AA}^3}$'

    @property
    def avg_volume(self) -> float:
        """Average unit cell volume in Ang^3."""
        if self.lattices is None:
            return self.structure.lattice.volume
        return np.mean([lat.volume for lat in self.lattices])

    def set_color_symbol(self, dict_or_string: dict | str) -> None:
        """
        Set the dictionary mapping chemical_symbol --> color
        used in the matplotlib plots.

        Args:
            dict_or_string: "VESTA", "Jmol"
        """
        if isinstance(dict_or_string, dict):
            self.color_symbol = dict_or_string
        else:
            self.color_symbol = get_color_symbol(style=dict_or_string)

        for symbol in self.structure.symbol_set:
            if symbol not in self.color_symbol:
                raise KeyError(f"Cannot find {symbol=} in color_symbol dictionary!")

    def get_it_ts(self, t0: float) -> tuple[int, np.ndarray]:
        """
        Return the index of time t0 in self.times and the array with the time values.
        """
        if t0 < self.times[0] or t0 > self.times[-1]:
            raise ValueError(f"Invalid {t0=}. It should be between {self.times[0]} and {self.times[-1]}")

        it0 = find_le(self.times, t0)
        return it0, self.times[it0:] - self.times[it0]

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose."""
        lines = []
        app = lines.append
        app(marquee("MD PARAMETERS", mark="="))
        app(pd.Series(self.get_params_dict()).to_string())
        if verbose:
            app(self.structure.spget_summary(verbose=verbose))
        app("\n")

        return "\n".join(lines)

    def iatoms_with_symbol(self, symbol: str, atom_inds=None) -> np.ndarray:
        """
        Array with the index of the atoms with the given chemical symbol.
        If atom_inds is not None, filter sites accordingly.
        """
        iatoms = [iat for iat in range(len(self.structure)) if self.structure[iat].specie.symbol == symbol]
        if atom_inds is not None:
            iatoms = [iat for iat in iatoms if iat in atom_inds]
        if not iatoms:
            raise ValueError(f"Empty list of iatoms indices for {symbol=} and {atom_inds=}")

        return np.array(iatoms)

    def _select_symbols(self, symbols) -> list[str]:
        if symbols == "all": return sorted(self.structure.symbol_set)
        return list_strings(symbols)

    def get_sqdt_iatom(self, iatom: int, it0: int = 0) -> np.array:
        """
        Compute the square displacement vs time for a given atomic index
        starting from time index it0.
        """
        return ((self.pos_atc[iatom,it0:] - self.pos_atc[iatom,it0]) ** 2).sum(axis=1)

    def get_sqdt_symbol(self, symbol: str, it0: int = 0, atom_inds=None) -> np.array:
        """
        Compute the square displacement vs time averaged over atoms with the same chemical symbol
        starting from time index it0. atoms_inds adds an additional filter on the site index.
        """
        for count, iatom in enumerate(self.iatoms_with_symbol(symbol, atom_inds=atom_inds)):
            if count == 0:
                sqdt = self.get_sqdt_iatom(iatom, it0=it0)
            else:
                sqdt += self.get_sqdt_iatom(iatom, it0=it0)

        sqdt /= (count + 1)
        return sqdt

    def get_dw_symbol(self, symbol, t0: float = 0.0, tmax=None, atom_inds=None):
        """
        Compute diffusion coefficient by performing a naive linear regression of the raw MQST.
        """
        it0, ts = self.get_it_ts(t0)
        it1 = -1
        if tmax is not None:
            it1, _ = self.get_it_ts(tmax)
            ts = ts[:it1]

        sqdt = self.get_sqdt_symbol(symbol, it0=it0, atom_inds=atom_inds)
        fit = linregress(ts[:it1], sqdt[:it1])
        naive_d = fit.slope * Ang2PsTocm2S / 6
        label = r"Linear fit: D={:.2E} cm$^2$/s, $r^2$={:.2f}".format(naive_d, fit.rvalue**2)
        return AttrDict(naive_d=naive_d, ts=ts, fit=fit, label=label)

    def get_msdtt0_symbol_tmax(self, symbol: str, tmax: float, atom_inds=None, nprocs=None) -> Msdtt0:
        r"""
        Calculates the MSD for every possible pair of time points using the formula:

            $$MSD(t,t_0) = \frac{1}{N} \sum_{i=1}^{N} (\vec{r}_i(t+t_0) - \vec{r}_i(t_0))^2$$

        where $N$ is the number of particles with the given symbol, and $\vec{r}_i(t)$ is the position vector.

        Args:
            symbols:
            tmax:
            atoms_ins
        """
        index_tmax, _ = self.get_it_ts(tmax)
        iatoms = self.iatoms_with_symbol(symbol, atom_inds=atom_inds)
        tac = self.pos_atc[iatoms].transpose(1, 0, 2).copy()
        arr_tt0 = msd_tt0_from_tac(tac, index_tmax, nprocs=nprocs)

        return Msdtt0(arr_tt0=arr_tt0, mda=self, index_tmax=index_tmax, symbol=symbol)

    @add_fig_kwargs
    def plot_sqdt_atoms(self, symbols="all", t0: float = 0.0, atom_inds=None,
                        ax=None, xy_log=None, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement of atoms vs time.

        Args:
            symbols: List of chemical symbols to consider.
            t0: Initial time in ps.
            atom_inds: List of atom indices to include. None to disable filtering.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
        """
        it0, ts = self.get_it_ts(t0)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            for iatom in self.iatoms_with_symbol(symbol, atom_inds=atom_inds):
                sd_t = self.get_sqdt_iatom(iatom, it0=it0)
                ax.plot(ts, sd_t, label=f"${symbol}_{{{iatom}}}$")

        ax.legend(fontsize=fontsize, loc="upper left")
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel(r'square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        #set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.latex_formula_n_temp}\n{self.latex_avg_volume}\n" +
                                   "sd(t, $t_0$ =" + str(int(self.times[it0])) + " ps)",
                                   loc="upper right", prop=dict(size=fontsize)))
        return fig

    @add_fig_kwargs
    def plot_sqdt_symbols(self, symbols, t0: float = 0.0, atom_inds=None, with_dw=0,
                          ax=None, xy_log=None, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement averaged over all atoms of the same specie vs time.

        Args:
            symbols: List of chemical symbols to consider. "all" for all symbols in structure.
            t0: Initial time in ps.
            atom_inds: List of atom indices to include. None to disable filtering.
            with_dw: If != 0 compute diffusion coefficient via least-squares fit in the time-interval [t0, with_dw].
                If with_dw < 0 e.g. -1 the time-interval is set to [t0, tmax].
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
        """
        it0, ts = self.get_it_ts(t0)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            ax.plot(ts, self.get_sqdt_symbol(symbol, it0=it0, atom_inds=atom_inds),
                    label=symbol + ' msd($t, t_0$ = ' + str(t0) + ' ps)',
                    color=self.color_symbol[symbol],
                    )

            if with_dw != 0:
                tmax = self.times[-1] if with_dw < 0 else with_dw
                dw = self.get_dw_symbol(symbol, t0=t0, tmax=with_dw, atom_inds=atom_inds)
                ax.plot(dw.ts, dw.fit.slope*dw.ts + dw.fit.intercept,
                        label=symbol + " " + dw.label, color=self.color_symbol[symbol],
                        )

        ax.legend(fontsize=fontsize, loc="upper left")
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel('mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        #set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.latex_formula_n_temp}\n{self.latex_avg_volume}",
                                   loc="upper right", prop=dict(size=fontsize)))
        return fig

    @add_fig_kwargs
    def plot_sqdt_symbols_tmax(self, symbols, tmax: float, atom_inds=None, nprocs=None,
                               ax=None, xy_log=None, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement averaged over all atoms of the same specie vs time.

        Args:
            symbols: List of chemical symbols to consider. "all" for all symbols in structure.
            tmax: Max time in ps.
            atom_inds: List of atom indices to include. None to disable filtering.
            nprocs: Number of procs to use.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
        """
        index_tmax, _ = self.get_it_ts(tmax)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            iatoms = self.iatoms_with_symbol(symbol, atom_inds=atom_inds)
            tac = self.pos_atc[iatoms].transpose(1, 0, 2).copy()
            msd_tt0 = msd_tt0_from_tac(tac, index_tmax, nprocs=nprocs)
            msd_t = np.mean(msd_tt0, axis=1)
            t_start = self.nt - index_tmax
            ts = self.times[t_start:] - self.times[t_start]

            ax.plot(ts, msd_t,
                    label=symbol + " <msd($t$, $t_0$)>$\{$t_0$\}$, $t$ = [0, " + str(int(self.times[index_tmax])) + " ps]",
                    color=self.color_symbol[symbol],
                    )

        set_axlims(ax, xlims, "x")
        ax.legend(fontsize=fontsize, loc="upper left")
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        #set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.latex_formula_n_temp}\n{self.latex_avg_volume}",
                                   loc="upper right", prop=dict(size=fontsize)))
        return fig

    @add_fig_kwargs
    def plot_lattices(self, what_list=("abc", "angles", "volume"), ax_list=None,
                      xy_log=None, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot lattice lengths/angles/volume as a function of time.

        Args:
            what_list: List of strings specifying the quantities to plot. Default all
            ax_list: List of axis or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
        """
        if self.lattices is None:
            cprint("MD simulation has been performed with fixed lattice!", color="red")
            return None

        what_list = list_strings(what_list)
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=1, ncols=len(what_list),
                                                sharex=True, sharey=False, squeeze=False)
        markers = ["o", "^", "v"]

        cnt = -1
        if "abc" in what_list:
            # plot lattice parameters.
            for i, label in enumerate(["a", "b", "c"]):
                cnt += 1
                ax = ax_list[cnt]
                ax.plot(self.times, [lattice.abc[i] for lattice in self.lattices],
                        label=label, marker=markers[i])
            ax.set_ylabel("abc (A)")

        if "angles" in what_list:
            # plot lattice angles.
            for i, label in enumerate(["alpha", "beta", "gamma"]):
                cnt += 1
                ax = ax_list[cnt]
                ax.plot(self.times, [lattice.angles[i] for lattice in self.lattices],
                        label=label, marker=markers[i])
            ax.set_ylabel(r"$\alpha\beta\gamma$ (degree)")

        if "volume" in what_list:
            # plot lattice volume.
            marker = "o"
            cnt += 1
            ax = ax_list[cnt]
            ax.plot(self.times, [lattice.volume for lattice in self.lattices],
                    label="Volume", marker=marker)
            ax.set_ylabel(r'$V\, (A^3)$')

        for ix, ax in enumerate(ax_list):
            set_logscale(ax, xy_log)
            set_axlims(ax, xlims, "x")
            if ix == len(ax_list) - 1:
                ax.set_xlabel('t (ps)', fontsize=fontsize)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig


@dataclasses.dataclass(kw_only=True)
class Msdtt0:
    r"""
    This object stores:

        $$MSD(t,t_0) = \frac{1}{N} \sum_{i=1}^{N} (\vec{r}_i(t+t_0) - \vec{r}_i(t_0))^2$$

    where $N$ is the number of particles of a particular chemical symbol and $\vec{r}_i(t)$ is the position vector.
    """
    index_tmax: int
    symbol: str
    arr_tt0: np.ndarray
    mda: MdAnalyzer

    @property
    def times(self) -> np.ndarray:
        return self.mda.times

    @property
    def temperature(self) -> float:
        return self.mda.temperature

    @lazy_property
    def msd_t(self) -> np.ndarray:
        """Average of MSD(t,t_0) over t0."""
        return np.mean(self.arr_tt0, axis=1)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level verbose."""
        lines = []
        app = lines.append
        r = self.get_linfit_results()
        app(r.label)
        return "\n".join(lines)

    def get_linfit_results(self):
        """
        Perform linear fit.
        """
        t_start = self.mda.nt - self.index_tmax
        ts = self.times[t_start:] - self.times[t_start]
        fit = linregress(ts, self.msd_t)
        naive_d = fit.slope * Ang2PsTocm2S / 6
        label = r"Linear fit: D={:.2E} cm$^2$/s, $r^2$={:.2f}".format(naive_d, fit.rvalue**2)
        return AttrDict(naive_d=naive_d, fit=fit, label=label)

    @add_fig_kwargs
    def plot(self, ax=None, xy_log=None, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot <msd($t, t_0$)>$ averaged over the initial time t0.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
        """
        index_tmax = self.index_tmax
        t_start = self.mda.nt - index_tmax
        ts = self.times[t_start:] - self.times[t_start]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(ts, self.msd_t,
                label=self.symbol + " <msd($t, t_0$)>$\{$t_0$\}$, t = [0, " + str(int(self.times[index_tmax])) + " ps]",
                color=self.mda.color_symbol[self.symbol],
                )

        # Linear fit.
        r = self.get_linfit_results()
        ax.plot(ts, r.fit.slope * ts + r.fit.intercept, label=r.label)

        set_axlims(ax, xlims, "x")
        ax.legend(fontsize=fontsize, loc="upper left")
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        #set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.mda.latex_formula_n_temp}\n{self.mda.latex_avg_volume}",
                                   loc="upper right", prop=dict(size=fontsize)))
        return fig

    def get_sigma_berend(self, t1: float, t2: float, nblock_step: int=1, tot_block: int=1000) -> SigmaBerend:
        """
        Args:
            t1:
            t2:
            nblock_step
            tot_block
        """
        # choose the time elapsed
        it1, _ = sfind_ind_val(self.times, t1)
        it2, _ = sfind_ind_val(self.times, t2)
        size_t = self.arr_tt0.shape[0]

        if it1 >= it2:
            raise ValueError(f"For input {t1=} and {t2=}, got {it1=} >= {it2=}")
        if it1 >= size_t:
            raise ValueError(f"For input {t1=}, got {it1=} >= {size_t=}")
        if it2 >= size_t:
            raise ValueError(f"For input {t2=}, got {it2=} >= {size_t=}")

        block_sizes1, sigmas1, delta_sigmas1 = sigma_berend(nblock_step, tot_block, self.arr_tt0[it1,:])
        block_sizes2, sigmas2, delta_sigmas2 = sigma_berend(nblock_step, tot_block, self.arr_tt0[it2,:])

        # Build instance from locals dict.
        mda = self.mda
        latex_formula = mda.latex_formula
        temperature = mda.temperature
        time1, time2 = mda.times[it1], mda.times[it2]
        data = locals()
        return SigmaBerend(**{k: data[k] for k in [field.name for field in dataclasses.fields(SigmaBerend)]})

    def get_diffusion_with_sigma(self,
                                 block_size1: int, block_size2: int,
                                 fit_time_start: float, fit_time_stop: float,
                                 sigma_berend: SigmaBerend) -> DiffusionData:
        """
        Compute diffusion coefficient with uncertainty.
        """
        sigmas1, block_sizes1 = sigma_berend.sigmas1, sigma_berend.block_sizes1
        sigmas2, block_sizes2 = sigma_berend.sigmas2, sigma_berend.block_sizes2
        it1, it2 = sigma_berend.it1, sigma_berend.it2

        mda = self.mda
        symbol = self.symbol
        temperature = mda.temperature
        latex_formula = mda.latex_formula
        engine = mda.engine
        avg_volume = mda.avg_volume
        ncarriers = len(mda.structure.indices_from_symbol(symbol))
        times = self.times

        ib1, size1 = sfind_ind_val(block_sizes1, block_size1)
        sig1 = sigmas1[ib1]
        print(f"{sig1=}")
        ib2, size2 = sfind_ind_val(block_sizes2, block_size2)
        sig2 = sigmas2[ib2]
        print(f"{sig2=}")

        # fit to a linear behaviour the errors
        mSigma = (sig2 - sig1) / (it2-it1)
        qSigma = sig1 - mSigma * it1
        mDataInBlock = (size2-size1) / (it2-it1)
        qDataInBlock = size1 - mDataInBlock*it1

        # and find error for anytime.
        msd_tt0 = self.arr_tt0
        err_msd = np.zeros(msd_tt0.shape[0], dtype=float)
        for t in range(msd_tt0.shape[0]):
            err_msd[t] = abs(mSigma * t + qSigma)

        # and find error for anytime
        dataScorrelated = np.zeros(msd_tt0.shape[0], dtype=int)
        for t in range(msd_tt0.shape[0]):
            dataScorrelated[t] = int(mDataInBlock * t + qDataInBlock)

        # average over the initial times.
        msd_t = np.mean(msd_tt0, axis=1)

        fit_istart, _ = sfind_ind_val(times, fit_time_start)
        fit_istop, _ = sfind_ind_val(times, fit_time_stop)

        msdSScorrelated, timeArrScorrelated, errMSDScorrelated = [], [], []
        counter, condition = fit_istart, True
        while condition:
            if counter >= fit_istop:
                condition = False
            else:
                index = dataScorrelated[counter]
                msdSScorrelated.append(msd_t[counter])
                timeArrScorrelated.append(times[counter])
                errMSDScorrelated.append(err_msd[counter])
                counter = counter + index

        msdSScorrelated = np.array(msdSScorrelated, dtype=float)
        timeArrScorrelated = np.array(timeArrScorrelated, dtype=float)
        errMSDScorrelated = np.array(errMSDScorrelated, dtype=float)

        best_angcoeff, quote, var_angcoeff, var_quote = linear_lsq_linefit(msdSScorrelated, timeArrScorrelated, 1/(errMSDScorrelated)**2)
        min_angcoeff = best_angcoeff - np.sqrt(var_angcoeff)
        max_angcoeff = best_angcoeff + np.sqrt(var_angcoeff)

        diffusion_coeff = best_angcoeff * Ang2PsTocm2S / 6
        err_diffusion_coeff = np.sqrt(var_angcoeff) * Ang2PsTocm2S / 6
        # TODO: Charge from oxidation state?
        conductivity = e2s / kbs * ncarriers * diffusion_coeff / avg_volume / temperature * 1.e09
        err_conductivity = e2s / kbs * ncarriers / avg_volume / temperature * err_diffusion_coeff * 1.e09

        print(f"{ncarriers=} for {symbol=}")
        print('{:.2E}'.format(best_angcoeff*Ang2PsTocm2S/6, 2))
        print(best_angcoeff)
        print(max_angcoeff)
        print(min_angcoeff)

        # Build instance from locals dict.
        data = locals()
        return DiffusionData(**{k: data[k] for k in [field.name for field in dataclasses.fields(DiffusionData)]})

        """
        ('timeStepJump = ' + str(timeStepJump) + '\n')
        ('time, cel and pos were cut before ' + str(timeArray[0]+tInitial)+ 'ps' + '\n')
        ('t elapsed max is ' + str(estimatedTMax)+ 'ps' + '\n')
        ('trajectory length is ' + str(timeArrayTmp[timeArrayTmp.shape[0]-1])+ 'ps' + '\n')
        ('error on msd(t) evaluated at ' +  str(estimatedFirstTElapsed) + 'ps' + ' and ' + str(estimatedSecondTElapsed) + 'ps' + '\n')
        ('evaluated n. of blocks at ' + str(estimatedFirstTElapsed) + 'ps' + ' is ' + str(block_size1)  + '\n')
        ('evaluated n. of blocks at ' + str(estimatedSecondTElapsed) + 'ps' + ' is ' + str(block_size2)  + '\n')
        ('number of decorrelated msd(t) data that we fit: ' + str(timeArrScorrelated.shape[0]) + '\n')
        """

    @add_fig_kwargs
    def plot_mat(self, cmap="jet", fontsize=8, ax=None, **kwargs) -> Figure:
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        im = ax.matshow(self.arr_tt0, cmap=cmap)
        fig.colorbar(im, ax=ax)
        ax.set_xlabel('t0 (ps)', fontsize=fontsize)
        ax.set_ylabel('t (ps)', fontsize=fontsize)
        ax.set_title(self.symbol + ": msd($t, t_0$)", fontsize=fontsize)
        return fig



class Msdtt0List(list):
    """
    A list of Msdtt0 objects.
    """
    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level verbose."""
        lines = []
        app = lines.append
        for i, msdtt0 in enumerate(self):
            app(msdtt0.to_string(verbose=verbose))
        return "\n".join(lines)

    @add_fig_kwargs
    def plot(self, sharex=True, sharey=True, fontsize=8, **kwargs) -> Figure:
        """
        Plot all Msdtt0 objects on a grid.
        """
        nrows, ncols = len(self), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=sharex, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()
        if len(self) % ncols != 0: ax_list[-1].axis("off")

        for ix, (msdtt0, ax) in enumerate(zip(self, ax_list)):
            msdtt0.plot(ax=ax, fontsize=fontsize, show=False, **kwargs)

        return fig



@dataclasses.dataclass(kw_only=True)
class SigmaBerend:
    """
    Stores the variance of correlated data as function of block number.
    """
    temperature: float
    latex_formula: str

    it1: int
    time1: float
    block_sizes1: np.ndarray
    sigmas1: np.ndarray
    delta_sigmas1: np.ndarray

    it2: int
    time2: float
    block_sizes2: np.ndarray
    sigmas2: np.ndarray
    delta_sigmas2: np.ndarray

    @add_fig_kwargs
    def plot(self, fontsize=8, ax_list=None, **kwargs) -> Figure:
        """
        Plot variance of correlated data as function of block number.
        """
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=1, ncols=2,
                                                sharex=True, sharey=True, squeeze=False)

        for ix, ax in enumerate(ax_list.ravel()):
            xs, ys, yerr, it, time = self.block_sizes1, self.sigmas1, self.delta_sigmas1, self.it1, self.time1
            if ix == 1:
                xs, ys, yerr, it, time = self.block_sizes2, self.sigmas2, self.delta_sigmas2, self.it2, self.time2

            ax.errorbar(xs, ys,
                        yerr=yerr, linestyle='-', # linewidth=0.5,
                        label="$\sigma(\mathrm{MSD}($" + '%2.1f' % time +" ps$))$ "+ '\n' +
                              self.latex_formula + ', '+ 'T = %4.0f' % self.temperature + 'K')
            ax.legend(fontsize=fontsize, loc="lower right")
            ax.set_xlabel('N. of data in block', fontsize=fontsize)
            ax.set_ylabel('$\sigma$ ($\AA^2$)', fontsize=fontsize)
            ax.grid(True)

        fig.suptitle("Variance of correlated data as function of block number")
        return fig


@dataclasses.dataclass(kw_only=True)
class DiffusionData(HasPickleIO):
    """
    Diffusion results for a given temperature.
    """
    diffusion_coeff: float
    err_diffusion_coeff: float
    conductivity: float
    err_conductivity: float
    temperature: float
    symbol: str
    latex_formula: str
    avg_volume: float
    ncarriers: int
    block_size1: int
    block_size2: int
    fit_time_start: float  # msd(t) fit starts at this time.
    fit_time_stop: float   # msd(t) fit ends at this time.
    min_angcoeff: float  # min angular coefficient of msd(t)
    max_angcoeff: float  # max angular coefficient of msd(t)
    best_angcoeff: float # best angular coefficient of msd(t)
    engine: str
    msd_t: np.ndarray
    err_msd: np.ndarray
    sigma_berend: SigmaBerend
    times: np.ndarray
    # TODO: Find better names
    quote: float
    timeArrScorrelated: np.ndarray
    msdSScorrelated: np.ndarray
    errMSDScorrelated: np.ndarray

    @add_fig_kwargs
    def plot(self, ax=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot MDS(t) with errors.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ts = self.times[:self.msd_t.shape[0]]

        ax.errorbar(ts, self.msd_t, yerr=self.err_msd, color='mediumblue', label=self.symbol)
        ax.errorbar(self.timeArrScorrelated, self.msdSScorrelated, yerr=self.errMSDScorrelated, linestyle='-')
        ax.errorbar(ts, self.best_angcoeff * ts + self.quote, linestyle='--')
        ax.errorbar(ts, self.min_angcoeff * ts + self.quote, linestyle='--')
        ax.errorbar(ts, self.max_angcoeff * ts + self.quote, linestyle='--')

        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel(r'$\mathrm{MSD}_\mathrm{tr}$ $\mathrm{(\AA}^2\mathrm{)}$', fontsize=fontsize)
        ax.add_artist(AnchoredText(
            'D$_{tr}$ = (' + str('{:.2E}'.format(self.diffusion_coeff)) + '\u00B1' +
            str('{:.2E}'.format(self.err_diffusion_coeff)) + ') cm$^2$/s',
            loc="upper left", prop=dict(size=fontsize)))
        ax.legend(fontsize=fontsize, loc="lower right")
        ax.add_artist(AnchoredText(f"{self.latex_formula}\nT = {self.temperature} K",
                                   loc="upper right", prop=dict(size=fontsize)))
        return fig


class DiffusionDataList(list):
    """
    A list of DiffusionData objects.
    """

    #@classmethod
    #def from_topdir(cls, topdir):
    #    return cls.from_files(filepaths)

    #@classmethod
    #def from_files(cls, filepaths):
    #    new = cls()
    #    for path in filepaths:
    #        new.append(DiffusionData.from_file(path))
    #    return new

    #def _nrows_ncols_nplots(self, size=None):
    #    size = size or len(self)
    #    nrows, ncols, nplots = 1, 1, size
    #    if nplots > 1:
    #        ncols = 2; nrows = nplots // ncols + nplots % ncols
    #    return nrows, ncols, nplots

    #def filter(self, filter_dict: dict) -> DiffusionDataList:
    #    """
    #    filter_dict = [{symbol: "Li"}
    #    """
    #    new = DiffusionDataList()
    #    for df_data in self:
    #        if all(getattr(df_data, a_name) == a_value for (a_name, a_value) in filter_dict.items()):
    #            new.append(df_data)
    #    return new

    #def write_csv(self, filename: PathLike) -> None:
    #    """
    #    Writes data to a file that can be easily plotted in other software.

    #    Args:
    #        filename: Supported formats are csv and dat.
    #            If the extension is csv, a csv file is written. Otherwise, a dat format is assumed.
    #    """
    #    delimiter = ", "
    #    with open(filename, "wt") as f:
    #        f.write("T,diffusion,err_diffusion,volume,symbol,composition\n")
    #        #for dt, msd, msdc, mscd in zip(
    #        #    self.dt, self.msd, self.msd_components, self.mscd
    #        #):
    #        #    f.write(delimiter.join([str(v) for v in [dt, msd, *list(msdc), mscd]]))
    #        #    f.write("\n")

    def get_dataframe(self, add_keys=None) -> pd.DataFrame:
        """
        Dataframe with diffusion results.

        Args:
            add_keys: optional list of attributes to add.
        """
        keys = [
            "temperature", "latex_formula", "symbol",
            "diffusion_coeff", "err_diffusion_coeff",
            "conductivity", "err_conductivity",
            "engine",
        ]
        if add_keys is not None:
            keys.extend(list_strings(add_keys))

        d = {k: np.array([getattr(data, k) for data in self]) for k in keys}
        df = pd.DataFrame(d).sort_values(["latex_formula", "engine", "temperature"])
        return df

    @add_fig_kwargs
    def plot_sigma_berend(self, **kwargs) -> Figure:
        """
        Plot variance of correlated data as function of block number
        for all objects stored in DiffusionDataList.
        """
        nrows, ncols = len(self), 2
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        for ix, (data, ax_list) in enumerate(zip(self, ax_mat)):
            data.sigma_berend.plot(ax_list=ax_list, show=False)

        return fig

    @add_fig_kwargs
    def plot(self, **kwargs) -> Figure:
        """
        Plot ...
        for all objects stored DiffusionDataList.
        """
        nrows, ncols = len(self), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()
        if len(self) % ncols != 0: ax_list[-1].axis("off")

        for i, (data, ax) in enumerate(zip(self, ax_list)):
            data.plot(ax=ax, show=False)

        return fig

    def get_arrhenius_nmtuples(self, df: pd.Dataframe, hue: str | None) -> list:
        """
        Return list of namedtuple objects.

        Args:
            df: |pandas-DataFrame|.
            hue: Variable defining how to group data. If None, no grouping is performed.
        """
        df = df.sort_values(["temperature"])

        if hue is None:
            if not df["temperature"].is_unique:
                raise ValueError("Found duplicated values in temperature column. Please specify hue")

            return [AttrDict(temps=df["temperature"].values,
                             d_coeffs=df["diffusion_coeff"].values,
                             d_errs=df["err_diffusion_coeff"].values,
                             )]

        nt_list = []
        for key, grp in df.groupby(hue):
            grp = grp.sort_values(["temperature"])
            if not grp["temperature"].is_unique:
                raise ValueError(f"Found duplicated values in temperature column for {hue=}")
            nt = AttrDict(temps=grp["temperature"].values,
                          d_coeffs=grp["diffusion_coeff"].values,
                          d_errs=grp["err_diffusion_coeff"].values,
                          )
            nt_list.append(nt)
        return nt_list

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot_sigma_berend(show=False)
        yield self.plot(show=False)

    def expose(self, exposer="mpl", **kwargs):
        from abipy.tools.plotting import Exposer
        with Exposer.as_exposer(exposer) as e:
            e(self.yield_figs(**kwargs))


class MultiMdAnalyzer(HasPickleIO):
    """
    High-level interface to analyze multiple MD trajectories
    """

    @classmethod
    def from_abiml_dirs(cls, directories: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
        """
        Build an instance from a list of directories produced by abiml.py md
        """
        p = pool_nprocs_pmode(len(directories), pmode=pmode)
        using_msg = f"Reading {len(directories)} abiml directories {p.using_msg}"
        args = [(dirpath, step_skip) for dirpath in directories]
        with p.pool_cls(p.nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_abiml_dir, args))

    @classmethod
    def from_lammps_dirs(cls, directories: list, step_skip=1, basename="in.lammps", pmode="processes") -> MdAnalyzer:
        """
        Build an instance from a list of directories containing LAMMPS results.
        """
        p = pool_nprocs_pmode(len(directories), pmode=pmode)
        using_msg = f"Reading {len(directories)} LAMMPS directories {p.using_msg}"
        args = [(dirpath, step_skip, basename) for dirpath in directories]
        with p.pool_cls(p.nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_lammps_dir, args))

    @classmethod
    def from_vaspruns(cls, vasprun_filepaths: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
        """
        Build an instance from a list of vasprun files.
        """
        p = pool_nprocs_pmode(len(vasprun_filepaths), pmode=pmode)
        using_msg = f"Reading {len(vasprun_filepaths)} vasprun files {p.using_msg}..."
        args = [(vrun, step_skip) for vrun in vasprun_filepaths]
        with p.pool_cls(p.nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_vaspruns, args))

    @classmethod
    def from_hist_files(cls, hist_filepaths: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
        """
        Build an instance from a list of ABINIT HIST.nc files.
        """
        p = pool_nprocs_pmode(len(hist_filepaths), pmode=pmode)
        using_msg = f"Reading {len(hist_filepaths)} HIST.nc files {p.using_msg}..."
        args = [(ncpath, step_skip) for ncpath in hist_filepaths]
        with p.pool_cls(p.nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_hist_file, args))

    #@classmethod
    #def from_qe_inputs(cls, filepaths: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
    #    """
    #    Build an instance from a list of QE input files.
    #    """
    #    p = pool_nprocs_pmode(len(filepaths), pmode=pmode)
    #    using_msg = f"Reading {len(filespaths)} QE input files {p.using_msg}..."
    #    args = [(path, step_skip) for path in filepaths]
    #    with p.pool_cls(p.nprocs) as pool, Timer(header=using_msg, footer="") as timer:
    #        return cls(pool.starmap(MdAnalyzer.from_qe_input, args))

    def __init__(self, mdas: list[MdAnalyzer], temp_colormap="jet"):
        """
        Args:
            mdas: List of MdAnalyzer
            colormap: matplotlib colormap for temperatures.
        """
        self.mdas = mdas

        # Sort analyzers according to temperature.
        if self.has_same_system():
            self.mdas = sorted(mdas, key=lambda x: x.temperature)

        self.set_temp_colormap(temp_colormap)

    def __iter__(self):
        return self.mdas.__iter__()

    def __len__(self) -> int:
        return len(self.mdas)

    def __getitem__(self, items):
        return self.mdas.__getitem__(items)

    def has_same_system(self) -> bool:
        return all(mda.latex_formula == self[0].latex_formula for mda in self)

    def set_temp_colormap(self, colormap) -> None:
        """Set the colormap for the list of temperatures."""
        import matplotlib.pyplot as plt
        self.temp_cmap = plt.get_cmap(colormap)

    def get_params_dataframe(self) -> pd.DataFrame:
        """Dataframe with the parameters of the different MdAnalyzers."""
        return pd.DataFrame([mda.get_params_dict() for mda in self])

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose."""
        lines = []
        app = lines.append
        app(marquee("MD PARAMETERS", mark="="))
        app(self.get_params_dataframe().to_string())
        app("\n")
        return "\n".join(lines)

    def iter_mdat(self):
        """Iterate over (MdAnalyzer, temperature)."""
        for itemp, mda in enumerate(self):
            yield mda, mda.temperature

    def iter_mdatc(self):
        """Iterate over (MdAnalyzer, temperature, color)."""
        for itemp, mda in enumerate(self):
            yield mda, mda.temperature, self.temp_cmap(float(itemp) / len(self))

    def get_msdtt0_symbol_tmax(self, symbol: str, tmax: float, atom_inds=None, nprocs=None) -> Msdtt0List:
        msdtt0_list = Msdtt0List()
        for mda in self:
            obj = mda.get_msdtt0_symbol_tmax(symbol, tmax, atom_inds=atom_inds, nprocs=nprocs)
            msdtt0_list.append(obj)

        return msdtt0_list

    #def color_itemp(self, itemp: int):
    #    return self.temp_cmap(float(itemp) / len(self))

    #def temps_colors(self) -> tuple[list, list]:
    #    return ([mda.temperature for mda in self],
    #            [self.color_itemp(itemp) for itemp in range(len(self))])

    def _nrows_ncols_nplots(self, size=None):
        size = size or len(self)
        nrows, ncols, nplots = 1, 1, size
        if nplots > 1:
            ncols = 2; nrows = nplots // ncols + nplots % ncols
        return nrows, ncols, nplots

    @add_fig_kwargs
    def plot_sqdt_symbols(self, symbols, t0: float = 0.0,
                          xy_log=None, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement averaged over all atoms of the same specie vs time
        for the different temperatures.

        Args:
            symbols: List of chemical symbols to consider. "all" for all symbols in structure.
            t0: Initial time in ps.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
        """
        symbols = self[0]._select_symbols(symbols)
        nrows, ncols, nplots = self._nrows_ncols_nplots(size=len(symbols))
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()
        if nplots % ncols != 0: ax_list[-1].axis("off")

        # Plot data.
        for itemp, (mda, temp, color) in enumerate(self.iter_mdatc()):
            it0, ts = mda.get_it_ts(t0)
            for ix, (ax, symbol) in enumerate(zip(ax_list, symbols)):
                ax.plot(ts, mda.get_sqdt_symbol(symbol, it0=it0),
                        label=f"T = {temp} K", color=color)

        # Decorate axes.
        for ix, (ax, symbol) in enumerate(zip(ax_list, symbols)):
            ax.set_title(symbol, fontsize=fontsize)
            set_axlims(ax, xlims, "x")
            ax.legend(fontsize=fontsize, loc="upper left")
            ax.set_xlabel('t (ps)', fontsize=fontsize)
            ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
            #set_ticks_fontsize(ax, fontsize)
            set_logscale(ax, xy_log)

        return fig


def sfind_ind_val(array, value, arr_is_sorted=False) -> tuple:
    if arr_is_sorted:
        # Use Log(N) bisection.
        ind = find_le(array, value)
        return ind, array[ind]

    array = np.asarray(array)
    ind = (np.abs(array - value)).argmin()
    return ind, array[ind]


def linear_lsq_linefit(x, z, weights):
    S00 = np.sum(weights)
    S10 = np.sum(z * weights)
    S01 = np.sum(x * weights)
    S20 = np.sum(z**2 * weights)
    S11 = np.sum((x*z) * weights)
    D = S00 * S20 - S10**2
    m = (S00*S11 - S10*S01) / D
    q = (S01*S20 - S11*S10) / D
    varM = S00 / D
    varQ = S20 / D
    return m, q, varM, varQ


def _func(size_t0: int, it: int, pos_tac: np.ndarray, msd_tt0: np.ndarray):
    msd_tt0[it,:] = np.mean(np.sum((pos_tac[it:it+size_t0,:,:] - pos_tac[:size_t0,:,:])**2, axis=2), axis=1)


def msd_tt0_from_tac(pos_tac: np.ndarray, size_t: int, nprocs=None) -> np.ndarray:
    r"""
    Calculates the MSD for every possible pair of time points in the input array, using the formula:

        $$MSD(t,t_0) = \frac{1}{N} \sum_{i=1}^{N} (\vec{r}_i(t+t_0) - \vec{r}_i(t_0))^2$$

    where $N$ is the number of particles, and $\vec{r}_i(t)$ is the position vector.

    Args:
        pos_tac
        size_t
        nprocs:
    """
    # Check if size_t is valid.
    n_time_points = pos_tac.shape[0]
    if size_t >= n_time_points:
        raise ValueError(f"{size_t=} must be less than {n_time_points}")

    size_t0 = pos_tac.shape[0] - size_t
    msd_tt0 = np.empty((size_t, size_t0), dtype=float)

    if nprocs == 1:
        for it in range(0, size_t):
            #for it0 in range(0, size_t0):
            #    msd_tt0[it,it0] = np.mean(np.sum((pos_tac[it+it0,:,:] - pos_tac[it0,:,:])**2, axis=1))
            msd_tt0[it,:] = np.mean(np.sum((pos_tac[it:it+size_t0,:,:] - pos_tac[:size_t0,:,:])**2, axis=2), axis=1)
    else:
        p = pool_nprocs_pmode(nprocs, pmode="threads")
        using_msg = f"Computing MSD(t,t_0) matrix of shape {(size_t, size_t0)} {p.using_msg}"
        args = [(size_t0, it, pos_tac, msd_tt0) for it in range(0, size_t)]
        with p.pool_cls(p.nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            pool.starmap(_func, args)

    return msd_tt0


def block_mean_var(data, data_mean, n_block) -> tuple[float, float]:
    """
    Perform the block mean and the block variance of data.
    """
    N = data.shape[0]
    n_inblock = int(N/n_block)
    sigma2 = 0
    for iblock in range(n_block):
        mean_inblock = 0
        for datavalue in data[n_inblock*iblock:(iblock+1)*n_inblock]:
            mean_inblock = mean_inblock + datavalue/n_inblock
        sigma2 = sigma2 + (mean_inblock - data_mean)**2/(n_block)

    sigma2 = sigma2 / (n_block-1)
    delta_sigma2 = np.sqrt(2./(n_block-1)**3)*sigma2

    return sigma2, delta_sigma2


def sigma_berend(nblock_step: int, tot_block: int, data: np.ndarray) -> tuple[float, float, float]:
    """
    Args:
        nblock_step:
        tot_block:
        data:

    Return: (data_in_block, sigma, delta_sigma)
    """
    Ndata = data.shape[0]
    mean = np.mean(data)
    # tot_block=int(data.shape[0]/nblock_step)
    sigma2 = np.zeros(tot_block, dtype=float)
    delta_sigma2 = np.zeros(tot_block, dtype=float)
    arr_nblock = np.zeros(tot_block, dtype=float)
    data_in_block = np.zeros(tot_block, dtype=float)
    counter = 1
    sigma2[0] = float("inf")
    delta_sigma2[0] = 0
    for nblock in range(1, nblock_step * tot_block, nblock_step):
        sigma2[counter], delta_sigma2[counter] = block_mean_var(data, mean, nblock+1)
        arr_nblock[counter] = nblock + 1
        data_in_block[counter] = int(Ndata/(nblock+1))
        counter = counter + 1

    sigma = np.sqrt(sigma2)
    delta_sigma = 0.5 * delta_sigma2 / sigma

    return data_in_block, sigma, delta_sigma


def _lin_fit(th_invt, log10d0, e_act):
    return log10d0 - e_act * th_invt


@dataclasses.dataclass(kw_only=True)
class ArrheniusEntry:
    """
    """
    key: str
    symbol: str
    composition: str
    temps: np.ndarray
    diffusions: np.ndarray
    err_diffusions: np.ndarray
    volumes: np.ndarray
    mpl_style: dict

    @classmethod
    def from_file(cls, filepath: PathLike, key, mpl_style) -> ArrheniusEntry:

        # Read data in CSV format. Assuming header with at least the following entries:
        #   temperature,diffusion,err_diffusion,volume,symbol,composition
        try:
            df = pd.read_csv(filepath, skipinitialspace=True)

            def get_unique(col):
                v0 = df[col].values[0]
                if np.any(v0 != df[col].values):
                    raise ValueError(f"All values for column: {col} should be unique while found:\n{df[col].values}")
                return v0

            symbol = get_unique("symbol")
            composition = Composition(get_unique("composition"))
            temps = df["temperature"].values
            volumes = df["volume"].values
            diffusions = df["diffusion"].values
            err_diffusions = np.zeros(len(temps))
            if "err_diffusion" in df.keys():
                err_diffusions = df["err_diffusion"].values

        except Exception as exc:
            raise RuntimeError(f"Exception while reading {filepath=}") from exc

        return cls(key=key,
                   symbol=symbol,
                   composition=composition,
                   temps=temps,
                   diffusions=diffusions,
                   err_diffusions=err_diffusions,
                   volumes=volumes,
                   mpl_style=mpl_style,
                   )

    #def __post_init__(self):
    #    self.latex_formula = latexify(self.formula)

    #@property
    #def latex_formula(self) -> str:
    #    return self._latex_formula

    #@latex_formula.setter
    #def latex_formula(self, value):
    #    """LaTeX formatted formula. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$."""
    #    self._latex_formula = latexify(value)

    def get_diffusion_data(self, fit_thinvt=None) -> AttrDict:
        """
        Fit diffusion(T) taking into account uncertainties.
        Return dict with activation energy in eV and fit parameters.
        """
        # NB: log10(x) = ln(x) log10(e) --> d_x log_10(x) = log10(e) / x
        th_invt = 1000 / self.temps
        log10 = np.log10(self.diffusions)
        err_log10 = np.log10(np.e) * self.err_diffusions / self.diffusions
        popt, pcov = optimize.curve_fit(_lin_fit, th_invt, log10, sigma=err_log10)
        e_act = popt[1] * 1000 * kBoltzEv * np.log(10)

        fit_log10 = None
        if fit_thinvt is not None:
            fit_log10 = _lin_fit(fit_thinvt, popt[0], popt[1])

        return AttrDict(
            th_invt=th_invt,
            log10=log10, err_log10=err_log10,
            fit_thinvt=fit_thinvt, fit_log10=fit_log10,
            e_act=e_act, popt=popt, pcov=pcov,
        )

    def get_conductivity_data(self, ncar: float, fit_thinvt=None) -> tuple[AttrDict, AttrDict]:
        """
        Compute conductivity sigma and T x sigma(T) assuming ncar carriers.
        Fit values taking into account uncertainties.

            sigma(T) =  Nq^2 D(T) / (VkT)
        """
        temps = self.temps
        diffusions = self.diffusions
        err_diffusions = self.err_diffusions
        volumes = self.volumes

        th_invt = 1000 / self.temps

        #if charge is None:
        #ncar = symb2oxi[symbol] * self.composition[self.symbol]
        #ncar = charge * self.composition[self.symbol]

        # NB: log10(x) = ln(x) log10(e) --> d_x log_10(x) = log10(e) / x
        conds = e2s/kbs * ncar * diffusions/volumes/temps * 1.e09
        err_conds = e2s/kbs * ncar * err_diffusions/volumes/temps * 1.e09
        log10_conds = np.log10(conds)
        err_log10_conds = np.log10(np.e) * err_conds / conds
        cond_popt, cond_pcov = optimize.curve_fit(_lin_fit, th_invt, log10_conds, sigma=err_log10_conds)

        fit_log10_cond = None
        if fit_thinvt is not None:
            fit_log10_cond = _lin_fit(fit_thinvt, cond_popt[0], cond_popt[1])

        # temp(K) * sigma(S/cm)
        t_conds = e2s/kbs * ncar * diffusions/volumes*1.e09
        err_tconds = e2s/kbs * ncar * err_diffusions/volumes*1.e09
        log10_tconds = np.log10(t_conds)
        err_log10_tconds = np.log10(np.e) * err_tconds / t_conds
        tcond_popt, tcond_pcov = optimize.curve_fit(_lin_fit, th_invt, log10_tconds, sigma=err_log10_tconds)

        fit_log10_tcond = None
        if fit_thinvt is not None:
            fit_log10_tcond = _lin_fit(fit_thinvt, tcond_popt[0], tcond_popt[1])

        return (
           AttrDict(
               th_invt=th_invt,
               log10=log10_conds, err_log10=err_log10_conds,
               fit_thinvt=fit_thinvt, fit_log10=fit_log10_cond,
               popt=cond_popt, cov=cond_pcov,
           ),
           AttrDict(
               th_invt=th_invt,
               log10=log10_tconds, err_log10=err_log10_tconds,
               fit_thinvt=fit_thinvt, fit_log10=fit_log10_tcond,
               popt=tcond_popt, cov=tcond_pcov,
           )
        )


class ArrheniusPlotter:
    """
    This object stores conductivities D(T) computed for different structures and/or
    different ML-potentials and allows one to produce Arrnehnius plots on the same figure.
    Internally, the results are indexed by a unique key that is be used as label in the matplotlib plot.
    The style for each key can be customized by setting the mpl_style dict
    with the options that will be passed to ax.plot.

    In the simplest case, one reads the data from external files in CSV format.

    Example:

        from abipy.dynamics.analyzer import ArrheniusPlotter

        key_path = {
            "matgl-MD":  "diffusion_cLLZO-matgl.csv",
            "m3gnet-MD": "diffusion_cLLZO-m3gnet.csv",
        }

        mpl_style_key = {
            "matgl-MD" : dict(c='blue'),
            "m3gnet-MD": dict(c='purple'),
        }

        plotter = ArrheniusPlotter()
        for key, path in key_path.items():
            plotter.add_entry_from_file(path, key, mpl_style=mpl_style_key[key])

        # temperature grid refined
        thinvt_arange = (0.6, 2.5, 0.01)
        xlims, ylims = (0.5, 2.0), (-7, -4)

        plotter.plot(thinvt_arange=thinvt_arange,
                     xlims=xlims, ylims=ylims, text='LLZO cubic', savefig=None)
    """

    def __init__(self, entries=None):
        self.entries = entries or []
        self.symb2oxi = common_oxidation_states()

    def __iter__(self):
        return self.entries.__iter__()

    def __len__(self) -> int:
        return len(self.entries)

    def __getitem__(self, items):
        return self.entries.__getitem__(items)

    def copy(self):
        """Deep copy of the object."""
        import copy
        return copy.deepcopy(self)

    def keys(self) -> list[str]:
        """List of keys. Must be unique"""
        return [entry.key for entry in self]

    def symbols(self) -> set[str]:
        """Return set with chemical symbols."""
        return set([entry.symbol for entry in self])

    def index_key(self, key: str) -> int:
        """Find the index of key. Raise KeyError if not found."""
        for i, entry in enumerate(self):
            if entry.key == key: return i
        raise KeyError(f"Cannot find {key=} in {self.keys=}")

    def pop_key(self, key: str) -> ArrheniusEntry:
        """Pop entry by key"""
        i = self.index_key(key)
        return self.entries.pop(i)

    def append(self, entry: ArrheniusEntry) -> None:
        """Append new entry."""
        if entry.key in self.keys():
            raise KeyError(f"{entry.key=} is already present.")
        self.entries.append(entry)

    def set_style(self, key: str, mpl_style: dict) -> None:
        """Set matplotlib style for key."""
        i = self.index_key(key)
        self.entries[i].mpl_style = mpl_style

    def get_min_max_temp(self) -> tuple[float, float]:
        """Compute the min and max temperature for all entries."""
        min_temp = min([e.temperatures.min() for e in self])
        max_temp = max([e.temperatures.max() for e in self])
        return min_temp, max_temp

    def add_entry_from_file(self, filepath: PathLike, key: str, mpl_style=None) -> None:
        """
        """
        mpl_style = mpl_style or {}
        self.append(ArrheniusEntry.from_file(filepath, key, mpl_style))

    @add_fig_kwargs
    def plot(self, thinvt_arange=None, what="diffusion", ncar=None, colormap="jet", with_t=True, text=None,
             ax=None, fontsize=8, xlims=None, ylims=None, **kwargs) -> Figure:
        """
        Arrhenius plot.

        Args:
            thinvt_arange: start, stop, step for 1000/T mesh. If None, the mesh is automatically computed.
            what: Selects the quantity to plot. Possibile values: "diffusion", "sigma", "tsigma".
            ncar: Number of carriers. Required if what is "sigma" or "tsigma".
            colormap: Colormap used to select the color if entry.mpl_style does not provide it.
            with_t: True to dd a twin axes with the value of T
            text:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: Similar to xlims but for the y-axis.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        # Build temperature grid for fit.
        if thinvt_arange is None:
            t_min, t_max = self.get_min_max_temp()
            xstart, xstop = 1000 / t_max, 1000 / t_min
            xstart, xstop = 0.8 * xstart, 1.2 * xstop
            thinvt_arange = (xstart, xstop, 0.01)

        fit_thinvt = np.arange(thinvt_arange[0], thinvt_arange[1], step=thinvt_arange[2])

        for ie, entry in enumerate(self):
            mpl_style = entry.mpl_style.copy()
            if "marker" not in mpl_style:
                mpl_style["marker"] = "."
            if "linestyle" not in mpl_style and "ls" not in mpl_style:
                mpl_style["linestyle"] = ""
            if "color" not in mpl_style and "c" not in mpl_style:
                mpl_style["color"] = cmap(ie / len(self))

            symbol, composition = entry.symbol, entry.composition

            diff = entry.get_diffusion_data(fit_thinvt=fit_thinvt)
            label = r"D$_{\mathrm{%s}}$ (%s): " % (symbol, entry.key)
            label += r"E$_\mathrm{a}$=" + str('{:.2F}'.format(diff.e_act)) + ' eV'

            if what == "diffusion":
                data = diff
            else:
                ncar_ = ncar
                if ncar_ is None:
                    if symbol not in self.symb2oxi:
                        raise ValueError(f"No entry for {symbol=} found in symb2oxi! Please add the oxistate manually!")
                    charge = self.symb2oxi[symbol]
                    ncar_ = charge * composition[symbol]
                    label += f", q={charge}"

                sigma_data, tsigma_data = entry.get_conductivity_data(ncar_, fit_thinvt=fit_thinvt)
                data = dict(sigma=sigma_data, tsigma=tsigma_data)[what]

            if data.err_log10 is not None:
                # Plot data with errors.
                lines = ax.errorbar(data.th_invt, data.log10, yerr=data.err_log10, label=label, capsize=5.0, **mpl_style)
            else:
                # Plot data without errors.
                lines = ax.plot(data.th_invt, data.log10, label=label, **mpl_style)

            ax.plot(fit_thinvt, data.fit_log10, linestyle='--', color=lines[0].get_color())

        ax.set_xlabel('1000 / T(K)', fontsize=18)
        ylabel = dict(
            diffusion=r"log$_{10}$ (D(cm$^2$/s))",
            sigma=r"log$_{10}$ [$\sigma$(S/cm)]",
            tsigma=r"log$_{10}$ [$\sigma$T(SK/cm)]",
        )[what]
        ax.set_ylabel(ylabel, fontsize=18)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        ax.legend(loc="lower left", fontsize=12)
        #set_ticks_fontsize(ax, fontsize=14)
        #from matplotlib.ticker import MultipleLocator
        #ax.yaxis.set_major_locator(MultipleLocator(1))
        #ax.yaxis.set_minor_locator(MultipleLocator(0.2))

        if text:
            ax.text(0.96, 0.85, text,
                    verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes,
                    color='black', fontsize=18, bbox=dict(facecolor='none', edgecolor='grey', pad=10))

        if with_t:
            # Add a twin axes and set its limits so it matches the first.
            ax_t = ax.twiny()
            ax_t.set_xlabel('T (K)', fontsize=18)
            ax_t.set_xlim(ax.get_xlim())
            # apply a function formatter
            import matplotlib.ticker as mticker
            formatter = mticker.FuncFormatter(lambda x, pos: '{:.0f}'.format(1000/x))
            ax_t.xaxis.set_major_formatter(formatter)

        return fig

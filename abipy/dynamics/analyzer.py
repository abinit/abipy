"""
Tools to analyze MD trajectories and compute diffusion coefficients.
"""
from __future__ import annotations

import dataclasses
import warnings
import json
import numpy as np
import pandas as pd
import pymatgen.core.units as units

from pathlib import Path
from matplotlib.offsetbox import AnchoredText
from monty.functools import lazy_property
from monty.bisect import find_le
from monty.string import list_strings, is_string, marquee
from monty.collections import AttrDict, dict2namedtuple
from pymatgen.util.string import latexify
from pymatgen.core.lattice import Lattice
from abipy.core.mixins import TextFile, NotebookWriter
from abipy.core.structure import Structure
from abipy.tools.typing import Figure
from abipy.tools.serialization import HasPickleIO
from abipy.tools.iotools import try_files
from abipy.tools.context_managers import Timer
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, get_color_symbol,
                                  set_ticks_fontsize, set_logscale, linear_fit_ax)


__author__ = "Giuliana Materzanini, Matteo Giantomassi"


Ang2PsTocm2S = 0.0001
e2s = 1.602188**2 # electron charge in Coulomb scaled by 10.d-19**2
kbs = 1.38066     # Boltzmann constant in Joule/K scaled by 10.d-23


def read_structure_postac_ucmats(traj_filepath: str, step_skip: int) -> tuple[Structure, np.ndarray, np.ndarray]:
    """
    Read all configurations from an ASE trajectory file.
    Returns (nsteps, natom, 3) array with the Cartesian coords and the initial structure.

    Args
        step_skip: Sampling frequency.
                time_step is multiplied by this number to get the real time between measurements.
    """
    from ase.io import read
    traj = read(traj_filepath, index=":")
    nsteps = len(traj)
    structure = Structure.as_structure(traj[0])

    pos_tac, ucmats = [], []
    for it in range(0, nsteps, step_skip):
        atoms = traj[it]
        pos_tac.append(atoms.positions)
        ucmats.append(atoms.cell.array)

    del traj
    pos_tac = np.array(pos_tac, dtype=float)
    ucmats = np.array(ucmats, dtype=float)

    return structure, pos_tac, ucmats


class MdAnalyzer(HasPickleIO):
    """
    High-level interface to read MD trajectories and metadata from external files,
    compute the MSQD and visualize the results.
    """

    @classmethod
    def from_abiml_dir(cls, directory, step_skip=1) -> MdAnalyzer:
        """
        Build an instance from a directory containing an ASE trajectory file and
        a JSON file with the MD parameters as produced by `abiml.py md`.
        """
        directory = Path(str(directory))
        structure, pos_tac, ucmats = read_structure_postac_ucmats(directory / "md.traj", step_skip)

        # Read metadata from the JSON file.
        with open(directory / "md.json", "rt") as fh:
            meta = json.load(fh)

        temperature = meta["temperature"]
        # Convert from fs to ps
        timestep = meta["timestep"] * 1e-3
        loginterval = meta["loginterval"]
        engine = meta["nn_name"]
        nsteps = len(pos_tac)
        times = (np.arange(0, nsteps) * timestep * loginterval)[::step_skip].copy()

        log_path = try_files([directory / "md.aselog", directory / "md.log"])
        from abipy.ml.aseml import AseMdLog
        with AseMdLog(log_path) as log:
            evp_df = log.df.copy()

        return cls(structure, temperature, times, pos_tac, ucmats, engine, evp_df=evp_df)

    @classmethod
    def from_hist_file(cls, hist_filepath: str, step_skip=1) -> MdAnalyzer:
        """
        Build an instance from a list of ABINIT HIST.nc files.
        """
        raise NotImplementedError()
        #from abipy.dynamics.hist import HistFile
        #with HistFile(hist_filepath) as hist:
        #    pos_tac = hist.r.read_value("xcart") * units.bohr_to_ang
        #times = (np.arange(0, nsteps) * timestep * loginterval)[::step_skip].copy()
        #temperature = None
        #evp_df = None
        #return cls(hist.structure.copy(), temperature, times, pos_tac, ucmats, "abinit", evp_df=evp_df)

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
    #def from_qe_dir(cls, directory: str, step_skip=1):

    #@classmethod
    #def from_qe_input(cls, filepath: str, step_skip=1):
    #    """
    #    Build an instance from a directory with CP results.
    #    """
    #    # Get structure from QE input.
    #    from pymatgen.io.pwscf import PWInput
    #    qe_inp = PWInput(filepath)
    #    times = None
    #    temperature = None
    #    # Get atomic positions from qe pos file.
    #    pos_tac = None
    #    engine = "qecp"
    #    engine = "qepw"
    #    return cls(qe_inp.structure, temperature, times, pos_tac, ucmats, engine)

    @classmethod
    def from_lammps_dir(cls, directory: str, step_skip=1) -> MdAnalyzer:
        """
        Build an instance from a directory containing LAMMPS results.
        """
        traj_filepath, log_filepath = None, None
        directory = Path(str(directory))
        for path in directory.listdir():
            if path.is_dir(): continue
            if path.suffix == ".lammpstrj":
                traj_filepath = path.absolute()

        if traj_filepath is None:
            raise ValueError(f"Cannot find file with .lammpstrj extension in {directory=})")
        if log_filepath is None:
            raise ValueError(f"Cannot find log file in {directory=}")

        return cls.from_lammpstrj(traj_filepath, log_filepath, step_skip=step_skip)

    @classmethod
    def from_lammpstrj(cls, traj_filepath: str, log_filepath: str, step_skip=1) -> MdAnalyzer:
        """
        Build an instance from a LAMMPS trajectory file and log file.
        """
        structure, pos_tac, ucmats = read_structure_postac_ucmats(traj_filepath, step_skip)
        temperature = None

        if log_filepath.endswith(".evp"):
            # Extract times from CP EVP file (Giuliana's way)
            from abipy.dynamics.cpx import EvpFile
            with EvpFile(log_filepath) as evp:
                times = evp.times.copy()
                evp_df = evp.df.copy()
        else:
            #from pymatgen.io.lammps.outputs import parse_lammps_log
            #def parse_lammps_log(filename: str = "log.lammps") -> list[pd.DataFrame]:
            evp_df = None
            raise NotImplementedError()

        return cls(structure, temperature, times, pos_tac, ucmats, "lammps", evp_df=evp_df)

    def __init__(self,
                structure: Structure,
                temperature: float,
                times: np.ndarray,
                cart_positions: np.ndarray,
                ucmats: np.ndarray,
                engine: str,
                pos_order: str="tac",
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

    def consistency_check(self):
        """Perform internal consistency check."""
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
        """Dictionary with important parameters."""
        attr_names = [
            "latex_formula", "temperature", "timestep", "nt", "max_time", "natom", "avg_volume", "engine",
        ]
        d = {aname: getattr(self, aname) for aname in attr_names}
        return d

    def deepcopy(self) -> MdAnalyzer:
        """Deep copy of the object."""
        import copy
        return copy.deepcopy(self)

    def resample_nsteps(self, start_nsteps: int, tmesh_nsteps) -> MdAnalyzer:
        return self.resample_time(start_time=start_nsteps * self.timestep, new_timestep=tmesh_nsteps * self.timestep)

    def resample_time(self, start_time: float, new_timestep: float) -> MdAnalyzer:
        """
        """
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

        if new_timestep < old_timestep:
            raise ValueError(f"Invalid {new_timestep=} should be >= {old_timestep}")

        istep = int(new_timestep / old_timestep)
        if istep != 1:
            new.pos_atc = new.pos_atc[:,::istep,:].copy()
            new.times = new.times[::istep] - new.times[0]
            if new.lattices is not None:
                new.lattices = new.lattices[::istep].copy()

        new.consistency_check()
        return new

    @property
    def timestep(self) -> float:
        """Timestep in ps."""
        return self.times[1] - self.times[0]

    @property
    def max_time(self) -> float:
        """Max. simulation time in ps."""
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

        return sqdt / (count + 1)

    def get_msq_tt0_symbol(self, symbol: str, tmax: float, atom_inds=None, nprocs=None) -> Msdtt0:
        """
        """
        index_tmax, _ = self.get_it_ts(tmax)
        iatoms = self.iatoms_with_symbol(symbol, atom_inds=atom_inds)
        tac = self.pos_atc[iatoms].transpose(1, 0, 2).copy()
        msd_tt0 = msd_tt0_from_tac(tac, index_tmax, nprocs=nprocs)

        return Msdtt0(msd_tt0=msd_tt0, mda=self, index_tmax=index_tmax, symbol=symbol)

    #def export_msdt(self, filename: str) -> None:
    #    """
    #    Writes MSD data to a file that can be easily plotted in other software.

    #    Args:
    #        filename: Supported formats are csv and dat.
    #            If the extension is csv, a csv file is written. Otherwise, a dat format is assumed.
    #    """
    #    fmt = "csv" if filename.lower().endswith(".csv") else "dat"
    #    delimiter = ", " if fmt == "csv" else " "
    #    with open(filename, "w") as f:
    #        if fmt == "dat": f.write("# ")
    #        f.write(delimiter.join(["t", "MSD", "MSD_a", "MSD_b", "MSD_c", "MSCD"]))
    #        f.write("\n")
    #        for dt, msd, msdc, mscd in zip(
    #            self.dt, self.msd, self.msd_components, self.mscd
    #        ):
    #            f.write(delimiter.join([str(v) for v in [dt, msd, *list(msdc), mscd]]))
    #            f.write("\n")

    #@add_fig_kwargs
    #def plot_evp(self, **kwargs) -> Figure:
    #    if self.evp_df is None:
    #        print("Cannot plot evp data as self.evp_df is None!")
    #        return None

    #    if self.engine == "lampps"
    #    #return self.evp.plot(**kwargs)

    @add_fig_kwargs
    def plot_sqdt_atoms(self, symbols="all", t0: float = 0.0, atom_inds=None,
                        ax=None, xy_log=None, fontsize=30, xlims=None, **kwargs) -> Figure:
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
                   or scalar e.g. ``left``. If left (right) is None, default values are used
        """
        it0, ts = self.get_it_ts(t0)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            for iatom in self.iatoms_with_symbol(symbol, atom_inds=atom_inds):
                sd_t = self.get_sqdt_iatom(iatom, it0=it0)
                ax.plot(ts, sd_t, label=f"${symbol}_{{{iatom}}}$")

        ax.legend(fontsize=5, loc=2)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel(r'square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.latex_formula_n_temp}\n{self.latex_avg_volume}\n" +
                                   'sd(t, $t_0$ =' + str(int(self.times[it0])) + ' ps)',
                                   loc=1, prop=dict(size=20)))
        return fig

    @add_fig_kwargs
    def plot_sqdt_symbols(self, symbols, t0: float = 0.0, atom_inds=None,
                          ax=None, xy_log=None, fontsize=30, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement averaged over all atoms of the same specie vs time.

        Args:
            symbols: List of chemical symbols to consider. "all" for all symbols in structure.
            t0: Initial time in ps.
            atom_inds: List of atom indices to include. None to disable filtering.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
        """
        it0, ts = self.get_it_ts(t0)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            ax.plot(ts, self.get_sqdt_symbol(symbol, it0=it0, atom_inds=atom_inds),
                    label=symbol + ' msd($t, t_0$ = ' + str(t0) + ' ps)',
                    color=self.color_symbol[symbol],
                    )

        ax.legend(fontsize=16, loc=2)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel('mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.latex_formula_n_temp}\n{self.latex_avg_volume}",
                                   loc=1, prop=dict(size=20)))
        return fig

    @add_fig_kwargs
    def plot_sqdt_symbols_tmax(self, symbols, tmax: float, atom_inds=None, nprocs=None,
                               ax=None, xy_log=None, fontsize=30, xlims=None, **kwargs) -> Figure:
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
                   or scalar e.g. ``left``. If left (right) is None, default values are used
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
                    label=symbol + ' <msd($t$, $t_0$)>$\{$t_0$\}$, $t$ = [0, ' + str(int(self.times[index_tmax])) + ' ps]',
                    color=self.color_symbol[symbol],
                    )

        set_axlims(ax, xlims, "x")
        ax.legend(fontsize=16, loc=2)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.latex_formula_n_temp}\n{self.latex_avg_volume}",
                                   loc=1, prop=dict(size=20)))
        return fig

    @add_fig_kwargs
    def plot_lattices(self, what_list=("abc", "angles", "volume"), ax_list=None,
                      xy_log=None, fontsize=8, xlims=None, **kwargs) -> Figure:
        """
        Plot lattice lengths/angles/volume as a function of time.
        """
        if self.lattices is None:
            print("MD simulation has been done with fixed lattice.")
            return None

        what_list = list_strings(what_list)
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=1, ncols=len(what_list),
                                                sharex=True, sharey=False, squeeze=False)
        markers = ["o", "^", "v"]

        if "abc" in what_list:
            # plot lattice parameters.
            for i, label in enumerate(["a", "b", "c"]):
                ax.plot(self.times, [lattice.abc[i] for lattice in self.lattices],
                        label=label, marker=markers[i])
            ax.set_ylabel("abc (A)")

        if "abc" in what_list:
            # plot lattice angles.
            for i, label in enumerate(["alpha", "beta", "gamma"]):
                ax.plot(self.times, [lattice.angles[i] for lattice in self.lattices],
                        label=label, marker=markers[i])
            ax.set_ylabel(r"$\alpha\beta\gamma$ (degree)")

        if "volume" in what_list:
            # plot lattice volume.
            marker = "o"
            ax.plot(self.times, [lattice.volume for lattice in self.lattices],
                    lable="Volume", marker=marker)
            ax.set_ylabel(r'$V\, (A^3)$')

        for ix, ax in enumerate(ax_list):
            set_logscale(ax, xy_log)
            set_axlims(ax, xlims, "x")
            if ix == len(ax_list) -1:
                ax.set_xlabel('t (ps)', fontsize=fontsize)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig


@dataclasses.dataclass(kw_only=True)
class Msdtt0:
    """
    """
    index_tmax: int
    symbol: str
    msd_tt0: np.ndarray
    mda: MdAnalyzer

    @property
    def times(self):
        return self.mda.times

    @property
    def temperature(self):
        return self.mda.temperature

    @add_fig_kwargs
    def plot(self, ax=None, xy_log=None, fontsize=30, xlims=None, **kwargs) -> Figure:
        """
        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
        """
        msd_t = np.mean(self.msd_tt0, axis=1)
        index_tmax = self.index_tmax
        t_start = self.mda.nt - index_tmax
        ts = self.times[t_start:] - self.times[t_start]

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ax.plot(ts, msd_t,
                label=self.symbol + ' <msd($t, t_0$)>$\{$t_0$\}$, t = [0, ' + str(int(self.times[index_tmax]))+ ' ps]',
                color=self.mda.color_symbol[self.symbol],
                )

        from scipy.stats import linregress
        fit = linregress(ts, msd_t)
        naive_d = fit.slope * Ang2PsTocm2S / 6
        label = r"Linear fit $\alpha={:.2f}$, $r^2$={:.2f}".format(fit.slope, fit.rvalue**2)
        ax.plot(ts, fit.slope*ts + fit.intercept, label=label)

        set_axlims(ax, xlims, "x")
        ax.legend(fontsize=16, loc=2)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(f"{self.mda.latex_formula_n_temp}\n{self.mda.latex_avg_volume}",
                                   loc=1, prop=dict(size=20)))
        return fig

    def get_sigma_berend(self, t1, t2, nblock_step=1, tot_block=1000) -> SigmaBerend:
        """
        """
        # choose the time elapsed
        it1, _ = sfind_ind_val(self.times, t1)
        it2, _ = sfind_ind_val(self.times, t2)
        size_t = self.msd_tt0.shape[0]

        if it1 >= it2:
            raise ValueError(f"For input {t1=} and {t2=}, got {it1=} >= {it2=}")
        if it1 >= size_t:
            raise ValueError(f"For input {t1=}, got {it1=} >= {size_t=}")
        if it2 >= size_t:
            raise ValueError(f"For input {t2=}, got {it2=} >= {size_t=}")

        block_sizes1, sigmas1, delta_sigmas1 = sigma_berend(nblock_step, tot_block, self.msd_tt0[it1,:])
        block_sizes2, sigmas2, delta_sigmas2 = sigma_berend(nblock_step, tot_block, self.msd_tt0[it2,:])

        # Build SigmaBerend instance from locals dict.
        mda = self.mda
        latex_formula = mda.latex_formula
        temperature = mda.temperature
        time1, time2 = mda.times[it1], mda.times[it2]
        data = locals()
        return SigmaBerend(**{k: data[k] for k in [field.name for field in dataclasses.fields(SigmaBerend)]})

    def get_diffusion_with_sigma(self,
                                 fit_time_start: float, fit_time_stop: float,
                                 block_size1: int, block_size2: int,
                                 sigma_berend: SigmaBerend,
                                 ax=None) -> DiffusionData:
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

        # and find error for anytime
        msd_tt0 = self.msd_tt0
        err_msd = np.zeros(msd_tt0.shape[0], dtype=float)
        for t in range(msd_tt0.shape[0]):
            err_msd[t] = abs(mSigma * t + qSigma)

        # and find error for anytime
        dataScorrelated = np.zeros(msd_tt0.shape[0], dtype=int)
        for t in range(msd_tt0.shape[0]):
            dataScorrelated[t] = int(mDataInBlock * t + qDataInBlock)

        # average over the initial times
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
        print('{:.2E}'.format(best_angcoeff*Ang2PsTocm2S/6,2))
        print(best_angcoeff)
        print(max_angcoeff)
        print(min_angcoeff)

        # Build DiffusionData instance from locals dict.
        data = locals()
        return DiffusionData(**{k: data[k] for k in [field.name for field in dataclasses.fields(DiffusionData)]})

        """
        outMSD_file.write('timeStepJump = ' + str(timeStepJump) + '\n')
        outMSD_file.write('time, cel and pos were cut before ' + str(timeArray[0]+tInitial)+ 'ps' + '\n')
        outMSD_file.write('t elapsed max is ' + str(estimatedTMax)+ 'ps' + '\n')
        outMSD_file.write('trajectory length is ' + str(timeArrayTmp[timeArrayTmp.shape[0]-1])+ 'ps' + '\n')
        outMSD_file.write('error on msd(t) evaluated at ' +  str(estimatedFirstTElapsed) + 'ps' + ' and ' + str(estimatedSecondTElapsed) + 'ps' + '\n')
        outMSD_file.write('evaluated n. of blocks at ' + str(estimatedFirstTElapsed) + 'ps' + ' is ' + str(block_size1)  + '\n')
        outMSD_file.write('evaluated n. of blocks at ' + str(estimatedSecondTElapsed) + 'ps' + ' is ' + str(block_size2)  + '\n')
        outMSD_file.write('number of decorrelated msd(t) data that we fit: ' + str(timeArrScorrelated.shape[0]) + '\n')
        outDiff_file.write(str(temperature) + ' ' + str(diffusion_coeff) + ' ' + str(err_diffusion_coeff) + ' ' +str(volumeArrayTmp[0])+ '\n')
        """

    @add_fig_kwargs
    def plot_mat(self, cmap="jet", fontsize=8, ax=None, **kwargs) -> Figure:
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        im = ax.matshow(self.msd_tt0, cmap=cmap)
        fig.colorbar(im, ax=ax)
        ax.set_xlabel('t0 (ps)', fontsize=fontsize)
        ax.set_ylabel('t (ps)', fontsize=fontsize)
        ax.set_title(self.symbol + ": msd($t, t_0$)", fontsize=fontsize)
        return fig



class Msdtt0List(list):
    """A list of Msdtt0 objects."""

    @add_fig_kwargs
    def plot(self, sharex=True, sharey=True, **kwargs) -> Figure:
        nrows, ncols = len(self), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=sharex, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()
        if len(self) % ncols != 0: ax_list[-1].axis("off")

        for ix, (obj, ax) in enumerate(zip(self, ax_list)):
            obj.plot(ax=ax, show=False)
        return fig



@dataclasses.dataclass(kw_only=True)
class SigmaBerend:
    """
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
    def plot(self, **kwargs) -> Figure:
        """
        Plot variance of correlated data as function of block number.
        """
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=1, ncols=2,
                                                sharex=True, sharey=True, squeeze=False)

        for ix, ax in enumerate(ax_list.ravel()):
            xs, ys, yerr, it, time = self.block_sizes1, self.sigmas1, self.delta_sigmas1, self.it1, self.time1
            if ix == 1:
                xs, ys, yerr, it, time = self.block_sizes2, self.sigmas2, self.delta_sigmas2, self.it2, self.time2

            ax.errorbar(xs, ys,
                        yerr=yerr, linestyle='-', linewidth=0.5,
                        label="$\sigma(\mathrm{MSD}($" + '%2.1f' % time +" ps$))$ "+ '\n' +
                               self.latex_formula + ', '+ 'T = %4.0f' % self.temperature + 'K')
            ax.legend(fontsize=16, loc=4)
            ax.set_xlabel('N. of data in block', fontsize=14)
            ax.set_ylabel('$\sigma$ ($\AA^2$)', fontsize=14)
            ax.grid(True)

        fig.suptitle("Variance of correlated data as function of block number")
        return fig


@dataclasses.dataclass(kw_only=True)
class DiffusionData:
    """
    """
    diffusion_coeff: float
    err_diffusion_coeff: float
    conductivity: float
    err_conductivity: float
    temperature: float
    symbol: str
    latex_formula: str
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

    #@classmethod
    #def from_file(cls, filepath):
    #    return cls

    #def json_dump(cls, filepath):

    @add_fig_kwargs
    def plot(self, ax=None, **kwargs) -> Figure:
        """
        Plot ...
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ts = self.times[:self.msd_t.shape[0]]

        ax.errorbar(ts, self.msd_t, yerr=self.err_msd, color='mediumblue', label=self.symbol)
        ax.errorbar(self.timeArrScorrelated, self.msdSScorrelated, yerr=self.errMSDScorrelated, linestyle='-')
        ax.errorbar(ts, self.best_angcoeff * ts + self.quote, linestyle='--')
        ax.errorbar(ts, self.min_angcoeff * ts + self.quote, linestyle='--')
        ax.errorbar(ts, self.max_angcoeff * ts + self.quote, linestyle='--')

        ax.set_xlabel('t (ps)', fontsize=18)
        ax.set_ylabel(r'$\mathrm{MSD}_\mathrm{tr}$ $\mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
        ax.add_artist(AnchoredText(
            'D$_{tr}$ = (' + str('{:.2E}'.format(self.diffusion_coeff)) + '\u00B1' +
            str('{:.2E}'.format(self.err_diffusion_coeff)) + ') cm$^2$/s',
            loc=2, prop=dict(size=14)))
        ax.legend(fontsize=12, loc=4)
        ax.add_artist(AnchoredText(f"{self.latex_formula}\nT = {self.temperature} K",
                                   loc=1, prop=dict(size=14)))
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

    #def append(self, data: DiffusionData) -> None
    #    super().append(data)

    #def _nrows_ncols_nplots(self, size=None):
    #    size = size or len(self)
    #    nrows, ncols, nplots = 1, 1, size
    #    if nplots > 1:
    #        ncols = 2; nrows = nplots // ncols + nplots % ncols
    #    return nrows, ncols, nplots

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
        for all items stored DiffusionDataList.
        """
        nrows, ncols = len(self), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()
        if len(self) % ncols != 0: ax_list[-1].axis("off")
        for ix, (data, ax) in enumerate(zip(self, ax_list)):
            data.sigma_berend.plot(ax=ax, show=False)

        return fig

    @add_fig_kwargs
    def plot(self, **kwargs) -> Figure:
        """
        Plot ...
        for all items stored DiffusionDataList.
        """
        nrows, ncols = len(self), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()
        if len(self) % ncols != 0: ax_list[-1].axis("off")
        for i, (data, ax) in enumerate(zip(self, ax_list)):
            data.plot(ax=ax, show=False)

        return fig

    #@add_fig_kwargs
    #def plot_arrhenius(self, ax=None, **kwargs) -> Figure:
    #    df = self.get_dataframe()
    #    ax, fig, plt = get_ax_fig_plt(ax=ax)
    #    return fig



class MultiMdAnalyzer(HasPickleIO):
    """
    High-level interface to analyze multiple MD trajectories
    """

    @classmethod
    def from_abiml_dirs(cls, directories: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
        """
        Build an instance from a list of directories produced by abiml.py md
        """
        nprocs, pool_cls, using_msg = nprocs_poolcls_msg(len(directories), pmode=pmode)
        using_msg = f"Reading {len(directories)} abiml directories {using_msg}"
        args = [(dirpath, step_skip) for dirpath in directories]
        with pool_cls(nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_abiml_dir, args))

    @classmethod
    def from_lammps_dirs(cls, directories: str, step_skip=1, pmode="processes") -> MdAnalyzer:
        """
        Build an instance from a list of directories containing LAMMPS results.
        """
        nprocs, pool_cls, using_msg = nprocs_poolcls_msg(len(directories), pmode=pmode)
        using_msg = f"Reading {len(directories)} LAMMPS directories {using_msg}"
        args = [(dirpath, step_skip) for dirpath in directories]
        with pool_cls(nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_lammps_dir, args))

    #@classmethod
    #def from_lammpstrj(cls, traj_filepath: str, log_filepath: str, step_skip=1) -> MdAnalyzer:
    #    """
    #    Build an instance from a LAMMPS trajectory.
    #    """
    #    nprocs, pool_cls, using_msg = nprocs_poolcls_msg(len(vasprun_filepaths), pmode=pmode)
    #    using_msg = f"Reading {len(vasprun_filespaths} vasprun files {using_msg}"
    #    args = [(vrun, step_skip) for vrun in vasprun_filepaths]
    #    with pool_cls(nprocs) as pool, Timer(header=using_msg, footer="") as timer:
    #        return cls(pool.starmap(MdAnalyzer.from_vaspruns, args))

    @classmethod
    def from_vaspruns(cls, vasprun_filepaths: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
        """
        Build an instance from a list of vasprun files.
        """
        nprocs, pool_cls, using_msg = nprocs_poolcls_msg(len(vasprun_filepaths), pmode=pmode)
        using_msg = f"Reading {len(vasprun_filespaths)} vasprun files {using_msg}..."
        args = [(vrun, step_skip) for vrun in vasprun_filepaths]
        with pool_cls(nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_vaspruns, args))

    @classmethod
    def from_hist_files(cls, hist_filepaths: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
        """
        Build an instance from a list of ABINIT HIST.nc files.
        """
        nprocs, pool_cls, using_msg = nprocs_poolcls_msg(len(hist_filepaths), pmode=pmode)
        using_msg = f"Reading {len(hist_filespaths)} HIST.nc files {using_msg}..."
        args = [(ncpath, step_skip) for ncpath in hist_filepaths]
        with pool_cls(nprocs) as pool, Timer(header=using_msg, footer="") as timer:
            return cls(pool.starmap(MdAnalyzer.from_hist_file, args))

    #@classmethod
    #def from_qe_files(cls, qe_filepaths: list, step_skip=1, pmode="processes") -> MultiMdAnalyzer:
    #    """Build an instance from a list of QE input files."""
    #    nprocs, pool_cls, using_msg = nprocs_poolcls_msg(len(hist_filepaths), pmode=pmode)
    #    using_msg = f"Reading {len(qe_filespaths)} QE files {using_msg}..."
    #    args = [(path, step_skip) for path in qe_filepaths]
    #    with pool_cls(nprocs) as pool, Timer(header=using_msg, footer="") as timer:
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

    def get_msq_tt0_symbol(self, symbol: str, tmax: float, atom_inds=None, nprocs=None) -> Msdtt0List:
        msdtt0_list = Msdtt0List()
        for mda in self:
            obj = mda.get_msq_tt0_symbol(symbol, tmax, atom_inds=atom_inds, nprocs=nprocs)
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

    #def get_arrhenius_data(self, symbol, params_list, nprocs=None) -> DiffusionDataList:
    #    """
    #    """
    #    if len(params_list) != len(self):
    #        raise ValueError(f"{len(params)=} != {len(self)=}")

    #    data_list = DiffusionDataList()
    #    for it, ((mda, temp), params) in enumerate(zip(self.iter_mdat(), params_list)):
    #        p = AttrDict(**params)
    #        msq_tt0 = mda.get_msq_tt0_symbol(symbol, p.tmax, nprocs=nprocs)
    #        sigma = msq_tt0.get_sigma_berend(t1=p.t1, t2=p.t2)
    #        data_t = msq_tt0.get_diffusion_with_sigma(p.fit_time_start, p.fit_time_stop, p.block_size1, p.block_size2, sigma)
    #        data_list.append(data_t)

    #    return data_list

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
                   or scalar e.g. ``left``. If left (right) is None, default values are used
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
                        label=f"T = {temp} K",
                        color=color,
                    )

        # Decorate axes.
        for ix, (ax, symbol) in enumerate(zip(ax_list, symbols)):
            ax.set_title(symbol, fontsize=fontsize)
            set_axlims(ax, xlims, "x")
            ax.legend(fontsize=fontsize, loc=2)
            ax.set_xlabel('t (ps)', fontsize=fontsize)
            ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
            set_ticks_fontsize(ax, fontsize)
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


def nprocs_poolcls_msg(nprocs: int | None, pmode: str) -> tuple:
    """

    Args:
        nprocs:
        pmode: "threads", "processes" or "seq"
    """
    from multiprocessing.pool import ThreadPool
    from multiprocessing import Pool
    import os
    max_nprocs = max(1, os.cpu_count())

    if pmode == "seq":
        nprocs, pool_cls = 1, ThreadPool

    elif pmode == "threads":
        pool_cls = ThreadPool
        if nprocs is not None:
            pool_cls = ThreadPool if nprocs > 0 else Pool
            nprocs = min(abs(nprocs), max_nprocs)

    elif pmode == "processes":
        pool_cls = Pool
        if nprocs is not None:
            pool_cls = Pool if nprocs > 0 else ThreadPool
            nprocs = min(abs(nprocs), max_nprocs)

    else:
        raise ValueError(f"Invalid value of {pmode=}, it should be in ['seq', 'threads', 'processes']")

    nprocs = nprocs or max_nprocs
    return nprocs, pool_cls, f"using {nprocs=} with {pmode=} and Pool class: {pool_cls.__name__} ..."


def _func(size_t0: int, it: int, pos_tac: np.ndarray, msd_tt0: np.ndarray):
    msd_tt0[it,:] = np.mean(np.sum((pos_tac[it:it+size_t0,:,:] - pos_tac[:size_t0,:,:])**2, axis=2), axis=1)


def msd_tt0_from_tac(pos_tac: np.ndarray, size_t: int, nprocs=None) -> np.ndarray:
    r"""
    Calculates the MSD for every possible pair of time points in the input array, using the formula:

        $$MSD(t,t_0) = \frac{1}{N} \sum_{i=1}^{N} (\vec{r}_i(t+t_0) - \vec{r}_i(t_0))^2$$

    where $N$ is the number of particles, and $\vec{r}_i(t)$ is the position vector.
    """
    # Check if size_t is valid.
    n_time_points = pos_tac.shape[0]
    if size_t >= n_time_points:
        raise ValueError(f"{size_t=} must be less than {n_time_points}")

    size_t0 = pos_tac.shape[0] - size_t
    msd_tt0 = np.empty((size_t, size_t0), dtype=float)

    #print(f"msd_tt0_from_tac: {size_t=}, {size_t0=} with {nprocs=}")
    if nprocs == 1:
        for it in range(0, size_t):
            #for it0 in range(0, size_t0):
            #    msd_tt0[it,it0] = np.mean(np.sum((pos_tac[it+it0,:,:] - pos_tac[it0,:,:])**2, axis=1))
            msd_tt0[it,:] = np.mean(np.sum((pos_tac[it:it+size_t0,:,:] - pos_tac[:size_t0,:,:])**2, axis=2), axis=1)
    else:
        nprocs, pool_cls, using_msg = nprocs_poolcls_msg(nprocs, pmode="threads")
        using_msg = f"Computing MSD(t,t_0) matrix of shape {(size_t, size_t0)} {using_msg}"
        args = [(size_t0, it, pos_tac, msd_tt0) for it in range(0, size_t)]
        with pool_cls(nprocs) as pool, Timer(header=using_msg, footer="") as timer:
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

    Return: data_in_block, sigma, delta_sigma
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

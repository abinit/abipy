"""
"""
from __future__ import annotations

#import dataclasses
#import functools
import warnings
import json
import numpy as np
import pandas as pd
import pymatgen.core.units as units

from pathlib import Path
from multiprocessing import Pool
from matplotlib.offsetbox import AnchoredText
from monty.functools import lazy_property
from monty.bisect import find_le
from monty.string import list_strings, is_string, marquee
from pymatgen.util.string import latexify
from abipy.core.mixins import TextFile, NotebookWriter
from abipy.core.structure import Structure
from abipy.tools.typing import Figure
from abipy.tools.serialization import HasPickleIO
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, get_color_symbol,
                                  set_ticks_fontsize, set_logscale)


__author__ = "Giuliana Materzanini, Matteo Giantomassi"


def _pos_tac_structure_from_traj(traj_filepath) -> tuple[np.ndarray, Structure]:
    """
    Read all configurations from an ASE trajectory file.
    Returns (nsteps, natom, 3) array with the Cartesian coords and the initial structure.
    """
    from ase.io import read
    traj = read(traj_filepath, index=":")
    nsteps, natoms = len(traj), len(traj[0])
    pos_tac = np.empty((nsteps, natoms, 3))
    for it, atoms in enumerate(traj):
        pos_tac[it] = atoms.positions

    structure = Structure.as_structure(traj[0])
    del traj
    return pos_tac, structure


class MdAnalyzer(HasPickleIO):
    """
    High-level interface to read MD trajectories and metadata from external files,
    compute the MSQD and visualize the results.
    """

    @classmethod
    def from_abiml_dir(cls, directory) -> MdAnalyzer:
        """
        Build an instance from a directory containing an ASE trajectory file and
        a JSON file with the MD parameters as produced by `abiml.py md`.
        """
        directory = Path(str(directory))
        pos_tac, structure = _pos_tac_structure_from_traj(directory / "md.traj")

        # Read metadata from the JSON file.
        with open(directory / "md.json", "rt") as fh:
            meta = json.load(fh)

        temperature = meta["temperature"]
        timestep = meta["timestep"]
        loginterval = meta["loginterval"]
        nsteps = len(pos_tac)
        times = np.arange(0, nsteps) * timestep * loginterval

        return cls(structure, temperature, times, pos_tac)

    #@classmethod
    #def from_hist_file(cls, hist_filepath: str) -> MdAnalyzer:
    #    """
    #    Build an instance from a list of ABINIT HIST.nc files.
    #    """
    #    from abipy.dynamics.hist import HistFile
    #    with HistFile(hist_filepath) as hist:
    #        pos_tac = self.r.read_value("xcart") * units.bohr_to_ang
    #    times = np.arange(0, nsteps) * timestep * loginterval
    #    temperature = None
    #    return cls(hist.structure.copy(), temperature, times, pos_tac)

    @classmethod
    def from_vaspruns(cls, filepaths) -> MdAnalyzer:
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
        times = np.arange(0, nsteps) * timestep * step_kip

        return cls(structure, temperature, times, pos_tac)

    #@classmethod
    #def from_cp_pos_and_input(cls, pos_filepath: str, inp_filepath: str):
    #    """
    #    Build an instance from a directory
    #    """
    #    # Get structure from QE input.
    #    from pymatgen.io.pwscf import PWInput
    #    qe_inp = PWInput(inp_filepath)
    #    # Get atomic positions from qe pos file.
    #    pos_tac = None
    #    times = None
    #    temperature =
    #    return cls(qe_inp.structure, temperature, times, pos_tac)

    @classmethod
    def from_lammpstrj(cls, traj_filepath: str, log_filepath: str) -> MdAnalyzer:
        """
        Build an instance from a LAMMPS trajectory.
        """
        pos_tac, structure = _pos_tac_structure_from_traj(traj_filepath)

        if log_filepath.endswith(".evp"):
            # Extract times from CP EVP file (Giuliana's way)
            from abipy.dynamics.cpx import EvpFile
            with EvpFile(log_filepath) as evp:
                times = evp.times.copy()
        else:
            raise NotImplementedError()

        temperature = None
        return cls(structure, temperature, times, pos_tac)

    def __init__(self,
                structure: Structure,
                temperature: float,
                times: np.ndarray,
                cart_positions: np.ndarray,
                lattices=None,
                pos_order: str="tac"):
        """
        Args:
            structure: Structure object (first geometry of the MD run).
            temperature: Temperature in Kelvin.
            times: Array with times in ps units.
            cart_positions: Cartesian positions in Ang. Default shape: (nt, natom, 3).
            lattices: Array of lattice matrix of every step. Used for NPT.
                For NVT-AIMD, the lattice at each time step is set to the lattice in the "structure" argument.
            pos_order: "tac" if cart_positions has shape (nt, natom, 3).
                       "atc" if cart_positions has shape (natom, nt, 3).
        """
        self.structure = structure
        self.times = times

        if pos_order == "tac":
           self.pos_atc = cart_positions.transpose(1, 0, 2).copy()
        elif pos_order == "atc":
           self.pos_atc = cart_positions
        else:
           raise ValueError(f"Invalid {pos_order=}")

        #self.lattices = lattices
        #if lattices is None:
        #    self.lattices = np.array([structure.lattice.matrix.tolist()])

        self.latex_formula = self.structure.latex_formula
        self.temperature = temperature

        self.set_color_symbol("VESTA")
        self.consistency_check()

    def consistency_check(self):

        # Consistencty check.
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

    def deepcopy(self) -> MdAnalyzer:
        """Deep copy of the object."""
        import copy
        return copy.deepcopy(self)

    def prune_nsteps(self, start_nsteps: int, tmesh_nsteps) -> None:
        return self.prune_time(start_time=100 * self.timestep, new_timestep=tmesh_nsteps * self.timestep)

    def prune_time(self, start_time: float, new_timestep: float) -> None:
        """
        """
        old_timestep = self.times[1] - self.times[0]
        if not (self.times[-1] > start_time > self.times[0]):
            raise ValueError(f"Invalid start_time should be between {self.times[0]} and {self.times[-1]})")

        it0 = int(start_time / old_timestep)
        print(it0)
        if it0 != 0:
            self.pos_atc = self.pos_atc[:,it0:,:]
            self.times = self.times[it0:] - self.times[it0]

        if new_timestep < old_timestep:
            raise ValueError(f"Invalid delta_time should be >= {old_timestep}")

        istep = int(new_timestep / old_timestep)
        if istep != 1:
            self.pos_atc = self.pos_atc[:,::istep,:].copy()
            self.times = self.times[::istep] - self.times[0]

        self.consistency_check()

    @property
    def timestep(self) -> float:
        """Timestep in ps."""
        return self.times[1] - self.times[0]

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
    def latex_formula(self) -> str:
        return self._latex_formula

    @latex_formula.setter
    def latex_formula(self, value):
        """LaTeX formatted formula. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$."""
        self._latex_formula = latexify(value)

    @lazy_property
    def avg_volume(self) -> float:
        """Average unit cell volume in Ang^3."""
        return self.structure.lattice.volume

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
        app(f"temperature = {self.temperature} K")
        app(f"timestep = {self.timestep} ps")
        app(f"max_time = {self.times[-1]} ps")
        app(f"trajectory_size = {self.nt}")
        #app(self.structure.to_string(verbose=verbose))
        #if verbose:
        #    app(self.structure.spget_summary(verbose=verbose))

        return "\n".join(lines)

    def iatoms_with_symbol(self, symbol: str, atom_inds=None) -> np.ndarray:
        """
        Array with the index of the atoms with the given chemical symbol.
        If atom_inds is not None, filter sites accordingly.
        """
        iatoms = [iat for iat in range(len(self.structure)) if self.structure[iat].specie.symbol == symbol]
        if atom_inds is not None:
            iatoms = [iat for iat in iatoms if iat in atom_inds]

        return np.array(iatoms)

    def select_symbols(self, symbols) -> list[str]:
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

    def export_msdt(self, filename: str) -> None:
        """
        Writes MSD data to a file that can be easily plotted in other software.

        Args:
            filename: Supported formats are csv and dat.
                If the extension is csv, a csv file is written. Otherwise, a dat format is assumed.
        """
        fmt = "csv" if filename.lower().endswith(".csv") else "dat"
        delimiter = ", " if fmt == "csv" else " "
        #with open(filename, "w") as f:
        #    if fmt == "dat": f.write("# ")
        #    f.write(delimiter.join(["t", "MSD", "MSD_a", "MSD_b", "MSD_c", "MSCD"]))
        #    f.write("\n")
        #    for dt, msd, msdc, mscd in zip(
        #        self.dt, self.msd, self.msd_components, self.mscd
        #    ):
        #        f.write(delimiter.join([str(v) for v in [dt, msd, *list(msdc), mscd]]))
        #        f.write("\n")

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

        for symbol in self.select_symbols(symbols):
            for iatom in self.iatoms_with_symbol(symbol, atom_inds=atom_inds):
                sd_t = self.get_sqdt_iatom(iatom, it0=it0)
                ax.plot(ts, sd_t, label=f"${symbol}_{{{iatom}}}$")

        ax.legend(fontsize=5, loc=2)
        ax.set_ylabel(r'square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(
                        self.latex_formula + '\n'+ 'T = ' + str(self.temperature) + ' K' + '\n' +
                        #'V' + '_' + 'ave = ' + str(self.avg_volume) + r'$\mathrm{{\AA}^3}$' + '\n' +
                        'sd(t, t0 =' + str(int(self.times[it0])) + ' ps)',
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

        for symbol in self.select_symbols(symbols):
            ax.plot(ts, self.get_sqdt_symbol(symbol, it0=it0, atom_inds=atom_inds),
                    label=symbol + ' msd(t, t0 = ' + str(t0) + ' ps)',
                    color=self.color_symbol[symbol],
                    )

        ax.legend(fontsize=16, loc=2)
        ax.set_ylabel('mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(
                self.latex_formula + '\n'+ 'T =' + str(self.temperature) + ' K', # +'\n', +
                #'V' + '_' + 'ave = ' + str(2186.87) + '$\mathrm{{\AA}^3}$',
                loc=1, prop=dict(size=20)))

        return fig

    def get_msq_tt0_symbol(self, symbol, tmax: float, atom_inds=None):
        index_tmax, _ = self.get_it_ts(tmax)
        print(f"{index_tmax=}")
        iatoms = self.iatoms_with_symbol(symbol, atom_inds=atom_inds)
        tac = self.pos_atc[iatoms].transpose(1, 0, 2).copy()
        msd_tt0 = timeMeanSquareDispAllinOne(tac, index_tmax)

        return MSDTT0(msd_tt0, self, index_tmax, symbol)


    @add_fig_kwargs
    def plot_sqdt_symbols_tmax(self, symbols, tmax: float, atom_inds=None,
                               ax=None, xy_log=None, fontsize=30, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement averaged over all atoms of the same specie vs time.

        Args:
            symbols: List of chemical symbols to consider. "all" for all symbols in structure.
            tmax: Max time in ps.
            atom_inds: List of atom indices to include. None to disable filtering.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
                "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
            fontsize: fontsize for legends and titles
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
        """
        index_tmax, _ = self.get_it_ts(tmax)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self.select_symbols(symbols):
            iatoms = self.iatoms_with_symbol(symbol, atom_inds=atom_inds)
            tac = self.pos_atc[iatoms].transpose(1, 0, 2).copy()
            #print(f"{tac.shape=}")
            #MSD = timeMeanSquareDispAllinOne(posDict['Li'], index_tmax)
            #ts = timeArray[indexTTotal-index_tmax:] - timeArray[indexTTotal-index_tmax]
            #TODO Finalize and optimize timeMeanSquareDispAllinOne
            MSD = timeMeanSquareDispAllinOne(tac, index_tmax)
            msdS = np.mean(MSD, axis=1)
            t_start = self.nt - index_tmax
            #print(f"{t_start=}, {len(self.times)=}")
            #print(f"{self.times[0]=}, {self.times[-1]=}")
            ts = self.times[t_start:] - self.times[t_start]

            ax.plot(ts, msdS,
                    label=symbol + ' <msd(t, t0)>$\{$t0$\}$, t = [0, ' + str(int(self.times[indexTMax]))+ ' ps]',
                    color=self.color_symbol[symbol],
                    )

        set_axlims(ax, xlims, "x")
        ax.legend(fontsize=16, loc=2)
        ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(
                         self.latex_formula + '\n' + 'T =' + str(self.temperature) + ' K', # +'\n' +
                         #'V' +'_'+'ave = ' + str(2186.87) + '$\mathrm{{\AA}^3}$',
                         loc=1, prop=dict(size=20)))

        return fig



class MSDTT0:

    def __init__(self, msd_tt0, mda, index_tmax, symbol):
        self.msd_tt0 = msd_tt0
        self.mda = mda
        self.index_tmax = index_tmax
        self.symbol = symbol

    #def __str__(self):
    #def to_string(self, verbose) -> str:

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
        msdS = np.mean(self.msd_tt0, axis=1)
        index_tmax = self.index_tmax
        t_start = self.mda.nt - index_tmax
        #print(f"{t_start=}, {len(self.times)=}")
        #print(f"{self.times[0]=}, {self.times[-1]=}")
        ts = self.times[t_start:] - self.times[t_start]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(ts, msdS,
                label=self.symbol + ' <msd(t, t0)>$\{$t0$\}$, t = [0, ' + str(int(self.times[index_tmax]))+ ' ps]',
                color=self.mda.color_symbol[self.symbol],
                )

        set_axlims(ax, xlims, "x")
        ax.legend(fontsize=16, loc=2)
        ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        set_ticks_fontsize(ax, fontsize)
        set_logscale(ax, xy_log)
        ax.add_artist(AnchoredText(
                         self.mda.latex_formula + '\n' + 'T =' + str(self.temperature) + ' K', # +'\n' +
                         #'V' +'_'+'ave = ' + str(self.avg_volume) + '$\mathrm{{\AA}^3}$',
                         loc=1, prop=dict(size=20)))

        return fig

    @add_fig_kwargs
    def plot_sigma(self, first_time, second_time, nblock_step=1, tot_block=1000,
                   fontsize=8, **kwargs) -> Figure:
        """
        """
        #print("In plot_sigma")
        # choose the time elapsed
        estimatedFirstTElapsed = first_time
        timeArray = self.times
        MSD = self.msd_tt0
        symbol = self.symbol
        temperature = self.temperature
        latex_formula = self.mda.latex_formula

        firstTElapsed = find_nearest(timeArray, estimatedFirstTElapsed)
        indexFirstTElapsed=int(np.nonzero(timeArray == firstTElapsed)[0])

        estimatedSecondTElapsed = second_time
        secondTElapsed = find_nearest(timeArray, estimatedSecondTElapsed)
        indexSecondTElapsed=int(np.nonzero(timeArray == secondTElapsed)[0])


        ax1, sigmaArrFirst, dataInBlockArrFirst = sigma_berend(nblock_step, tot_block,
                                                               MSD[indexFirstTElapsed,:], timeArray[indexFirstTElapsed],
                                                               self.temperature, self.mda.latex_formula)



        estimatedDataInBlockFirst = 200
        dataInBlockFirst = find_nearest(dataInBlockArrFirst, estimatedDataInBlockFirst)
        indexDataInBlockFirst=int(np.nonzero(dataInBlockArrFirst == dataInBlockFirst)[0])
        sigmaFirst = sigmaArrFirst[indexDataInBlockFirst]
        print(sigmaFirst)

        ax2, sigmaArrSecond, dataInBlockArrSecond = sigma_berend(nblock_step, tot_block,
                                                                 MSD[indexSecondTElapsed,:], timeArray[indexSecondTElapsed],
                                                                 self.temperature, self.mda.latex_formula)

        estimatedDataInBlockSecond = 200
        dataInBlockSecond = find_nearest(dataInBlockArrSecond, estimatedDataInBlockSecond)
        indexDataInBlockSecond=int(np.nonzero(dataInBlockArrSecond == dataInBlockSecond)[0])
        sigmaSecond = sigmaArrSecond[indexDataInBlockSecond]
        print(sigmaSecond)

        # fit to a linear behaviour the errors
        mSigma = (sigmaSecond-sigmaFirst)/(indexSecondTElapsed-indexFirstTElapsed)
        qSigma = sigmaFirst - mSigma*indexFirstTElapsed
        mDataInBlock = (dataInBlockSecond-dataInBlockFirst)/(indexSecondTElapsed-indexFirstTElapsed)
        qDataInBlock = dataInBlockFirst - mDataInBlock*indexFirstTElapsed
        #print(qDataInBlock,mDataInBlock )

        # and find error for anytime
        errMSD = np.zeros(MSD.shape[0],dtype=float)
        for t in range(MSD.shape[0]):
            errMSD[t] = abs(mSigma * t + qSigma)

        # and find error for anytime
        dataScorrelated = np.zeros(MSD.shape[0],dtype=int)
        for t in range(MSD.shape[0]):
            dataScorrelated[t] = int(mDataInBlock * t + qDataInBlock)

        errMSD_Li = errMSD
        msdS=np.mean(MSD,axis=1)  # msdS= msd_small, as it is the average over the initial times

        #print(timeArray.shape[0])
        #print(msdS.shape[0])
        timeArrayHalf = timeArray[:msdS.shape[0]]
        msdSScorrelatedList = []
        timeArrScorrelatedList = []
        errMSDScorrelatedList = []
        #timeIndexMin = int(1000/timeStepJump)
        # counter= timeIndexMin

        estimatedStartFitTime=50
        # TO BE DECIDED!!!!!!!
        estimatedEndFitTime=150
        #estimatedEndFitTime=timeArrayHalf[timeArrayHalf.shape[0]-1]


        timeArrStartFit = find_nearest(timeArray, estimatedStartFitTime)
        indextimeArrStartFit=int(np.nonzero(timeArray==timeArrStartFit)[0])

        timeArrEndFit = find_nearest(timeArray, estimatedEndFitTime)
        indextimeArrEndFit=int(np.nonzero(timeArray==timeArrEndFit)[0])

        counter= indextimeArrStartFit
        # print (indextimeArrStartFit)
        condition=True

        while condition:
            if counter >= indextimeArrEndFit:
                condition=False
            else:
                index=dataScorrelated[counter]
                msdSScorrelatedList.append(msdS[counter])
                timeArrScorrelatedList.append(timeArray[counter])
                errMSDScorrelatedList.append(errMSD[counter])
                counter = counter + index
        msdSScorrelated=np.array(msdSScorrelatedList,dtype=float)
        timeArrScorrelated=np.array(timeArrScorrelatedList,dtype=float)
        errMSDScorrelated=np.array(errMSDScorrelatedList,dtype=float)

        msdS_Li = msdS

        angCoeffMSD, quoteMSD, varangCoeffMSD,varquoteMSD = linear_lsq_linefit(msdSScorrelated,timeArrScorrelated,1/(errMSDScorrelated)**2)

        nCarriers = len(self.mda.structure.indices_from_symbol(self.symbol))
        print(f"{nCarriers=} for {self.symbol=}")

        volAve = 2186.87
        ax1, fig, plt = get_ax_fig_plt(ax=None)
        # plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
        # ax1.set_ylabel(r'$\langle\mathrm{MSD}_\mathrm{tr}^\mathrm{Li}\rangle \mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
        ax1.set_ylabel(r'$\mathrm{MSD}_\mathrm{tr}$ $\mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
        ax1.set_xlabel('t (ps)', fontsize=18)

        print(e2s)
        DiffCoeff=angCoeffMSD*Ang2PsTocm2S/6
        ErrDiffCoeff=np.sqrt(varangCoeffMSD)*Ang2PsTocm2S/6
        Conduct=e2s/kbs*nCarriers*DiffCoeff/volAve/temperature * 1.e09
        ErrConduct=e2s/kbs*nCarriers/volAve/temperature*ErrDiffCoeff * 1.e09

        print('{:.2E}'.format(angCoeffMSD*Ang2PsTocm2S/6,2))

        anchored_text = AnchoredText('D$_{tr}$ = (' + str('{:.2E}'.format(DiffCoeff)) + '\u00B1' + str('{:.2E}'.format(ErrDiffCoeff)) + ') cm$^2$/s', loc=2, prop=dict(size=14))

        ax1.add_artist(anchored_text)

        ax1.errorbar(timeArrayHalf,msdS_Li,yerr=errMSD_Li, color='mediumblue', label = 'Li')#,label='Li (D ~ 1.5e-05cm$^2$/s)')

        ax1.errorbar(timeArrScorrelated,msdSScorrelated,yerr=errMSDScorrelated, linestyle='-')
        ax1.errorbar(timeArrayHalf,angCoeffMSD*timeArrayHalf+quoteMSD, linestyle='--')
        ax1.errorbar(timeArrayHalf,(angCoeffMSD-np.sqrt(varangCoeffMSD))*timeArrayHalf+quoteMSD, linestyle='--')
        ax1.errorbar(timeArrayHalf,(angCoeffMSD+np.sqrt(varangCoeffMSD))*timeArrayHalf+quoteMSD, linestyle='--')
        print (angCoeffMSD)
        print (angCoeffMSD+np.sqrt(varangCoeffMSD))
        print (angCoeffMSD-np.sqrt(varangCoeffMSD))
        ax1.legend(fontsize=12, loc=4)
        anchored_text_1 = AnchoredText(latex_formula +'\n' +
                                '192 atoms' +'\n'+ 'T =' + str(temperature) + ' K' +'\n'+ '500 ps', loc=1, prop=dict(size=14))
        ax1.add_artist(anchored_text_1)
        #ax1.set_ylim([0,300])

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_ylabel(r'$\mathrm{MSD}_\mathrm{tr}$ $\mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
        ax1.set_xlabel('t (ps)', fontsize=18)
        #ax1.set_title('Li msd averaged over the initial times: ' + str(latex_formula) + ', ' + str(temperature) + 'K-NVT' )
        #ax1.set_title('Li msd averaged over the initial times: ' + str(latex_formula) + ', NVE(' + '929' + 'K)')


        #print(e2s)
        DiffCoeff=angCoeffMSD*Ang2PsTocm2S/6
        ErrDiffCoeff=np.sqrt(varangCoeffMSD)*Ang2PsTocm2S/6
        Conduct=e2s/kbs*nCarriers*DiffCoeff/volAve/temperature * 1.e09
        ErrConduct=e2s/kbs*nCarriers/volAve/temperature*ErrDiffCoeff * 1.e09

        print('{:.2E}'.format(angCoeffMSD*Ang2PsTocm2S/6,2))

        anchored_text = AnchoredText('D$_{tr}$ = (' + str('{:.2E}'.format(DiffCoeff)) + '\u00B1' + str('{:.2E}'.format(ErrDiffCoeff)) + ') cm$^2$/s', loc=2, prop=dict(size=18))
        ax1.add_artist(anchored_text)

        ax1.errorbar(timeArrayHalf,msdS_Li,yerr=errMSD_Li, color='mediumblue', label = 'Li')#,label='Li (D ~ 1.5e-05cm$^2$/s)')

        ax1.errorbar(timeArrScorrelated,msdSScorrelated,yerr=errMSDScorrelated, linestyle='-')
        ax1.errorbar(timeArrayHalf,angCoeffMSD*timeArrayHalf+quoteMSD, linestyle='--')
        ax1.errorbar(timeArrayHalf,(angCoeffMSD-np.sqrt(varangCoeffMSD))*timeArrayHalf+quoteMSD, linestyle='--')
        ax1.errorbar(timeArrayHalf,(angCoeffMSD+np.sqrt(varangCoeffMSD))*timeArrayHalf+quoteMSD, linestyle='--')
        print (angCoeffMSD)
        print (angCoeffMSD+np.sqrt(varangCoeffMSD))
        print (angCoeffMSD-np.sqrt(varangCoeffMSD))
        # plt.xlabel('time (ps)')
        # plt.ylabel('<MSD>(A2)')
        # plotMSD = plt.errorbar(timeArray[:msdS.shape[0]],msdS,yerr=errMSD)
        # plotMSD = plt.errorbar(timeArrScorrelated,msdSScorrelated,yerr=errMSDScorrelated)
        ax1.legend(fontsize=12)
        #plt.savefig(fileMsdPngName,dpi=300, bbox_inches='tight', pad_inches=0)

        """
        outMSD_file=open(fileMsdOutName, "w+")
        outDiff_file=open(fileDiffusionName, "a")
        outCond_file=open(fileConductivityName, "a")
        #outMSD_file.write('IONS CORRELATED'+ '\n')
        outMSD_file.write('timeStepJump = ' + str(timeStepJump) + '\n')
        outMSD_file.write('temperature = ' + str(temperature) + '\n')
        #outMSD_file.write('input pos file = ' + str(filePosName) + '\n')
        outMSD_file.write('time, cel and pos were cut before ' + str(timeArray[0]+tInitial)+ 'ps' + '\n')
        outMSD_file.write('t elapsed max is ' + str(estimatedTMax)+ 'ps' + '\n')
        outMSD_file.write('trajectory length is ' + str(timeArrayTmp[timeArrayTmp.shape[0]-1])+ 'ps' + '\n')
        outMSD_file.write('error on msd(t) evaluated at ' +  str(estimatedFirstTElapsed) + 'ps' + ' and ' + str(estimatedSecondTElapsed) + 'ps' + '\n')
        outMSD_file.write('evaluated n. of blocks at ' + str(estimatedFirstTElapsed) + 'ps' + ' is ' + str(estimatedDataInBlockFirst)  + '\n')
        outMSD_file.write('evaluated n. of blocks at ' + str(estimatedSecondTElapsed) + 'ps' + ' is ' + str(estimatedDataInBlockSecond)  + '\n')
        outMSD_file.write('msd(t) fit starts at ' + str(estimatedStartFitTime)+ 'ps' + '\n')
        outMSD_file.write('msd(t) fit ends at ' + str(estimatedEndFitTime)+ 'ps' + '\n')
        outMSD_file.write('number of decorrelated msd(t) data that we fit: ' + str(timeArrScorrelated.shape[0]) + '\n')
        outMSD_file.write('min angular coefficient of msd(t): ' + str(angCoeffMSD-np.sqrt(varangCoeffMSD)) + '\n')
        outMSD_file.write('max angular coefficient of msd(t): ' + str(angCoeffMSD+np.sqrt(varangCoeffMSD)) + '\n')
        outMSD_file.write('best angular coefficient of msd(t): ' + str(angCoeffMSD) + '\n')
        outMSD_file.write('diffusion coefficient = ' + str(DiffCoeff) + '+-' + str(ErrDiffCoeff) + '\n')
        outMSD_file.write('conductivity = ' + str(Conduct) + '+-' + str(ErrConduct) + '\n')
        outMSD_file.write(str(temperature) + ' ' + str(DiffCoeff) + ' ' + str(ErrDiffCoeff) + '\n')
        outMSD_file.write(str(temperature) + ' ' + str(Conduct) + ' ' + str(ErrConduct) + '\n')
        outDiff_file.write(str(temperature) + ' ' + str(DiffCoeff) + ' ' + str(ErrDiffCoeff) + ' ' +str(volumeArrayTmp[0])+ '\n')
        outCond_file.write(str(temperature) + ' ' + str(Conduct) + ' ' + str(ErrConduct) + '\n')
        outMSD_file.close()
        outDiff_file.close()
        outCond_file.close()

        #
        outMSDT_file=open(fileMsdTOutName, "w+")
        for i in range(indexTMax):
            outMSDT_file.write(str(timeArray[i])+ ' '  + str(msdS[i]) + '\n')
        outMSDT_file.close()
        """





class MultiMdAnalyzer(HasPickleIO):
    """
    High-level interface to analyze MD trajectories computed
    for the same system at different temperatures.
    """

    @classmethod
    def from_abiml_dirs(cls, directories) -> MultiMdAnalyzer:
        """Build instance from a list of directories produced by abiml.py md"""
        return cls([MdAnalyzer.from_abiml_dir(d) for d in directories])

    def __init__(self, mdas: list[MdAnalyzer], colormap="jet"):
        """
        Args:
            mdas: List of MdAnalyzer
            colormap: matplotlib colormap per temperatures.
        """
        # Sort analyzers according to temperature.
        self.mdas = sorted(mdas, key=lambda x: x.temperature)

        # Consistency check.
        #for i in range(1, len(self)):
        #    if np.any((self.mdas[i].times - self.mdas[i-1].times).abs() > 1e-3):
        #        raise ValueError(f"Found different time meshes for indices: {i} and {i-1}")

        self.set_colormap(colormap)

    def __iter__(self):
        return self.mdas.__iter__()

    def __len__(self) -> int:
        return len(self.mdas)

    def __getitem__(self, items):
        return self.mdas.__getitem__(items)

    #@property
    #def ntemps(self) -> int:
    #    """Number of temperatures."""
    #    return len(self)

    def set_colormap(self, colormap) -> None:
        """Set the colormap for the list of temperatures."""
        import matplotlib.pyplot as plt
        self.cmap = plt.get_cmap(colormap)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose."""
        lines = []
        app = lines.append
        #lines = [mda.to_string(verbose=verbose) for mda in self]

        return "\n".join(lines)

    def iter_atc(self):
        """Iterate over (MdAnalyzer, temperature, color)."""
        for itemp, mda in enumerate(self):
            yield mda, mda.temperature, self.cmap(float(itemp) / len(self))

    #def color_itemp(self, itemp: int):
    #    return self.cmap(float(itemp) / len(self))

    #def temps_colors(self) -> tuple[list, list]:
    #    return ([mda.temperature for mda in self],
    #            [self.color_itemp(itemp) for itemp in range(len(self))])

    def _nrows_ncols_nplots(self, size=None):
        size = size or len(self)
        nrows, ncols, nplots = 1, 1, size
        if nplots > 1:
            ncols = 2; nrows = nplots // ncols + nplots % ncols
        return nrows, ncols, nplots

    #def get_arrenhius_data(self):
    #    return

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
        symbols = self[0].select_symbols(symbols)
        nrows, ncols, nplots = self._nrows_ncols_nplots(size=len(symbols))
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols, sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()
        if nplots % ncols != 0: ax_list[-1].axis("off")

        # Plot data.
        for itemp, (mda, temp, color) in enumerate(self.iter_atc()):
            it0, ts = mda.get_it_ts(t0)
            for ix, (ax, symbol) in enumerate(zip(ax_list, symbols)):
                ax.plot(ts, mda.get_sqdt_symbol(symbol, it0=it0),
                        #label=f" msd(t, t0 = {t0}  ps), T = {temp} K",
                        label=f"T = {temp} K",
                        color=color,
                    )

        # Decorate axes.
        for ix, (ax, symbol) in enumerate(zip(ax_list, symbols)):
            ax.set_title(symbol, fontsize=fontsize)
            set_axlims(ax, xlims, "x")
            ax.legend(fontsize=fontsize, loc=2)
            ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
            ax.set_xlabel('t (ps)', fontsize=fontsize)
            set_ticks_fontsize(ax, fontsize)
            set_logscale(ax, xy_log)

        return fig

    @add_fig_kwargs
    def plot_arrhenius(self, ax=None, **kwargs) -> Figure:
        """
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        return fig





Ang2PsTocm2S=0.0001
bohr2a = 0.529177
e2s = 1.602188**2 # electron charge in Coulomb scaled by 10.d-19**2
kbs = 1.38066     # Boltzmann constant in Joule/K scaled by 10.d-23


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


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


#def read_cel_from_file_cel(f_name: str, nstep_cel: int) -> tuple[np.ndarray, int]:
#    """
#    """
#    with open(str(f_name), 'rt') as file_cel:
#        cel_every_step = []
#        file_cel_lines = file_cel.readlines()
#        nlinetot = len(file_cel_lines)
#        nt = int(nlinetot/nstep_cel)
#        celArrayB = np.zeros((nt,nstep_cel-1,3), dtype=float)
#        #
#        for it in range(nt):
#            cel_every_step = []
#            for line in file_cel_lines[(nstep_cel*it)+1:nstep_cel*(it+1)]:
#                y = line.split()
#                y = np.array(y, dtype=float)
#                cel_every_step.append(y)
#            celArrayB[it,:,:] = np.array(cel_every_step, dtype=float)
#        celArray = np.copy(celArrayB*bohr2a)
#        return celArray, nt


def timeMeanSquareDispAllinOne(pos_tac: np.ndarray, site_t: int) -> np.ndarray:
    r"""
    Calculates the MSD for every possible pair of time points in the input array, using the formula:

        $$MSD(t,t_0) = \frac{1}{N} \sum_{i=1}^{N} (\vec{r}_i(t+t_0) - \vec{r}_i(t_0))^2$$

    where $N$ is the number of particles, $\vec{r}_i(t)$ is the position vector
    """
    # Check if site_t is valid.
    n_time_points = pos_tac.shape[0]
    if site_t >= n_time_points:
        raise ValueError(f"site_t must be less than {n_time_points}")
    size_t0 = pos_tac.shape[0] - site_t

    print(f"timeMeanSquareDispAllinOne: {site_t=}, {size_t0=}")
    msd_tt0 = np.zeros((site_t, size_t0), dtype=float)

    for it in range(0, site_t):
        for it0 in range(0, size_t0):
            msd_tt0[it,it0] = np.mean(np.sum((pos_tac[it+it0,:,:] - pos_tac[it0,:,:])**2, axis=1))

    #msd_tt0 = np.mean((pos_tac[site_t:, :, :] - pos_tac[:n_time_points - site_t, :, :])**2, axis=(-1, -2))
    print("done", msd_tt0.shape)

    return msd_tt0


def block_mean_var(data, data_mean, n_block):
    """
    Perform the block mean and the block variance of data
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


def sigma_berend(nblock_step, tot_block, data, timeElapsedConsidered, temperature, material, ax=None):
    """
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
    for nblock in range(1, nblock_step * tot_block,nblock_step):
        sigma2[counter], delta_sigma2[counter] = block_mean_var(data, mean, nblock+1)
        arr_nblock[counter] = nblock + 1
        data_in_block[counter] = int(Ndata/(nblock+1))
        counter = counter + 1

    ax, fig, plt = get_ax_fig_plt(ax=ax)
    ax.set_ylabel('$\sigma$ ($\AA^2$)', fontsize=14)
    ax.set_xlabel('N. of data in block', fontsize=14)
    #ax.set_title('Variance of correlated data as function of block number.')
    ax.grid(True)
    from matplotlib.ticker import MultipleLocator
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    set_ticks_fontsize(ax, 14)
    #ax.xaxis.grid(True, which='minor')
    sigma = np.sqrt(sigma2)
    delta_sigma = 0.5 * delta_sigma2 / sigma
    ax.errorbar(data_in_block, sigma, yerr=delta_sigma,
                linestyle='-', linewidth=0.5,
                label="$\sigma(\mathrm{MSD}($"+'%2.1f' % timeElapsedConsidered+" ps$))$ "+ '\n' + material + ', '+ '%4.0f' % temperature + 'K')
    ax.legend(fontsize=16, loc=4)

    return ax, sigma, data_in_block


#if __name__ == "__main__":
#    main()

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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
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
                                  set_ticks_fontsize)


def _pos_tac_structure_from_traj(traj_filepath) -> tuple[np.ndarray, Structure]:
    """
    Read all configurations from an ASE trajectory file.
    Returns (nsteps, natom, 3) array with cartesian coords and the initial structure.
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
    """

    @classmethod
    def from_abiml_dir(cls, directory) -> MdAnalyzer:
        """
        Build an instance from an ASE traj file and a JSON file with MD parameters.
        """
        directory = Path(str(directory))
        pos_tac, structure = _pos_tac_structure_from_traj(directory / "md.traj")

        # Read parameters from the JSON file produced by `abiml.py md`
        with open(directory / "md.json", "rt") as fh:
            meta = json.load(fh)

        temperature = meta["temperature"]
        timestep = meta["timestep"]
        loginterval = meta["loginterval"]

        nsteps = len(pos_tac)
        times = np.arange(0, nsteps) * timestep * loginterval
        new = cls(structure, times, pos_tac)
        return new

    #@classmethod
    #def from_hist_file(cls, hist_filepath: str) -> MdAnalyzer:
    #    """
    #    Build an instance from a list of ABINIT HIST.nc files.
    #    """
    #    from abipy.dynamics.hist import HistFile
    #    with HistFile(hist_filepath) as hist:
    #        pos_tac = self.r.read_value("xcart") * units.bohr_to_ang
    #    times = None
    #    temperature = None
    #    new = cls(hist.structure.copy(), times, pos_tac)
    #    return new

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

        # Extract cartesian positions.
        pos_tac = []
        for i, strc in enumerate(s):
            if i == 0: structure = strc
            pos_tac.append(strc.coords)

        nsteps, natom = i + 1, len(structure)
        pos_tac = np.reshape(pos_tac, (nsteps, natom, 3))

        times = np.arange(0, nsteps) * timestep * step_kip
        new = cls(structure, times, pos_tac)
        #new.temperature = temperature
        return new

    #@classmethod
    #def from_cp_pos_and_input(cls, qepos_filepath: str, qeinp_filepath: str):
    #    """
    #    Build an instance from a directory
    #    """
    #    # Get structure from QE input.
    #    from pymatgen.io.pwscf import PWInput
    #    qe_inp = PWInput(qeinp_filepath)
    #    # Get atomic positions from qe pos file.
    #    pos_tac = None
    #    times = None
    #    temperature =
    #    new = cls(qe_inp.structure, times, pos_tac)
    #    return new

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
                times = evp.df["tps(ps)"].values.copy()
        else:
            raise NotImplementedError

        temperature = None
        new = cls(structure, times, pos_tac)
        return new

    def __init__(self,
                structure: Structure,
                times: np.ndarray,
                cart_positions: np.ndarray,
                lattices=None,
                pos_order: str="tac"):
        """
        Args:
            structure: Structure object (first geometry of the MD run).
            times: List of times in ps units.
            cart_positions: Cartesian positions in Ang.
            lattices: array of lattice matrix of every step. Used for NPT.
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
        self.temperature = 0.0

        # Consistencty check.
        if self.pos_atc.shape != (self.natom, self.nt, 3):
            raise ValueError(f"Invalid shape {self.pos_atc.shape=}, expecting: {(self.natom, self.nt, 3)}")
        if len(self.times) != self.nt:
            raise ValueError(f"{len(self.times)=} != {self.nt=}")

        # Check times mesh
        ierr = 0
        for it in range(self.nt-1):
            dt = self.times[it+1] - self.times[it]
            if abs(dt - self.timestep) > 1e-3:
                ierr += 1
                if ierr < 10: print(f"{dt=} != {self.timestep=}")
        if ierr:
            raise ValueError(f"Time-mesh is not linear. There are {ierr} points with wrong timestep")

        self.set_color_symbol("VESTA")

    @lazy_property
    def timestep(self) -> float:
        """Timestep in ps."""
        return self.times[1] - self.times[0]

    @lazy_property
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
        self._temperature = float(value)

    @property
    def latex_formula(self) -> str:
        return self._latex_formula

    @latex_formula.setter
    def latex_formula(self, value):
        """LaTeX formatted formula. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$."""
        self._latex_formula = latexify(value)

    def set_color_symbol(self, dict_or_string) -> None:
        """
        Set the mapping chemical_symbol --> color used in the matplotlib plots.

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

    def get_it0_ts(self, t0: float) -> tuple[int, np.ndarray]:
        """
        Return the index of time t0 in self.times and the array with the time values.
        """
        it0 = find_le(self.times, t0)
        return it0, self.times[it0:] - self.times[it0]

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose."""
        lines = []
        app = lines.append
        app(self.structure.to_string(verbose=verbose))

        return "\n".join(lines)

    def iatoms_with_symbol(self, symbol) -> np.ndarray:
        """Array with the index of the atoms with chemical symbol."""
        return np.array([iatom for iatom in range(len(self.structure)) if self.structure[iatom].specie.symbol == symbol])

    def _select_symbols(self, symbols) -> list[str]:
        if symbols == "all": return sorted(self.structure.symbol_set)
        return list_strings(symbols)

    def get_sqdt_iatom(self, iatom: int, it0: int = 0) -> np.array:
        """
        Compute the square displacement vs time for a given atomic index
        starting from time index it0
        """
        return ((self.pos_atc[iatom,it0:] - self.pos_atc[iatom,it0]) ** 2).sum(axis=1)

    def get_sqdt_symbol(self, symbol: str, it0: int = 0) -> np.array:
        """
        Compute the square displacement vs time averaged over atoms with the same chemical symbol
        starting from time index it0
        """
        for count, iatom in enumerate(self.iatoms_with_symbol(symbol)):
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
    def plot_sqdt_atoms(self, symbols="all", t0: float = 0.0,
                        ax=None, fontsize=30, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement of each atom vs time.

        Args:
            symbols:
            t0: Initial time in ps.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
        """
        it0, ts = self.get_it0_ts(t0)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            for iatom in self.iatoms_with_symbol(symbol):
                sd_t = self.get_sqdt_iatom(iatom, it0=it0)
                ax.plot(ts, sd_t, label=f"${symbol}_{{iatom}}$")

        ax.legend(fontsize=5, loc=2)
        ax.set_ylabel(r'square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        set_ticks_fontsize(ax, fontsize)
        anchored_text = AnchoredText(self.latex_formula + '\n'+ 'T = ' + str(self.temperature) + ' K' + '\n' +
                                     #'V' + '_' + 'ave = ' + str(2186.87) + r'$\mathrm{{\AA}^3}$' + '\n' +
                                     'sd(t, t0 =' + str(int(self.times[it0])) + ' ps)',
                                     loc=1, prop=dict(size=20))
        ax.add_artist(anchored_text)

        return fig

    @add_fig_kwargs
    def plot_sqdt_symbols(self, symbols, t0: float = 0.0,
                          ax=None, fontsize=30, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement averaged over all atoms of the same specie vs time.

        Args:
            symbols: List of chemical symbols.
            t0: Initial time in ps.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
        """
        it0, ts = self.get_it0_ts(t0)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            ax.plot(ts, self.get_sqdt_symbol(symbol, it0=it0),
                    label=symbol + ' msd(t, t0 = ' + str(t0) + ' ps)',
                    color=self.color_symbol[symbol],
                    )

        ax.legend(fontsize=16, loc=2)
        ax.set_ylabel('mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        set_axlims(ax, xlims, "x")
        set_ticks_fontsize(ax, fontsize)
        anchored_text = AnchoredText(self.latex_formula + '\n'+ 'T =' + str(self.temperature) + ' K', # +'\n', +
                                     #'V' + '_' + 'ave = ' + str(2186.87) + '$\mathrm{{\AA}^3}$',
                                     loc=1, prop=dict(size=20))
        ax.add_artist(anchored_text)

        return fig

    @add_fig_kwargs
    def plot_sqdt_symbols_tmax(self, symbols, tmax: float,
                               ax=None, fontsize=30, xlims=None, **kwargs) -> Figure:
        """
        Plot the square displacement averaged over all atoms of the same specie vs time.

        Args:
            symbols: List of chemical symbols.
            tmax: Max time in ps.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
        """
        #TODO Finalize and optimize timeMeanSquareDispAllinOne
        indexTMax, _ = self.get_it0_ts(tmax)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for symbol in self._select_symbols(symbols):
            iatoms = self.iatoms_with_symbol(symbol)
            tac = self.pos_atc[iatoms].transpose(1, 0, 2).copy()
            #print(f"{tac.shape=}")
            MSD = timeMeanSquareDispAllinOne(tac, indexTMax)
            msdS = np.mean(MSD, axis=1)
            t_start = self.nt - indexTMax
            ts = self.times[t_start:] - self.times[t_start]
            #MSD = timeMeanSquareDispAllinOne(posDict['Li'], indexTMax)
            #ts = timeArray[indexTTotal-indexTMax:] - timeArray[indexTTotal-indexTMax]

            ax.plot(ts, msdS,
                    label=symbol + ' <msd(t, t0)>$\{$t0$\}$, t = [0, ' + str(int(self.times[indexTMax]))+ ' ps]',
                    color=self.color_symbol[symbol],
                    )

        set_axlims(ax, xlims, "x")
        ax.legend(fontsize=16, loc=2)
        ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=fontsize)
        ax.set_xlabel('t (ps)', fontsize=fontsize)
        set_ticks_fontsize(ax, fontsize)
        anchored_text = AnchoredText(self.latex_formula + '\n' + 'T =' + str(self.temperature) + ' K', # +'\n' +
                                     #'V' +'_'+'ave = ' + str(2186.87) + '$\mathrm{{\AA}^3}$',
                                     loc=1, prop=dict(size=20))
        ax.add_artist(anchored_text)

        return fig


__author__ = "Giuliana Materzanini"

Ang2PsTocm2S=0.0001
bohr2a = 0.529177
e2s = 1.602188**2 # electron charge in Coulomb scaled by 10.d-19**2
kbs = 1.38066     # Boltzmann constant in Joule/K scaled by 10.d-23


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def linearLSQLineFit(x, z, weights):
    S00 = np.sum(weights)
    S10 = np.sum(z*weights)
    S01 = np.sum(x*weights)
    S20 = np.sum(z**2*weights)
    S11 = np.sum((x*z)*weights)
    D = S00*S20-S10**2
    m = (S00*S11-S10*S01)/D
    q = (S01*S20-S11*S10)/D
    varM = S00/D
    varQ = S20/D
    return m, q, varM, varQ


def cumulativeSpecies(speciesDict: dict):
    """
    """
    counter=0
    listTmp=[]
    while counter < len(list(speciesDict.keys())):
        for name, speciesList in speciesDict.items():
            if speciesList[0]==counter:
                listTmp.append(speciesList[1])
                counter=counter+1
    listCumul=[]
    speciesCumulStartDict={}
    speciesCumulEndDict={}
    counter=0
    for i in range(len(listTmp)):
        listCumul.append(listTmp[i]+counter)
        counter=counter+listTmp[i]
        for name, speciesList in speciesDict.items():
             if speciesList[1]==listTmp[i]:
                speciesCumulEndDict[name]=listCumul[i]
                if i==0:
                    speciesCumulStartDict[name]=0
                else:
                    speciesCumulStartDict[name]=listCumul[i-1]

    return speciesCumulStartDict, speciesCumulEndDict


def read_from_file_evp(f_name: str) -> np.ndarray:
    """
    """
    with open(str(f_name), 'rt') as file_evp:
        every_step = []
        file_evp_lines = file_evp.readlines()
        nlinetot = len(file_evp_lines)
        nt = nlinetot
        for it in range(0,nlinetot):
            y = file_evp_lines[it].split()
            y = np.array(y, dtype=float)
            every_step.append(y)
        timeArray = np.array(every_step, dtype=float)
        return timeArray


def read_pos_from_file_pos(f_name, nstep_pos) -> tuple[np.ndarray, int]:
    """
    """
    with open(str(f_name), 'rt') as file_pos:
        pos_every_step = []
        file_pos_lines = file_pos.readlines()
        nlinetot = len(file_pos_lines)
        nt = int(len(file_pos_lines)/nstep_pos)
        posArrayB = np.zeros((nt,nstep_pos-1,3), dtype=float)

        for it in range(nt):
            pos_every_step = []
            for line in file_pos_lines[(nstep_pos*it)+1:nstep_pos*(it+1)]:
                y = line.split()
                y = np.array(y, dtype=float)
                pos_every_step.append(y)
            posArrayB[it,:,:] = np.array(pos_every_step, dtype=float)

        posArray = np.copy(posArrayB*bohr2a)
        return posArray, nt


def read_cel_from_file_cel(f_name: str, nstep_cel: int) -> tuple[np.ndarray, int]:
    """
    """
    with open(str(f_name), 'rt') as file_cel:
        cel_every_step = []
        file_cel_lines = file_cel.readlines()
        nlinetot = len(file_cel_lines)
        nt = int(nlinetot/nstep_cel)
        celArrayB = np.zeros((nt,nstep_cel-1,3), dtype=float)
        #
        for it in range(nt):
            cel_every_step = []
            for line in file_cel_lines[(nstep_cel*it)+1:nstep_cel*(it+1)]:
                y = line.split()
                y = np.array(y, dtype=float)
                cel_every_step.append(y)
            celArrayB[it,:,:] = np.array(cel_every_step, dtype=float)
        celArray = np.copy(celArrayB*bohr2a)
        return celArray, nt


def read_propTris_from_file_prop(f_name, nstep_cel, nThrow, nThrow1, nCol) -> tuple[np.ndarray, int]:
    """"""
    with open(str(f_name), 'rt') as file_prop:
        prop_every_step = []
        file_prop_lines = file_prop.readlines()
        nlinetot = len(file_prop_lines)
        nstep_tot = nstep_cel+nThrow+nThrow1
        nt = int(len(file_prop_lines)/(nstep_tot))
        propArray = np.zeros((nt,nstep_cel,nCol), dtype=float)

        for it in range(nt):
            prop_every_step = []
            for line in file_prop_lines[nThrow+(nstep_tot*it):nThrow+(nstep_tot*it)+nstep_cel]:
                y = line.split()
                y = np.array(y, dtype=object)
                prop_every_step.append(y)
            propArray[it,:] = np.array(prop_every_step, dtype=float)

        return propArray, nt


def read_propTrisObj_from_file_prop(f_name, nstep_cel, nThrow, nThrow1, nCol) -> tuple[np.ndarray, int]:
    """"""
    with open(f_name, 'rt') as file_prop:
        prop_every_step = []
        file_prop_lines = file_prop.readlines()
        nlinetot = len(file_prop_lines)
        nstep_tot = nstep_cel+nThrow+nThrow1
        nt = int(len(file_prop_lines)/(nstep_tot))
        propArray = np.zeros((nt,nstep_cel,nCol), dtype=object)

        for it in range(nt):
            prop_every_step = []
            for line in file_prop_lines[nThrow+(nstep_tot*it):nThrow+(nstep_tot*it)+nstep_cel]:
                y = line.split()
                y = np.array(y, dtype=object)
                prop_every_step.append(y)
            propArray[it,:] = np.array(prop_every_step, dtype=object)

        return propArray, nt


def fromPos2DictPos(speciesDict: dict, posArray: np.ndarray) -> dict:
    """
    Args:
        speciesDict:
        posArray:
    """
    dictPos = {}
    speciesCumulStartDict, speciesCumulEndDict = cumulativeSpecies(speciesDict)
    for species in speciesCumulStartDict:
        dictPos[species] = np.copy(posArray[:,speciesCumulStartDict[species]:speciesCumulEndDict[species],:])

    return dictPos


#def timeMeanSquareDispOneAtomOneDimTrial(posT: np.ndarray) -> np.ndarray:
#    """"""
#    indexTmax = int(posT.shape[0]/2)
#    msd = np.zeros(indexTmax, dtype=float)
#    for t in range(indexTmax):
#        summT0 = 0
#        for t0 in range(indexTmax):
#            summT0 = summT0 + (posT[t+t0]-posT[t0])**2/indexTmax
#        msd[t] = summT0
#
#    return msd


def timeMeanSquareDispOneAtomOneDim(posT: np.ndarray) -> np.ndarray:
    """
    Calculate MSD of each atom depending on t and t0
    """
    indexTmax = int(posT.shape[0]/2)
    sizeT = int(indexTmax)
    sizeT0 = int(indexTmax)
    msd = np.zeros((sizeT, sizeT0), dtype=float)
    for it in range(0, sizeT):
        tStar = it
        for it0 in range(0, sizeT0):
            t0Star = it0
            msd[it,it0] = (posT[tStar+t0Star] - posT[t0Star])**2

    return msd


def timeMeanSquareDispAllinOne(pos_tac: np.ndarray, max_index: int) -> np.ndarray:
    # Check if max_index is valid.
    n_time_points = pos_tac.shape[0]
    if max_index >= n_time_points:
        raise ValueError(f"max_index must be less than {n_time_points}")
    size_t0 = pos_tac.shape[0] - max_index

    print(f"timeMeanSquareDispAllinOne: {max_index=}, {size_t0=}")
    msd_tt0 = np.zeros((max_index, size_t0), dtype=float)

    for it in range(0, max_index):
        for it0 in range(0, size_t0):
            msd_tt0[it,it0] = np.mean(np.sum((pos_tac[it+it0,:,:] - pos_tac[it0,:,:])**2, axis=1))

    #msd_tt0 = np.mean((pos_tac[max_index:, :, :] - pos_tac[:n_time_points - max_index, :, :])**2, axis=(-1, -2))
    print("done", msd_tt0.shape)

    return msd_tt0


def SquareDispOneAtom(pos_tac: np.ndarray, it0: int = 0) -> np.ndarray:
    """
    Calculate SD of each atom
    """
    nt = pos_tac.shape[0] - it0
    natoms = pos_tac.shape[1]
    sd = np.zeros((natoms, nt), dtype=float)
    for ia in range(natoms):
        for it in range(nt):
            sd[ia,it] = np.sum((pos_tac[it0+it,ia,:] - pos_tac[it0,ia,:])**2)

    return sd


def MeanSquareDisp(pos_tac: np.ndarray, it0: int = 0) -> np.ndarray:
    """
    """
    nt = pos_tac.shape[0] - it0
    natoms = pos_tac.shape[1]
    sd = np.zeros((natoms, nt), dtype=float)
    msd = np.zeros(nt, dtype=float)
    for ia in range(natoms):
        for it in range(nt):
            sd[ia,t] = np.sum((pos_tac[it0+it,ia,:] - pos_tac[it0,ia,:])**2)
    msd = np.mean(sd, axis = 0)

    return msd


def timeMeanSquareDispAllinOneIonCorrel(pos_tac: np.ndarray, indexTmax: int) -> np.ndarray:
    """"""
    sizeT = int(indexTmax)
    sizeT0 = int((pos_tac.shape[0] - indexTmax))
    msd = np.zeros((sizeT, sizeT0), dtype=float)
    for it in range(0, sizeT):
        tStar = it
        for t0 in range(0, sizeT0):
            t0Star = t0
            rDiff = pos_tac[tStar+t0Star,:,:] - pos_tac[t0Star,:,:]
            msd[it,t0] = np.sum(np.dot(rDiff,rDiff.T)) / pos_tac.shape[1]

    return msd


#def timeMeanSquareDispOneAtom(pos3DT: np.ndarray) -> np.ndarray:
#    """"""
#    msdListDim = []
#    for ic in range(3):
#        msdListDim.append(timeMeanSquareDispOneAtomOneDim(pos3DT[:,ic]))
#
#    msd = np.zeros((msdListDim[0].shape[0], msdListDim[0].shape[1]), dtype=float)
#    for msdDim in msdListDim:
#        msd = msd + msdDim
#
#    return msd


#def timeMeanSquareDisp(pos_tac: np.ndarray) -> np.ndarray:
#    """"""
#    msdListAtoms = []
#    for ia in range(pos_tac.shape[1]):
#        msdListAtoms.append(timeMeanSquareDispOneAtom(pos_tac[:,ia,:]))
#
#    msd = np.zeros((msdListAtoms[0].shape[0], msdListAtoms[1].shape[0]), dtype=float)
#    for msdAtom in msdListAtoms:
#        msd = msd + msdAtom / len(msdListAtoms)
#
#    return msd


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


def sigma_berend(nblock_step, tot_block, data, timeElapsedConsidered, temperature, material):
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

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel('$\sigma$ ($\AA^2$)', fontsize=14)
    ax.set_xlabel('N. of data in block', fontsize=14)
    #ax.set_title('Variance of correlated data as function of block number.')
    ax.grid(True)
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


def main():
    import matplotlib.pyplot as plt
    material = 'c-LLZO'
    materiall = 'LLZO'
    temperature = 600
    path = ""
    path0 = ""

    fileEvpName = path + materiall + '.evp'  # this is the log file but only the energy lines!
    #filePosName = path + materiall + '.pos'
    fileDumpName = path + 'dump.lammpstrj'

    fileMsdOutName = path + 'diff_coeff_' + material + '_' + str(temperature) + 'K.out'
    fileMsdTOutName = path + 'msd_' + material + '_' + str(temperature) + 'K.out'
    fileMsdPngName = path + 'msd-ave_' + material + '_' + str(temperature) + 'K.png'
    fileSigmaMSDPngName1 = path + 'sigma_' + material + '_' + str(temperature) + '_1.png'
    fileSigmaMSDPngName2 = path + 'sigma_' + material + '_' + str(temperature) + '_1.png'
    fileSdAllPngName = path + 'sd_' + material + '_' + str(temperature) + '.png'
    fileDiffusionName = path0 + 'diffusion_' + material + '.dat'
    fileConductivityName = path0 + 'conductivity_' + material + '.dat'
    fileVolumeOutName = path + 'volume_' + material + '_' + str(temperature) + 'K.out'
    fileSdName = path + 'sd_' + material + '_' + str(temperature) + '.dat'

    # natoms = 192
    speciesDict={'Li':[0,56],'La':[1,24],'Zr':[2,16],'O':[3,96]}
    nCarriers=speciesDict['Li'][1]
    natoms = speciesDict['Li'][1] + speciesDict['La'][1] + speciesDict['Zr'][1] + speciesDict['O'][1]
    print('material = ', material)
    print('temperature = ', temperature,'K')
    print('nCarriers = ', nCarriers)
    print('natoms = ', natoms)

    nblockStep = 1
    timeStepJump = 100

    tot_array_Tmp = read_from_file_evp(fileEvpName)

    mdArrayTmp = np.copy(tot_array_Tmp[:,0])
    timeArrayTmp = np.copy(tot_array_Tmp[:,1])
    tIonsArrayTmp = np.copy(tot_array_Tmp[:,2])
    eCPArrayTmp = np.copy(tot_array_Tmp[:,3])
    ekinIonsArrayTmp = np.copy(tot_array_Tmp[:,4])
    totEneArrayTmp = np.copy(tot_array_Tmp[:,5])
    pressArrayTmp = np.copy(tot_array_Tmp[:,6])
    volumeArrayTmp = np.copy(tot_array_Tmp[:,7])
    densityArrayTmp = np.copy(tot_array_Tmp[:,8])

    eConsArrayTmp = eCPArrayTmp + ekinIonsArrayTmp
    #eRestArrayTmp = eContArrayTmp - eConsArrayTmp - ekincArrayTmp
    nstepsTot = timeArrayTmp.shape[0]
    #print(timeArrayTmp[0])
    #print(nstepsTot)

    # stepAarray are the steps printed out in .evp file
    stepArrayTmp = np.arange(0, nstepsTot, dtype=int)
    #print(stepArrayTmp)
    #print(timeArrayTmp.shape[0])

    # estimate indexTInitial and cut timeArrayTmp from tInitial
    estimatedTForEnergy = 0.2
    tEne = find_nearest(timeArrayTmp, estimatedTForEnergy)
    indexTEne = int(np.nonzero(timeArrayTmp == tEne)[0])
    #print(indexTEne, tEne)
    #print(timeArrayTmp.shape[0])

    # Plot ionic temperature.
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.plot(timeArrayTmp[indexTEne:], tIonsArrayTmp[indexTEne:], label="$Temp_{ions}$", c='tab:red')
    ax.set_xlabel('time (ps)')
    ax.set_ylabel('temperature (K)')
    ax.text(0.35, 0.9, material,
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, color='black', fontsize=15)
    plt.legend(loc=1, bbox_to_anchor=(1, 1))
    plt.show()

    # Plot other energies.
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.plot(timeArrayTmp[indexTEne:], eCPArrayTmp[indexTEne:], label="$E_{CP}$", c='tab:cyan')
    ax.plot(timeArrayTmp[indexTEne:], eConsArrayTmp[indexTEne:], label="$E_{cons}$", c='tab:red')
    #ax.plot(timeArrayTmp[indexTEne:], totEneArrayTmp[indexTEne:], label="$E_{cont}$", c='black')
    ax.set_xlabel('time (ps)')
    ax.set_ylabel('energy (Hartree)')
    ax.text(0.8, 0.9, material,
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, color='black', fontsize=15)
    plt.show()

    estimatedTInitial = 1.0
    tInitial = find_nearest(timeArrayTmp, estimatedTInitial)
    indexTInitial = int(np.nonzero(timeArrayTmp == tInitial)[0])
    #print(indexTInitial, timeArrayTmp.shape[0])
    #print(timeArrayTmp.shape[0])

    # RE-READ timeArrayTmp so that you can run this cell at any time
    timeArrayTmp = read_from_file_evp(fileEvpName)[:,1]
    posArrayTotTmp, ntt = read_propTris_from_file_prop(fileDumpName, natoms, 9, 0, 5)
    posArrayTmp = np.zeros((ntt, natoms, 3), dtype=float)
    for i in range(ntt):
        for j in range(natoms):
            posArrayTmp[i,j,0] = np.copy(posArrayTotTmp[i,j,2])
            posArrayTmp[i,j,1] = np.copy(posArrayTotTmp[i,j,3])
            posArrayTmp[i,j,2] = np.copy(posArrayTotTmp[i,j,4])

    nstepsTot = timeArrayTmp.shape[0]
    stepArrayTmp =np.arange(0, nstepsTot, dtype=int)

    timeStepJump = 10
    #print('timeArrayTot = ', timeArrayTmp)
    #print('stepArrayTot = ', stepArrayTmp)
    posArrayTmp = posArrayTmp[indexTInitial:,:,:]
    timeArrayTmp = timeArrayTmp[indexTInitial:]
    stepArrayTmp = stepArrayTmp[indexTInitial:]

    # mask for timeStepJump
    timeArray = np.zeros(int(timeArrayTmp.shape[0]/timeStepJump), dtype=float)
    stepArray = np.zeros(int(stepArrayTmp.shape[0]/timeStepJump), dtype=int)
    posArray = np.zeros((int(posArrayTmp.shape[0]/timeStepJump), posArrayTmp.shape[1], posArrayTmp.shape[2]), dtype=float)
    for i in range(int(timeArrayTmp.shape[0]/timeStepJump)):
        t = i * timeStepJump
        timeArray[i] = timeArrayTmp[t]
        stepArray[i] = stepArrayTmp[t]
        posArray[i,:,:] = np.copy(posArrayTmp[t,:,:])

    #print('time0 = ',timeArray[0])
    #print('step0 = ',stepArray[0])
    timeArray = timeArray - timeArray[0]
    stepArray = stepArray - stepArray[0]
    posDict = fromPos2DictPos(speciesDict, posArray)
    #print('timeArray = ', timeArray)
    #print('stepArray = ', stepArray)

    # not important now
    estimatedTForStr1 = 5.
    tStr1 = find_nearest(timeArrayTmp, estimatedTForStr1)
    indexTStr1 = int(np.nonzero(timeArrayTmp == tStr1)[0])
    #print(indexTStr1, tStr1)
    #print(timeArrayTmp.shape[0])

    # not important now
    estimatedTForStr2 = 10.
    tStr2 = find_nearest(timeArrayTmp, estimatedTForStr2)
    indexTStr2 = int(np.nonzero(timeArrayTmp == tStr2)[0])
    #print(indexTStr2,tStr2)
    #print(timeArrayTmp.shape[0])

    sd = SquareDispOneAtom(posDict['Li'])
    with open(fileSdName, "w+") as fileSd:
        for t in range(timeArray.shape[0]):
            for i in range(sd.shape[0]):
                fileSd.write(str(timeArray[t]) + ' ' + str(sd[i,t]))
            fileSd.write('\n')

    species2Plot = 'Li'
    it0 = 2000
    sd = SquareDispOneAtom(posDict[species2Plot], it0=it0)
    #print(sd.shape[1])

    fig = plt.figure(figsize=(6.4*2, 4.8*2))
    ax = fig.add_subplot(111)
    for speciesN in range(sd.shape[0]):
        ax.plot(timeArray[it0:] - timeArray[it0], sd[speciesN,:], label=species2Plot + '_' + str(speciesN+1))
        #ax.plot(stepArray,sd[speciesN,:],label=species2Plot+'_'+str(speciesN+1))

    ax.legend(fontsize=5, loc=2)
    ax.set_ylabel('square displacement ($\mathrm{{\AA}^2}$)', fontsize=30)
    ax.set_xlabel('t (ps)', fontsize=30)
    ax.set_xlim(0, 500)
    set_ticks_fontsize(ax, 30)
    anchored_text = AnchoredText(material +'\n'+ 'T =' + str(temperature) + ' K' +'\n' +
                                 'V' +'_'+'ave = ' + str(2186.87) + '$\mathrm{{\AA}^3}$' +'\n'+
                                 'sd(t, t0 = ' + str(int(timeArray[it0]))+' ps)', loc=1, prop=dict(size=20))
    ax.add_artist(anchored_text)
    plt.show()

    it0 = 2000
    meansdLi = MeanSquareDisp(posDict[species2Plot], it0=it0)
    meansdLa = MeanSquareDisp(posDict['La'], it0=it0)
    meansdZr = MeanSquareDisp(posDict['Zr'], it0=it0)
    meansdO = MeanSquareDisp(posDict['O'], it0=it0)

    fig = plt.figure(figsize=(6.4*2, 4.8*2))
    ax = fig.add_subplot(111)
    ax.plot(timeArray[it0:] - timeArray[it0], meansdLi,
             label=species2Plot+' msd(t, t0 = ' + str(int(timeArray[it0])) + ' ps)',
             color='tab:blue')
    ax.plot(timeArray[it0:] - timeArray[it0],meansdLa,
             label='La msd(t, t0 = ' + str(int(timeArray[it0])) + ' ps)',
             color='black')
    ax.plot(timeArray[it0:] - timeArray[it0],meansdZr,
             label='Zr msd(t, t0 = ' + str(int(timeArray[it0])) + ' ps)',
             color='green')
    ax.plot(timeArray[it0:] - timeArray[it0],meansdO,
             label='O msd(t, t0 = ' + str(int(timeArray[it0]))+' ps)',
             color='red')
    #ax.set_xlim(0, 500)
    ax.legend(fontsize=16, loc=2)
    ax.set_ylabel('mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=30)
    ax.set_xlabel('t (ps)', fontsize=30)
    set_ticks_fontsize(ax, 30)
    anchored_text = AnchoredText(material +'\n'+ 'T =' + str(temperature) + ' K' +'\n' +
                                 'V' +'_'+'ave = ' + str(2186.87) + '$\mathrm{{\AA}^3}$', loc=1, prop=dict(size=20))
    ax.add_artist(anchored_text)
    plt.show()

    # estimate indexTMax
    estimatedTMax = 150
    tMax = find_nearest(timeArray, estimatedTMax)
    indexTMax = int(np.nonzero(timeArray == tMax)[0])
    #print(tMax)

    MSD = timeMeanSquareDispAllinOne(posDict['Li'], indexTMax)
    #MSD = timeMeanSquareDispAllinOneIonCorrel(posDict['Li'], indexTMax)

    indexTTotal = timeArray.shape[0]
    #print(indexTTotal)
    msdS = np.mean(MSD, axis=1)

    fig = plt.figure(figsize=(6.4*2, 4.8*2))
    ax = fig.add_subplot(111)
    ax.plot(timeArray[indexTTotal-indexTMax:] - timeArray[indexTTotal-indexTMax],msdS,
            label=species2Plot+' <msd(t, t0)>$\{$t0$\}$, t = [0, ' + str(int(timeArray[indexTMax]))+ ' ps]',
            color='tab:blue')
    ax.set_xlim(0, 500)
    ax.legend(fontsize=16, loc=2)
    ax.set_ylabel('average mean square displacement ($\mathrm{{\AA}^2}$)', fontsize=30)
    ax.set_xlabel('t (ps)', fontsize=30)
    set_ticks_fontsize(ax, 30)
    anchored_text = AnchoredText(material +'\n'+ 'T =' + str(temperature) + ' K' +'\n' +
                                 'V' +'_'+'ave = ' + str(2186.87) + '$\mathrm{{\AA}^3}$',
                                 loc=1, prop=dict(size=20))
    ax.add_artist(anchored_text)
    #plt.savefig(fileSdAllPngName,dpi=300, bbox_inches='tight', pad_inches=0.01)
    plt.show()

    # choose the time elapsed
    timeInitial = np.zeros(timeArray.shape[0], dtype=float)
    estimatedFirstTElapsed = 20
    firstTElapsed = find_nearest(timeArray, estimatedFirstTElapsed)
    indexFirstTElapsed = int(np.nonzero(timeArray == firstTElapsed)[0])
    estimatedSecondTElapsed = 50
    secondTElapsed = find_nearest(timeArray, estimatedSecondTElapsed)
    indexSecondTElapsed = int(np.nonzero(timeArray == secondTElapsed)[0])

    ax, sigmaArrFirst, dataInBlockArrFirst = sigma_berend(1, 1000, MSD[indexFirstTElapsed,:], timeArray[indexFirstTElapsed],
                                                          temperature, material)
    #plt.savefig(fileSigmaMSDPngName1,dpi=300, bbox_inches='tight', pad_inches=0)

    estimatedDataInBlockFirst = 200
    dataInBlockFirst = find_nearest(dataInBlockArrFirst, estimatedDataInBlockFirst)
    indexDataInBlockFirst = int(np.nonzero(dataInBlockArrFirst == dataInBlockFirst)[0])
    sigmaFirst = sigmaArrFirst[indexDataInBlockFirst]
    #print(sigmaFirst)

    ax2, sigmaArrSecond, dataInBlockArrSecond = sigma_berend(1, 1000, MSD[indexSecondTElapsed,:], timeArray[indexSecondTElapsed],
                                                            temperature, material)
    #plt.savefig(fileSigmaMSDPngName2,dpi=300, bbox_inches='tight', pad_inches=0)

    estimatedDataInBlockSecond = 200
    dataInBlockSecond = find_nearest(dataInBlockArrSecond, estimatedDataInBlockSecond)
    indexDataInBlockSecond = int(np.nonzero(dataInBlockArrSecond == dataInBlockSecond)[0])
    sigmaSecond = sigmaArrSecond[indexDataInBlockSecond]
    #print(sigmaSecond)

    # fit to a linear behaviour the errors
    mSigma = (sigmaSecond-sigmaFirst) / (indexSecondTElapsed-indexFirstTElapsed)
    qSigma = sigmaFirst - mSigma*indexFirstTElapsed
    mDataInBlock = (dataInBlockSecond-dataInBlockFirst) / (indexSecondTElapsed-indexFirstTElapsed)
    qDataInBlock = dataInBlockFirst - mDataInBlock*indexFirstTElapsed
    #print(qDataInBlock,mDataInBlock )

    # and find error for anytime
    errMSD = np.zeros(MSD.shape[0], dtype=float)
    for t in range(MSD.shape[0]):
        errMSD[t] = abs(mSigma * t + qSigma)

    # and find error for anytime
    dataScorrelated = np.zeros(MSD.shape[0], dtype=int)
    for t in range(MSD.shape[0]):
        dataScorrelated[t] = int(mDataInBlock * t + qDataInBlock)
        #print(t,dataScorrelated[t])

    errMSD_Li = errMSD
    #errMSD_S = errMSD
    #errMSD_Ge = errMSD
    #errMSD_P = errMSD

    msdS = np.mean(MSD, axis=1)  # msdS= msd_small, as it is the average over the initial times

    print(timeArray.shape[0])
    print(msdS.shape[0])
    timeArrayHalf = timeArray[:msdS.shape[0]]
    msdSScorrelatedList = []
    timeArrScorrelatedList = []
    errMSDScorrelatedList = []
    timeIndexMin = int(1000/timeStepJump)
    # counter= timeIndexMin

    estimatedStartFitTime = 50
    # TO BE DECIDED!!!!!!!
    estimatedEndFitTime = 150
    #estimatedEndFitTime=timeArrayHalf[timeArrayHalf.shape[0]-1]

    timeArrStartFit = find_nearest(timeArray, estimatedStartFitTime)
    indextimeArrStartFit = int(np.nonzero(timeArray==timeArrStartFit)[0])

    timeArrEndFit = find_nearest(timeArray, estimatedEndFitTime)
    indextimeArrEndFit = int(np.nonzero(timeArray==timeArrEndFit)[0])

    counter= indextimeArrStartFit
    #print (indextimeArrStartFit)
    condition = True
    while condition:
        if counter >= indextimeArrEndFit:
            condition = False
        else:
            index=dataScorrelated[counter]
            msdSScorrelatedList.append(msdS[counter])
            timeArrScorrelatedList.append(timeArray[counter])
            errMSDScorrelatedList.append(errMSD[counter])
            counter = counter + index

    msdSScorrelated = np.array(msdSScorrelatedList, dtype=float)
    timeArrScorrelated = np.array(timeArrScorrelatedList, dtype=float)
    errMSDScorrelated = np.array(errMSDScorrelatedList, dtype=float)

    msdS_Li = msdS
    #msdS_S = msdS
    #msdS_Ge = msdS
    #msdS_P = msdS

    angCoeffMSD, quoteMSD, varangCoeffMSD, varquoteMSD = linearLSQLineFit(msdSScorrelated, timeArrScorrelated, 1/(errMSDScorrelated)**2)

    volAve = 2186.87
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    #ax.set_ylabel(r'$\langle\mathrm{MSD}_\mathrm{tr}^\mathrm{Li}\rangle \mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{MSD}_\mathrm{tr}$ $\mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
    ax.set_xlabel('t (ps)', fontsize=18)
    set_ticks_fontsize(ax, 14)
    #ax.set_title('Li msd averaged over the initial times: ' + str(material) + ', ' + str(temperature) + 'K-NVT' )
    #ax.set_title('Li msd averaged over the initial times: ' + str(material) + ', NVE(' + '929' + 'K)')

    #print(e2s)
    DiffCoeff = angCoeffMSD*Ang2PsTocm2S/6
    ErrDiffCoeff = np.sqrt(varangCoeffMSD)*Ang2PsTocm2S/6
    Conduct = e2s/kbs*nCarriers*DiffCoeff/volAve/temperature * 1.e09
    ErrConduct = e2s/kbs*nCarriers/volAve/temperature*ErrDiffCoeff * 1.e09
    #print('{:.2E}'.format(angCoeffMSD*Ang2PsTocm2S/6,2))

    anchored_text = AnchoredText(
        'D$_{tr}$ = (' + str('{:.2E}'.format(DiffCoeff)) + '\u00B1' + str('{:.2E}'.format(ErrDiffCoeff)) + ') cm$^2$/s',
        loc=2, prop=dict(size=14))

    ax.add_artist(anchored_text)
    ax.errorbar(timeArrayHalf, msdS_Li, yerr=errMSD_Li, color='mediumblue', label = 'Li')#,label='Li (D ~ 1.5e-05cm$^2$/s)')
    #ax.errorbar(timeArrayHalf, msdS_P,yerr=errMSD_S, color='darkviolet',label='P')
    #ax.errorbar(timeArrayHalf, msdS_Ge,yerr=errMSD_S, color='violet', label='Ge')
    #ax.errorbar(timeArrayHalf, msdS_S,yerr=errMSD_S, color='navajowhite', label='S')

    ax.errorbar(timeArrScorrelated, msdSScorrelated, yerr=errMSDScorrelated, linestyle='-')
    ax.errorbar(timeArrayHalf, angCoeffMSD*timeArrayHalf+quoteMSD, linestyle='--')
    ax.errorbar(timeArrayHalf, (angCoeffMSD-np.sqrt(varangCoeffMSD))*timeArrayHalf+quoteMSD, linestyle='--')
    ax.errorbar(timeArrayHalf, (angCoeffMSD+np.sqrt(varangCoeffMSD))*timeArrayHalf+quoteMSD, linestyle='--')
    #print(angCoeffMSD)
    #print(angCoeffMSD + np.sqrt(varangCoeffMSD))
    #print(angCoeffMSD - np.sqrt(varangCoeffMSD))

    #plt.xlabel('time (ps)')
    #plt.ylabel('<MSD>(A2)')
    #plotMSD = plt.errorbar(timeArray[:msdS.shape[0]],msdS,yerr=errMSD)
    #plotMSD = plt.errorbar(timeArrScorrelated,msdSScorrelated,yerr=errMSDScorrelated)
    ax.legend(fontsize=12, loc=4)
    #plt.savefig(fileMsdPngName,dpi=300, bbox_inches='tight', pad_inches=0)
    anchored_text_1 = AnchoredText(material +'\n' +
                                  '192 atoms' +'\n'+ 'T =' + str(temperature) + ' K' +'\n'+ '500 ps', loc=1, prop=dict(size=14))
    ax.add_artist(anchored_text_1)
    ax.set_ylim([0, 300])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    # ax.set_ylabel(r'$\langle\mathrm{MSD}_\mathrm{tr}^\mathrm{Li}\rangle \mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{MSD}_\mathrm{tr}$ $\mathrm{(\AA}^2\mathrm{)}$', fontsize=18)
    ax.set_xlabel('t (ps)', fontsize=18)
    set_ticks_fontsize(ax, 14)

    # ax.set_title('Li msd averaged over the initial times: ' + str(material) + ', ' + str(temperature) + 'K-NVT' )
    # ax.set_title('Li msd averaged over the initial times: ' + str(material) + ', NVE(' + '929' + 'K)')
    #print(e2s)
    DiffCoeff = angCoeffMSD*Ang2PsTocm2S/6
    ErrDiffCoeff = np.sqrt(varangCoeffMSD)*Ang2PsTocm2S/6
    Conduct = e2s/kbs*nCarriers*DiffCoeff/volAve/temperature * 1.e09
    ErrConduct = e2s/kbs*nCarriers/volAve/temperature*ErrDiffCoeff * 1.e09

    #print('{:.2E}'.format(angCoeffMSD*Ang2PsTocm2S/6,2))
    anchored_text = AnchoredText(
        'D$_{tr}$ = (' + str('{:.2E}'.format(DiffCoeff)) + '\u00B1' + str('{:.2E}'.format(ErrDiffCoeff)) + ') cm$^2$/s',
        loc=2, prop=dict(size=18))
    ax.add_artist(anchored_text)

    # ax.text(50, 100, 'D = ' + str('{:.2E}'.format(DiffCoeff)) + '\u00B1' + str('{:.2E}'.format(ErrDiffCoeff)) + ' cm2/s' + '\n' + '(fit starts at '+str(estimatedStartFitTime)+' ps)', withdash=True)
    ax.errorbar(timeArrayHalf, msdS_Li, yerr=errMSD_Li, color='mediumblue', label='Li') #,label='Li (D ~ 1.5e-05cm$^2$/s)')

    #ax.errorbar(timeArrayHalf,msdS_P,yerr=errMSD_S, color='darkviolet',label='P')
    #ax.errorbar(timeArrayHalf,msdS_Ge,yerr=errMSD_S, color='violet',label='Ge')
    #ax.errorbar(timeArrayHalf,msdS_S,yerr=errMSD_S, color='navajowhite',label='S')
    ax.errorbar(timeArrScorrelated, msdSScorrelated, yerr=errMSDScorrelated, linestyle='-')
    ax.errorbar(timeArrayHalf, angCoeffMSD*timeArrayHalf+quoteMSD, linestyle='--')
    ax.errorbar(timeArrayHalf, (angCoeffMSD-np.sqrt(varangCoeffMSD)) * timeArrayHalf+quoteMSD, linestyle='--')
    ax.errorbar(timeArrayHalf, (angCoeffMSD+np.sqrt(varangCoeffMSD)) * timeArrayHalf+quoteMSD, linestyle='--')
    #print(angCoeffMSD)
    #print(angCoeffMSD + np.sqrt(varangCoeffMSD))
    #print(angCoeffMSD - np.sqrt(varangCoeffMSD))
    #plt.xlabel('time (ps)')
    #plt.ylabel('<MSD>(A2)')
    #plotMSD = plt.errorbar(timeArray[:msdS.shape[0]],msdS,yerr=errMSD)
    #plotMSD = plt.errorbar(timeArrScorrelated,msdSScorrelated,yerr=errMSDScorrelated)
    ax.legend(fontsize=12)
    #plt.savefig(fileMsdPngName,dpi=300, bbox_inches='tight', pad_inches=0)

    # to resize png borders
    #plt.savefig(fileMsdPngName,dpi=300)
    #img = Image.open(fileMsdPngName)
    #img_with_border = ImageOps.expand(img,border=300,fill='black')
    #img_with_border.save(fileMsdPngExpandedName)
    # ImageOps.expand(Image.open(fileMsdPngName),border=300,fill='black').save('fileMsdPngExpandedName', dpi=300)

    outMSD_file = open(fileMsdOutName, "w+")
    outDiff_file = open(fileDiffusionName, "a")
    outCond_file = open(fileConductivityName, "a")
    #outMSD_file.write('IONS CORRELATED'+ '\n')
    outMSD_file.write('timeStepJump = ' + str(timeStepJump) + '\n')
    outMSD_file.write('temperature = ' + str(temperature) + '\n')
    #outMSD_file.write('input pos file = ' + str(filePosName) + '\n')
    outMSD_file.write('time, cel and pos were cut before ' + str(timeArray[0]+tInitial)+ 'ps' + '\n')
    outMSD_file.write('t elapsed max is ' + str(estimatedTMax)+ 'ps' + '\n')
    outMSD_file.write('trajectory length is ' + str(timeArrayTmp[timeArrayTmp.shape[0]-1])+ 'ps' + '\n')
    outMSD_file.write('error on msd(t) evaluated at ' + str(estimatedFirstTElapsed) + 'ps' + ' and ' + str(estimatedSecondTElapsed) + 'ps' + '\n')
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
    outMSD_file.close()

    outDiff_file.write(str(temperature) + ' ' + str(DiffCoeff) + ' ' + str(ErrDiffCoeff) + ' ' +str(volumeArrayTmp[0])+ '\n')
    outDiff_file.close()

    outCond_file.write(str(temperature) + ' ' + str(Conduct) + ' ' + str(ErrConduct) + '\n')
    outCond_file.close()

    with open(fileMsdTOutName, "w+") as outMSDT_file:
        for i in range(indexTMax):
            outMSDT_file.write(str(timeArray[i])+ ' '  + str(msdS[i]) + '\n')
    #print(DiffCoeff,ErrDiffCoeff)


if __name__ == "__main__":
    main()

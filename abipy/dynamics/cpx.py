"""
Parser for the output files produced by CP code.

See: https://www.quantum-espresso.org/Doc/INPUT_CP.html

      prefix.pos : atomic positions
      prefix.vel : atomic velocities
      prefix.for : atomic forces
      prefix.cel : cell parameters
      prefix.str : stress tensors
      prefix.evp : energies
      prefix.hrs : Hirshfeld effective volumes (ts-vdw)
      prefix.eig : eigen values
      prefix.nos : Nose-Hoover variables
      prefix.spr : spread of Wannier orbitals
      prefix.wfc : center of Wannier orbitals
      prefix.ncg : number of Poisson CG steps (PBE0)
      prefix_ndw.save/ : write restart folder
      prefix_ndr.save/ : read restart folder
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import dataclasses
import abipy.core.abinit_units as abu

from pathlib import Path
from monty.functools import lazy_property
from monty.bisect import find_le
from abipy.core.mixins import TextFile, NotebookWriter
from abipy.tools.typing import PathLike, Figure
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    get_ax3d_fig_plt, rotate_ticklabels, set_visible, plot_unit_cell, set_ax_xylabels, get_figs_plotly,
    get_fig_plotly, add_plotly_fig_kwargs, PlotlyRowColDesc, plotly_klabels, plotly_set_lims)


def parse_file_with_blocks(filepath: PathLike, block_len: int) -> tuple[int, np.ndarray, list]:
    """
    Parse a QE file whose format is:

          20    0.00241888
         0.31882936800000E+01     0.14832370390000E+02     0.12288296100000E+01
         0.78323146900000E+01     0.67870403900000E+01     0.12288296100000E+01
         0.20744346700000E+01     0.59953799200000E+01     0.47375825000000E+01
          20    0.00241888
         0.31889894893205E+01     0.14833320999242E+02     0.12271083047604E+01
         0.78345550025202E+01     0.67881458009071E+01     0.12273445304341E+01
         0.20762325213327E+01     0.59955571558384E+01     0.47335647385293E+01

    i.e. a header followed by block_len lines

    Args:
        filepath: Filename
        block_len: Number of lines in block.

    Return: (number_of_steps, array_to_be_reshaped, list_of_headers)
    """
    nsteps, data, headers = 0, [], []
    with open(filepath, "rt") as fh:
        for il, line in enumerate(fh):
            toks = line.split()
            if il % (block_len + 1) == 0:
                #assert len(toks) == 2
                headers.append([int(toks[0]), float(toks[1])])
                nsteps += 1
            else:
                data.append([float(t) for t in toks])

    return nsteps, np.array(data, dtype=float), headers


def parse_file_with_header(filepath: PathLike) -> pd.DataFrame:
    """
    Parse a file with an header followed by csv-like columns e.g. an .evp file
    """
    from io import StringIO
    sbuf = StringIO()
    with open(filepath, "rt") as fh:
        for il, line in enumerate(fh):
            if il == 0:
                if not line.startswith("#"):
                    raise ValueError(f"Line number {il=} in {filepath=} should start with '#', but got {line=}")
                names = line[1:].split()
                continue

            if line.startswith("#"):
                # Found other header, likely due to restart.
                # Check that column names agree with the first header.
                new_names = line[1:].split()
                if new_names != names:
                    raise ValueError("Line number {il=} in {filepath=} has different names\n: {new_names=}\n{names=}")
                continue

            sbuf.write(line)

    sbuf.seek(0)
    df = pd.read_csv(sbuf, delim_whitespace=True, names=names)
    sbuf.close()
    return df


@dataclasses.dataclass
class Key:
    name: str
    info: str = "No info available"
    color: str = "b"


class EvpFile(TextFile, NotebookWriter):
    """
    The evp file is a file that contains the electronic and ionic energies and pressures
    for each time step of a CP simulation.
    The evp file has the following format:

    # nfi tps(ps) ekinc Tcell(K) Tion(K) etot enthal econs
          0            0          300   -2099.0164    7.4066063   -2091.6098    19363.353    2188.2406    5.0978969
          10         0.01    557.43004   -2105.3429    13.762216   -2091.5807    16917.626    2188.2406    5.0978969
    ...

    NB: Energies are in Hartree and not in Ry.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EvpFile
    """

    COLS_DICT = {k.name: k for k in [
         Key(name="nfi"    , color="b", info="number of iterations"),
         Key(name="ekinc"  , color="b", info="electron fake kinetic energy"),
         Key(name="temphc" , color="b", info="temperature due to “cell” kinetic energy in K"),
         Key(name="tempp"  , color="b", info="temperature due to ionic displacements within the cell"),
         Key(name="etot"   , color="b", info="dft total energy (without ionic and cell kinetic energy)"),
         Key(name="enthal" , color="b", info="etot+external_pressure*volume"),
         Key(name="econs"  , color="b", info="etot+kinetic_energy_due_to_ions_moving"),
         Key(name="econt"  , color="b", info="econs+ekinc+(thermostat_total_energy)"),
         Key(name="time(ps)", color="b", info="time"),
         Key(name="tps(ps)", color="b", info="time"),
         #Volume        Pressure(GPa)        EXX               EVDW
    ]}

    @lazy_property
    def time_key(self) -> str:
        key1 = "time(ps)"
        key2 = "tps(ps)"
        if key1 in self.df.keys(): return key1
        if key2 in self.df.keys(): return key2
        raise ValueError(f"Cannot find time key in {self.df.keys()}")

    @lazy_property
    def df(self) -> pd.DataFrame:
        """Dataframe with the results."""
        # TODO: Handle restart
        df = parse_file_with_header(self.filepath)
        return df

    @lazy_property
    def times(self) -> np.ndarray:
        """Array with times in ps units."""
        return np.array(self.df["tps(ps)"].values, dtype=float)

    def to_string(self, verbose=0) -> str:
        """String representation with verbosity level verbose."""
        lines = []; app = lines.append
        app(self.df.describe(percentiles=None).to_string())
        if verbose:
            for name, col in self.COLS_DICT.items():
                app(f"{name}: {col.info}")
        else:
            app("Use verbose to print the meaning of each column")
        return "\n".join(lines)

    @add_fig_kwargs
    def plot(self, **kwargs) -> Figure:
        """
        """
        ax_mat = self.df.plot.line(x=self.time_key, subplots=True,
                                   y=[k for k in self.df.keys() if k not in ("nfi", self.time_key)])
        return ax_mat.ravel()[0].get_figure()

        """
        # plot(self, ax_mat=None, t0=0.0, **kwargs) -> Figure:
        times = self.df["tps(ps)"]
        it0 = find_le(times, t0) if t0 != 0.0 else 0
        nrows, ncols = (1, 1)
        ax_mat, fig, plt = get_axarray_fig_plt(ax_array=ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        ax = ax_mat[0,0]

        yvals = self.df["ekinc"]
        # Plot ionic temperature.
        ax.plot(times[it0:], yvals[it0:], label="$Temp_{ions}$", c='tab:red')
        ax.set_xlabel('time (ps)')
        ax.set_ylabel('temperature (K)')
        material = "foobar"
        ax.text(0.35, 0.9, material,
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes, color='black', fontsize=15)
        ax.legend(loc=1, bbox_to_anchor=(1, 1))
        #for irow in range(nrows):
        #    for iax, (ax, data, label) in enumerate(zip(ax_mat[irow], select_irow[irow], label_irow[irow])):
        return fig
        """

    @add_fig_kwargs
    def plot_hist(self, **kwargs) -> Figure:
        """Plot one histogram per column."""
        ax_mat = self.df.hist(column=[k for k in self.df.keys() if k not in ("nfi", self.time_key)])
        return ax_mat.ravel()[0].get_figure()

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)
        yield self.plot_hist(show=False)

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("evp = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("evp.plot();"),
            nbv.new_code_cell("evp.plot_hist();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


def traj_to_qepos(traj_filepath: PathLike, pos_filepath: PathLike) -> None:
    """
    Convert an ASE trajectory file to QE POS file.

    Args:
        traj: Name of ASE trajectory file or Trajectory instance.
        pos_filepath: Name of output file.
    """
    from ase.io import read
    traj = read(traj_filepath, index=":")

    nsteps, natoms = len(traj), len(traj[0])
    pos_tac = np.empty((nsteps, natoms, 3))
    for it, atoms in enumerate(traj):
        pos_tac[it] = atoms.positions

    with open(str(pos_filepath), "wt") as fh:
        for it in range(nsteps):
            fh.write(str(it) + '\n')
            for ia in range(natoms):
                fh.write(str(pos_tac[it,ia,0]) + ' ' + str(pos_tac[it,ia,1]) + ' ' + str(pos_tac[it,ia,2]) + '\n')


class Qe2Extxyz:
    """
    Convert QE/CP output files into ASE extended xyz format.

    Example:

        from abipy.dynamics.cpx import Qe2Extxyz
        converter = Qe2Extxyz.from_input("cp.in", code="cp")
        converter.write_xyz("extended.xyz", take_every=1)
    """

    @classmethod
    def from_input(cls, qe_input, code, prefix=None):
        """
        Build object from QE/CP input assuming all the other output files
        are located in the same directory with the given prefix.
        """
        qe_input = Path(str(qe_input)).absolute()
        directory = qe_input.cwd()
        from ase.io.espresso import read_fortran_namelist
        with open(qe_input, "rt") as fh:
            sections, card_lines = read_fortran_namelist(fh)

        # https://www.quantum-espresso.org/Doc/INPUT_CP.html
        control = sections["control"]
        calculation = control["calculation"]
        #outdir = control["outdir"]

        default_prefix = "cp" if code == "cp" else "pwscf"
        if prefix is None:
            prefix = control.get("prefix", default_prefix)

        pos_filepath = directory / f"{prefix}.pos"
        for_filepath = directory / f"{prefix}.for"
        cell_filepath = directory / f"{prefix}.cel"
        evp_filepath = directory / f"{prefix}.evp"
        str_filepath = directory / f"{prefix}.str"
        if not str_filepath.exists():
            # File with stresses is optional.
            str_filepath = None

        return cls(code, qe_input, pos_filepath, for_filepath, cell_filepath, evp_filepath, str_filepath=str_filepath)

    def __init__(self, code, qe_input, pos_filepath, for_filepath, cell_filepath, evp_filepath, str_filepath=None):
        """
        Args:
            code: Either "pwscf" or "cp".
            qe_input: QE/CP input file
            pos_filepath: File with cartesian positions.
            for_filepath: File with cartesian forces.
            cell_filepath: File with lattice vectors.
            evp_filepath: File with energies
        """
        print(f"""
Reading symbols from: {qe_input=}
Reading positions from: {pos_filepath=}
Reading forces from: {for_filepath=}
Reading cell from: {cell_filepath=}
Reading energies from: {evp_filepath=}
Reading stresses from: {str_filepath=}
        """)

        # Parse input file to get initial_atoms and symbols.
        from ase.io.espresso import read_espresso_in
        with open(qe_input, "rt") as fh:
            self.initial_atoms = read_espresso_in(fh)
            natom = len(self.initial_atoms)
            #print("initial_atoms:", self.initial_atoms)

        # NB: in QE/CP lengths are in Bohr
        # but CP uses Hartree for energies whereas PW uses Rydberg.
        e_fact = 1.0 if code == "cp" else 0.5

        # Parse lattice vectors NB: in the cel file, lattice vectors are stored as column-vectors.
        # so he have to transpose the data as ASE uses row-vectors.
        cell_nsteps, cells, cell_headers = parse_file_with_blocks(cell_filepath, 3)
        cells = np.reshape(cells, (cell_nsteps, 3, 3)) * abu.Bohr_Ang
        self.cells = np.transpose(cells, (0, 2, 1)).copy()

        # Get energies from the evp file and convert from Ha to eV.
        with EvpFile(evp_filepath) as evp:
            self.energies = evp.df["etot"].values * (e_fact * abu.Ha_to_eV)
            if len(self.energies) != cell_nsteps:
                raise RuntimeError(f"{len(energies)=} != {cell_nsteps=}")

        # Parse Cartesian positions (Bohr --> Ang).
        pos_nsteps, pos_cart, pos_headers = parse_file_with_blocks(pos_filepath, natom)
        self.pos_cart = np.reshape(pos_cart, (pos_nsteps, natom, 3)) * abu.Bohr_Ang
        if pos_nsteps != cell_nsteps:
            raise RuntimeError(f"{pos_nsteps=} != {cell_nsteps=}")

        # Parse Cartesian forces (Ha/Bohr --> eV/Ang).
        for_nsteps, forces_list, for_headers = parse_file_with_blocks(for_filepath, natom)
        self.forces_step = np.reshape(forces_list, (for_nsteps, natom, 3)) * (e_fact * abu.Ha_to_eV / abu.Bohr_Ang)
        if for_nsteps != cell_nsteps:
            raise RuntimeError(f"{for_nsteps=} != {cell_nsteps=}")

        # Optionally, parse stress. In ASE, stress is eV/Ang^3.
        self.stresses = None
        if str_filepath is not None:
            str_nsteps, stresses, str_headers = parse_file_with_blocks(str_filepath, 3)
            self.stresses = np.reshape(stresses, (str_nsteps, 3, 3)) * (e_fact * abu.Ha_to_eV / abu.Bohr_Ang**3)
            if str_nsteps != cell_nsteps:
                raise RuntimeError(f"{str_nsteps=} != {cell_nsteps=}")

        self.nsteps = cell_nsteps

    def write_xyz(self, xyz_filepath, take_every=1, pbc=True) -> None:
        """
        Write results in ASE extended xyz format.

        Args:
            xyz_filepath: Name of the XYZ file.
            take_every: Used to downsample the trajectory.
            pbc: Values of pbc.
        """
        print(f"Writing results in extended xyz format in: {xyz_filepath}")
        from ase.io import write
        with open(xyz_filepath, "wt") as fh:
            for istep, atoms in enumerate(self.yield_atoms(pbc)):
                if istep % take_every != 0: continue
                write(fh, atoms, format='extxyz', append=True)

    def yield_atoms(self, pbc):
        """Yields ASE atoms along the trajectory."""
        from ase.atoms import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator
        for istep in range(self.nsteps):
            atoms = Atoms(symbols=self.initial_atoms.symbols,
                          positions=self.pos_cart[istep],
                          cell=self.cells[istep],
                          pbc=pbc,
                          )

            # Attach calculator with results.
            atoms.calc = SinglePointCalculator(atoms,
                                               energy=self.energies[istep],
                                               free_energy=self.energies[istep],
                                               forces=self.forces_step[istep],
                                               stress=self.stresses[istep] if self.stresses is not None else None,
                                               )
            yield atoms


def downsample_xyz(input_xyz, take_every, output_xyz, skip_head=None, verbose=1) -> int:
    """
    Downsample an XYZ file. Return the number of configurations written in the new file.

    Args:
        input_xyz: Input XYZ file.
        take_every: Extract configurations from input_xyz every `take_every` step.
        output_xyz: Output XYZ file.
        skip_head: If not None, skip the first skip_head configurations.
        verbose: Verbosity level.
    """
    from ase.io import write, iread
    count = 0
    with open(output_xyz, "wt") as out_xyz:
        for istep, atoms in enumerate(iread(input_xyz)):
            if istep % take_every != 0: continue
            if skip_head is not None and (istep + 1) <= skip_head: continue
            count += 1
            write(out_xyz, atoms, format='extxyz', append=True)

    if verbose: print(f"Wrote {count=} configurations to {output_xyz=} with {take_every=} and {skip_head=}")
    return count


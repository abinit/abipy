"""
Output files produced by CP code.
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

from monty.functools import lazy_property
from monty.bisect import find_le
from abipy.core.mixins import TextFile, NotebookWriter
from abipy.tools.typing import Figure
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    get_ax3d_fig_plt, rotate_ticklabels, set_visible, plot_unit_cell, set_ax_xylabels, get_figs_plotly,
    get_fig_plotly, add_plotly_fig_kwargs, PlotlyRowColDesc, plotly_klabels, plotly_set_lims)


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

    # nfi tps(ps) ekinc T cell(K) Tion(K) etot enthal econs
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
         Key(name="tps(ps)", color="b", info="time"),
         #Volume        Pressure(GPa)        EXX               EVDW
    ]}

    @lazy_property
    def df(self) -> pd.DataFrame:
        """Dataframe with the results."""
        # Read colums from first row.
        header = pd.read_csv(self.filepath, nrows=0).columns[0]
        names = header[1:].split()
        # Build dataframe from the remaining rows and set column names
        df = pd.read_csv(self.filepath, delim_whitespace=True, skiprows=1, names=names)
        return df

    @lazy_property
    def times(self):
        """Array with times in ps units."""
        return self.df["tps(ps)"].values

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
        ax_mat = self.df.plot.line(x='tps(ps)', subplots=True,
                                   y=[k for k in self.df.keys() if k not in ("nfi", "tps(ps)")])
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
        ax_mat = self.df.hist(column=[k for k in self.df.keys() if k not in ("nfi", "tps(ps)")])
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

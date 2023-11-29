# coding: utf-8
"""Classes to analyse LRUJ results."""
from __future__ import annotations

import dataclasses
import numpy as np
import pandas as pd
from pathlib import Path
#from typing import Any
#from monty.string import is_string, list_strings, marquee
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.iotools import ETSF_Reader
from abipy.tools.iotools import yaml_safe_load
from abipy.core.structure import Structure
from abipy.tools.typing import Figure, PathLike
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    get_ax3d_fig_plt, rotate_ticklabels, set_visible, plot_unit_cell, set_ax_xylabels, get_figs_plotly)



#class LrujFile(AbinitNcFile, Has_Header, Has_Structure): #, Has_ElectronBands, NotebookWriter):
#    """
#    File containing the results of a ground-state calculation.
#
#    Usage example:
#
#    .. code-block:: python
#
#        with GsrFile("foo_GSR.nc") as gsr:
#            print("energy: ", gsr.energy)
#            gsr.ebands.plot()
#
#    .. rubric:: Inheritance Diagram
#    .. inheritance-diagram:: LrujFile
#    """
#
#    @classmethod
#    def from_file(cls, filepath: str) -> GsrFile:
#        """Initialize the object from a netcdf_ file."""
#        return cls(filepath)
#
#    def __init__(self, filepath: str):
#        super().__init__(filepath)
#        self.r = r = EtsfReader(filepath)



@dataclasses.dataclass(kw_only=True)
class LrujResults:
    """
    This object stores the results produced by lruj.
    """
    npert: int
    maxdeg: int
    pawujat: int
    macro_uj: int
    dmatpuopt: int
    diem: float
    alphas: np.ndarray
    occ_unscr: np.ndarray
    occ_scr: np.ndarray
    chi0_coefficients: dict
    chi_coefficients: dict

    @classmethod
    def from_file(cls, filepath: PathLike):
        """
        Extract results from the main ouput file produced by lruj.
        """
        with open(filepath, "rt") as fh:
            lines = [line.lstrip() for line in fh]

        # Extract the Yaml document with the chi/chi0 coefficients
        in_doc = False
        yaml_lines = []
        for i, line in enumerate(lines):

            if line.startswith("--- !LRUJ_Abipy_Plots"):
                in_doc = True
                continue

            if in_doc and line.startswith("..."):
                data = yaml_safe_load("".join(yaml_lines))
                print("data:", data)
                break

            if in_doc:
                yaml_lines.append(line)

        chi0_coefficients = {}
        chi_coefficients = {}
        for k, v in data.items():
            magic = "chi0_coefficients_degree"
            if k.startswith(magic):
                degree = int(k.replace(magic, ""))
                chi0_coefficients[degree] = v
            magic = "chi_coefficients_degree"
            if k.startswith(magic):
                degree = int(k.replace(magic, ""))
                chi_coefficients[degree] = v

        #print(f"{chi0_coefficients=}")
        #print(f"{chi_coefficients=}")

        def find(header, dtype=None):
            for i, line in enumerate(lines):
                if line.startswith(header):
                    after = line.replace(header, "", 1).strip()
                    if dtype: after = dtype(after)
                    return i, after
            raise ValueError(f"Cannot find {header=} in {filepath=}")

        _, npert = find("Number of perturbations detected:", dtype=int)
        _, maxdeg = find("Maximum degree of polynomials analyzed:", dtype=int)
        _, pawujat = find("Index of perturbed atom:", dtype=int)
        _, macro_uj = find("Value of macro_uj:", dtype=int)
        _, dmatpuopt = find("Value of dmatpuopt:", dtype=int)
        _, diem = find("Mixing constant factored out of Chi0:", dtype=float)

        # Parse the section with perturbations and occupations.
        """
        Perturbations           Occupations
        -------------- -----------------------------
          alpha [eV]     Unscreened      Screened
        -------------- -----------------------------
         0.0000000000   8.6380182458   8.6380182458
        -0.1500000676   8.6964981922   8.6520722003

        """
        i, _ = find("alpha [eV]     Unscreened      Screened")
        i += 2
        vals = []
        for ipert in range(npert):
            vals.append([float(t) for t in lines[i+ipert].split()])
        vals = np.reshape(vals, (npert, 3))
        alphas, occ_unscr, occ_scr = vals[:,0], vals[:,1], vals[:,2]
        """
                                                                               RMS Errors
                                                                 ---------------------------------------
         Regression   Chi0 [eV^-1]   Chi [eV^-1]      U [eV]    | Chi0 [eV^-1]  Chi [eV^-1]     U [eV]
        --------------------------------------------------------|---------------------------------------
         Linear:        -0.8594082    -0.0949434     9.3689952  |    0.0023925    0.0000878    0.1139297
         Quadratic:     -0.8574665    -0.0955791     9.2963099  |    0.0023777    0.0000129    0.0681722
         Cubic:         -0.8007858    -0.0952726     9.2474220  |    0.0001546    0.0000015    0.0200543
        *************************************************************************************************
        """
        i, _ = find("Regression   Chi0 [eV^-1]")
        i += 2
        keys = ["Chi0", "Chi", "U", "rms_Chi0", "rms_Chi", "rms_U"]
        dict_list = []
        for irow in range(maxdeg):
            l = lines[i+irow].replace("|", " ")
            tokens = l.split()
            d = dict(zip(keys, [float(t) for t in tokens[1:]]))
            d["degree"] = irow + 1
            dict_list.append(d)

        fit_df = pd.DataFrame(dict_list)
        #print("fit_df:\n", fit_df)

        # Build instance from locals dict.
        data = locals()
        return cls(**{k: data[k] for k in [field.name for field in dataclasses.fields(cls)]})

    @add_fig_kwargs
    def plot(self, ax=None, fontsize=12, **kwargs) -> Figure:
        """
        Plot

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Plot raw data.
        ax.scatter(self.alphas, self.occ_unscr, marker="o", c="b", label="Unscreened")
        ax.scatter(self.alphas, self.occ_scr, marker="x", c="r", label="Screened")
        ax.axvline(x=0.0, color="k", linestyle="--", lw=0.5)

        # Plot regression fit (only linear part)
        xstart, xstop = self.alphas.min(), self.alphas.max()
        xstart, xstop = 0.9 * xstart, 1.1 * xstop
        xs = np.arange(xstart, xstop, step=0.01)
        #for ideg in range(maxdeg):
        #    ax.plot(xs, lin_coeff * xs, color=, label="Degree {ideg}")

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        ax.grid(True)
        ax.set_xlabel(r"$\alpha$ (eV)")
        #ylabel = r"$N(n^{\up} + n^{\down})" is U else r"$N(n^{\up} - n^{\down})"
        #ax.set_ylabel(ylabel)

        return fig


class LrujAnalyzer:
    """
    Analyzes multiple sets of LRUJ files.
    """
    def __init__(self, manager=None, verbose=0):
        self.ncfiles_of_key = {}
        self.results_of_key = {}
        self.manager = manager
        self.verbose = verbose

    #def __str__(self):
    #    return self.to_string()

    #def to_string(self, verbose: int = 0) -> str:
    #    lines = []
    #    app = lines.append
    #    return "\n".join(lines)

    #@classmethod
    #def from_dir(cls, key: str, directory: PathLike) -> None:
    #    new = cls()
    #    new.scan_dir(key, directory)
    #    return new

    #def walk_dir(self, key: str, top: PathLike) -> None:
    #    nc_paths = None
    #    self.add_ncpaths(key, nc_paths)

    def add_ncpaths(self, key: str, nc_paths: list[str]) -> None:
        self.ncfiles_of_key[key] = nc_paths
        self.results_of_key[key] = None
        self.run()

    def run(self, workdir=None) -> None:
        """
        Invoke lruj
        """
        from abipy.flowtk.wrappers import Lruj
        from abipy.core.globals import get_workdir
        workdir = Path(get_workdir(workdir))
        for key, nc_paths in self.ncfiles_of_key.items():
            if self.results_of_key[key] is not None: continue
            wdir = workdir / f"run_{key}"
            wdir.mkdir()
            Lruj(manager=self.manager, verbose=self.verbose).run(nc_paths, workdir=wdir)
            self.results_of_key[key] = LrujResults.from_file(wdir / "lruj.stdout")

    @add_fig_kwargs
    def plot(self, **kwargs) -> Figure:
        """
        Plot

        Args:
        """
        keys = self.results_of_key.keys()

        # Build grid of plots.
        num_plots, ncols, nrows = len(keys), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for ix, (key, ax) in enumerate(zip(keys, ax_list)):
            res = self.results_of_key[key]
            res.plot(ax=ax, show=False, title=key)

        return fig


#class LrujInputGenerator
#
#    def __init__(self, gs_input, elements, site_inds=None, dmatpuopt=3):
#        self.dmatpuopt = dmatpuopt
#        self.gs_input = gs_input.copy()
#        self.gs_input.set_vars(
#            chkprim=0,            # Will complain otherwise with AFII magnetic state
#            chksymbreak=0,        # Don't check for symmetry breaking
#            # DFT+U
#            usepawu=1,            # Alert Abinit to use of DFT+U
#            lpawu=[2, 2, 1],      # Subspaces of interest: Ni 3d, O 2p
#            upawu="*0.0",         # Raw (non-corrected) XC functional required to establish U(J)
#            jpawu="*0.0",
#            dmatpuopt=self.dmatpuopt, # PAW projector for density matrix
#        )
#
#    def gen_inputs(self, macro_uj, alphas, supercells):
#        #supercells = np.respahe(supercells, (-1,
#        for supercell in supercells:
#            for alpha in alphas:
#                if alpha == 0.0: continue
#                gs_input.new_with_vars(
#                    macro_uj=macro_uj,
#                    pawujat=1,              # Label of perturbed atom
#                    pawujv1=f"{alpha} eV",  # List of perturbation strengths
#                    prtwf=0,
#                    prtebands 0             # Don't print ebands
#                )

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


#===============================================================================================================
#===============================================================================================================
@dataclasses.dataclass(kw_only=True)
class LrujResults:
    """
    This object stores the results produced by lruj.
    """
    natom: int
    npert: int
    ndata: int
    pawujat: int
    macro_uj: int
    diem_token: str
    diem: float
    alphas: np.ndarray
    occ_unscr: np.ndarray
    occ_scr: np.ndarray
    chi0_coefficients: dict
    chi_coefficients: dict
    maxdeg: int
    dmatpuopt: int
    pert_name: str
    parname: str
    metric: str
    fit_df: pd.DataFrame

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
                break

            if in_doc:
                yaml_lines.append(line)

        natom = data['natom']
        ndata = data['ndata']
        pawujat = data['pawujat']
        macro_uj = data['macro_uj']
        diem_token = data['diem_token']
        diem = data['diem']
        npert = ndata - 1
        if macro_uj==4:
          pert_name = r"$\beta$"
          metric = r"M $(n^{\uparrow} - n^{\downarrow})$"
          parname = 'J'
        else:
          pert_name = r"$\alpha$"
          metric = r"N $(n^{\uparrow} + n^{\downarrow})$"
          parname = 'U'

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

        #print(f"{chi_coefficients=}")

        def find(header, dtype=None):
            for i, line in enumerate(lines):
                if line.startswith(header):
                    after = line.replace(header, "", 1).strip()
                    if dtype: after = dtype(after)
                    return i, after
            raise ValueError(f"Cannot find {header=} in {filepath=}")

        _, maxdeg = find("Maximum degree of polynomials analyzed:", dtype=int)
        _, dmatpuopt = find("Value of dmatpuopt:", dtype=int)

        # Parse the section with perturbations and occupations.
        """
        Perturbations           Occupations
        -------------- -----------------------------
          alpha [eV]     Unscreened      Screened
        -------------- -----------------------------
         0.0000000000   8.6380182458   8.6380182458
        -0.1500000676   8.6964981922   8.6520722003

       -OR-

        Perturbations         Magnetizations
        --------------- -----------------------------
           beta [eV]     Unscreened      Screened
        --------------- -----------------------------

        """
        i, _ = find("Perturbations",dtype=None)
        i += 4
        vals = []
        for ipert in range(ndata):
            vals.append([float(t) for t in lines[i+ipert].split()])
        vals = np.reshape(vals, (ndata, 3))
        alphas, occ_unscr, occ_scr = vals[:,0], vals[:,1], vals[:,2]
        """
                                                                               RMS Errors
                                                                 ---------------------------------------
         Regression   Chi0 [eV^-1]   Chi [eV^-1]      J [eV]    | Chi0 [eV^-1]  Chi [eV^-1]     J [eV]
        --------------------------------------------------------|---------------------------------------
         Linear:        -0.8594082    -0.0949434     9.3689952  |    0.0023925    0.0000878    0.1139297
         Quadratic:     -0.8574665    -0.0955791     9.2963099  |    0.0023777    0.0000129    0.0681722
         Cubic:         -0.8007858    -0.0952726     9.2474220  |    0.0001546    0.0000015    0.0200543
        *************************************************************************************************
        """
        i, _ = find("Regression   Chi0 [eV^-1]")
        i += 2
        keys = ["Chi0", "Chi", "HP", "rms_Chi0", "rms_Chi", "rms_HP"]
        dict_list = []
        for irow in range(maxdeg):
            l = lines[i+irow].replace("|", " ")
            tokens = l.split()
            d = dict(zip(keys, [float(t) for t in tokens[-6:]]))
            d["degree"] = irow + 1
            dict_list.append(d)

        fit_df = pd.DataFrame(dict_list)

        # Build instance from locals dict.
        _data = locals()
        return cls(**{k: _data[k] for k in [field.name for field in dataclasses.fields(cls)]})

#===============================================================================================================
#===============================================================================================================
    @add_fig_kwargs
    def plot(self, ax=None, degrees="all", inset=True, insetdegree=1, insetlocale="lower left",
             ptcolor0='k', ptcolor='k', gradcolor1='#3575D5',gradcolor2='#FDAE7B',
             ptitle="default", fontsize=12, **kwargs) -> Figure:
        """
        Plot

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            degrees: List of degrees to plot. If None, no polynomials are plotted.
            ptcolor0: Color of unscreened response data points (default: black)
            ptcolor: Color of screened response data points (default: black)
            gradcolor1: Hex code of linear regression color (default: Blue #3575D5)
            gradcolor2: Hex code of color of regression of maximum degree in list <degrees> (default: Salmon #FDAE7B)
            ptitle: Plot title (default: "Linear Response for atom <pawujat>)
            inset: Plots inset with LR parameter for polynomial fit of degree <insetdegree> (default: True)
            insetdegree: Degree of polynomial fit information printed in the inset (default: 1)
            insetlocale: Position of inset in the plot. Standard matplotlob locations. (default: "lower left")
        """
        import seaborn as sns
        sns.set()
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Plot data
        yshift = self.occ_unscr[np.where(self.alphas == 0.0000)] * (self.diem - 1.0)
        data0 = 1.0/self.diem * (self.occ_unscr + yshift)
        ax.scatter(self.alphas, data0, s=70, color=ptcolor0, facecolors='none', linewidths=2, label="Unscreened")
        ax.scatter(self.alphas, self.occ_scr, s=70, color=ptcolor, label="Screened")
        ax.axvline(x=0.0, color="white", linestyle="--", lw=0.5)

        # Generate mesh for polynomial functions
        xstart, xstop = 1.1 * self.alphas.min(), 1.1 * self.alphas.max()
        xs = np.arange(xstart, xstop, step=0.01)

        # Prepare colors and coefficients for polynomials the use wants to plot
        if degrees == "all":
          degrees = self.chi0_coefficients.keys()

        def hex_to_RGB(hex_str):
          return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

        def get_color_gradient(c1, c2, n):
          assert n > 1
          c1_rgb = np.array(hex_to_RGB(c1))/255
          c2_rgb = np.array(hex_to_RGB(c2))/255
          mix_pcts = [x/(n-1) for x in range(n)]
          rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
          return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]

        linecolors=get_color_gradient(gradcolor1,gradcolor2,len(degrees))

        # Plot polynomial functions
        def poly0(coeffs):
          return lambda x: 1.0/self.diem * (sum((coeff*x**i for i,coeff in enumerate(coeffs))) + yshift)

        def poly(coeffs):
          return lambda x: sum((coeff*x**i for i,coeff in enumerate(coeffs)))

        for degree in degrees:
          polynomial0=poly0(self.chi0_coefficients[degree])
          polynomial=poly(self.chi_coefficients[degree])
          if degree == 1:
            Labelstring='Linear'
          elif degree == 2:
            Labelstring='Quadratic'
          elif degree == 3:
            Labelstring='Cubic'
          else:
            Labelstring=' '.join(['Degree',str(degree)])

          if insetdegree==degree:
            deginfo = ' '.join(['Parameters for',Labelstring,'fit'])
            insetcolor = linecolors[degree-1]

          ax.plot(xs,polynomial0(xs),color=linecolors[degree-1],linewidth=2.0,linestyle='dashed')
          ax.plot(xs,polynomial(xs),color=linecolors[degree-1],linewidth=2.0,label=Labelstring)

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        if ptitle=="default":
          ptitle=' '.join(["Linear Response",self.parname,"on atom",str(self.pawujat)])
        plt.title(ptitle)
        ax.grid(True)
        ax.set_xlabel(' '.join([self.pert_name,"(eV)"]))
        ax.set_ylabel(self.metric)

        # Generate inset with numerical information on LR
        if inset:
          from matplotlib.offsetbox import AnchoredText
          select =  self.fit_df["degree"] == insetdegree
          def dfvalue(keywo):
            tablevalue = self.fit_df[select][keywo]
            return "%.4f" % tablevalue.values[0]

          X0str = ' '.join([r"$\chi_0$","=",dfvalue("Chi0"),r"$\pm$",dfvalue("rms_Chi0"),r"[eV]$^{-1}$"])
          Xstr = ' '.join([r"$\chi$","=",dfvalue("Chi"),r"$\pm$",dfvalue("rms_Chi"),r"[eV]$^{-1}$"])
          HPstr = ' '.join([self.parname,"=",dfvalue("HP"),r"$\pm$",dfvalue("rms_HP"),r"[eV]"])
          insettxt = '\n'.join([deginfo,X0str,Xstr,HPstr])
          parambox = AnchoredText(insettxt,loc=insetlocale)
          parambox.patch.set_linewidth(2)
          parambox.patch.set_edgecolor(insetcolor)
          parambox.patch.set_facecolor('white')
          ax.add_artist(parambox)

        return fig



#===============================================================================================================
#===============================================================================================================
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

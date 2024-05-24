# coding: utf-8
"""
Objects to analyze and visualize the results of GWR calculations.
"""
from __future__ import annotations

import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from collections.abc import Iterable
from typing import Union, Any
from monty.functools import lazy_property
from monty.string import list_strings, marquee
from monty.termcolor import cprint
from abipy.core.structure import Structure
from abipy.core.kpoints import Kpoint, KpointList, Kpath, IrredZone, has_timrev_from_kptopt
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.iotools import ETSF_Reader
from abipy.tools import duck
from abipy.tools.typing import Figure, KptSelect
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, Marker,
    set_axlims, set_ax_xylabels, set_visible, rotate_ticklabels, set_grid_legend, hspan_ax_line, Exposer)
from abipy.tools.numtools import data_from_cplx_mode
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.electrons.gw import SelfEnergy, QPState, QPList #, SigresFile, SigresRobot
from abipy.abio.robots import Robot
from abipy.abio.enums import GWR_TASK

__all__ = [
    "GwrFile",
    "GwrRobot",
    "TchimVsSus",
]


class _MyQpkindsList(list):
    """Returned by find_qpkinds."""


@dataclasses.dataclass(kw_only=True)
class MinimaxMesh:
    """
    The minimax mesh reported in the GWR file.
    """
    ntau: int              # Number of points.
    tau_mesh: np.ndarray   # tau points.
    tau_wgs: np.ndarray    # tau weights for integration.
    iw_mesh: np.ndarray    # omega points along the imag. axis.
    iw_wgs: np.ndarray     # omega weights for integration.
    cosft_wt: np.ndarray   # weights for cosine transform (tau --> omega).
    cosft_tw: np.ndarray   # weights for cosine transform (omega --> tau).
    sinft_wt: np.ndarray   # weights for sine transform (tau --> omega).

    min_transition_energy_eV: float   # Minimum transition energy.
    max_transition_energy_eV: float   # Maximum transition energy.
    eratio: float
    ft_max_err_t2w_cos: float
    ft_max_err_w2t_cos: float
    ft_max_err_t2w_sin: float
    cosft_duality_error: float
    regterm: float                    # Regularization term.

    @classmethod
    def from_ncreader(cls, reader: ETSF_Reader) -> MinimaxMesh:
        """
        Build minimax mesh from a netcdf reader.
        """
        d = {}
        for field in dataclasses.fields(cls):
            name = field.name
            if name == "ntau":
                value = reader.read_dimvalue(name)
            else:
                value = reader.read_value(name)
                if field.type == "float":
                    value = float(value)
            d[name] = value

        return cls(**d)

    @add_fig_kwargs
    def plot_ft_weights(self, other: MinimaxMesh, self_name="self", other_name="other",
                        with_sinft=False, fontsize=6, **kwargs):
        """
        Plot the Fourier transformt weights of two minimax meshes (self and other)
        """
        if self.ntau != other.ntau:
            raise ValueError("Cannot compare minimax meshes with different ntau")

        import matplotlib.pyplot as plt
        nrows, ncols = (4, 2) if with_sinft else (3, 2)
        fig, ax_mat = plt.subplots(nrows=nrows, ncols=ncols,
                                   sharex=False, sharey=False, squeeze=False,
                                   figsize=(12, 8),
                                   #subplot_kw={'xticks': [], 'yticks': []},
        )

        I_mat = np.eye(self.ntau)
        select_irow = {
            0: [(self.cosft_wt @ self.cosft_tw) - I_mat,
                (other.cosft_wt @ other.cosft_tw) - I_mat], # , other.cosft_wt @ other.cosft_tw],
            1: [self.cosft_wt, other.cosft_wt], # self.cosft_wt - other.cosft_wt],
            2: [self.cosft_tw, other.cosft_tw], # self.cosft_tw - other.cosft_tw],
            3: [self.sinft_tw, other.sinft_tw], # self.sinft_tw - other.sinft_tw],
        }

        label_irow = {
            0: [f"(cosft_wt @ cosft_tw) - I ({self_name})", f"(cosft_wt @ cosft_tw) - I ({other_name})"],
            1: [f"cosft_wt ({self_name})", f"cosft_wt ({other_name})"],
            2: [f"cosft_tw ({self_name})", f"cosft_tw ({other_name})"],
            3: [f"sinft_tw ({self_name})", f"sinft_tw ({other_name})"],
        }

        for irow in range(nrows):
            for iax, (ax, data, label) in enumerate(zip(ax_mat[irow], select_irow[irow], label_irow[irow])):
                im = ax.matshow(data, cmap='seismic')
                fig.colorbar(im, ax=ax)
                ax.set_title(label, fontsize=fontsize)

        return fig


class GwrFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This object provides an high-level interface to the GWR.nc file produced by the GWR code.
    """

    # Markers used for up/down bands.
    marker_spin = {0: "^", 1: "v"}

    color_spin = {0: "k", 1: "r"}

    @classmethod
    def from_file(cls, filepath: str) -> GwrFile:
        """Initialize an instance from filepath."""
        return cls(filepath)

    def __init__(self, filepath: str):
        """Read data from the netcdf filepath."""
        super().__init__(filepath)
        self.r = self.reader = GwrReader(self.filepath)

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.r.structure

    @property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| with the KS energies."""
        return self.r.bands

    @property
    def sigma_kpoints(self) -> KpointList:
        """The k-points where QP corrections have been calculated."""
        return self.r.sigma_kpoints

    @property
    def nkcalc(self) -> int:
        """Number of k-points in sigma_kpoints"""
        return self.r.nkcalc

    @lazy_property
    def ks_dirgaps(self) -> np.ndarray:
       """KS direct gaps in eV. Shape: [nsppol, nkcalc]"""
       return self.r.read_value("ks_gaps") * abu.Ha_eV

    @lazy_property
    def qpz0_dirgaps(self) -> np.ndarray:
       """
       QP direct gaps in eV computed with the Z factor at the KS energy
       Shape: [nsppol, nkcalc]
       """
       return self.r.read_value("qpz_gaps") * abu.Ha_eV

    @lazy_property
    def minimax_mesh(self) -> MinimaxMesh:
        """Object storing minimax mesh and weights."""
        return MinimaxMesh.from_ncreader(self.r)

    @lazy_property
    def qplist_spin(self) -> tuple[QPList]:
        """Tuple of QPList objects indexed by spin."""
        return self.r.read_allqps()

    def find_qpkinds(self, qp_kpoints) -> _MyQpkindsList:
        """
        Find kpoints for QP corrections from user input.
        Return list of (kpt, ikcalc) tuples where kpt is a |Kpoint| and
        ikcalc is the index in the nkcalc array..
        """
        if isinstance(qp_kpoints, _MyQpkindsList):
            return qp_kpoints

        if isinstance(qp_kpoints, Kpoint):
            qp_kpoints = [qp_kpoints]

        if qp_kpoints is None or (duck.is_string(qp_kpoints) and qp_kpoints == "all"):
            # qp_kpoints in (None, "all")
            items = self.sigma_kpoints, list(range(self.nkcalc))

        elif duck.is_intlike(qp_kpoints):
            # qp_kpoints = 1
            ikc = int(qp_kpoints)
            items = [self.sigma_kpoints[ikc]], [ikc]

        elif isinstance(qp_kpoints, Iterable):
            # either [0, 1] or [[0, 0.5, 0]]
            # note possible ambiguity with [0, 0, 0] that is treated as list of integers.
            if duck.is_intlike(qp_kpoints[0]):
                ik_list = duck.list_ints(qp_kpoints)
                items = [self.sigma_kpoints[ikc] for ikc in ik_list], ik_list
            else:
                ik_list = [self.r.kpt2ikcalc(kpt) for kpt in qp_kpoints]
                qp_kpoints = [self.sigma_kpoints[ikc] for ikc in ik_list]
                items = qp_kpoints, ik_list
        else:
            raise TypeError("Don't know how to interpret `%s`" % (type(qp_kpoints)))

        # Check indices
        errors = []
        eapp = errors.append
        for ikc in items[1]:
            if ikc >= self.nkcalc:
                eapp("K-point index %d >= nkcalc %d, check input qp_kpoints" % (ikc, self.nkcalc))
        if errors:
            raise ValueError("\n".join(errors))

        return _MyQpkindsList(zip(items[0], items[1]))

    @lazy_property
    def params(self) -> dict:
        """
        dict with parameters that might be subject to convergence studies e.g ecuteps.
        """
        minimax_mesh = self.minimax_mesh
        r = self.r
        return dict(
            gwr_ntau=r.read_dimvalue("ntau"),
            nband=self.ebands.nband,
            ecuteps=r.read_value("ecuteps"),
            ecutsigx=r.read_value("ecutsigx"),
            ecut=r.read_value("ecut"),
            gwr_boxcutmin=r.read_value("gwr_boxcutmin"),
            nkpt=self.ebands.nkpt,
            symchi=r.read_value("symchi"),
            symsigma=r.read_value("symsigma"),
            regterm=minimax_mesh.regterm,
        )

    def close(self) -> None:
        """Close the netcdf file."""
        self.r.close()

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level ``verbose``."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(title="KS Electron Bands", with_structure=False))
        app("")

        # GWR section.

        app(marquee("GWR parameters", mark="="))
        app("gwr_task: %s" % self.r.gwr_task)

        if self.r.gwr_task == GWR_TASK.RPA_ENERGY:
            pass
            #d = self.get_rpa_ene_dict()

        else:
            app("Number of k-points in Sigma_{nk}: %d" % (len(self.r.sigma_kpoints)))
            app("Number of bands included in e-e self-energy sum: %d" % (self.nband))
            keys = self.params.keys() if verbose else ["ecuteps", "ecutsigx", "ecut", "gwr_boxcutmin"]
            for k in keys:
                app("%s: %s" % (k, self.params[k]))
            if verbose:
                app(marquee("k-points in Sigma_nk", mark="="))
                app(self.r.sigma_kpoints.to_string(title=None, pre_string="\t"))
            app("")

            app(marquee("QP direct gaps in eV", mark="="))
            app(str(self.get_dirgaps_dataframe(with_params=False)))
            app("")

        #app(str(self.minimax_mesh))

        return "\n".join(lines)

    #def get_qpgap(self, spin, kpoint, with_ksgap=False):
    #    """
    #    Return the QP gap in eV at the given (spin, kpoint)
    #    """
    #    ik = self.reader.kpt2fileindex(kpoint)
    #    if not with_ksgap:
    #        return self.qpgaps[spin, ik]
    #    else:
    #        return self.qpgaps[spin, ik], self.ksgaps[spin, ik]

    def get_dirgaps_dataframe(self, with_params: bool=True, with_geo: bool=False) -> pd.DataFrame:
        """
        Return a pandas DataFrame with the QP direct gaps in eV.

        Args:
            with_params: True if GWR parameters should be included.
            with_geo: True if geometry info should be included.
        """
        d = {}
        d["kpoint"] = [k.frac_coords for k in self.sigma_kpoints] * self.nsppol
        d["kname"] = [k.name for k in self.sigma_kpoints] * self.nsppol
        d["ks_dirgaps"] = self.ks_dirgaps.ravel()
        d["qpz0_dirgaps"] = self.qpz0_dirgaps.ravel()
        #d["qp_pade_dirgaps"] = self.qp_pade_dirgaps.ravel()
        d["spin"] = [0] * len(self.sigma_kpoints)
        if self.nsppol == 2: d["spin"].extend([1] * len(self.sigma_kpoints))

        if with_params:
            for k, v in self.params.items():
                 d[k] = [v] * len(self.sigma_kpoints) * self.nsppol

        if with_geo:
            d.update(**self.structure.get_dict4pandas(with_spglib=True))

        return pd.DataFrame(d)

    def get_dataframe_sk(self, spin: int, kpoint: KptSelect, index=None,
                         ignore_imag=False, with_params=True, with_geo=False) -> pd.Dataframe:
        """
        Returns |pandas-DataFrame| with the QP results for the given (spin, k-point).

        Args:
            spin:
            kpoint:
            index:
            ignore_imag: Only real part is returned if ``ignore_imag``.
            with_params: True if GWR parameters should be included.
            with_geo: True if geometry info should be included.
        """
        rows, bands = [], []

        if with_geo:
            geo_dict = self.structure.get_dict4pandas(with_spglib=True)

        qp_list = self.r.read_qplist_sk(spin, kpoint, ignore_imag=ignore_imag)
        for qp in qp_list:
            bands.append(qp.band)
            d = qp.as_dict()

            # Add other entries that may be useful when comparing different calculations.
            if with_params:
                d.update(self.params)
            if with_geo:
                d.update(**geo_dict)

            rows.append(d)

        index = len(bands) * [index] if index is not None else bands
        return pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))

    def _get_include_bands(self, include_bands: Any, spin: int) -> Union[set, None]:
        """
        Helper function to return include_bands for given spin.
        """
        if include_bands is None:
            return None

        if duck.is_listlike(include_bands) or isinstance(include_bands, set):
            return set(include_bands)

        if duck.is_string(include_bands):
            if include_bands == "all":
                return None
            if include_bands == "gap":
                lumo_band = self.ebands.lumos[spin].band
                homo_band = self.ebands.homos[spin].band
                return set(range(homo_band, lumo_band + 1))

        raise ValueError(f"Invalid value for include_bands: {include_bands}")

    def get_rpa_ene_dict(self) -> dict:
        """
        Read RPA energies computed for different ecut_chi (all in Ha)
        """
        # nctkarr_t("ecut_chi", "dp", "ncut")
        # nctkarr_t("ec_rpa_ecut", "dp", "ncut")
        d = {}
        for key in ("ecut_chi", "ec_rpa_ecut", "ec_mp2_ecut"):
            d[key] = self.r.read_value(key)

        # Extrapolate value for ecut --> oo.
        xs = d["ecut_chi"] ** (-3/2.0)
        for key in ("ec_rpa_ecut", "ec_mp2_ecut"):
            ys = d[key]
            coef = np.polyfit(xs, ys, 1)
            # poly1d_fn is a function which takes in x and returns an estimate for y
            poly1d_fn = np.poly1d(coef)
            extrap_value = poly1d_fn(0.0)
            #print(f"{extrap_value=}")
            d[key + "_inf"] = extrap_value

        #print(d)
        return d

    @add_fig_kwargs
    def plot_sigma_imag_axis(self, kpoint: KptSelect,
                             spin=0, include_bands="gap", with_tau=True,
                             fontsize=8, ax_mat=None, **kwargs) -> Figure:
        """
        Plot Sigma_nk(iw) along the imaginary axis for given k-point, spin and list of bands.

        Args:
            kpoint: k-point in self-energy. Accepts |Kpoint|, vector or index.
            spin: Spin index.
            include_bands: List of bands to include. None means all.
            with_tau:
            fontsize: Legend and title fontsize.
            ax_mat:

        Returns: |matplotlib-Figure|
        """
        nrows, ncols = (2, 2) if with_tau else (1, 2)
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=not with_tau,
                                               sharey=False, squeeze=False)
        ax_mat = np.array(ax_mat)

        # Read Sigma_nk in sigma_of_band
        ikcalc, kpoint = self.r.get_ikcalc_kpoint(kpoint)
        include_bands = self._get_include_bands(include_bands, spin)
        sigma_of_band = self.r.read_sigma_bdict_sikcalc(spin, ikcalc, include_bands)

        # Plot Sigmac_nk(iw)
        re_ax, im_ax = ax_list = ax_mat[0]
        style = dict(marker="o")
        for band, sigma in sigma_of_band.items():
            sigma.c_iw.plot_ax(re_ax, cplx_mode="re", label=f"Re band: {band}", **style)
            sigma.c_iw.plot_ax(im_ax, cplx_mode="im", label=f"Im band: {band}", **style)

        re_ax.set_ylabel(r"$\Re{\Sigma_c}(i\omega)$ (eV)")
        im_ax.set_ylabel(r"$\Im{\Sigma_c}(i\omega)$ (eV)")
        set_grid_legend(ax_list, fontsize, xlabel=r"$i\omega$ (eV)")

        if with_tau:
            # Plot Sigmac_nk(itau)
            re_ax, im_ax = ax_list = ax_mat[1]
            style = dict(marker="o")
            for band, sigma in sigma_of_band.items():
                sigma.c_tau.plot_ax(re_ax, cplx_mode="re", label=f"Re band: {band}", **style)
                sigma.c_tau.plot_ax(im_ax, cplx_mode="im", label=f"Im band: {band}", **style)

            re_ax.set_ylabel(r"$\Re{\Sigma_c}(i\tau)$ (eV)")
            im_ax.set_ylabel(r"$\Im{\Sigma_c}(i\tau)$ (eV)")
            set_grid_legend(ax_list, fontsize, xlabel=r"$i\tau$ (a.u.)")

        fig.suptitle(r"$\Sigma_{nk}$" +  f" at k-point: {kpoint}, spin: {spin}", fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_sigma_real_axis(self, kpoint: KptSelect, spin=0, include_bands="gap",
                             fontsize=8, ax_mat=None, **kwargs) -> Figure:
        """
        Plot Sigma_nk(w) along the real-axis for given k-point, spin and set of bands.

        Args:
            kpoint: k-point in self-energy. Accepts |Kpoint|, vector or index.
            spin: Spin index.
            include_bands: List of bands to include. None means all.
            fontsize: Legend and title fontsize.
            ax_mat:
        """
        nrows, ncols = 1, 2
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        ax_mat = np.array(ax_mat)

        ikcalc, kpoint = self.r.get_ikcalc_kpoint(kpoint)
        include_bands = self._get_include_bands(include_bands, spin)
        sigma_of_band = self.r.read_sigma_bdict_sikcalc(spin, ikcalc, include_bands)
        nwr = self.r.nwr

        ax_list = re_ax, im_ax = ax_mat[0]
        # Show KS gap as filled area.
        self.ebands.add_fundgap_span(ax_list, spin)

        for band, sigma in sigma_of_band.items():
            ib = band - self.r.min_bstart
            e0 = self.r.e0_kcalc[spin, ikcalc, ib]

            ys = sigma.xc.values.real
            l = sigma.xc.plot_ax(re_ax, cplx_mode="re", label=f"Re band: {band}")
            point_style = dict(color=l[0].get_color(), marker="^", markersize=5.0)
            re_ax.plot(e0, ys[nwr//2], **point_style)

            ys = sigma.xc.values.imag
            l = sigma.xc.plot_ax(im_ax, cplx_mode="im", label=f"Im band: {band}")
            point_style = dict(color=l[0].get_color(), marker="^", markersize=5.0)
            im_ax.plot(e0, ys[nwr//2], **point_style)

        re_ax.set_ylabel(r"$\Re{\Sigma_{xc}(\omega)}$ (eV)")
        im_ax.set_ylabel(r"$\Im{\Sigma_{xc}(\omega)}$ (eV)")
        set_grid_legend(ax_list, fontsize=fontsize, xlabel=r"$\omega$ (eV)")

        fig.suptitle(r"$\Sigma_{nk}(\omega)$" +  f" at k-point: {kpoint}", fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_qps_vs_e0(self, with_fields="all", exclude_fields=None, e0="fermie",
                       xlims=None, sharey=False, ax_list=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot QP results stored in the GWR file as function of the KS energy.

        Args:
            with_fields: The names of the qp attributes to plot as function of eKS.
                Accepts: List of strings or string with tokens separated by blanks.
                See :class:`QPState` for the list of available fields.
            exclude_fields: Similar to ``with_fields`` but excludes fields
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            ax_list: List of |matplotlib-Axes| for plot. If None, new figure is produced.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            sharey: True if y-axis should be shared.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        with_fields = QPState.get_fields_for_plot(with_fields, exclude_fields)

        # Because qplist does not have the fermi level.
        fermie = self.ebands.get_e0(e0) if e0 is not None else None
        for spin in range(self.nsppol):
            fig = self.qplist_spin[spin].plot_qps_vs_e0(
                with_fields=with_fields, exclude_fields=exclude_fields, fermie=fermie,
                xlims=xlims, sharey=sharey, ax_list=ax_list, fontsize=fontsize,
                marker=self.marker_spin[spin], show=False, **kwargs)
            ax_list = fig.axes

        return fig

    @add_fig_kwargs
    def plot_spectral_functions(self, include_bands=None, ax_list=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the spectral function A_{nk}(w) for all k-points, bands and
        spins available in the GWR file.

        Args:
            include_bands: List of bands to include. None means all.
            ax_list:
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        # Build grid of plots.
        nrows, ncols = len(self.sigma_kpoints), 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        for ikcalc, (kcalc, ax) in enumerate(zip(self.sigma_kpoints, ax_list)):
            for spin in range(self.nsppol):
                include_bands_ks = self._get_include_bands(include_bands, spin)
                sigma_of_band = self.r.read_sigma_bdict_sikcalc(spin, ikcalc, include_bands_ks)

                for band, sigma in sigma_of_band.items():
                    label = r"$A(\omega)$: band: %d, spin: %d" % (band, spin)
                    l = sigma.plot_ax(ax, what="a", label=label, fontsize=fontsize)
                    # Show position of the KS energy as vertical line.
                    ib = band - self.r.min_bstart
                    ax.axvline(self.r.e0_kcalc[spin, ikcalc, ib],
                               lw=1, color=l[0].get_color(), ls="--")

            # Show KS gap as filled area.
            self.ebands.add_fundgap_span(ax, spin)

            ax.set_xlabel(r"$\omega$ (eV)")
            ax.set_ylabel(r"$A(\omega)$ (1/eV)")
            ax.set_title("k-point: %s" % repr(kcalc), fontsize=fontsize)

        return fig

    #@add_fig_kwargs
    #def plot_sparsity(self, what, origin="lower", **kwargs):
    #   x_mat
    #   ax.spy(mat, precision=0.1, markersize=5, origin=origin)

    #@add_fig_kwargs
    #def plot_mat(self, what, origin="lower", **kwargs):
    #   x_mat
    #   ax.spy(mat, precision=0.1, markersize=5, origin=origin)

    #def get_panel(self, **kwargs):
    #    """
    #    Build panel with widgets to interact with the GWR.nc either in a notebook or in panel app.
    #    """
    #    from abipy.panels.gwr import GwrFilePanel
    #    return GwrFilePanel(self).get_panel(**kwargs)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        verbose = kwargs.pop("verbose", 0)

        include_bands = "all" if verbose else "gaps"
        yield self.plot_spectral_functions(include_bands=include_bands, show=False)

        # TODO
        for spin in range(self.nsppol):
            for ik, kpoint in enumerate(self.sigma_kpoints):
                kws = dict(spin=spin, include_bands = "gaps", show=False)
                yield self.plot_sigma_imag_axis(ik, **kws)
                yield self.plot_sigma_real_axis(ik, **kws)

    def write_notebook(self, nbpath=None, title=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=title)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class GwrReader(ETSF_Reader):
    r"""
    This object provides methods to read data from the GWR.nc file.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GwrReader
    """

    def __init__(self, path: str):
        super().__init__(path)

        # Save important quantities needed to simplify the API.
        self.bands = ElectronBands.from_file(path)
        self.structure = self.read_structure()

        self.nsppol = self.bands.nsppol
        self.nwr = self.read_dimvalue("nwr")
        self.nkcalc = self.read_dimvalue("nkcalc")
        self.smat_bsize1 = self.read_dimvalue("smat_bsize1")
        self.smat_bsize2 = self.read_dimvalue("smat_bsize2")
        char_list = self.read_value("gwr_task")
        for i, ch in enumerate(reversed(char_list)):
            if not ch: break
        i = len(char_list) + i
        self.gwr_task = "".join(c.decode("utf-8") for c in char_list[:i])
        self.sig_diago = bool(self.read_value("sig_diago"))

        # The k-points where QP corrections have been calculated.
        kcalc_red_coords = self.read_value("kcalc")
        self.sigma_kpoints = KpointList(self.structure.reciprocal_lattice, kcalc_red_coords)
        # Find k-point name
        for kpoint in self.sigma_kpoints:
            kpoint.set_name(self.structure.findname_in_hsym_stars(kpoint))

        # Read mapping kcalc --> IBZ and convert to C indexing.
        # nctkarr_t("kcalc2ibz", "int", "nkcalc, six") &
        self.kcalc2ibz = self.read_variable("kcalc2ibz")[0,:] - 1

        # Note conversion between Fortran and python convention.
        self.bstart_sk = self.read_value("bstart_ks") - 1
        self.bstop_sk = self.read_value("bstop_ks")
        # min band index for GW corrections over spins and k-points
        self.min_bstart = np.min(self.bstart_sk)

    @lazy_property
    def iw_mesh(self) -> np.ndarray:
        """Frequency mesh in eV along the imaginary axis."""
        return self.read_value("iw_mesh") * abu.Ha_eV

    @lazy_property
    def tau_mesh(self) -> np.ndarray:
        """Tau mesh in a.u."""
        return self.read_value("tau_mesh")

    @lazy_property
    def wr_step(self) -> float:
        """Frequency-mesh along the real axis in eV."""
        return self.read_value("wr_step") * abu.Ha_eV

    @lazy_property
    def e0_kcalc(self) -> np.ndarray:
        # nctkarr_t("e0_kcalc", "dp", "smat_bsize1, nkcalc, nsppol"), &
        return self.read_value("e0_kcalc") * abu.Ha_eV

    def get_ikcalc_kpoint(self, kpoint) -> tuple(int, Kpoint):
        """
        Return the ikcalc index and the Kpoint
        """
        ikcalc = self.kpt2ikcalc(kpoint)
        kpoint = self.sigma_kpoints[ikcalc]
        return ikcalc, kpoint

    def kpt2ikcalc(self, kpoint) -> int:
        """
        Returns the index of the k-point in the sigma_kpoints array.
        Used to access data in the arrays that are dimensioned [0:nkcalc]
        """
        if duck.is_intlike(kpoint):
            return int(kpoint)
        else:
            return self.sigma_kpoints.index(kpoint)

    def get_wr_mesh(self, e0: float) -> np.ndarray:
        """
        The frequency mesh in eV is linear and centered around KS e0.
        """
        nwr = self.nwr
        return np.linspace(start=e0 - self.wr_step * (nwr // 2),
                           stop=e0 + self.wr_step * (nwr // 2),  num=nwr)

    def read_sigee_skb(self, spin, kpoint, band) -> SelfEnergy:
        """"
        Read self-energy for (spin, kpoint, band).
        """
        ikcalc, kpoint = self.get_ikcalc_kpoint(kpoint)
        ib = band - self.min_bstart
        ib2 = 0 if self.sig_diago else ib

        e0 = self.e0_kcalc[spin, ikcalc, ib]
        wmesh = self.get_wr_mesh(e0)

        # nctkarr_t("sigxc_rw_diag", "dp", "two, nwr, smat_bsize1, nkcalc, nsppol"), &
        sigxc_values = self.read_variable("sigxc_rw_diag")[spin,ikcalc,ib,:,:] * abu.Ha_eV
        sigxc_values = sigxc_values[:,0] + 1j *sigxc_values[:,1]

        # nctkarr_t("spfunc_diag", "dp", "nwr, smat_bsize1, nkcalc, nsppol") &
        spf_values = self.read_variable("spfunc_diag")[spin,ikcalc,ib,:] / abu.Ha_eV

        # nctkarr_t("sigc_iw_mat", "dp", "two, ntau, smat_bsize1, smat_bsize2, nkcalc, nsppol"), &
        sigc_iw = self.read_value("sigc_iw_mat", cmode="c") * abu.Ha_eV
        c_iw_values = sigc_iw[spin,ikcalc,ib2,ib]

        # nctkarr_t("sigc_it_mat", "dp", "two, two, ntau, smat_bsize1, smat_bsize2, nkcalc, nsppol"), &
        sigc_tau = self.read_value("sigc_it_mat", cmode="c") * abu.Ha_eV
        c_tau_pm = sigc_tau[spin,ikcalc,ib2,ib]
        tau_mp_mesh = np.concatenate((-self.tau_mesh[::-1], self.tau_mesh))
        c_tau_mp_values = np.concatenate((c_tau_pm[::-1,1], c_tau_pm[:,0]))

        return SelfEnergy(spin, kpoint, band, wmesh, sigxc_values, spf_values,
                          iw_mesh=self.iw_mesh, c_iw_values=c_iw_values,
                          tau_mp_mesh=tau_mp_mesh, c_tau_mp_values=c_tau_mp_values)

    def read_sigma_bdict_sikcalc(self, spin: int, ikcalc: int, include_bands) -> dict[int, SelfEnergy]:
        """
        Return dict of self-energy objects for given (spin, ikcalc) indexed by the band index.
        """
        sigma_of_band = {}
        for band in range(self.bstart_sk[spin, ikcalc], self.bstop_sk[spin, ikcalc]):
            if include_bands and band not in include_bands: continue
            sigma_of_band[band] = self.read_sigee_skb(spin, ikcalc, band)
        return sigma_of_band

    def read_allqps(self, ignore_imag=False) -> tuple[QPList]:
        """
        Return list with ``nsppol`` items. Each item is a :class:`QPList` with the QP results

        Args:
            ignore_imag: Only real part is returned if ``ignore_imag``.
        """
        qps_spin = self.nsppol * [None]

        for spin in range(self.nsppol):
            qps = []
            for kpoint in self.sigma_kpoints:
                ikcalc = self.kpt2ikcalc(kpoint)
                qps.extend(self.read_qplist_sk(spin, kpoint, ignore_imag=ignore_imag))

            qps_spin[spin] = QPList(qps)

        return tuple(qps_spin)

    def read_qplist_sk(self, spin, kpoint, ignore_imag=False) -> QPList:
        """
        Read and return a QPList object for the given spin, kpoint.

        Args:
            ignore_imag: Only real part is returned if ``ignore_imag``.
        """
        ikcalc, kpoint = self.get_ikcalc_kpoint(kpoint)

        def ri(a):
            return np.real(a) if ignore_imag else a

        # TODO: Finalize the implementation.
        #sigxme = sigx_mat
        sigxme = 0.0
        #self._sigxme[spin, ikcalc, ib],
        qp_list = QPList()
        for band in range(self.bstart_sk[spin, ikcalc], self.bstop_sk[spin, ikcalc]):
            ib = band - self.min_bstart

            qpe = self.read_variable("qpz_ene")[spin, ikcalc, ib] * abu.Ha_meV
            qpe = qpe[0] + 1j*qpe[1]

            ze0 = self.read_variable("ze0_kcalc")[spin, ikcalc, ib]
            ze0 = ze0[0] + 1j*ze0[1]

            # TODO Finalize the implementation
            qp_list.append(QPState(
                spin=spin,
                kpoint=kpoint,
                band=band,
                e0=self.e0_kcalc[spin, ikcalc, ib],
                qpe=ri(qpe),
                qpe_diago=0.0,
                #vxcme=self._vxcme[spin, ikcalc, ib],
                vxcme=0.0,
                sigxme=sigxme,
                #sigcmee0=ri(self._sigcmee0[spin, ikcalc, ib]),
                sigcmee0=0.0,
                vUme=0.0,
                ze0=ri(ze0),
            ))
        return qp_list


class GwrRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple GWR.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GwrRobot
    """
    # Try to have API similar to SigresRobot
    EXT = "GWR"

    # matplotlib option to fill convergence window.
    HATCH = "/"

    def __init__(self, *args):
        super().__init__(*args)
        if len(self.abifiles) in (0, 1): return

        # Check dimensions and self-energy states and issue warning.
        warns = []; wapp = warns.append
        nc0 : GwrFile = self.abifiles[0]
        same_nsppol, same_nkcalc = True, True
        if any(nc.nsppol != nc0.nsppol for nc in self.abifiles):
            same_nsppol = False
            wapp("Comparing ncfiles with different values of nsppol.")
        if any(nc.r.nkcalc != nc0.r.nkcalc for nc in self.abifiles):
            same_nkcalc = False
            wapp("Comparing ncfiles with different number of k-points in self-energy. Doh!")

        if same_nsppol and same_nkcalc:
            # FIXME
            # Different values of bstart_ks are difficult to handle
            # Because the high-level API assumes an absolute global index
            # Should decide how to treat this case: either raise or interpret band as an absolute band index.
            if any(np.any(nc.r.bstart_sk != nc0.r.bstart_sk) for nc in self.abifiles):
                wapp("Comparing ncfiles with different values of bstart_sk")
            if any(np.any(nc.r.bstop_sk != nc0.r.bstop_sk) for nc in self.abifiles):
                wapp("Comparing ncfiles with different values of bstop_sk")

        if warns:
            for w in warns:
                cprint(w, color="yellow")

    def _check_dims_and_params(self) -> None:
        """
        Test nsppol, sigma_kpoints.
        """
        if not len(self.abifiles) > 1:
            return

        nc0: GwrFile = self.abifiles[0]
        errors = []
        eapp = errors.append

        if any(nc.nsppol != nc0.nsppol for nc in self.abifiles[1:]):
            eapp("Files with different values of `nsppol`")

        if any(nc.nkcalc != nc0.nkcalc for nc in self.abifiles[1:]):
            cprint("Files with different values of `nkcalc`", color="yellow")

        for nc in self.abifiles[1:]:
            for k0, k1 in zip(nc0.sigma_kpoints, nc.sigma_kpoints):
                if k0 != k1:
                    cprint("Files with different values of `sigma_kpoints`\n"+
                           "Specify the kpoint via reduced coordinates and not via the index", "yellow")
                    break

        if errors:
            raise ValueError("Cannot compare multiple GWR.nc files. Reason:\n %s" % "\n".join(errors))

    def get_dataframe_sk(self, spin, kpoint, with_params=True, ignore_imag=False) -> pd.DataFrame:
        """
        Return |pandas-Dataframe| with QP results for this spin, k-point

        Args:
            spin: Spin index
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            with_params: True to add convergence parameters.
            ignore_imag: only real part is returned if ``ignore_imag``.
        """
        df_list = []; app = df_list.append
        for label, ncfile in self.items():
            df = ncfile.get_dataframe_sk(spin, kpoint, index=None,
                                         with_params=with_params, ignore_imag=ignore_imag)
            app(df)

        return pd.concat(df_list)

    def get_dirgaps_dataframe(self, sortby="kname", with_params=True) -> pd.DataFrame:
        """
        Returns |pandas-DataFrame| with QP direct gaps for all the files treated by the GWR robot.

        Args:
            sortby: Name to sort by.
            with_params: False to exclude calculation parameters from the dataframe.
        """
        with_geo = self.has_different_structures()

        df_list = []; app = df_list.append
        for _, ncfile in self.items():
            app(ncfile.get_dirgaps_dataframe(with_params=with_params, with_geo=with_geo))

        df = pd.concat(df_list)
        if sortby and sortby in df: df = df.sort_values(sortby)
        return df

    def get_dataframe(self, sortby="kname", with_params=True, ignore_imag=False) -> pd.DataFrame:
        """
        Return |pandas-Dataframe| with QP results for all k-points, bands and spins
        present in the files treated by the GWR robot.

        Args:
            sortby: Name to sort by.
            with_params: True to add parameters.
            ignore_imag: only real part is returned if ``ignore_imag``.
        """
        df_list = []; app = df_list.append
        for _, ncfile in self.items():
            for spin in range(ncfile.nsppol):
                for ikc, _ in enumerate(ncfile.sigma_kpoints):
                    app(ncfile.get_dataframe_sk(spin, ikc, with_params=with_params,
                                                ignore_imag=ignore_imag))

        df = pd.concat(df_list)
        if sortby and sortby in df: df = df.sort_values(sortby)
        return df

    def get_rpa_ene_dataframe(self, with_params=True) -> pd.DataFrame:
        """
        Return |pandas-Dataframe| with RPA energies for all the files
        treated by the GWR robot.

        Args:
            with_params: True to add parameters.
        """
        keys = ["ec_rpa_ecut_inf", "ec_mp2_ecut_inf"]
        dict_list = []
        for _, ncfile in self.items():
            d = ncfile.get_rpa_ene_dict()
            d = {k: d[k] for k in keys}
            if with_params:
                d.update(ncfile.params)
            dict_list.append(d)

        df = pd.DataFrame(dict_list)
        #if sortby and sortby in df: df = df.sort_values(sortby)
        return df

    @add_fig_kwargs
    def plot_selfenergy_conv(self, spin, kpoint, band, axis="wreal", sortby=None, hue=None,
                             colormap="viridis", xlims=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the convergence of the e-e self-energy wrt to the ``sortby`` parameter.
        Values can be optionally grouped by `hue`.

        Args:
            spin: Spin index.
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band: Band index.
            axis: "wreal": to plot Sigma(w) and A(w) along the real axis.
                  "wimag": to plot Sigma(iw)
                  "tau": to plot Sigma(itau)) along the imag axis.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            colormap: matplotlib color map.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colormap)

        # Make sure nsppol and sigma_kpoints are consistent.
        self._check_dims_and_params()
        ebands0 = self.abifiles[0].ebands

        if hue is None:
            # Build grid depends on axis.
            nrows = {"wreal": 3, "wimag": 2, "tau": 2}[axis]
            ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=1,
                                                    sharex=True, sharey=False, squeeze=False)
            ax_list = np.array(ax_list).ravel()

            lnp_list = self.sortby(sortby)
            for ix, (nclabel, ncfile, param) in enumerate(lnp_list):
                label = "%s: %s" % (self._get_label(sortby), param)
                kws = dict(label=label or nclabel, color=cmap(ix / len(lnp_list)))
                sigma = ncfile.r.read_sigee_skb(spin, kpoint, band)

                if axis == "wreal":
                    # Plot Sigma(w) along the real axis.
                    sigma.plot_reima_rw(ax_list, **kws)
                elif axis == "wimag":
                    # Plot Sigma(iw) along the imaginary axis.
                    sigma.plot_reimc_iw(ax_list, **kws)
                elif axis == "tau":
                    # Plot Sigma(itau) along the imaginary axis.
                    sigma.plot_reimc_tau(ax_list, **kws)

            if axis == "wreal": ebands0.add_fundgap_span(ax_list, spin)
            set_grid_legend(ax_list, fontsize)

        else:
            # group_and_sortby and build (3, ngroups) subplots
            groups = self.group_and_sortby(hue, sortby)
            nrows = {"wreal": 3, "wimag": 2, "tau": 2}[axis]
            ncols = len(groups)
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                   sharex=True, sharey=False, squeeze=False)
            for ig, g in enumerate(groups):
                subtitle = "%s: %s" % (self._get_label(hue), g.hvalue)
                ax_list = ax_mat[:,ig]
                ax_list[0].set_title(subtitle, fontsize=fontsize)

                for ix, (nclabel, ncfile, param) in enumerate(g):
                    kws = dict(label="%s: %s" % (self._get_label(sortby), param), color=cmap(ix / len(g)))
                    sigma = ncfile.r.read_sigee_skb(spin, kpoint, band)

                    if axis == "wreal":
                        lines = sigma.plot_reima_rw(ax_list, **kws)
                        if ix == 0:
                            # Show position of KS energy as vertical line.
                            ikcalc, _ = ncfile.r.get_ikcalc_kpoint(kpoint)
                            ib = band - ncfile.r.min_bstart
                            for ax, l in zip(ax_list, lines):
                                ax.axvline(ncfile.r.e0_kcalc[spin, ikcalc, ib], lw=1, color=l[0].get_color(), ls="--")

                    elif axis == "wimag":
                        # Plot Sigma_c(iw) along the imaginary axis.
                        sigma.plot_reimc_iw(ax_list, **kws)

                    elif axis == "tau":
                        # Plot Sigma(itau) along the imaginary axis.
                        sigma.plot_reimc_tau(ax_list, **kws)

                if axis == "wreal": ebands0.add_fundgap_span(ax_list, spin)
                set_grid_legend(ax_list, fontsize)
                for ax in ax_list:
                    set_axlims(ax, xlims, "x")

            if ig != 0:
                set_visible(ax_mat[:, ig], False, "ylabel")

        _, kpoint = self.abifiles[0].r.get_ikcalc_kpoint(kpoint)
        fig.suptitle(r"$\Sigma_{nk}$" + f" at k-point: {kpoint}, band: {band}, spin: {spin}",
                     fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_qpgaps_convergence(self, qp_kpoints="all", qp_type="qpz0", sortby=None, hue=None, abs_conv=0.01,
                                plot_qpmks=True, fontsize=8, **kwargs) -> Figure:
        """
        Plot the convergence of the direct QP gaps for all the k-points and spins treated by the GWR robot.

        Args:
            qp_kpoints: List of k-points in self-energy. Accept integers (list or scalars), list of vectors,
                or "all" to plot all k-points.
            qp_type: "qpz0" for linear qp equation with Z factor computed at KS e0,
                     "otms" for on-the-mass-shell values.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            abs_conv: If not None, show absolute convergence window.
            plot_qpmks: If False, plot QP_gap, KS_gap else (QP_gap - KS_gap)
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        # Make sure nsppol and sigma_kpoints are the same.
        self._check_dims_and_params()

        nc0 = self.abifiles[0]
        nsppol = nc0.nsppol
        qpkinds = nc0.find_qpkinds(qp_kpoints)
        if len(qpkinds) > 10:
            cprint("More that 10 k-points in file. Only 10 k-points will be shown. Specify kpt index expliclty", "yellow")
            qpkinds = qpkinds[:10]

        # Build grid with (nkpt, 1) plots.
        nrows, ncols = len(qpkinds), 1
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        if hue is None:
            labels, ncfiles, xs = self.sortby(sortby, unpack=True)
        else:
            groups = self.group_and_sortby(hue, sortby)

        #if qp_type not in {"qpz0", "otms"}:
        if qp_type not in {"qpz0", }:
            raise ValueError("Invalid qp_type: %s" % qp_type)

        name = "QP dirgap" if not plot_qpmks else "QP - KS dirgap"
        name = "%s (%s)" % (name, qp_type.upper())

        for ix, ((kpt, ikc), ax_row) in enumerate(zip(qpkinds, ax_mat)):
            ax = ax_row[0]
            for spin in range(nsppol):
                ax.set_title("%s k:%s" % (name, repr(kpt)), fontsize=fontsize)

                # Extract QP dirgap for [spin, ikcalc]
                if hue is None:
                    if qp_type == "qpz0": yvals = [ncfile.qpz0_dirgaps[spin, ikc] for ncfile in ncfiles]
                    #if qp_type == "otms": yvals = [ncfile.qp_dirgaps_otms_t[spin, ikc, itemp] for ncfile in ncfiles]
                    if plot_qpmks:
                        yvals = np.array(yvals) - np.array([ncfile.ks_dirgaps[spin, ikc] for ncfile in ncfiles])

                    lines = self.plot_xvals_or_xstr_ax(ax, xs, yvals, fontsize, marker=nc0.marker_spin[spin],
                                                       **kwargs)
                    hspan_ax_line(ax, lines[0], abs_conv, self.HATCH)

                else:
                    for g in groups:
                        if qp_type == "qpz0": yvals = [ncfile.qpz0_dirgaps[spin, ikc] for ncfile in g.abifiles]
                        #if qp_type == "otms": yvals = [ncfile.qp_dirgaps_otms_t[spin, ikc, itemp] for ncfile in g.abifiles]
                        if plot_qpmks:
                            yvals = np.array(yvals) - np.array([ncfile.ks_dirgaps[spin, ikc] for ncfile in g.abifiles])

                        label = "%s: %s" % (self._get_label(hue), g.hvalue)
                        lines = ax.plot(g.xvalues, yvals, marker=nc0.marker_spin[spin], label=label)

                        hspan_ax_line(ax, lines[0], abs_conv, self.HATCH)

            ax.grid(True)
            if ix == len(qpkinds) - 1:
                ax.set_ylabel("%s (eV)" % name)
                ax.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax, 15)
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_qpdata_conv_skb(self, spin, kpoint, band,
                             sortby=None, hue=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the convergence of the QP results for given (spin, kpoint, band).

        Args:
            spin: Spin index.
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band: Band index.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        # Make sure that nsppol and sigma_kpoints are consistent.
        self._check_dims_and_params()

        # TODO: Add more quantities DW, Fan(0)
        # TODO: Decide how to treat complex quantities, avoid annoying ComplexWarning
        # TODO: Format for g.hvalue
        # Treat fundamental gaps
        # Quantities to plot.
        what_list = ["re_qpe", "imag_qpe", "ze0"]

        # Build grid plot.
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        nc0: GwrFile = self.abifiles[0]
        ikc = nc0.kpt2ikcalc(kpoint)
        kpoint = nc0.sigma_kpoints[ikc]

        # Sort and read QP data.
        if hue is None:
            labels, ncfiles, params = self.sortby(sortby, unpack=True)
            qplist = [ncfile.reader.read_qp(spin, kpoint, band) for ncfile in ncfiles]
        else:
            groups = self.group_and_sortby(hue, sortby)
            qplist_group = []
            for g in groups:
                lst = [ncfile.reader.read_qp(spin, kpoint, band) for ncfile in g.abifiles]
                qplist_group.append(lst)

        for ix, (ax, what) in enumerate(zip(ax_list, what_list)):
            if hue is None:
                # Extract QP data.
                #yvals = [getattr(qp, what)[itemp] for qp in qplist]
                if not duck.is_string(params[0]):
                    ax.plot(params, yvals, marker=nc0.marker_spin[spin])
                else:
                    # Must handle list of strings in a different way.
                    xn = range(len(params))
                    ax.plot(xn, yvals, marker=nc0.marker_spin[spin])
                    ax.set_xticks(xn)
                    ax.set_xticklabels(params, fontsize=fontsize)
            else:
                for g, qplist in zip(groups, qplist_group):
                    # Extract QP data.
                    yvals = [getattr(qp, what)[itemp] for qp in qplist]
                    label = "%s: %s" % (self._get_label(hue), g.hvalue)
                    ax.plot(g.xvalues, yvals, marker=nc0.marker_spin[spin],
                            label=label if ix == 0 else None)

            ax.grid(True)
            ax.set_ylabel(what)
            if ix == len(what_list) - 1:
                ax.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax, 15)
            if ix == 0 and hue is not None:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        if "title" not in kwargs:
            title = "QP results spin: %s, k:%s, band: %s, T = %.1f K" % (
                    spin, repr(kpoint), band, nc0.tmesh[itemp])
            fig.suptitle(title, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_qpfield_vs_e0(self, field, reim="real", function=lambda x: x, sortby=None, hue=None,
                           fontsize=8, colormap="jet", e0="fermie", **kwargs) -> Figure:
        """
        For each file in the GWR robot, plot one of the attributes of :class:`QpTempStat
        as a function of the KS energy.

        Args:
            field (str): String defining the attribute to plot.
            reim: Plot the real or imaginary part
            function: Apply a function to the results before plotting
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it is assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            colormap: matplotlib color map.
            fontsize: legend and label fontsize.
            e0: Option used to define the zero of energy in the band structure plot.

        .. note::

            For the meaning of the other arguments, see other robot methods.

        Returns: |matplotlib-Figure|
        """
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colormap)

        if hue is None:
            lnp_list = self.sortby(sortby)
            ax_list = None
            for i, (label, ncfile, param) in enumerate(lnp_list):
                if sortby is not None:
                    label = "%s: %s" % (self._get_label(sortby), param)
                fig = ncfile.plot_qps_vs_e0(with_fields=list_strings(field),
                    reim=reim, function=function, e0=e0, ax_list=ax_list,
                    color=cmap(i / len(lnp_list)), fontsize=fontsize,
                    label=label, show=False)
                ax_list = fig.axes
        else:
            # group_and_sortby and build (ngroups,) subplots
            groups = self.group_and_sortby(hue, sortby)
            nrows, ncols = 1, len(groups)
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                   sharex=True, sharey=True, squeeze=False)
            for ig, g in enumerate(groups):
                subtitle = "%s: %s" % (self._get_label(hue), g.hvalue)
                ax_mat[0, ig].set_title(subtitle, fontsize=fontsize)
                for i, (nclabel, ncfile, param) in enumerate(g):
                    fig = ncfile.plot_qps_vs_e0(with_fields=list_strings(field),
                        reim=reim, function=function,
                        e0=e0, ax_list=ax_mat[:, ig], color=cmap(i / len(g)), fontsize=fontsize,
                        label="%s: %s" % (self._get_label(sortby), param), show=False)

                if ig != 0:
                    set_visible(ax_mat[:, ig], False, "ylabel")

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        verbose = kwargs.pop("verbose", 0)
        yield self.plot_qpgaps_convergence(qp_kpoints="all", show=False)

        # Visualize the convergence of the self-energy for
        # all the k-points and the most important bands.
        nc0: GwrFile = self.abifiles[0]

        for spin in range(nc0.nsppol):
            for ikcalc in range(nc0.nkcalc):
                ik_ibz = nc0.r.kcalc2ibz[ikcalc]
                band_v = nc0.ebands.homo_sk(spin, ik_ibz).band
                band_c = nc0.ebands.lumo_sk(spin, ik_ibz).band
                for band in range(band_v, band_c + 1):
                    for axis in ("wreal", "wimag", "tau"):
                        yield self.plot_selfenergy_conv(spin, ikcalc, band, axis=axis, show=False)

    def write_notebook(self, nbpath=None, title=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=title)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            nbv.new_code_cell("robot = abilab.SigEPhRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.get_params_dataframe()"),
            nbv.new_code_cell("# data = robot.get_dataframe()\ndata"),
            nbv.new_code_cell("robot.plot_qpgaps_convergence(itemp=0, sortby=None, hue=None);"),
            nbv.new_code_cell("""\
nc0 = robot.abifiles[0]
for spin in range(nc0.nsppol):
    for ikc, sigma_kpoint in enumerate(nc0.sigma_kpoints):
        for band in range(nc0.r.bstart_sk[spin, ikc], nc0.bstop_sk[spin, ikc]):
            robot.plot_qpdata_conv_skb(spin, sigma_kpoint, band, itemp=0, sortby=None, hue=None);"""),

            nbv.new_code_cell("""\
#nc0 = robot.abifiles[0]
#for spin in range(nc0.nsppol):
#    for ikc, sigma_kpoint in enumerate(nc0.sigma_kpoints):
#        for band in range(nc0.r.bstart_sk[spin, ikc], nc0.bstop_sk[spin, ikc]):
#           robot.plot_selfenergy_conv(spin, sigma_kpoint, band, itemp=0, sortby=None);"),"""),
        ])

        # Mixins.
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)



#class GwrSigresComparator:
#    """
#    This object provides a high-level API to compare the results stored in
#    a GWR.nc and a SIGRES.nc file (usually produced with the AC method).
#    """
#
#    def __init__(self, gwr_filepaths, sigres_filepaths):
#        """
#        Args:
#            gwr_filepaths: Path(s) of the GWR.nc file.
#            sigres_filepaths: Path(s) of the SIGRES.nc file.
#        """
#        self.gwr_robot = GwrRobot.from_files(gwr_filepaths)
#        self.sigres_robot = SigresRobot.from_files(sigres_filepaths)
#
#    def __enter__(self):
#        return self
#
#    def __exit__(self, exc_type, exc_val, exc_tb):
#        """
#        Activated at the end of the with statement. It automatically closes the files.
#        """
#        self.gwr_robot.close()
#        self.sigres_robot.close()
#
#    def __str__(self) -> str:
#        return self.to_string()
#
#    def to_string(self, verbose: int = 0) -> str:
#        """String representation with verbosity level ``verbose``."""
#        lines = []; app = lines.append
#        app(self.gwr_robot.to_string(verbose=verbose))
#        app("")
#        app(self.sigres_robot.to_string(verbose=verbose))
#        return "\n".join(lines)


def _find_g(gg, gvecs):
    for ig, g_sus in enumerate(gvecs):
        if all(gg == g_sus): return ig
    raise ValueError(f"Cannot find g-vector: {gg}")



class TchimVsSus:
    """
    Object used to compare the polarizability computed in the GWR code with the one
    produced by the legacy algorithm based on the Adler-Wiser expression.

    Example:

        gpairs = [
            ((0, 0, 0), (1, 0, 0)),
            ((1, 0, 0), (0, 1, 0)),
        ]
        qpoint_list = [
            [0, 0, 0],
            [0.5, 0.5, 0],
        ]

        with TchimVsSus("runo_DS3_TCHIM.nc", "AW_CD/runo_DS3_SUS.nc") as o
            o.expose_qpoints_gpairs(qpoint_list, gpairs, exposer="mpl")
    """
    def __init__(self, tchim_filepath: str, sus_filepath: str):
        """
        Args:
            tchim_filepath: TCHIM filename.
            sus_filepath: SUS filename.
        """
        from abipy.electrons.scr import SusFile
        self.sus_file = SusFile(sus_filepath)
        self.tchi_reader = ETSF_Reader(tchim_filepath)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Activated at the end of the with statement. It automatically closes the files.
        """
        self.sus_file.close()
        self.tchi_reader.close()

    def find_tchim_qpoint(self, qpoint) -> tuple[int, np.ndarray]:
        """
        Find the index of the q-point in the IBZ. Return index and q-points.
        """
        tchi_qpoints = self.tchi_reader.read_variable("qibz")

        if duck.is_intlike(qpoint):
            return qpoint, tchi_qpoints

        for tchi_iq, tchi_qq in enumerate(tchi_qpoints):
            if np.all(abs(tchi_qq - qpoint) < 1e-6):
                return tchi_iq, tchi_qpoints

        raise ValueError(f"Cannot find q-point: {qpoint} in TCHIM file")

    @add_fig_kwargs
    def plot_qpoint_gpairs(self, qpoint, gpairs,
                           fontsize=8, spins=(0, 0), with_title=True, **kwargs) -> Figure:
        """
        Plot the Fourier components of the polarizability for given q-point and list of (g, g') pairs.

        Args:
            qpoint: q-point or q-point index.
            gpairs: List of (g,g') pairs
        """
        gpairs = np.array(gpairs)

        sus_reader = self.sus_file.reader
        _, sus_iq = sus_reader.find_kpoint_fileindex(qpoint)

        # int reduced_coordinates_plane_waves_dielectric_function(number_of_qpoints_dielectric_function,
        # number_of_coefficients_dielectric_function, number_of_reduced_dimensions) ;
        # SUSC file uses Gamma-centered G-spere so read at iq = 0
        susc_iwmesh = sus_reader.read_value("frequencies_dielectric_function", cmode="c").imag
        susc_gvecs = sus_reader.read_variable("reduced_coordinates_plane_waves_dielectric_function")[0]

        tchi_iq, _ = self.find_tchim_qpoint(qpoint)
        tchi_iwmesh = self.tchi_reader.read_value("iw_mesh")
        tchi_gvecs = self.tchi_reader.read_variable("gvecs")[tchi_iq]

        # Find indices of gpairs in the two files
        sus_inds = [(_find_g(g1, susc_gvecs), _find_g(g2, susc_gvecs)) for g1, g2 in gpairs]
        chi_inds = [(_find_g(g1, tchi_gvecs), _find_g(g2, tchi_gvecs)) for g1, g2 in gpairs]

        if sus_reader.path.endswith("SUS.nc"):
            sus_var_name = "polarizability"
        elif sus_reader.path.endswith("SCR.nc"):
            sus_var_name = "inverse_dielectric_function"
        else:
            raise ValueError(f"Cannot detect varname for file: {sus_reader.path}")

        # Build grid of plots.
        num_plots, ncols, nrows = 2 * len(gpairs), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        for i, ((g1, g2), (sus_ig1, sus_ig2), (chi_ig1, chi_ig2)) in enumerate(zip(gpairs, sus_inds, chi_inds)):
            ax_re, ax_im = ax_mat[i]

            # number_of_qpoints_dielectric_function, number_of_frequencies_dielectric_function,
            # number_of_spins, number_of_spins, number_of_coefficients_dielectric_function,
            # number_of_coefficients_dielectric_function, complex)
            sus_data = sus_reader.read_variable(sus_var_name)[sus_iq, :, 0, 0, sus_ig2, sus_ig1]
            sus_data = sus_data[:, 0] + 1j * sus_data[:, 1]

            # nctkarr_t("mats", "dp", "two, mpw, mpw, ntau, nqibz, nsppol")
            tchi_data = self.tchi_reader.read_variable("mats")[0, tchi_iq, :, chi_ig2, chi_ig1, :]
            tchi_data = tchi_data[:, 0] + 1j * tchi_data[:, 1]

            style = dict(markersize=4)
            ax_re.plot(susc_iwmesh[1:], sus_data[1:].real, ls="--", marker="o", label="AW Re", **style)
            ax_re.plot(tchi_iwmesh, tchi_data[0:].real, ls=":", marker="^", label="minimax Re", **style)
            ax_re.grid(True)

            ax_im.plot(susc_iwmesh[1:], sus_data[1:].imag, ls="--", marker="o", label="AW Im", **style)
            ax_im.plot(tchi_iwmesh, tchi_data[0:].imag, ls=":", marker="^", label="minimax Im", **style)
            ax_im.grid(True)

            tex_g1, tex_g2, tex_qpt = r"${\bf{G}}_1$", r"${\bf{G}}_2$", r"${\bf{q}}$"
            ax_re.set_title(f"{tex_g1}: {g1}, {tex_g2}: {g2}, {tex_qpt}: {qpoint}", fontsize=fontsize)
            ax_im.set_title(f"{tex_g1}: {g1}, {tex_g2}: {g2}, {tex_qpt}: {qpoint}", fontsize=fontsize)

            ax_re.legend(loc="best", fontsize=fontsize, shadow=True)
            ax_im.legend(loc="best", fontsize=fontsize, shadow=True)

            #if i == 0:
            ax_re.set_ylabel(r'$\Re{\chi^0_{\bf{G}_1 \bf{G}_2}(i\omega)}$')
            ax_im.set_ylabel(r'$\Im{\chi^0_{\bf{G}_1 \bf{G}_2}(i\omega)}$')

            if i == len(gpairs) - 1:
                ax_re.set_xlabel(r'$i \omega$ (Ha)')
                ax_im.set_xlabel(r'$i \omega$ (Ha)')

            # The two files may store chi on different meshes.
            # Here we select the min value for comparison purposes.
            w_max = min(susc_iwmesh[-1], tchi_iwmesh[-1]) + 5
            xlims = [0, w_max]
            xlims = [-0.2, 7]

            set_axlims(ax_re, xlims, "x")
            set_axlims(ax_im, xlims, "x")

        if with_title:
            tex1 = r"$\chi^0(i\omega)$"
            tex2 = r"$\chi^0(i\tau) \rightarrow\chi^0(i\omega)$"
            fig.suptitle(f"Comparison between Adler-Wiser {tex1} and minimax {tex2}\nLeft: real part. Right: imag part",
                         fontsize=fontsize)

        return fig

    def expose_qpoints_gpairs(self, qpoint_list, gpairs, exposer="mpl"):
        """
        Plot the Fourier components of the polarizability for a list of q-points and,
        for each q-point, a list of (g, g') pairs.

        Args:
            qpoint_list: List of q-points to consider.
            gpairs: List of (g,g') pairs.
            exposer: "mpl" for matplotlib, "panel" for web interface.
        """
        with Exposer.as_exposer(exposer) as e:
            for qpoint in qpoint_list:
                e(self.plot_qpoint_gpairs(qpoint, gpairs, show=False))
            return e

    @add_fig_kwargs
    def plot_matdiff(self, qpoint, iw_index: int, npwq=101, with_susmat=False,
                     cmap="jet", fontsize=8, spins=(0, 0), **kwargs) -> Figure:
        """
        Plot G-G' matrix with the absolute value of chi_AW and chi_GWR for given q-point and omega index.

        Args:
            qpoint: q-point or q-point index.
            iw_index: Frequency index.
            with_susmat: True to add additional subplot with absolute value of susmat_{g,g'}.
        """
        sus_reader = self.sus_file.reader
        _, sus_iq = sus_reader.find_kpoint_fileindex(qpoint)

        # int reduced_coordinates_plane_waves_dielectric_function(number_of_qpoints_dielectric_function,
        # number_of_coefficients_dielectric_function, number_of_reduced_dimensions) ;
        # SUSC file uses Gamma-centered G-spere so read at iq = 0
        susc_iwmesh = sus_reader.read_value("frequencies_dielectric_function", cmode="c").imag
        susc_gvecs = sus_reader.read_variable("reduced_coordinates_plane_waves_dielectric_function")[0]

        tchi_iq, qibz = self.find_tchim_qpoint(qpoint)
        tchi_iwmesh = self.tchi_reader.read_value("iw_mesh")
        tchi_gvecs = self.tchi_reader.read_variable("gvecs")[tchi_iq]
        #FIXME: Requires new format
        npwq = self.tchi_reader.read_variable("chinpw_qibz")[tchi_iq]
        #npwq = min(npwq, len(tchi_gvecs))
        tchi_gvecs = tchi_gvecs[:npwq]

        # Find indices of tchi_gvecs in susc_gvecs to remap matrix
        ig_tchi2sus = np.array([_find_g(g1, susc_gvecs) for g1 in tchi_gvecs])

        if sus_reader.path.endswith("SUS.nc"):
            sus_var_name = "polarizability"
        elif sus_reader.path.endswith("SCR.nc"):
            sus_var_name = "inverse_dielectric_function"
        else:
            raise ValueError(f"Cannot detect varname for file: {sus_reader.path}")

        # Read full matrix from file, extract smaller submatrix using sus_index, chi_inds
        #
        # number_of_qpoints_dielectric_function, number_of_frequencies_dielectric_function,
        # number_of_spins, number_of_spins, number_of_coefficients_dielectric_function,
        # number_of_coefficients_dielectric_function, complex)
        sus_data = sus_reader.read_variable(sus_var_name)[sus_iq, iw_index, 0, 0]
        sus_data = (sus_data[:,:, 0] + 1j * sus_data[:,:,1]).T.copy()

        # TODO: Very inefficient for large npwq.
        sus_mat = np.zeros((npwq, npwq), dtype=complex)
        for ig2 in range(npwq):
            ig2_sus = ig_tchi2sus[ig2]
            for ig1 in range(npwq):
                ig1_sus = ig_tchi2sus[ig1]
                sus_mat[ig1,ig2] = sus_data[ig1_sus,ig2_sus]
        sus_data = sus_mat

        # nctkarr_t("mats", "dp", "two, mpw, mpw, ntau, nqibz, nsppol")
        tchi_data = self.tchi_reader.read_variable("mats")[0, tchi_iq, iw_index]
        tchi_data = (tchi_data[:,:,0] + 1j * tchi_data[:,:,1]).T.copy()
        tchi_data = tchi_data[:npwq,:npwq]

        # compute and visualize diff mat
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=1, ncols=2 if with_susmat else 1,
                                                squeeze=False)
        ax_list = ax_list.ravel()

        import matplotlib.colors as colors
        mat = data_from_cplx_mode("abs", tchi_data - sus_data)
        ax = ax_list[0]
        im = ax.matshow(mat, cmap=cmap, norm=colors.LogNorm(vmin=mat.min(), vmax=mat.max()))
        fig.colorbar(im, ax=ax)
        qq = qibz[tchi_iq]
        qq_str = "[%.2f, %.2f, %.2f]" % (qq[0], qq[1], qq[2])
        q_info = r"at ${\bf{q}}: %s$" % qq_str
        ww_ev = susc_iwmesh[iw_index] * abu.Ha_eV
        w_info = r", $i\omega$: %.2f eV" % ww_ev
        label = r"$|\chi^0_{AW}({\bf{G}}_1,{\bf{G}}_2) - \chi^0_{GWR}({\bf{G}}_1,{\bf{G}}_2))|$ " + q_info + w_info
        ax.set_title(label, fontsize=10)

        if with_susmat:
            mat = data_from_cplx_mode("abs", sus_data)
            ax = ax_list[1]
            im = ax.matshow(mat, cmap=cmap, norm=colors.LogNorm(vmin=mat.min(), vmax=mat.max()))
            fig.colorbar(im, ax=ax)
            label = r"$|\chi^0_{AW}(\bf{G}_1,\bf{G}_2)|$"
            ax.set_title(label, fontsize=10)

        return fig

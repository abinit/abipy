# coding: utf-8
"""
Objects to analyze the TCHIM file produced by the GWR code.
"""
from __future__ import annotations

import itertools
import numpy as np
#import pandas as pd
import abipy.core.abinit_units as abu

from typing import Any
from monty.functools import lazy_property
from monty.string import list_strings, marquee
#from monty.termcolor import cprint
from abipy.core.structure import Structure
#from abipy.core.kpoints import Kpoint, KpointList
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands # , NotebookWriter
from abipy.electrons.ebands import ElectronBands, RobotWithEbands, ElectronsReader
from abipy.tools import duck
from abipy.tools.typing import Figure, KptSelect, GvecSelect
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, Marker,
    set_axlims, set_ax_xylabels, set_visible, rotate_ticklabels, set_grid_legend, hspan_ax_line, Exposer)
from abipy.tools.numtools import data_from_cplx_mode


__all__ = [
    "TchimFile",
    "TchimVsSus",
]


def _find_g(gg, gvecs) -> int:
    """Return the index of gg vector in gvecs."""
    for ig, g_sus in enumerate(gvecs):
        if all(gg == g_sus): return ig
    raise ValueError(f"Cannot find g-vector: {gg}")


class TchimFile(AbinitNcFile, Has_Structure, Has_ElectronBands):
    """
    This file stores the e-ph matrix elements produced by the EPH code of Abinit
    and provides methods to analyze and plot results.
    The same object can be used to post-process _WCIMW.nc files with
    the correlated part of W along the imaginary axis.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: TchimFile
    """
    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = TchimReader(filepath)

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| with the KS energies."""
        return self.r.read_ebands()

    def close(self) -> None:
        """Close the file."""
        self.r.close()

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        #app(f"nsppol: {self.r.nsppol}")
        #app(f"nqibz: {self.r.nqibz}")
        #app(f"gstore_with_vk: {self.r.with_vk}")
        #app(f"gstore_kptopt: {self.r.kptopt}")
        #app(f"gstore_qptopt: {self.r.qptopt}")

        return "\n".join(lines)

    @lazy_property
    def params(self) -> dict:
        """
        dict with parameters that might be subject to convergence studies e.g ecuteps.
        """
        r = self.r
        return dict(
            gwr_ntau=r.read_dimvalue("ntau"),
            nband=self.ebands.nband,
            #ecuteps=r.read_value("ecuteps"),
            #ecutsigx=r.read_value("ecutsigx"),
            #ecut=r.read_value("ecut"),
            #gwr_boxcutmin=r.read_value("gwr_boxcutmin"),
            #nkpt=self.ebands.nkpt,
            #symchi=r.read_value("symchi"),
            #symsigma=r.read_value("symsigma"),
        )

    @add_fig_kwargs
    def plot_mats_ggp(self,
                      g1: GvecSelect,
                      g2: GvecSelect,
                      spin: int = 0,
                      wt_space: str = "omega",
                      method: str = "scan",
                      verbose: int = 0,
                      fontsize: int = 6,
                      **kwargs) -> Figure:
        """
        Plot the matrix elements for the given (g1, g2) and all the q-points in the IBZ  along the imaginary axis

        Args:
            g1, g2: g-vector or index.
            spin: Spin index.
            wt_space:
            method:
            verbose: Verbosity level.
            fontsize:
        """
        qibz = self.r.qibz
        nqibz = len(qibz)

        # ncols 2 for Re/Im
        nrows, ncols = nqibz, 2
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        style = dict(markersize=4, ls="--", marker="o")
        fit_style = dict(markersize=4, ls="-.", marker="x")

        xs = self.r.iw_mesh * abu.Ha_eV if wt_space == "omega" else self.r.tau_mesh

        for iq_ibz, qpoint in enumerate(qibz):
            gvec1, gvec2, qpoint, cvals = self.r.read_mat_ggp_at_qpt(iq_ibz, g1, g2, spin, wt_space)
            # Fit data.
            fit_xs = xs
            fit_ys = Fit(xs, cvals, self.r.tau_wgs, self.r.iw_wgs, wt_space, verbose).eval(fit_xs, method)

            # Plot results.
            re_ax, im_ax = ax_mat[iq_ibz, :]
            title = f"qpt={qpoint}, g1={gvec1}, g2={gvec2}, {spin=}"

            re_ax.plot(xs, cvals.real, label="ab-initio Re", **style)
            re_ax.plot(fit_xs, fit_ys.real, label="Fit Re", **fit_style)
            set_grid_legend(re_ax, fontsize, title=title)

            im_ax.plot(xs, cvals.imag, label="ab-initio Im", **style)
            im_ax.plot(fit_xs, fit_ys.imag, label="Fit Im", **fit_style)
            set_grid_legend(im_ax, fontsize, title=title)

            if iq_ibz == len(qibz) -1:
                xlabel = r"$i\omega$ (eV)" if wt_space == "omega" else r"$i\tau$ (a.u.)"
                re_ax.set_xlabel(xlabel)
                im_ax.set_xlabel(xlabel)

        return fig

    @add_fig_kwargs
    def plot_mats_all_ggp(self,
                          spin: int = 0,
                          take_every: int = 10,
                          wt_space: str = "omega",
                          method: str = "scan",
                          verbose: int = 0,
                          fontsize: int = 6,
                          **kwargs) -> Figure:
        """
        Plot the matrix elements for the given (g1, g2) and all the q-points in the IBZ  along the imaginary axis

        Args:
            spin: Spin index.
            wt_space:
            method:
            verbose: Verbosity level.
            fontsize:
        """
        qibz = self.r.qibz
        nqibz = len(qibz)

        # ncols 2 for Re/Im
        nrows, ncols = nqibz, 2
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        def get_error(x_ref, y, mode, epsilon=1e-10):
            """Compute "distance" between x and y according to mode."""
            if mode == "percent_error":
                # Compute relative difference based on reference value with a small constant epsilon.
                diff = np.abs(x_ref - y)
                return 100 * diff / np.maximum(np.abs(x_ref), epsilon)
                #return 100 * diff / np.maximum((np.abs(x_ref) + np.abs(y) / 2), epsilon)

            elif mode == "diff":
                return x_ref - y

            elif mode == "abs_diff":
                return np.abs(x_ref - y)

            raise ValueError(f"Invalid {mode=}")

        mode = "diff"
        mode = "abs_diff"
        mode = "percent_error"

        xs = self.r.iw_mesh * abu.Ha_eV if wt_space == "omega" else self.r.tau_mesh
        fit_xs = xs
        style = dict(markersize=2, ls="--", marker="o")

        for iq_ibz, qpoint in enumerate(qibz):
            gvecs, qpoint, cmat_g1g2 = self.r.read_mat_all_ggp_at_qpt(iq_ibz, spin, wt_space)
            npw_q = len(gvecs)
            re_ax, im_ax = ax_mat[iq_ibz, :]

            for count, (ig1, ig2) in enumerate(itertools.product(range(npw_q), range(npw_q))):
                if count % take_every != 0:
                    continue

                cvals = cmat_g1g2[ig1, ig2]
                # Fit data.
                fit_ys = Fit(xs, cvals, self.r.tau_wgs, self.r.iw_wgs, wt_space, verbose).eval(fit_xs, method)

                # Plot results.
                ys = get_error(cvals.real, fit_ys.real, mode)
                re_ax.plot(xs, ys, **style)
                ys = get_error(cvals.imag, fit_ys.imag, mode)
                im_ax.plot(xs, ys, **style)

            set_grid_legend(re_ax, fontsize, title=f"Real part qpt={qpoint}, {spin=}")
            set_grid_legend(im_ax, fontsize, title=f"Imag part qpt={qpoint}, {spin=}")

            if iq_ibz == len(qibz) -1:
                xlabel = r"$i\omega$ (eV)" if wt_space == "omega" else r"$i\tau$ (a.u.)"
                re_ax.set_xlabel(xlabel)
                im_ax.set_xlabel(xlabel)

        return fig


class Fit:
    """
    Fit complex values along the imaginary axis, either in frequency or time domain.
    """

    def __init__(self, xs, ys, tau_wgs, iw_wgs, wt_space, verbose):
        """
        Args:
            xs:
            ys:
            tau_wgs:
            iw_wgs:
            wt_space:
            verbose: Verbosity level.
        """
        self.xs, self.ys = xs, ys
        self.tau_wgs, self.iw_wgs = tau_wgs, iw_wgs
        self.wt_space = wt_space
        self.verbose = verbose

        if wt_space not in ("omega", "tau"):
            raise ValueError(f"Invalid {wt_space=}")

    def eval(self, fit_xs: np.ndarray, method: str) -> np.ndarray:
        """
        Fit data and evaluate the fit on `fit_xs` points.
        """
        if self.wt_space == "omega":
            return self._eval_omega(fit_xs, method)
        if self.wt_space == "tau":
            return self._eval_tau(fit_xs, method)

        raise ValueError(f"Invalida {self.wt_space=}")

    def _minimize_loss(self, fit_with_second_index, method):
        """
        Fit the signal by changing the second point, and select
        the one which minimizes the loss function.
        """
        # Select weights for the loss function.
        weights = self.iw_wgs if self.wt_space == "omega" else self.tau_wgs

        # Compute loss functions for all possible values of second_index.
        losses = []
        for second_index in range(1, len(self.xs)):
            ys_fit = fit_with_second_index(second_index, self.xs)
            loss = np.sum(weights * np.abs(ys_fit - self.ys)**2)
            losses.append((second_index, loss))

        # Fin min of losses.
        return min(losses, key=lambda t: t[1])

    def _eval_omega(self, fit_xs: np.ndarray, method: str) -> np.ndarray | None:
        """
        Fit values in imaginary frequency using A/(B^2 + omega^2).
        """
        def _fit_with_second_index(second_index, xs) -> np.ndarray:
            # First and second_index data points
            w0, y0 = self.xs[0], self.ys[0]
            wn, yn = self.xs[second_index], self.ys[second_index]

            # TODO: Check equations
            # Solve for B^2 using real part (assuming the same B^2 for real and imaginary parts)
            b2 = (y0.real * w0**2 - yn.real * wn**2) / (yn.real - y0.real)
            if b2 < 0.0:
                a, b = 0.0, 0.0
                return a / (b**2 + xs**2)
                #b2 = 0

            b = np.sqrt(b2)
            #print(f"{b2=}, {b=}")
            # Solve for Re[A] and Im[A]
            a_real = y0.real * (b2 + w0**2)
            a_imag = y0.imag * (b2 + w0**2)
            # Generate the fitted signal
            a = a_real + 1j * a_imag
            return a / (b**2 + xs**2)

        second_index, _ = self._minimize_loss(_fit_with_second_index, method)
        return _fit_with_second_index(second_index, fit_xs)

    def _eval_tau(self, fit_xs: np.ndarray, method: str) -> np.ndarray:
        """
        Fit values in imaginary time using B exp^{-a t} with a > 0.
        """
        def _fit_with_second_index(second_index, xs) -> np.ndarray:
            # First and second_index data points
            w0, y0 = self.xs[0], self.ys[0]
            wn, yn = self.xs[second_index], self.ys[second_index]
            # TODO: Check equations
            # Compute b and a
            b = y0
            a = -np.log(yn / y0) / (wn - w0)
            # Generate the fitted signal
            return b * np.exp(-a * (xs - w0))

        second_index, _ = self._minimize_loss(_fit_with_second_index, method)
        return _fit_with_second_index(second_index, fit_xs)


class TchimReader(ElectronsReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: TchimReader
    """
    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("number_of_spins")
        self.nqibz = self.read_dimvalue("nqibz")
        self.ntau = self.read_dimvalue("ntau")

        # Read important variables.
        self.iw_mesh = self.read_value("iw_mesh")
        self.iw_wgs = self.read_value("iw_wgs")
        self.tau_mesh = self.read_value("tau_mesh")
        self.tau_wgs = self.read_value("tau_wgs")
        self.qibz = self.read_value("qibz")

    def find_qpoint(self, qpoint, tol=1e-6):
        """
        Find the index of the q-point in the IBZ. Return index and q-point.
        """
        if duck.is_intlike(qpoint):
            iq_ibz = qpoint
            return iq_ibz, self.qibz[iq_ibz]

        for iq_ibz, qq in enumerate(self.qibz):
            if np.all(abs(qq - qpoint) < tol):
                return iq_ibz, qq

        raise ValueError(f"Cannot find {qpoint=}, {type(qpoint)=}")

    @staticmethod
    def _find_ig_g(gg, gvecs):
        """Return the idex of gg in gvecs and the g-vector."""
        if duck.is_intlike(gg):
            return gg, gvecs[gg]

        ig = _find_g(gg, gvecs)
        return ig, gvecs[ig]

    def read_mat_ggp_at_qpt(self,
                            qpoint: KptSelect,
                            g1: GvecSelect,
                            g2: GvecSelect,
                            spin: int,
                            wt_space: str):
        """
        Args:
            qpoint:
            g1, g2: g-vector or index.
            spin: Spin index.
            wt_space:
        """
        # Fortran arrays on disk.
        # nctkarr_t("chinpw_qibz", "int", "nqibz")
        # nctkarr_t("gvecs", "int", "three, mpw, nqibz")
        # nctkarr_t("mats_w", "dp", "two, mpw, mpw, ntau, nqibz, nsppol")
        # nctkarr_t("mats_tau", "dp", "two, mpw, mpw, ntau, nqibz, nsppol")
        iq_ibz, qpoint = self.find_qpoint(qpoint)
        npw_q = self.read_variable("chinpw_qibz")[iq_ibz]
        gvecs_q = self.read_variable("gvecs")[iq_ibz, 0:npw_q]
        ig1, g1 = self._find_ig_g(g1, gvecs_q)
        ig2, g2 = self._find_ig_g(g2, gvecs_q)

        varname = {"omega": "mats_w", "tau": "mats_tau"}[wt_space]
        if varname not in self.rootgrp.variables:
            raise ValueError(f"Cannot find {varname=} in rootgrp.variables")

        mat = self.read_variable(varname)[spin, iq_ibz, :, ig2, ig1]  # Note exchange of g1,g2
        mat = mat[..., 0] + 1j*mat[..., 1]

        return g1, g2, qpoint, mat

    def read_mat_all_ggp_at_qpt(self,
                                qpoint: KptSelect,
                                spin: int,
                                wt_space: str,
                                ):
        """
        Args:
            qpoint:
            spin: Spin index.
            wt_space:
        """
        iq_ibz, qpoint = self.find_qpoint(qpoint)
        npw_q = self.read_variable("chinpw_qibz")[iq_ibz]
        gvecs_q = self.read_variable("gvecs")[iq_ibz, 0:npw_q]

        varname = {"omega": "mats_w", "tau": "mats_tau"}[wt_space]
        if varname not in self.rootgrp.variables:
            raise ValueError(f"Cannot find {varname=} in rootgrp.variables")

        mat = self.read_variable(varname)[spin, iq_ibz]
        mat = mat[..., 0] + 1j*mat[..., 1]
        # (ntau, g2, g1) --> (g1, g2, ntau)
        mat = mat.transpose((2, 1, 0)).copy()

        return gvecs_q, qpoint, mat


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

    def find_tchim_qpoint(self, qpoint, tol=1e-6) -> tuple[int, np.ndarray]:
        """
        Find the index of the q-point in the IBZ. Return index and q-points.
        """
        tchi_qpoints = self.tchi_reader.read_variable("qibz")

        if duck.is_intlike(qpoint):
            return qpoint, tchi_qpoints

        for tchi_iq, tchi_qq in enumerate(tchi_qpoints):
            if np.all(abs(tchi_qq - qpoint) < tol):
                return tchi_iq, tchi_qpoints

        raise ValueError(f"Cannot find {qpoint=} in TCHIM file")

    @add_fig_kwargs
    def plot_qpoint_gpairs(self,
                           qpoint: KptSelect,
                           gpairs,
                           fontsize=8,
                           spins=(0, 0),
                           with_title: bool = True,
                           **kwargs) -> Figure:
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
            re_ax, im_ax = ax_mat[i]

            # number_of_qpoints_dielectric_function, number_of_frequencies_dielectric_function,
            # number_of_spins, number_of_spins, number_of_coefficients_dielectric_function,
            # number_of_coefficients_dielectric_function, complex)
            sus_data = sus_reader.read_variable(sus_var_name)[sus_iq, :, 0, 0, sus_ig2, sus_ig1]
            sus_data = sus_data[:, 0] + 1j * sus_data[:, 1]

            # nctkarr_t("mats_w", "dp", "two, mpw, mpw, ntau, nqibz, nsppol")
            tchi_data = self.tchi_reader.read_variable("mats_w")[0, tchi_iq, :, chi_ig2, chi_ig1, :]
            tchi_data = tchi_data[:, 0] + 1j * tchi_data[:, 1]

            style = dict(markersize=4)
            re_ax.plot(susc_iwmesh[1:], sus_data[1:].real, ls="--", marker="o", label="AW Re", **style)
            re_ax.plot(tchi_iwmesh, tchi_data[0:].real, ls=":", marker="^", label="minimax Re", **style)
            re_ax.grid(True)

            im_ax.plot(susc_iwmesh[1:], sus_data[1:].imag, ls="--", marker="o", label="AW Im", **style)
            im_ax.plot(tchi_iwmesh, tchi_data[0:].imag, ls=":", marker="^", label="minimax Im", **style)
            im_ax.grid(True)

            tex_g1, tex_g2, tex_qpt = r"${\bf{G}}_1$", r"${\bf{G}}_2$", r"${\bf{q}}$"
            re_ax.set_title(f"{tex_g1}: {g1}, {tex_g2}: {g2}, {tex_qpt}: {qpoint}", fontsize=fontsize)
            im_ax.set_title(f"{tex_g1}: {g1}, {tex_g2}: {g2}, {tex_qpt}: {qpoint}", fontsize=fontsize)

            re_ax.legend(loc="best", fontsize=fontsize, shadow=True)
            im_ax.legend(loc="best", fontsize=fontsize, shadow=True)

            #if i == 0:
            re_ax.set_ylabel(r'$\Re{\chi^0_{\bf{G}_1 \bf{G}_2}(i\omega)}$')
            im_ax.set_ylabel(r'$\Im{\chi^0_{\bf{G}_1 \bf{G}_2}(i\omega)}$')

            if i == len(gpairs) - 1:
                re_ax.set_xlabel(r'$i \omega$ (Ha)')
                im_ax.set_xlabel(r'$i \omega$ (Ha)')

            # The two files may store chi on different meshes.
            # Here we select the min value for comparison purposes.
            w_max = min(susc_iwmesh[-1], tchi_iwmesh[-1]) + 5
            xlims = [0, w_max]
            xlims = [-0.2, 7]

            set_axlims(re_ax, xlims, "x")
            set_axlims(im_ax, xlims, "x")

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
    def plot_mat_diff(self,
                      qpoint,
                      iw_index: int,
                      npwq: int = 101,
                      with_susmat=False,
                      cmap="jet",
                      fontsize=8,
                      spins=(0, 0),
                      **kwargs) -> Figure:
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

        # nctkarr_t("mats_w", "dp", "two, mpw, mpw, ntau, nqibz, nsppol")
        tchi_data = self.tchi_reader.read_variable("mats_w")[0, tchi_iq, iw_index]
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

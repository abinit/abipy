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
from typing import Any
from monty.collections import dict2namedtuple
from monty.functools import lazy_property
from monty.string import list_strings, marquee
from monty.termcolor import cprint
from abipy.core.func1d import Function1D
from abipy.core.structure import Structure
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import Kpoint, KpointList, Kpath, IrredZone, has_timrev_from_kptopt
from abipy.iotools import ETSF_Reader
from abipy.tools import duck
from abipy.tools.typing import Figure, KptSelect
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, Marker,
    set_axlims, set_ax_xylabels, set_visible, rotate_ticklabels, set_grid_legend, hspan_ax_line, Exposer)
from abipy.abio.robots import Robot
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.electrons.gw import SelfEnergy, QPState, QPList
from abipy.abio.enums import GWR_TASK
from abipy.tools.pade import pade, dpade, SigmaPade


__all__ = [
    "GwrFile",
    "GwrRobot",
]


class _MyQpkindsList(list):
    """Returned by find_qpkinds."""


@dataclasses.dataclass(kw_only=True)
class MinimaxMesh:
    """
    The minimax mesh stored in the GWR file.
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
        """Build a minimax mesh instance from a netcdf reader."""
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

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level `verbose`."""
        lines = []
        app = lines.append
        app(super().__str__())
        return "\n".join(lines)

    #def ft_w2t_even(self, f_w: np.ndarray) -> np.ndarray:
    #def ft_t2w_even(self, f_t: np.ndarray) -> np.ndarray:

    def get_ft_mptau(self, f_mpt: np.ndarray) -> np.ndarray:
        """
        Transform a signal in imaginary time to imaginary frequency with inhomogeneous FT.

        Args:
            f_mpt: numpy array of shape [2, ntau] or [2*ntau] where
            [0,:] corresponds to negative taus and [1,:] to positive taus.
        """
        if (ntau := len(f_mpt)) % 2 != 0:
            raise ValueError("Expecting an even number of points but got {ntau=}")
        f_mpt = np.reshape(f_mpt, (2, ntau // 2))

        # f(t) = E(t) + O(t) = (f(t) + f(-t)) / 2  + (f(t) - f(-t)) / 2
        even_t = (f_mpt[1] + f_mpt[0]) / 2.0
        odd_t = (f_mpt[1] - f_mpt[0]) / 2.0
        return self.cosft_wt @ even_t + 1j * self.sinft_wt @ odd_t

    @add_fig_kwargs
    def plot_ft_weights(self,
                        other: MinimaxMesh,
                        self_name: str ="self",
                        other_name: str ="other",
                        with_sinft: bool= False,
                        fontsize: int = 6,
                        **kwargs) -> Figure:
        """
        Plot the Fourier transform weights of two minimax meshes (self and other)

        Args:
            other:
            self_name:
            other_name:
            with_sinft:
            fontsize:
        """
        if self.ntau != other.ntau:
            raise ValueError(f"Cannot compare minimax meshes with different ntau: {self.ntau=}, {other.ntau=}")

        import matplotlib.pyplot as plt
        nrows, ncols = (4, 2) if with_sinft else (3, 2)
        fig, ax_mat = plt.subplots(nrows=nrows, ncols=ncols,
                                   sharex=False, sharey=False, squeeze=False,
                                   figsize=(12, 8),
                                  )

        I_mat = np.eye(self.ntau)
        select_irow = {
            0: [(self.cosft_wt @ self.cosft_tw) - I_mat,
                (other.cosft_wt @ other.cosft_tw) - I_mat], # , other.cosft_wt @ other.cosft_tw],
            1: [self.cosft_wt, other.cosft_wt], # self.cosft_wt - other.cosft_wt],
            2: [self.cosft_tw, other.cosft_tw], # self.cosft_tw - other.cosft_tw],
            #3: [self.sinft_tw, other.sinft_tw], # self.sinft_tw - other.sinft_tw],
        }

        label_irow = {
            0: [f"(cosft_wt @ cosft_tw) - I ({self_name})", f"(cosft_wt @ cosft_tw) - I ({other_name})"],
            1: [f"cosft_wt ({self_name})", f"cosft_wt ({other_name})"],
            2: [f"cosft_tw ({self_name})", f"cosft_tw ({other_name})"],
            #3: [f"sinft_tw ({self_name})", f"sinft_tw ({other_name})"],
        }

        for irow in range(nrows):
            for iax, (ax, data, label) in enumerate(zip(ax_mat[irow], select_irow[irow], label_irow[irow])):
                im = ax.matshow(data, cmap='seismic')
                fig.colorbar(im, ax=ax)
                ax.set_title(label, fontsize=fontsize)

        return fig


@dataclasses.dataclass(kw_only=True)
class PadeData:
    """Container for the Pade' results."""
    w_vals: np.ndarray
    sigxc_w: np.ndarray
    aw: np.ndarray
    e0: float
    ze0: float


@dataclasses.dataclass(kw_only=True)
class SigmaTauFit:
    """Stores the fit for Sigma(i tau)"""

    tau_mesh: np.ndarray
    values: np.ndarray
    a_mtau: complex   # Coefficient
    beta_mtau: float  # exp(beta tau)
    a_ptau: complex
    beta_ptau: float  # exp(-beta tau)

    def eval_omega(self, ws: np.ndarray) -> np.ndarray:
        r"""
        Compute the Fourier transform of the piecewise function.

            e^{-a t}, & t > 0 \quad (a > 0) \\
            e^{b t}, & t < 0 \quad (b > 0),

        that is:

            F(\omega) = \frac{1}{b - i \omega} + \frac{1}{a + i \omega}.
        """
        # TODO: Enforce time-ordering with zcut.
        cvals = np.empty(len(ws), dtype=complex)
        zcut = 1j * 0.1
        for iw, ww in enumerate(ws):
            cvals[iw] = self.a_mtau / (self.beta_mtau - ww - zcut) + self.a_ptau / (self.beta_ptau + ww + zcut)
        return cvals


class GwrSelfEnergy(SelfEnergy):
    """
    Extends SelfEnergy by composing it with a MinimaxMesh instance.
    """

    PADE_METHODS = [
      "abinit_pade",  # Pade results produced by ABINIT
      "abipy_pade",   # Pade results produced by AbiPy (should be equal to abinit_pade).
      "tau_fit",      # Fit Sigma_c(tau) and apply pade to the difference.
    ]

    def __init__(self, *args, **kwargs):
        mx_mesh = kwargs.pop("mx_mesh")
        super().__init__(*args, **kwargs)
        self.mx_mesh = mx_mesh

    def tau_fit(self, first, last, xs) -> tuple[np.ndarray, complex, float]:
        """
        Performs the exponential fit in imaginary time.
        """
        mp_taus = self.c_tau.mesh
        vals_mptaus = self.c_tau.values
        w0, f0 = mp_taus[first], vals_mptaus[first]
        wn, fn = mp_taus[last], vals_mptaus[last]
        # NB: take the real part of the log to avoid oscillatory behaviour in the exp.
        a = -np.log(fn / f0) / (wn - w0)
        a = a.real
        # If something goes wrong, disable the fit.
        # Note that the sign of a depends whether as we working with positive or negative tau.
        if wn >= 0 and a <= 1e-12: f0 = 0.0j
        if wn < 0 and a >= -1e-12: f0 = 0.0j
        #print(f"{f0=}")
        return f0 * np.exp(-a * (xs - w0)), f0, a

    def _minimize_loss_tau(self, tau_fit, zone: str):
        """
        Compute loss functions for all possible values of second_index.
        """
        mp_taus, vals_mptaus = self.c_tau.mesh, self.c_tau.values
        ntau = len(mp_taus) // 2
        if len(mp_taus) % 2 != 0:
            raise ValueError("Expecting an even number of points but got {len(mp_taus)=}")

        # Note weighted sum with minimax weights tau_wgs.
        losses = []
        if zone == "+":
            xs, ys = mp_taus[ntau:], vals_mptaus[ntau:]
            first = ntau
            for last in range(first+1, len(mp_taus)):
                ys_fit, alpha, beta = tau_fit(first, last, xs)
                losses.append((last, np.sum(self.mx_mesh.tau_wgs * np.abs(ys_fit - ys)**2), ys_fit, alpha, beta))

        elif zone == "-":
            xs, ys = mp_taus[:ntau], vals_mptaus[:ntau]
            first = ntau - 1
            for last in range(first):
                ys_fit, alpha, beta = tau_fit(first, last, xs)
                losses.append((last, np.sum(self.mx_mesh.tau_wgs * np.abs(ys_fit - ys)**2), ys_fit, alpha, beta))

        else:
            raise ValueError(f"Invalid {zone=} should be in (-, +)")

        # Find min of losses.
        min_loss = min(losses, key=lambda t: t[1])
        return dict2namedtuple(imin=min_loss[0], loss=min_loss[1], values=min_loss[2],
                               alpha=min_loss[3], beta=min_loss[4])

    def get_exp_tau_fit(self) -> SigmaTauFit:
        """
        Fit Sigma_c(i tau) values for tau < 0 and tau > 0 with two different exponentials.
        """
        fit_m = self._minimize_loss_tau(self.tau_fit, "-")
        fit_p = self._minimize_loss_tau(self.tau_fit, "+")

        return SigmaTauFit(tau_mesh=self.c_tau.mesh,
                           values=np.concatenate((fit_m.values, fit_p.values), axis=0),
                           a_mtau=fit_m.alpha,
                           beta_mtau=fit_m.beta,
                           a_ptau=fit_p.alpha,
                           beta_ptau=fit_p.beta)

    def get_pade_data(self, w_vals: np.ndarray, e0: float, pade_method: str) -> PadeData:
        """
        Use the Pade' method to get the self-energy on the real axis,
        the renormalization factor ze0 at the KS energy and the spectral function A(w).

        Args:
            w_vals: Energy values in eV. Ignored if pade_method == "abinit_pade".
            e0: bare KS energy used to compute ze0 (eV units)
            pade_method: String defining the Pade' algorithm.
        """
        if pade_method == "abinit_pade":
            # Get data from the file produced by ABINIT.
            sigc_w = self.xc.values - self.x_val
            aw = self.aw.values
            ze0 = self.ze0

        elif pade_method == "abipy_pade":
            # Compute the AC using the python version implemented in Abipy.
            # It should give the same results as abinit_pade as these functions
            # have been translated from Fortran to python.
            if not self.has_c_iw:
                raise ValueError("GwrSelfEnergy instance does not have data on the imaginary axis!")

            # if z_eval in 2 or 3 quadrant, avoid the branch cut in the complex plane using Sigma(-iw) = Sigma(iw)*.
            # See also sigma_pade_eval in m_dyson_solver.F90
            zs, f_zs = 1j * self.c_iw.mesh, self.c_iw.values
            spade = SigmaPade(zs, f_zs)
            sigc_w, dsigc_dw = spade.eval(w_vals)
            #ze0 = sigc_w = spade.eval_dz(e0)

            # FIXME
            ze0 = dpade(zs, f_zs, e0)
            aw = f_zs

            from pprint import pprint as p
            print("zs:\n", p(zs.tolist()))
            print("f_zs:\n", p(f_zs.tolist()))
            print("w_vals\n", p(w_vals.tolist()[:4]))
            print("sigc_w:\n", p(sigc_w.tolist()))
            import sys
            sys.exit(0)

        elif pade_method == "tau_fit":
            # Fit values for tau < 0, tau > 0 with two exponentials with real argument
            # and remove them from the signal.
            # Use Pade' to perform the AC of the difference and add analytic transform.
            # Trasform F(i tau) --> F(i ww) using minimax_mesh and inhomogeneous FT.
            # Ordering is first negative then positive tau values.
            fit = self.get_exp_tau_fit()
            diff_tau = self.c_tau.values - fit.values
            diff_iw = self.mx_mesh.get_ft_mptau(diff_tau)
            zs = 1j * self.c_iw.mesh
            spade = SigmaPade(zs, diff_iw)
            diff_sigc_w, diff_dsigc_dw = spade.eval(w_vals)

            sigc_w = diff_sigc_w + fit.eval_omega(w_vals)
            sigc_w = fit.eval_omega(w_vals)

            # FIXME
            aw = sigc_w
            ze0 = 0
            #raise NotImplementedError()

        else:
            raise ValueError(f"Invalid {pade_method=}, should be in {self.PADE_METHODS=}")

        # Add the static exchange part.
        sigxc_w = sigc_w + self.x_val
        return PadeData(w_vals=w_vals, sigxc_w=sigxc_w, aw=aw, e0=e0, ze0=ze0)

    @add_fig_kwargs
    def plot_pade(self,
                  pade_methods: list[str],
                  wmesh: None | np.ndarray = None,
                  ref_data: PadeData | None = None,
                  ax_mat=None,
                  fontsize: int = 8,
                  **kwargs) -> Figure:
        """
        Args:
            pade_methods: string or list of strings defining the Pade' algorithm.
            wmesh: Real frequency mesh. If None the internal mesh is used.
            ref_data: Reference data to compare with.
            ax_mat: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and title fontsize.
        """
        nrows, ncols = (3, 1)
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols)
        ax_list = np.array(ax_mat).flatten()
        ax_re, ax_im, ax_aw = ax_list

        if wmesh is None:
            wmesh = self.wmesh

        # FIXME
        e0 = 0.0
        for pade_method in list_strings(pade_methods):
           pdata = self.get_pade_data(wmesh, e0, pade_method)
           #print(pdata)
           ax_re.plot(pdata.w_vals, pdata.sigxc_w.real, label=pade_method)
           ax_im.plot(pdata.w_vals, pdata.sigxc_w.imag, label=pade_method)
           #ax_aw.plot(pdata.w_vals, pdata.aw, label=pade_method)

        if ref_data is not None:
           ax_re.plot(ref_data.w_vals, ref_data.sigxc_w.real, label="Ref")
           ax_im.plot(ref_data.w_vals, ref_data.sifxc_w.imag, lable="Ref")
           #ax_aw.plot(ref_data.w_vals, ref_data.aw, label="Ref")

        for ax in ax_list:
            set_grid_legend(ax, fontsize) #, xlabel="Iteration") #, ylabel=ylabel)

        ax_re.set_ylabel(r"$\Re\Sigma(\omega)$ (eV)")
        ax_im.set_ylabel(r"$\Im\Sigma(\omega)$ (eV)")
        ax_aw.set_ylabel(r"$A(\omega)$")

        return fig

    @add_fig_kwargs
    def plot_tau_fit(self, ax_list=None, fontsize: int = 8, **kwargs) -> Figure:
        """
        Plot Sig(i tau) and the exponential fit.

        Args:
            ax_list: List of |matplotlib-Axes| for plot. If None, new figure is produced.
            fontsize: Legend and title fontsize.
        """
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=2, ncols=1)
        ax_list = ax_re, ax_im = np.ravel(ax_list)
        fit = self.get_exp_tau_fit()

        # Plot data.
        self.plot_reimc_tau(ax_list, marker="o")
        ax_re.plot(fit.tau_mesh, fit.values.real)
        ax_im.plot(fit.tau_mesh, fit.values.imag)

        for ax in ax_list:
            set_grid_legend(ax, fontsize=fontsize)

        return fig


class GwrFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This object provides a high-level interface to the GWR.nc file produced by the GWR code.
    This file stores the QP energies, the self-energy along the imaginary/real frequency axis
    as well as metadata useful for performing convergence studies.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GwrFile
    """

    # Markers used for up/down bands.
    marker_spin = {0: "^", 1: "v"}

    #color_spin = {0: "k", 1: "r"}

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
        return self.r.ebands

    @property
    def sigma_kpoints(self) -> KpointList:
        """The k-points where the QP corrections have been calculated."""
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
       QP direct gaps in eV computed with the renormalization Z factor at the KS energy
       Shape: [nsppol, nkcalc]
       """
       return self.r.read_value("qpz_gaps") * abu.Ha_eV

    @lazy_property
    def minimax_mesh(self) -> MinimaxMesh:
        """Object storing the minimax mesh and weights."""
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
            raise TypeError(f"Don't know how to interpret {type(qp_kpoints)}")

        # Check indices
        errors = []
        eapp = errors.append
        for ikc in items[1]:
            if ikc >= self.nkcalc:
                eapp(f"K-point index {ikc} >= {self.nkcalc=}. Please check input qp_kpoints")

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
        app(f"gwr_task: {self.r.gwr_task}")

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

        if verbose:
            app(self.minimax_mesh.to_string(verbose=verbose))

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

    def get_dataframe_sk(self,
                         spin: int,
                         kpoint: KptSelect,
                         index=None,
                         ignore_imag: bool = False,
                         with_params: bool = True,
                         with_geo: bool = False) -> pd.Dataframe:
        """
        Returns a |pandas-DataFrame| with the QP results for the given (spin, k-point).

        Args:
            spin: Spin index
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
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

    def _get_include_bands(self, include_bands: Any, spin: int) -> set | None:
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

    def interpolate(self,
                    lpratio: int = 5,
                    ks_ebands_kpath: ElectronBands | None = None,
                    ks_ebands_kmesh: ElectronBands | None = None,
                    ks_degatol: float = 1e-4,
                    vertices_names: list[tuple] | None = None,
                    line_density: int = 20,
                    filter_params: list | None = None,
                    only_corrections: bool = False,
                    verbose: int = 0):
        """
        Interpolate the QP corrections in k-space on a k-path and, optionally, on a k-mesh
        using the star-functions method.

        Args:
            lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
                The default should be OK in many systems, larger values may be required for accurate derivatives.
            ks_ebands_kpath: KS |ElectronBands| on a k-path. If present,
                the routine interpolates the QP corrections and apply them on top of the KS band structure
                This is the recommended option because QP corrections are usually smoother than the
                QP energies and therefore easier to interpolate. If None, the QP energies are interpolated
                along the path defined by ``vertices_names`` and ``line_density``.
            ks_ebands_kmesh: KS |ElectronBands| on a homogeneous k-mesh. If present, the routine
                interpolates the corrections on the k-mesh (used to compute QP the DOS)
            ks_degatol: Energy tolerance in eV. Used when either ``ks_ebands_kpath`` or ``ks_ebands_kmesh`` are given.
                KS energies are assumed to be degenerate if they differ by less than this value.
                The interpolator may break band degeneracies (the error is usually smaller if QP corrections
                are interpolated instead of QP energies). This problem can be partly solved by averaging
                the interpolated values over the set of KS degenerate states.
                A negative value disables this ad-hoc symmetrization.
            vertices_names: Used to specify the k-path for the interpolated QP band structure
                when ``ks_ebands_kpath`` is None.
                It is a list of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            line_density: Number of points in the smallest segment of the k-path. Used with ``vertices_names``.
            filter_params: List with parameters used to filter high-frequency components (Eq 9 of PhysRevB.61.1639)
                First item gives rcut, second item sigma. Ignored if None.
            only_corrections: If True, the output contains the interpolated QP corrections instead of the QP energies.
                Available only if ks_ebands_kpath and/or ks_ebands_kmesh are used.
            verbose: Verbosity level.

        Returns:

            :class:`namedtuple` with the following attributes::

                * qp_ebands_kpath: |ElectronBands| with the QP energies interpolated along the k-path.
                * qp_ebands_kmesh: |ElectronBands| with the QP energies interpolated on the k-mesh.
                    None if ``ks_ebands_kmesh`` is not passed.
                * ks_ebands_kpath: |ElectronBands| with the KS energies interpolated along the k-path.
                * ks_ebands_kmesh: |ElectronBands| with the KS energies on the k-mesh..
                    None if ``ks_ebands_kmesh`` is not passed.
                * interpolator: |SkwInterpolator| object.
        """
        # TODO: Consistency check.
        errlines = []
        eapp = errlines.append
        if len(self.sigma_kpoints) != len(self.ebands.kpoints):
            eapp("QP energies should be computed for all k-points in the IBZ but nkibz != nkptgw")
        if len(self.sigma_kpoints) == 1:
            eapp("QP Interpolation requires nkptgw > 1.")
        #if (np.any(self.bstop_sk[0, 0] != self.bstop_sk):
        #    cprint("Highest bdgw band is not constant over k-points. QP Bands will be interpolated up to...")
        #if (np.any(self.bstart_sk[0, 0] != self.bstart_sk):
        #if (np.any(self.bstart_sk[0, 0] != 0):
        if errlines:
            raise ValueError("\n".join(errlines))

        # Get symmetries from abinit spacegroup (read from file).
        abispg = self.structure.abi_spacegroup
        fm_symrel = [s for (s, afm) in zip(abispg.symrel, abispg.symafm) if afm == 1]

        if ks_ebands_kpath is None:
            # Generate k-points for interpolation. Will interpolate all bands available in the sigres file.
            bstart, bstop = 0, -1
            if vertices_names is None:
                vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]
            kpath = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)
            kfrac_coords, knames = kpath.frac_coords, kpath.names

        else:
            # Use list of k-points from ks_ebands_kpath.
            ks_ebands_kpath = ElectronBands.as_ebands(ks_ebands_kpath)
            kfrac_coords = [k.frac_coords for k in ks_ebands_kpath.kpoints]
            knames = [k.name for k in ks_ebands_kpath.kpoints]

            # Find the band range for the interpolation.
            bstart, bstop = 0, ks_ebands_kpath.nband
            bstop = min(bstop, self.r.min_bstop)
            if ks_ebands_kpath.nband < self.r.min_bstop:
                cprint("Number of bands in KS band structure smaller than the number of bands in GW corrections", "red")
                cprint("Highest GW bands will be ignored", "red")

            if not ks_ebands_kpath.kpoints.is_path:
                cprint("Energies in ks_ebands_kpath should be along a k-path!", "red")

        # Interpolate QP energies if ks_ebands_kpath is None else interpolate QP corrections
        # and re-apply them on top of the KS band structure.
        gw_kcoords = [k.frac_coords for k in self.sigma_kpoints]

        # Read GW energies from file (real part) and compute corrections if ks_ebands_kpath.
        # This is the section in which the fileoformat (SIGRES.nc, GWR.nc) enters into play...

        # nctkarr_t("ze0_kcalc", "dp", "two, smat_bsize1, nkcalc, nsppol"), &
        # nctkarr_t("qpz_ene", "dp", "two, smat_bsize1, nkcalc, nsppol"), &
        # nctkarr_t("qp_pade", "dp", "two, smat_bsize1, nkcalc, nsppol"), &

        # where
        #   smat_bsize1 = gwr%b2gw - gwr%b1gw + 1
        #   smat_bsize2 = merge(1, gwr%b2gw - gwr%b1gw + 1, gwr%sig_diago)

        # Read QP energies
        varname = "qpz_ene"
        egw_rarr = self.r.read_value(varname, cmode="c").real * abu.Ha_eV
        if ks_ebands_kpath is not None:
            if ks_ebands_kpath.structure != self.structure:
                cprint("sigres.structure and ks_ebands_kpath.structures differ. Check your files!", "red")
            # Compute QP corrections
            egw_rarr -= (self.r.read_value("e0_kcalc") * abu.Ha_eV)

        # Note there's no guarantee that the sigma_kpoints and the corrections have the same k-point index.
        # Be careful because the order of the k-points and the band range stored in the SIGRES file may differ ...
        qpdata = np.empty(egw_rarr.shape)
        kcalc2ibz = self.r.read_value("kcalc2ibz")
        kpt2ibz = kcalc2ibz[0,:] - 1

        for ikcalc, gwk in enumerate(self.sigma_kpoints):
            #ik_ibz = self.r.kpt2ibz(gwk)
            ik_ibz = kpt2ibz[ikcalc]
            #print(f"{ikcalc=}, {ik_ibz=}")
            for spin in range(self.nsppol):
                qpdata[spin, ik_ibz, :] = egw_rarr[spin, ik_ibz, :]

        # Build interpolator for QP corrections.
        from abipy.core.skw import SkwInterpolator
        cell = (self.structure.lattice.matrix, self.structure.frac_coords, self.structure.atomic_numbers)
        qpdata = qpdata[:, :, bstart:bstop]
        has_timrev = has_timrev_from_kptopt(self.r.read_value("kptopt"))

        skw = SkwInterpolator(lpratio, gw_kcoords, qpdata, self.ebands.fermie, self.ebands.nelect,
                              cell, fm_symrel, has_timrev,
                              filter_params=filter_params, verbose=verbose)

        if ks_ebands_kpath is None:
            # Interpolate QP energies.
            eigens_kpath = skw.interp_kpts(kfrac_coords).eigens
        else:
            # Interpolate QP energies corrections and add them to KS.
            ref_eigens = ks_ebands_kpath.eigens[:, :, bstart:bstop]
            qp_corrs = skw.interp_kpts_and_enforce_degs(kfrac_coords, ref_eigens, atol=ks_degatol).eigens
            eigens_kpath = qp_corrs if only_corrections else ref_eigens + qp_corrs

        # Build new ebands object with k-path.
        kpts_kpath = Kpath(self.structure.reciprocal_lattice, kfrac_coords, weights=None, names=knames)
        occfacts_kpath = np.zeros(eigens_kpath.shape)

        # Finding the new Fermi level of the interpolated bands is not trivial, in particular if metals
        # because one should first interpolate the QP bands on a mesh. Here I align the QP bands
        # at the HOMO of the KS bands.
        homos = ks_ebands_kpath.homos if ks_ebands_kpath is not None else self.ebands.homos
        qp_fermie = max([eigens_kpath[e.spin, e.kidx, e.band] for e in homos])
        #qp_fermie = self.ebands.fermie; qp_fermie = 0.0

        qp_ebands_kpath = ElectronBands(self.structure, kpts_kpath, eigens_kpath, qp_fermie, occfacts_kpath,
                                        self.ebands.nelect, self.ebands.nspinor, self.ebands.nspden,
                                        smearing=self.ebands.smearing)

        qp_ebands_kmesh = None
        if ks_ebands_kmesh is not None:
            # Interpolate QP corrections on the same k-mesh as the one used in the KS run.
            ks_ebands_kmesh = ElectronBands.as_ebands(ks_ebands_kmesh)
            if bstop > ks_ebands_kmesh.nband:
                raise ValueError("Not enough bands in ks_ebands_kmesh, found %s, minimum expected %d\n" % (
                    ks_ebands_kmesh.nband, bstop))
            if ks_ebands_kpath.structure != self.structure:
                cprint("sigres.structure and ks_ebands_kpath.structures differ. Check your files!", "red")
            if not ks_ebands_kmesh.kpoints.is_ibz:
                cprint("Energies in ks_ebands_kmesh should be given in the IBZ", "red")

            # K-points and weights for DOS are taken from ks_ebands_kmesh.
            dos_kcoords = [k.frac_coords for k in ks_ebands_kmesh.kpoints]
            dos_weights = [k.weight for k in ks_ebands_kmesh.kpoints]

            # Interpolate QP corrections from bstart to bstop.
            ref_eigens = ks_ebands_kmesh.eigens[:, :, bstart:bstop]
            qp_corrs = skw.interp_kpts_and_enforce_degs(dos_kcoords, ref_eigens, atol=ks_degatol).eigens
            eigens_kmesh = qp_corrs if only_corrections else ref_eigens + qp_corrs

            # Build new ebands object with k-mesh.
            kpts_kmesh = IrredZone(self.structure.reciprocal_lattice, dos_kcoords, weights=dos_weights,
                                   names=None, ksampling=ks_ebands_kmesh.kpoints.ksampling)
            occfacts_kmesh = np.zeros(eigens_kmesh.shape)
            qp_ebands_kmesh = ElectronBands(self.structure, kpts_kmesh, eigens_kmesh, qp_fermie, occfacts_kmesh,
                                            self.ebands.nelect, self.ebands.nspinor, self.ebands.nspden,
                                            smearing=self.ebands.smearing)

        return dict2namedtuple(qp_ebands_kpath=qp_ebands_kpath,
                               qp_ebands_kmesh=qp_ebands_kmesh,
                               ks_ebands_kpath=ks_ebands_kpath,
                               ks_ebands_kmesh=ks_ebands_kmesh,
                               interpolator=skw,
                               )
    @add_fig_kwargs
    def plot_sigma_imag_axis(self,
                             kpoint: KptSelect,
                             spin: int = 0,
                             include_bands="gap",
                             with_tau: bool = True,
                             fontsize: int = 8,
                             ax_mat=None,
                             **kwargs) -> Figure:
        """
        Plot Sigma_nk(iw) along the imaginary axis for given k-point, spin and list of bands.

        Args:
            kpoint: k-point in self-energy. Accepts |Kpoint|, vector or index.
            spin: Spin index.
            include_bands: List of bands to include. None means all.
            with_tau:
            fontsize: Legend and title fontsize.
            ax_mat:
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
    def plot_sigma_real_axis(self,
                             kpoint: KptSelect,
                             spin: int = 0,
                             include_bands="gap",
                             fontsize: int = 8,
                             ax_mat=None,
                             **kwargs) -> Figure:
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
    def plot_qps_vs_e0(self,
                       with_fields="all",
                       exclude_fields=None,
                       e0="fermie",
                       xlims=None,
                       sharey=False,
                       ax_list=None,
                       fontsize=8,
                       **kwargs) -> Figure:
        """
        Plot the QP results stored in the GWR file as function of the KS energy.

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
    def plot_all_spectral_functions(self,
                                    include_bands=None,
                                    ax_mat=None,
                                    fontsize=8,
                                    **kwargs) -> Figure:
        """
        Plot the spectral function A_{nk}(w) for all k-points, bands and
        spins available in the GWR file.

        Args:
            include_bands: List of bands to include. None means all.
            ax_mat:
            fontsize: Legend and title fontsize.
        """
        # Build grid of plots.
        nrows, ncols = len(self.sigma_kpoints), self.nsppol
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        for ikcalc, kcalc in enumerate(self.sigma_kpoints):
            for spin in range(self.nsppol):
                ax = ax_mat[ikcalc, spin]
                self.plot_spectral_function(ikcalc, spin=spin, include_bands=include_bands,
                                            ax=ax, fontsize=fontsize, show=False)

        return fig

    @add_fig_kwargs
    def plot_spectral_function(self,
                               kpoint: KptSelect,
                               spin: int = 0,
                               include_bands=None,
                               ax=None,
                               fontsize=8,
                               **kwargs) -> Figure:
        """
        Plot the spectral function A_{nk}(w) for the given k-point, spin and bands.

        Args:
            include_bands: List of bands to include. None means all.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and title fontsize.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ikcalc, kpoint = self.r.get_ikcalc_kpoint(kpoint)

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
            ax.set_title("k-point: %s" % repr(kpoint), fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_tau_fit_sk(self,
                        spin: int,
                        kpoint: KptSelect,
                        fontsize: int = 8,
                        **kwargs) -> Figure:
        """
        Plot the ab-initio results and the fit in imaginary-time
        for all bands at the given kpoint and spin index.

        Args
            spin: Spin index.
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            fontsize: Legend and title fontsize.
        """
        ikcalc, kpoint = self.r.get_ikcalc_kpoint(kpoint)
        band_range = range(self.r.bstart_sk[spin, ikcalc], self.r.bstop_sk[spin, ikcalc])
        nrows, ncols = len(band_range), 2
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols)
        for ib, band in enumerate(band_range):
            ax_list = ax_mat[ib]
            sigma = self.r.read_sigee_skb(spin=spin, kpoint=kpoint, band=band)
            sigma.plot_tau_fit(ax_list=ax_list, show=False)

        fig.suptitle(r"$\Sigma_{nk}$" + f" at k-point: {kpoint}, spin: {spin}", fontsize=fontsize)

        return fig

    #@add_fig_kwargs
    #def plot_sig_mat(self, what, origin="lower", **kwargs):
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

        """
        print("In hacked yield_figs")
        sigma = self.r.read_sigee_skb(spin=0, kpoint=0, band=4)
        return sigma.plot_pade(["abinit_pade", "tau_fit"], show=False)
        #return sigma.plot_pade(["abinit_pade", "abipy_pade"], show=False)
        return None

        for ik in range(len(self.sigma_kpoints)):
            yield self.plot_tau_fit_sk(spin=0, kpoint=ik, show=False)
        return None
        """

        #include_bands = "all" if verbose else "gaps"
        #yield self.plot_spectral_functions(include_bands=include_bands, show=False)

        # TODO
        for spin in range(self.nsppol):
            for ik, kpoint in enumerate(self.sigma_kpoints):
                kws = dict(spin=spin, include_bands="gap", show=False)
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
        self.ebands = ElectronBands.from_file(path)
        self.structure = self.read_structure()

        self.nsppol = self.ebands.nsppol
        self.nwr = self.read_dimvalue("nwr")
        self.nkcalc = self.read_dimvalue("nkcalc")
        self.smat_bsize1 = self.read_dimvalue("smat_bsize1")
        self.smat_bsize2 = self.read_dimvalue("smat_bsize2")
        self.gwr_task = self.read_string("gwr_task")
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
        self.min_bstop = np.min(self.bstop_sk)

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

    def get_ikcalc_kpoint(self, kpoint: KptSelect) -> tuple[int, Kpoint]:
        """
        Return the ikcalc index and the Kpoint
        """
        ikcalc = self.kpt2ikcalc(kpoint)
        kpoint = self.sigma_kpoints[ikcalc]
        return ikcalc, kpoint

    def kpt2ikcalc(self, kpoint: KptSelect) -> int:
        """
        Return the index of the k-point in the sigma_kpoints array.
        Used to access data in the arrays that are dimensioned as [0:nkcalc].
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

    def read_sigee_skb(self, spin: int, kpoint: KptSelect, band: int) -> GwrSelfEnergy:
        """"
        Read self-energy for the given (spin, kpoint, band).
        """
        ikcalc, kpoint = self.get_ikcalc_kpoint(kpoint)
        ib = band - self.min_bstart
        ib2 = 0 if self.sig_diago else ib
        e0 = self.e0_kcalc[spin, ikcalc, ib]
        wmesh = self.get_wr_mesh(e0)

        # nctkarr_t("sigxc_rw_diag", "dp", "two, nwr, smat_bsize1, nkcalc, nsppol"), &
        xc_vals = self.read_variable("sigxc_rw_diag")[spin,ikcalc,ib,:,:] * abu.Ha_eV
        xc_vals = xc_vals[:,0] + 1j *xc_vals[:,1]

        # nctkarr_t("spfunc_diag", "dp", "nwr, smat_bsize1, nkcalc, nsppol") &
        aw_vals = self.read_variable("spfunc_diag")[spin,ikcalc,ib,:] / abu.Ha_eV

        # nctkarr_t("sigc_iw_mat", "dp", "two, ntau, smat_bsize1, smat_bsize2, nkcalc, nsppol"), &
        sigc_iw = self.read_value("sigc_iw_mat", cmode="c") * abu.Ha_eV
        c_iw_values = sigc_iw[spin, ikcalc, ib2, ib]

        # nctkarr_t("sigc_it_mat", "dp", "two, two, ntau, smat_bsize1, smat_bsize2, nkcalc, nsppol")
        sigc_tau = self.read_value("sigc_it_mat", cmode="c") * abu.Ha_eV
        c_tau_pm = sigc_tau[spin,ikcalc,ib2,ib]
        tau_mp_mesh = np.concatenate((-self.tau_mesh[::-1], self.tau_mesh))
        c_tau_mp_values = np.concatenate((c_tau_pm[::-1,1], c_tau_pm[:,0]))

        # nctkarr_t("sigx_mat", "dp", "two, smat_bsize1, smat_bsize2, nkcalc, nsppol")
        x_val = self.read_variable("sigx_mat")[spin, ikcalc, ib2, ib, 0] * abu.Ha_eV
        mx_mesh = MinimaxMesh.from_ncreader(self)
        # nctkarr_t("ze0_kcalc", "dp", "two, smat_bsize1, nkcalc, nsppol")
        ze0 = self.read_variable("ze0_kcalc")[spin, ikcalc, ib]
        ze0 = ze0[0] + 1j*ze0[1]

        return GwrSelfEnergy(spin, kpoint, band, wmesh, xc_vals, x_val, ze0, aw_vals,
                             iw_mesh=self.iw_mesh, c_iw_values=c_iw_values,
                             tau_mp_mesh=tau_mp_mesh, c_tau_mp_values=c_tau_mp_values,
                             mx_mesh=mx_mesh)

    def read_sigma_bdict_sikcalc(self, spin: int, ikcalc: int, include_bands: bool) -> dict[int, GwrSelfEnergy]:
        """
        Return dictionary of self-energy objects for the given (spin, ikcalc) indexed by the band index.
        """
        sigma_of_band = {}
        for band in range(self.bstart_sk[spin, ikcalc], self.bstop_sk[spin, ikcalc]):
            if include_bands and band not in include_bands: continue
            sigma_of_band[band] = self.read_sigee_skb(spin, ikcalc, band)

        return sigma_of_band

    def read_allqps(self, ignore_imag: bool = False) -> tuple[QPList]:
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

    def read_qplist_sk(self, spin: int, kpoint: KptSelect, band: int = None, ignore_imag: bool = False) -> QPList:
        """
        Read and return a QPList object for the given spin, kpoint.

        Args:
            spin: Spin index
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band: band index. If None all bands are considered.
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
        for sigma_band in range(self.bstart_sk[spin, ikcalc], self.bstop_sk[spin, ikcalc]):
            if band is not None and sigma_band != band: continue
            ib = sigma_band - self.min_bstart

            qpe = self.read_variable("qpz_ene")[spin, ikcalc, ib] * abu.Ha_meV
            qpe = qpe[0] + 1j*qpe[1]

            ze0 = self.read_variable("ze0_kcalc")[spin, ikcalc, ib]
            ze0 = ze0[0] + 1j*ze0[1]

            # TODO Finalize the implementation
            qp_list.append(QPState(
                spin=spin,
                kpoint=kpoint,
                band=sigma_band,
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

    def get_dataframe_sk(self,
                         spin: int,
                         kpoint: KptSelect,
                         with_params: bool = True,
                         ignore_imag: bool = False) -> pd.DataFrame:
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
    def plot_selfenergy_conv(self,
                             spin: int,
                             kpoint: KptSelect,
                             band: int,
                             axis: str = "wreal",
                             sortby=None,
                             hue=None,
                             colormap="viridis",
                             xlims=None,
                             fontsize: int = 8,
                             **kwargs) -> Figure:
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
    def plot_qpgaps_convergence(self,
                                qp_kpoints="all",
                                qp_type="qpz0",
                                sortby=None,
                                hue=None,
                                abs_conv=0.01,
                                plot_qpmks=True,
                                fontsize=8,
                                **kwargs) -> Figure:
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
    def plot_qpdata_conv_skb(self,
                             spin: int,
                             kpoint: KptSelect,
                             band: int,
                             sortby=None,
                             hue=None,
                             fontsize=8,
                             **kwargs) -> Figure:
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
        #ikc = nc0.kpt2ikcalc(kpoint)
        ikc, kpoint = nc0.r.get_ikcalc_kpoint(kpoint)
        kpoint = nc0.sigma_kpoints[ikc]

        # Sort and read QP data.
        if hue is None:
            labels, ncfiles, params = self.sortby(sortby, unpack=True)
            qplist = [ncfile.r.read_qp(spin, kpoint, band) for ncfile in ncfiles]
        else:
            groups = self.group_and_sortby(hue, sortby)
            qplist_group = []
            for g in groups:
                lst = [ncfile.r.read_qp(spin, kpoint, band) for ncfile in g.abifiles]
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
    def plot_qpfield_vs_e0(self,
                           field: str,
                           reim: str = "real",
                           function=lambda x: x,
                           sortby=None,
                           hue=None,
                           fontsize: int = 8,
                           colormap="jet",
                           e0="fermie",
                           **kwargs) -> Figure:
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

        # Visualize the convergence of the self-energy for all the k-points and the most important bands.
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

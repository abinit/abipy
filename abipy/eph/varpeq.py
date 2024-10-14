"""
This module contains objects for postprocessing polaron calculations
using the results stored in the VARPEQ.nc file.

For a theoretical introduction see ...
"""
from __future__ import annotations

import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from collections import defaultdict
from monty.string import marquee #, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.func1d import Function1D
from abipy.core.structure import Structure
from abipy.core.kpoints import kpoints_indices, kmesh_from_mpdivs, map_grid2ibz, IrredZone
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.typing import PathLike
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles, Marker, set_grid_legend)
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.dfpt.phonons import PhononBands
from abipy.dfpt.ddb import DdbFile
from abipy.tools.typing import Figure
from abipy.tools.numtools import BzRegularGridInterpolator, gaussian
from abipy.iotools import bxsf_write
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader


#TODO Finalize the implementation. Look at Pedro's implementation for GFR
#from abipy.electrons.effmass_analyzer import EffMassAnalyzer
#
#class FrohlichAnalyzer:
#
#    def __init__(self, gsr_kpath, ddb, verbose = 0, **anaddb_kwargs):
#        """
#        """
#        ebands_kpath = ElectronBands.as_ebands(gsr_kpath)
#        self.ddb = DdbFile.as_ddb(ddb)
#        self.verbose = verbose
#
#        r = self.ddb.anaget_epsinf_and_becs(verbose=verbose, return_input=True, **anaddb_kwargs)
#        self.epsinf, self.becs = r.epsinf, r.becs
#
#        self.diel_gen, diel_inp = self.ddb.anaget_dielectric_tensor_generator(verbose=verbose, return_input=True, **anaddb_kwargs)
#        self.eps0_tensor = self.diel_gen.tensor_at_frequency(0.0)
#
#        # Spherical average of eps_inf and eps_0 (real part only)
#        einf_savg, e0_savg = self.epsinf.trace() / 3, self.eps0_tensor.real.trace() / 3
#
#        self.kappa = 1 / (1/einf_savg - 1/e0_savg)
#
#        self.emana = EffMassAnalyzer(ebands_kpath)
#
#    def __str__(self) -> str:
#        return self.to_string()
#
#    def to_string(self, verbose: int = 0) -> str:
#        """String representation with verbosity level verbose"""
#        lines = []
#        app = lines.append
#
#        app("epsilon_infinity in Cartesian coordinates:")
#        app(str(self.epsinf))
#        app("BECS:")
#        app(str(self.becs))
#        app("eps0 tensor:")
#        app(str(self.eps0_tensor))
#        app(f"kappa = {self.kappa}")
#
#        return "\n".join(lines)
#
#    def analyze_band_edges(self):
#        self.emana.select_cbm()
#        self.emana.summarize()
#        self.emana.plot_emass()
#        self.emana.select_vbm()
#        self.emana.summarize()
#        self.emana.plot_emass()
#        #self.emana.select_band_edges()


@dataclasses.dataclass(kw_only=True)
class Entry:
    name: str
    latex: str
    info: str
    utype: str
    #color: str


# NB: All quantities are in atomic units!
_ALL_ENTRIES = [
    Entry(name="E_pol", latex=r'$E_{pol}$', info="Formation energy", utype="energy"),
    Entry(name="E_el", latex=r'$E_{el}$', info="Electronic part", utype="energy"),
    Entry(name="E_ph", latex=r'$E_{ph}$', info="Phonon part", utype="energy"),
    Entry(name="elph", latex=r'$E_{elph}$', info="e-ph term", utype="energy"),
    Entry(name="epsilon", latex=r"$\varepsilon$", info="Polaron eigenvalue", utype="energy"),
    Entry(name="gsr", latex=r"$|\nabla|$", info="||gradient||", utype="gradient"),
]

# Convert to dictionary: name --> Entry
_ALL_ENTRIES = {e.name: e for e in _ALL_ENTRIES}


class VarpeqFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file stores the results of a VARPEQ calculations: SCF cycle, A_nk, B_qnu
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        from abipy.eph.varpeq import VarpeqFile
        with VarpeqFile("out_VARPEQ.nc") as varpeq:
            print(varpeq)
            varpeq.plot_scf_cycle()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqFile
    """

    @classmethod
    def from_file(cls, filepath: PathLike) -> VarpeqFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = VarpeqReader(filepath)

    @lazy_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    def close(self) -> None:
        """Close the file."""
        self.r.close()

    @lazy_property
    def polaron_spin(self) -> list[Polaron]:
        """List of polaron objects, one for each spin (if any)."""
        return [Polaron.from_varpeq(self, spin) for spin in range(self.r.nsppol)]

    @lazy_property
    def params(self) -> dict:
        """dict with the convergence parameters, e.g. ``nbsum``."""
        r = self.r

        ksampling = self.ebands.kpoints.ksampling
        ngkpt, shifts = ksampling.mpdivs, ksampling.shifts
        nkbz = np.prod(ngkpt)

        od = dict([
            ("nkbz", nkbz),
            ("ngkpt", ngkpt),
            ("invsc_size", 1.0 / (nkbz * ((abu.Ang_Bohr * self.structure.lattice.volume) ** (1/3)))),
            ("frohl_ntheta", r.frohl_ntheta),
        ])
        return od

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

        app("")
        app("VARPEQ parameters:")
        app(f"varpeq_pkind: {self.r.varpeq_pkind}")
        #app(f"gstore_cplex: {self.r.cplex}")
        #app(f"gstore_kzone: {self.r.kzone}")
        #app(f"gstore_kfilter: {self.r.kfilter}")
        #app(f"gstore_kptopt: {self.r.kptopt}")
        #app(f"gstore_qptopt: {self.r.qptopt}")

        for spin in range(self.nsppol):
            polaron = self.polaron_spin[spin]
            df = polaron.get_final_results_df()
            app(f"Last SCF iteration. Energies in eV units")
            app(str(df))
            app("")

        return "\n".join(lines)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        for polaron in self.polaron_spin:
            yield polaron.plot_scf_cycle(show=False)

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("varpeq = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(varpeq)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


@dataclasses.dataclass(kw_only=True)
class Polaron:
    """
    This object stores the polaron coefficients A_kn, B_qnu for a given spin and all the nstate polaron states.
    Provides methods to plot |A_nk|^2 or |B_qnu|^2 together with the band structures (fatbands-like plots).
    """
    spin: int          # Spin index.
    nstates: int       # Number of polaronic states.
    nb: int            # Number of bands in A_kn,
    nk: int            # Number of k-points in A_kn, (including filtering if any)
    nq: int            # Number of q-points in B_qnu (including filtering if any)
    bstart: int        # First band starts at bstart
    bstop: int         # Last band (python convention)
    varpeq: VarpeqFile

    @classmethod
    def from_varpeq(cls, varpeq: VarpeqFile, spin: int) -> Polaron:
        """
        Build an istance from a VarpeqFile and the spin index.
        """
        r = varpeq.r
        nstates, nk, nq, nb = r.nstates, r.nk_spin[spin], r.nq_spin[spin], r.nb_spin[spin]
        bstart, bstop = r.brange_spin[spin]

        data = locals()
        return cls(**{k: data[k] for k in [field.name for field in dataclasses.fields(Polaron)]})

    @property
    def structure(self) -> Structure:
        """Crystalline structure."""
        return self.varpeq.structure

    @property
    def ebands(self) -> ElectronBands:
        """Electron bands."""
        return self.varpeq.ebands

    @lazy_property
    def kpoints(self) -> np.ndarray:
        """Reduced coordinates of the k-points."""
        return self.varpeq.r.read_value("kpts_spin")[self.spin, :self.nk]

    @lazy_property
    def qpoints(self) -> np.ndarray:
        """Reduced coordinates of the q-points."""
        return self.varpeq.r.read_value("qpts_spin")[self.spin, :self.nq]

    @lazy_property
    def a_kn(self) -> np.ndarray:
        """A_{pnk} coefficients for this spin."""
        return self.varpeq.r.read_value("a_spin", cmode="c")[self.spin, :self.nstates, :self.nk, :self.nb]

    @lazy_property
    def b_qnu(self) -> np.ndarray:
        """B_{pqnu} coefficients for this spin."""
        return self.varpeq.r.read_value("b_spin", cmode="c")[self.spin, :self.nstates, :self.nq]

    @lazy_property
    def scf_df_state(self) -> list[pd.DataFrame]:
        """
        List of dataframes with the SCF iterations. One dataframe for each polaron.
        NB: Energies are converted from Ha to eV.
        """
        # int nstep2cv_spin(nsppol, nstates) ;
        #   Number of steps to convergence at each state for each spin
        # double scf_hist_spin(nsppol, nstates, nstep, six) ;
        #   SCF optimization history at each state for each spin
        # int cvflag_spin(nsppol, nstates) ;
        #   0 --> calculation is not converged
        #   1 --> calculation is converged
        spin = self.spin
        r = self.varpeq.r
        nstep2cv = r.read_variable("nstep2cv_spin")[spin]
        scf_hist = r.read_variable("scf_hist_spin")[spin]
        cvflag = r.read_variable("cvflag_spin")[spin]

        def ufact_k(k):
            """Convert energies to eV"""
            entry = _ALL_ENTRIES[k]
            if entry.utype == "energy":
                return abu.Ha_eV
            if entry.utype == "gradient":
                return 1.0
            raise ValueError(f"Don't know how to convert {entry=}")

        # Build list of dataframe.
        df_list = []
        for pstate in range(self.nstates):
            n = nstep2cv[pstate]
            dct = {k: scf_hist[pstate, :n, i] * ufact_k(k) for i, k in enumerate(_ALL_ENTRIES)}
            df = pd.DataFrame(dct)
            # Add metadata to the attrs dictionary
            df.attrs["converged"] = bool(cvflag[pstate])
            df_list.append(df)

        return df_list

    def get_final_results_df(self, with_params: bool=False) -> pd.DataFrame:
        """
        Return daframe with the last iteration for all polaronic states.
        NB: Energies are in eV.
        """
        row_list = []
        for pstate in range(self.nstates):
            df = self.scf_df_state[pstate]
            row = {"pstate": pstate, "spin": self.spin}
            row.update(df.iloc[-1].to_dict())
            row["converged"] = df.attrs["converged"]
            if with_params:
                row.update(self.varpeq.params)
            row_list.append(row)

        return pd.DataFrame(row_list)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int=0) -> str:
        """
        String representation with verbosiy level verbose.
        """
        lines = []; app = lines.append

        app(marquee(f"Ank for spin: {self.spin}", mark="="))
        app(f"nstates: {self.nstates}")
        app(f"nk: {self.nk}")
        app(f"nq: {self.nq}")
        app(f"nb: {self.nb}")
        app(f"bstart: {self.bstart}, bstop: {self.bstop}")
        ksampling = self.ebands.kpoints.ksampling
        ngkpt, shifts = ksampling.mpdivs, ksampling.shifts
        app(f"ksampling: {str(ksampling)}")
        ngqpt = self.varpeq.r.ngqpt
        app(f"q-mesh: {ngqpt}")
        if verbose:
            for pstate in range(self.nstates):
                app("For %d: 1/N_k sum_nk |A_nk|^2: %f" % (pstate, np.sum(np.abs(self.a_kn[pstate]) ** 2) / self.nk))

        return "\n".join(lines)

    @lazy_property
    def ngkpt_and_shifts(self) -> tuple:
        """
        Return k-mesh divisions and shifts.
        """
        ksampling = self.ebands.kpoints.ksampling
        ngkpt, shifts = ksampling.mpdivs, ksampling.shifts

        if ngkpt is None:
            raise ValueError("Non diagonal k-meshes are not supported!")
        if len(shifts) > 1:
            raise ValueError("Multiple k-shifts are not supported!")

        return ngkpt, shifts

    def get_title(self, with_gaps: bool=True) -> str:
        """
        Return string with title for matplotlib plots.
        """
        varpeq = self.varpeq
        pre = "" if varpeq.ebands.nsppol == 1 else f"spin={self.spin}"
        if not with_gaps:
            return f"{pre}{varpeq.r.varpeq_pkind} polaron"
        else:
            gaps_string = varpeq.ebands.get_gaps_string()
            return f"{pre}{varpeq.r.varpeq_pkind} polaron, {gaps_string}"

    def insert_a_inbox(self, fill_value=None) -> tuple:
        """
        Return a_data, ngkpt, shifts where a_data is a
        (nstates, nb, nkx, nky, nkz)) array with A_{pnk} with p the polaron index.
        """
        # Need to know the shape of the k-mesh.
        ngkpt, shifts = self.ngkpt_and_shifts
        k_indices = kpoints_indices(self.kpoints, ngkpt)
        nx, ny, nz = ngkpt

        shape = (self.nstates, self.nb, nx, ny, nz)
        a_data = np.empty(shape, dtype=complex) if fill_value is None else np.full(shape, fill_value, dtype=complex)

        for ip in range(self.nstates):
            for ib in range(self.nb):
                for a_cplx, k_inds in zip(self.a_kn[ip,:,ib], k_indices, strict=True):
                    ix, iy, iz = k_inds
                    a_data[ip, ib, ix, iy, iz] = a_cplx

        return a_data, ngkpt, shifts

    def insert_b_inbox(self, fill_value=None) -> tuple:
        """
        Return b_data, ngqpt, shifts where b_data is a
        (nstates, nb, nqx, nqy, nqz)) array with B_{pqnu} with p the polaron index.
        """
        # Need to know the shape of the q-mesh (always Gamma-centered)
        ngqpt, shifts = self.varpeq.r.ngqpt, [0, 0, 0]
        q_indices = kpoints_indices(self.qpoints, ngqpt)

        natom3 = 3 * len(self.structure)
        nx, ny, nz = ngqpt

        shape = (self.nstates, natom3, nx, ny, nz)
        b_data = np.empty(shape, dtype=complex) if fill_value is None else np.full(shape, fill_value, dtype=complex)

        for ip in range(self.nstates):
            for nu in range(natom3):
                for b_cplx, q_inds in zip(self.b_qnu[ip,:,nu], q_indices, strict=True):
                    ix, iy, iz = q_inds
                    b_data[ip, nu, ix, iy, iz] = b_cplx

        return b_data, ngqpt, shifts

    def get_a2_interpolator_state(self) -> BzRegularGridInterpolator:
        """
        Build and return an interpolator for |A_nk|^2 for each polaronic state.
        """
        a_data, ngkpt, shifts = self.insert_a_inbox()

        return [BzRegularGridInterpolator(self.structure, shifts, np.abs(a_data[pstate])**2, method="linear")
                for pstate in range(self.nstates)]

    def get_b2_interpolator_state(self) -> BzRegularGridInterpolator:
        """
        Build and return an interpolator for |B_qnu|^2 for each polaronic state.
        """
        b_data, ngqpt, shifts = self.insert_b_inbox()

        return [BzRegularGridInterpolator(self.structure, shifts, np.abs(b_data[pstate])**2, method="linear")
                for pstate in range(self.nstates)]

    def write_a2_bxsf(self, filepath: PathLike, fill_value=0.0) -> None:
        r"""
        Export \sum_n |A_{pnk}|^2 in BXSF format suitable for visualization with xcrysden (use ``xcrysden --bxsf FILE``).
        Requires gamma-centered k-mesh.

        Args:
            filepath: BXSF filename.
        """
        # NB: the kmesh must be gamma-centered, multiple shifts are not supported.
        # Init values with 0. This is relevant only if kfiltering is being used.
        a_data, ngkpt, shifts = self.insert_a_inbox(fill_value=fill_value)

        # Compute \sum_n A^2_{pnk}.
        a2_data = np.sum(np.abs(a_data)**2, axis=1)
        fermie = a2_data.mean()

        bxsf_write(filepath, self.structure, 1, self.nstates, ngkpt, a2_data, fermie, unit="Ha")

    def write_b2_bxsf(self, filepath: PathLike, fill_value=0.0) -> None:
        r"""
        Export \sum_{\nu} |B_{q\nu}|^2 in BXSF format suitable for visualization with xcrysden (use ``xcrysden --bxsf FILE``).

        Args:
            filepath: BXSF filename.
        """
        # NB: qmesh must be gamma-centered, multiple shifts are not supported.
        # Init values with 0. This is relevant only if kfiltering is being used.
        b_data, ngqpt, shifts = self.insert_b_inbox(fill_value=fill_value)

        # Compute \sum_{\nu} B^2_{pq\nu}.
        b2_data = np.sum((np.abs(b_data)**2), axis=1)
        fermie = b2_data.mean()

        bxsf_write(filepath, self.structure, 1, self.nstates, ngqpt, b2_data, fermie, unit="Ha")

    #@add_fig_kwargs
    #def plot_bz_sampling(self, what="kpoints", fold=False,
    #                     ax=None, pmg_path=True, with_labels=True, **kwargs) -> Figure:
    #    """
    #    Plots a 3D representation of the Brillouin zone with the sampling.

    #    Args:
    #        what: "kpoints" or "qpoints"
    #        fold: whether the points should be folded inside the first Brillouin Zone.
    #            Defaults to False.
    #    """
    #    bz_points = dict(kpoints=self.kpoints, qpoints=self.qpoints)[what]
    #    kws = dict(ax=ax, pmg_path=pmg_path, with_labels=with_labels, fold=fold, kpoints=bz_points)

    #    return self.structure.plot_bz(show=False, **kws)

    @add_fig_kwargs
    def plot_scf_cycle(self, ax_mat=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the SCF cycle.

        Args:
            ax_max: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
        """
        # Build grid of plots.
        nrows, ncols = self.nstates, 2
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        for pstate in range(self.nstates):
            df = self.scf_df_state[pstate]
            niter = len(df)
            xs = np.arange(1, niter + 1)

            for iax, ax in enumerate(ax_mat[pstate]):
                # Create a twin Axes sharing the x-axis
                #grad_ax = ax
                grad_ax = ax.twinx()

                for ilab, (name, entry) in enumerate(_ALL_ENTRIES.items()):
                    # Convert energies to Hartree. Keep gradient as it is.
                    ys = df[name].to_numpy()

                    _ax, energy_like = ax, True
                    if entry.utype == "gradient":
                        _ax, energy_like = grad_ax, False

                    if iax == 0:
                        # Plot values linear scale.
                        _ax.plot(xs, ys, label=entry.latex)
                    else:
                        # Plot deltas in logscale.
                        _ax.plot(xs, np.abs(ys - ys[-1]), label=entry.latex)
                        _ax.set_yscale("log")

                _ax.set_xlim(1, niter)

                if energy_like:
                    ylabel = "Energy (eV)" if iax == 0 else r"$|\Delta E|$ (eV)"
                else:
                    ylabel = r"$|\nabla|$" if iax == 0 else r"$|\Delta |\nabla|$"

                set_grid_legend(_ax, fontsize, xlabel="Iteration") #, ylabel=ylabel)
                _ax.set_ylabel(ylabel)

        fig.suptitle(self.get_title(with_gaps=True))
        fig.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_ank_with_ebands(self, ebands_kpath,
                             ebands_kmesh=None, lpratio: int=5, method="gaussian", step: float=0.05, width: float=0.1,
                             nksmall: int=20, normalize: bool=False, with_title=True,
                             ax_mat=None, ylims=None, scale=10, marker_color="gold", fontsize=12, **kwargs) -> Figure:
        """
        Plot electron bands with markers with size proportional to |A_nk|^2.

        Args:
            ebands_kpath: ElectronBands or Abipy file providing an electronic band structure along a k-path.
            ebands_kmesh: ElectronBands or Abipy file providing an electronic band structure with k in the IBZ.
            nksmall:
            normalize: Rescale the two DOS to plot them on the same scale.
            lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
                The default should be OK in many systems, larger values may be required for accurate derivatives.
            method: Integration scheme for DOS
            step: Energy step (eV) of the linear mesh for DOS computation.
            width: Standard deviation (eV) of the gaussian for DOS computation.
            with_title: True to add title with chemical formula and gaps.
            ax_mat: List of |matplotlib-Axes| or None if a new figure should be created.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            scale: Scaling factor for |A_nk|^2.
            marker_color: Color for markers.
            fontsize: fontsize for legends and titles
        """
        nrows, ncols = self.nstates, 2
        gridspec_kw = {'width_ratios': [2, 1]}
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=True, squeeze=False, gridspec_kw=gridspec_kw)
        # Get interpolators for A_nk
        a2_interp_state = self.get_a2_interpolator_state()

        # DEBUG SECTION
        #ref_kn = np.abs(self.a_kn) ** 2
        #for ik, kpoint in enumerate(self.kpoints):
        #    interp = a2_interp_state[0].eval_kpoint(kpoint)
        #    print("MAX (A2 ref - A2 interp) at qpoint", kpoint)
        #    print((np.abs(ref_kn[ik] - interp)).max())

        df = self.get_final_results_df()

        ebands_kpath = ElectronBands.as_ebands(ebands_kpath)
        ymin, ymax = +np.inf, -np.inf
        for pstate in range(self.nstates):
            x, y, s = [], [], []
            for ik, kpoint in enumerate(ebands_kpath.kpoints):
                enes_n = ebands_kpath.eigens[self.spin, ik, self.bstart:self.bstop]
                for e, a2 in zip(enes_n, a2_interp_state[pstate].eval_kpoint(kpoint), strict=True):
                    x.append(ik); y.append(e); s.append(scale * a2)
                    ymin, ymax = min(ymin, e), max(ymax, e)

            points = Marker(x, y, s, color=marker_color, edgecolors='gray', alpha=0.8, label=r'$|A_{n\mathbf{k}}|^2$')
            ax = ax_mat[pstate, 0]
            ebands_kpath.plot(ax=ax, points=points, show=False)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

            # Add energy and convergence status
            with_info = True
            if with_info:
                data = (df[df["pstate"] == pstate]).to_dict(orient="list")
                e_pol_ev, converged = float(data["E_pol"][0]), bool(data["converged"][0])
                title = f"Formation energy: {e_pol_ev:.3f} eV, {converged=}"
                ax.set_title(title, fontsize=8)

            if pstate != self.nstates - 1:
                set_visible(ax, False, *["legend", "xlabel"])

        vertices_names = [(k.frac_coords, k.name) for k in ebands_kpath.kpoints]

        if ebands_kmesh is None:
            print(f"Computing ebands_kmesh with star-function interpolation and {nksmall=} ...")
            edos_ngkpt = self.structure.calc_ngkpt(nksmall)
            r = self.ebands.interpolate(lpratio=lpratio, vertices_names=vertices_names, kmesh=edos_ngkpt)
            ebands_kmesh = r.ebands_kmesh

        # Get electronic DOS from ebands_kmesh.
        edos_kws = dict(method=method, step=step, width=width)
        edos = ebands_kmesh.get_edos(**edos_kws)
        edos_mesh = edos.spin_dos[self.spin].mesh
        e0 = self.ebands.fermie

        ##################
        # Compute A_nk DOS
        ##################
        # NB: A_nk does not necessarily have the symmetry of the lattice so we have to loop over the full BZ.
        # Here we get the mapping BZ --> IBZ needed to obtain the KS eigenvalues e_nk from the IBZ for the DOS.
        kmesh = ebands_kmesh.get_bz2ibz_bz_points()

        for pstate in range(self.nstates):
            # Compute A^2(E) DOS with A_nk in the full BZ
            ank_dos = np.zeros(len(edos_mesh))
            for ik_ibz, kpoint in zip(kmesh.bz2ibz, kmesh.bz_kpoints):
                enes_n = ebands_kmesh.eigens[self.spin, ik_ibz, self.bstart:self.bstop]
                a2_n = a2_interp_state[pstate].eval_kpoint(kpoint)
                for e, a2 in zip(enes_n, a2_n):
                    ank_dos += a2 * gaussian(edos_mesh, width, center=e-e0)
            ank_dos /= np.product(kmesh.ngkpt)
            ank_dos = Function1D(edos_mesh, ank_dos)
            print(f"For {pstate=}, A^2(E) integrates to:", ank_dos.integral_value, " Ideally, it should be 1.")

            ax = ax_mat[pstate, 1]
            edos_opts = {"color": "black",} if self.spin == 0 else {"color": "red"}
            edos.plot_ax(ax, e0, spin=self.spin, normalize=normalize, exchange_xy=True, label="eDOS(E)", **edos_opts)
            ank_dos.plot_ax(ax, exchange_xy=True, normalize=normalize, label=r"$A^2$(E)", color=marker_color)

            # Computes A2(E) using only k-points in the IBZ. This is just for testing.
            # A2_IBZ(E) should be equal to A2(E) only if A_nk fullfills the lattice symmetries. See notes above.
            with_ibz_a2dos = True
            if with_ibz_a2dos:
                ank_dos = np.zeros(len(edos_mesh))
                for ik_ibz, kpoint in enumerate(ebands_kmesh.kpoints):
                    weight = kpoint.weight
                    enes_n = ebands_kmesh.eigens[self.spin, ik_ibz, self.bstart:self.bstop]
                    for e, a2 in zip(enes_n, a2_interp_state[pstate].eval_kpoint(kpoint), strict=True):
                        ank_dos += weight * a2 * gaussian(edos_mesh, width, center=e-e0)
                ank_dos = Function1D(edos_mesh, ank_dos)
                print(f"For {pstate=}, A2_IBZ(E) integrates to:", ank_dos.integral_value, " Ideally, it should be 1.")
                ank_dos.plot_ax(ax, exchange_xy=True, normalize=normalize, label=r"$A^2_{IBZ}$(E)", color=marker_color, ls="--")

            set_grid_legend(ax, fontsize, xlabel="Arb. unit")

            if pstate != self.nstates - 1:
                set_visible(ax, False, *["legend", "xlabel"])

        if ylims is None:
            # Automatic ylims.
            ymin -= 0.1 * abs(ymin)
            ymax += 0.1 * abs(ymax)
            ylims = [ymin - e0, ymax - e0]

        for ax in ax_mat.ravel():
            set_axlims(ax, ylims, "y")

        if with_title:
            fig.suptitle(self.get_title(with_gaps=True))

        return fig

    @add_fig_kwargs
    def plot_bqnu_with_ddb(self, ddb, with_phdos=True, anaddb_kwargs=None, **kwargs) -> Figure:
        """
        High-level interface to plot phonon energies with markers with size proportional to |B_qnu|^2.
        Similar to plot_bqnu_with_phbands but this function receives in input a DdbFile or a
        path to a ddb file and automates the computation of the phonon bands by invoking anaddb.

        Args:
            ddb: DdbFile or path to file.
            with_phdos: True if phonon DOS should be computed and plotter.
            anaddb_kwargs: Optional arguments passed to anaddb.
            kwargs: Optional arguments passed to plot_bqnu_with_phbands.
        """
        ddb = DdbFile.as_ddb(ddb)
        anaddb_kwargs = {} if anaddb_kwargs is None else anaddb_kwargs

        with ddb.anaget_phbst_and_phdos_files(**anaddb_kwargs) as g:
            phbst_file, phdos_file = g[0], g[1]
            phbands_qpath = phbst_file.phbands
            return self.plot_bqnu_with_phbands(phbands_qpath,
                                               phdos_file=phdos_file if with_phdos else None,
                                               ddb=ddb, **kwargs)

    @add_fig_kwargs
    def plot_bqnu_with_phbands(self, phbands_qpath,
                               phdos_file=None, ddb=None, width=0.001, normalize: bool=True,
                               verbose=0, anaddb_kwargs=None, with_title=True,
                               ax_mat=None, scale=10, marker_color="gold", fontsize=12, **kwargs) -> Figure:
        """
        Plot phonon energies with markers with size proportional to |B_qnu|^2.

        Args:
            phbands_qpath: PhononBands or Abipy file providing a phonon band structure.
            phdos_file:
            ddb: DdbFile or path to file.
            width: Standard deviation (eV) of the gaussian.
            normalize: Rescale the two DOS to plot them on the same scale.
            verbose:
            anaddb_kwargs: Optional arguments passed to anaddb.
            with_title: True to add title with chemical formula and gaps.
            ax_mat: List of |matplotlib-Axes| or None if a new figure should be created.
            scale: Scaling factor for |B_qnu|^2.
            marker_color: Color for markers.
            fontsize: fontsize for legends and titles.
        """
        with_phdos = phdos_file is not None  and ddb is not None
        nrows, ncols, gridspec_kw = self.nstates, 1, None
        if with_phdos:
            ncols, gridspec_kw = 2, {'width_ratios': [2, 1]}

        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=True, squeeze=False, gridspec_kw=gridspec_kw)

        phbands_qpath = PhononBands.as_phbands(phbands_qpath)

        # Get interpolators for B_qnu
        b2_interp_state = self.get_b2_interpolator_state()

        for pstate in range(self.nstates):
            x, y, s = [], [], []
            for iq, qpoint in enumerate(phbands_qpath.qpoints):
                omegas_nu = phbands_qpath.phfreqs[iq,:]
                for w, b2 in zip(omegas_nu, b2_interp_state[pstate].eval_kpoint(qpoint), strict=True):
                    x.append(iq); y.append(w); s.append(scale * b2)

            ax = ax_mat[pstate, 0]
            points = Marker(x, y, s, color=marker_color, edgecolors='gray', alpha=0.8, label=r'$|B_{\nu\mathbf{q}}|^2$')
            phbands_qpath.plot(ax=ax, points=points, show=False)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

            if pstate != self.nstates - 1:
                set_visible(ax, False, *["legend", "xlabel"])

        if not with_phdos:
            if with_title: fig.suptitle(self.get_title(with_gaps=True))
            return fig

        ####################
        # Compute B_qnu DOS
        ####################
        # NB: The B_qnu do not necessarily have the symmetry of the lattice so we have to loop over the full BZ.
        # Add phdos and |B_qn| dos. Mesh is given in eV, values are in states/eV.
        phdos = phdos_file.phdos
        phdos_ngqpt = np.diagonal(phdos_file.qptrlatt) # Use same q-mesh as phdos
        phdos_shifts = [0.0, 0.0, 0.0]
        phdos_nqbz = np.product(phdos_ngqpt)
        wmesh = phdos.mesh

        # Here we get the mapping BZ --> IBZ needed to obtain the ph frequencies omega_qnu from the IBZ for the DOS.
        bz_qpoints = kmesh_from_mpdivs(phdos_ngqpt, phdos_shifts)
        #bz2ibz = map_grid2ibz(self.structure, ibz_qpoints, phdos_ngqpt, has_timrev=True)

        #ibz = IrredZone.from_ngkpt(self.structure, phdos_ngqpt, phdos_shifts, kptopt=1)

        # Call anaddb (again) to get phonons on the nqpt mesh.
        anaddb_kwargs = {} if anaddb_kwargs is None else anaddb_kwargs
        phbands_qmesh = ddb.anaget_phmodes_at_qpoints(qpoints=bz_qpoints, ifcflag=1, verbose=verbose, **anaddb_kwargs)
        if len(phbands_qmesh.qpoints) != np.product(phdos_ngqpt):
            raise RuntimeError(f"{len(phbands_qmesh.qpoints)=} != {np.product(phdos_ngqpt)=}")
        #print(phbands_qmesh.qpoints)

        for pstate in range(self.nstates):
            # TODO New version using the BZ. Requires new VARPEQ.nc file with all symmetries
            # The B_qnu do not necessarily have the symmetry of the lattice so we have to loop over the full BZ.
            # Get mapping BZ --> IBZ needed to obtain the KS eigenvalues e_nk from the IBZ for the DOS
            """
            bqnu_dos = np.zeros(len(wmesh))
            for iq_ibz, qpoint in zip(bz2ibz, bz_qpoints):
                freqs_nu = phbands_qmesh.phfreqs[iq_ibz]
                for w, b2 in zip(freqs_nu, b2_interp_state[pstate].eval_kpoint(qpoint), strict=True)
                    bqnu_dos += b2 gaussian(wmesh, width, center=w)
            bqnu_dos /= np.product(phdos_ngqpt)
            """

            # Compute B2(E) using only q-points in the IBZ. This is just for testing.
            # B2_IBZ(E) should be equal to B2(E) only if B_qnu fullfill the lattice symmetries. See notes above.
            #with_ibz_b2dos = True
            #if with_ibz_b2dos:
            bqnu_dos = np.zeros(len(wmesh))
            for iq, qpoint in enumerate(phbands_qmesh.qpoints):
                weight = qpoint.weight
                #print(weight, 1.0/phdos_nqbz)
                freqs_nu = phbands_qmesh.phfreqs[iq]
                for w, b2 in zip(freqs_nu, b2_interp_state[pstate].eval_kpoint(qpoint), strict=True):
                    bqnu_dos += weight * b2 * gaussian(wmesh, width, center=w)

            bqnu_dos = Function1D(wmesh, bqnu_dos)

            ax = ax_mat[pstate, 1]
            phdos.plot_ax(ax, exchange_xy=True, normalize=normalize, label="phDOS(E)", color="black")
            bqnu_dos.plot_ax(ax, exchange_xy=True, normalize=normalize, label=r"$B^2$(E)", color=marker_color)
            set_grid_legend(ax, fontsize, xlabel="Arb. unit")

            if pstate != self.nstates - 1:
                set_visible(ax, False, *["legend", "xlabel"])

        if with_title:
            fig.suptitle(self.get_title(with_gaps=True))

        return fig


class VarpeqReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqReader
    """

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        # Netcdf Variables

        # int eph_task ;
        # int nkbz ;
        # int nqbz ;
        # int frohl_ntheta ;
        # double varpeq_tolgrs ;
        # double e_frohl ;
        # char varpeq_pkind(fnlen) ;
        # char varpeq_aseed(fnlen) ;
        # int ngkpt(three) ;
        # int gstore_ngqpt(three) ;
        # int nk_spin(nsppol) ;
        # int nq_spin(nsppol) ;
        # int nb_spin(nsppol) ;
        # int brange_spin(nsppol, two) ;
        # int cvflag_spin(nsppol, nstates) ;
        # int nstep2cv_spin(nsppol, nstates) ;
        # double scf_hist_spin(nsppol, nstates, nstep, six) ;
        # double kpts_spin(nsppol, max_nk, three) ;
        # double qpts_spin(nsppol, max_nq, three) ;
        # double cb_min_spin(nsppol) ;
        # double vb_max_spin(nsppol) ;
        # double a_spin(nsppol, nstates, max_nk, max_nb, two) ;
        # double b_spin(nsppol, nstates, max_nq, natom3, two) ;

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("nsppol")
        self.nstates = self.read_dimvalue("nstates")
        self.nk_spin = self.read_value("nk_spin")
        self.nq_spin = self.read_value("nq_spin")
        #self.nkbz = self.read_dimvalue("nkbz")
        #self.nqbz = self.read_dimvalue("nqbz")
        self.varpeq_pkind = self.read_string("varpeq_pkind")
        #self.varpeq_aseed = self.read_string("varpeq_aseed")
        self.ngqpt = self.read_value("gstore_ngqpt")
        self.frohl_ntheta = self.read_value("frohl_ntheta")

        # Read important variables.
        #self.completed = self.read_value("gstore_completed")
        #self.done_spin_qbz = self.read_value("gstore_done_qbz_spin")
        #self.qptopt = self.read_value("gstore_qptopt")
        #self.kptopt = self.read_value("kptopt")
        #self.kzone = self.read_string("gstore_kzone")
        #self.qzone = self.read_string("gstore_qzone")
        #self.kfilter = self.read_string("gstore_kfilter")
        #self.gmode = self.read_string("gstore_gmode")

        # Note conversion Fortran --> C for the isym index.
        self.brange_spin = self.read_value("brange_spin")
        self.brange_spin[:,0] -= 1
        self.nb_spin = self.brange_spin[:,1] - self.brange_spin[:,0]

        #self.erange_spin = self.read_value("gstore_erange_spin")
        # Total number of k/q points for each spin after filtering (if any)
        #self.glob_spin_nq = self.read_value("gstore_glob_nq_spin")
        #self.glob_nk_spin = self.read_value("gstore_glob_nk_spin")


class VarpeqRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple VARPEQ.nc files.

    Usage example:

    .. code-block:: python

        robot = VarpeqRobot.from_files([
            "out1_VARPEQ.nc",
            "out2_VARPEQ.nc",
            ])

        print(robot)
        df = robot.get_final_results_df()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqRobot
    """

    EXT = "VARPEQ"

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append
        df = self.get_final_results_df()
        lines.append(str(df))

        return "\n".join(lines)

    def get_final_results_df(self, spin=None, sortby=None, with_params: bool=True) -> pd.DataFrame:
        """
        Return dataframe with the last iteration for all polaronic states.
        NB: Energies are in eV.

        Args:
            spin:
            sortby: Name to sort by.
            with_params:
        """
        df_list = []
        for abifile in self.abifiles:
            if spin is None:
                for polaron in abifile.polaron_spin:
                    df_list.append(polaron.get_final_results_df(with_params=with_params))
            else:
                polaron = abifile.polaron_spin[spin]
                df_list.append(polaron.get_final_results_df(with_params=with_params))

        df = pd.concat(df_list)
        if sortby and sortby in df: df = df.sort_values(sortby)
        return df

    @add_fig_kwargs
    def plot_kconv(self, colormap="jet", fontsize=12, **kwargs) -> Figure:
        """
        Plot the convergence of the results wrt to the k-point sampling.

        Args:
            colormap: Color map. Have a look at the colormaps here and decide which one you like:
            fontsize: fontsize for legends and titles
        """
        nsppol = self.getattr_alleq("nsppol")

        # Build grid of plots.
        nrows, ncols = len(_ALL_ENTRIES), nsppol
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        cmap = plt.get_cmap(colormap)
        for spin in range(nsppol):
            df = self.get_final_results_df(spin=spin, sortby=None)
            xs = df["invsc_size"]
            xvals = np.linspace(0.0, 1.1 * xs.max(), 100)

            for ix, ylabel in enumerate(_ALL_ENTRIES):
                ax = ax_mat[ix, spin]
                ys = df[ylabel]

                # Plot ab-initio points.
                ax.scatter(xs, ys, color="red", marker="o")

                # Plot fit using the first nn points.
                for nn in range(1, len(xs)):
                    color = cmap((nn - 1) / len(xs))
                    p = np.poly1d(np.polyfit(xs[:nn+1], ys[:nn+1], deg=1))
                    ax.plot(xvals, p(xvals), color=color, ls="--")

                xlabel = "Inverse supercell size (Bohr$^-1$)" if ix == len(_ALL_ENTRIES) - 1 else None
                set_grid_legend(ax, fontsize, xlabel=xlabel, ylabel=f"{ylabel} (eV)", legend=False)
                ax.tick_params(axis='x', color='black', labelsize='20', pad=5, length=5, width=2)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        #yield self.plot_scf_cycle(show=False)
        yield self.plot_kconv()

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            nbv.new_code_cell("robot = abilab.VarpeqRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

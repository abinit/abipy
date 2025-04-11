"""
This module contains objects for post-processing polaron calculations
using the results stored in the VPQ.nc file.

For a theoretical introduction see ...
"""
from __future__ import annotations

import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from monty.string import marquee
from monty.functools import lazy_property
from scipy.interpolate import interp1d
#from monty.termcolor import cprint
from abipy.core.func1d import Function1D
from abipy.core.structure import Structure
from abipy.core.kpoints import kpoints_indices, kmesh_from_mpdivs, map_grid2ibz
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.typing import PathLike
from abipy.tools.plotting import (add_fig_kwargs, get_axarray_fig_plt, set_axlims, set_visible,
    Marker, set_grid_legend)
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
    name: str   # Entry name
    latex: str  # Latex label
    info: str   # Description string
    utype: str  # Unit type
    #color: str


# NB: All quantities are in atomic units!
_ALL_ENTRIES = [
    Entry(name="E_pol", latex=r'$E_{pol}$', info="Formation energy", utype="energy"),
    Entry(name="E_el", latex=r'$E_{el}$', info="Electronic part", utype="energy"),
    Entry(name="E_ph", latex=r'$E_{ph}$', info="Phonon part", utype="energy"),
    Entry(name="elph", latex=r'$E_{elph}$', info="e-ph term", utype="energy"),
    Entry(name="epsilon", latex=r"$\varepsilon$", info="Polaron eigenvalue", utype="energy"),
    Entry(name="grs", latex=r"$|\nabla|$", info="||gradient||", utype="gradient"),
]

# Convert to dictionary: name --> Entry
_ALL_ENTRIES = {e.name: e for e in _ALL_ENTRIES}


class VpqFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file stores the results of a VPQ calculations: SCF cycle, A_nk, B_qnu coefficients
    It also provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        from abipy.eph.vpq import VpqFile
        with VpqFile("out_VPQ.nc") as vpq:
            print(vpq)
            for polaron in vpq.polaron_spin:
                print(polaron)
                polaron.plot_scf_cycle()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VpqFile
    """

    @classmethod
    def from_file(cls, filepath: PathLike) -> VpqFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = VpqReader(filepath)

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
        """List of Polaron objects, one for each spin (if any)."""
        return [Polaron.from_vpq(self, spin) for spin in range(self.r.nsppol)]

    @lazy_property
    def params(self) -> dict:
        """dict with the convergence parameters, e.g. ``nbsum``."""
        r = self.r

        ksampling = self.ebands.kpoints.ksampling
        ngkpt, shifts = ksampling.mpdivs, ksampling.shifts
        nkbz = np.prod(ngkpt)

        avg_g = r.read_value("vpq_avg_g")
        e_frohl = r.read_value("e_frohl") # in Ha

        d = dict(
            avg_g = bool(avg_g),
            e_frohl = e_frohl * abu.Ha_eV,
            ngkpt=tuple(ngkpt),
            inv_k = 1. / np.cbrt(nkbz),
            invsc_linsize = 1. / np.cbrt(nkbz * self.structure.lattice.volume),
        )

        return d

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))

        app("")
        app("VPQ parameters:")
        app(f"vpq_pkind: {self.r.vpq_pkind}")
        #app(f"gstore_cplex: {self.r.cplex}")
        #app(f"gstore_kzone: {self.r.kzone}")
        #app(f"gstore_kfilter: {self.r.kfilter}")
        #app(f"gstore_kptopt: {self.r.kptopt}")
        #app(f"gstore_qptopt: {self.r.qptopt}")

        for spin in range(self.nsppol):
            polaron = self.polaron_spin[spin]
            df = polaron.get_final_results_df()
            app("Last SCF iteration. Energies in eV units")
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
    nstates: int       # Number of polaronic states for this spin.
    nb: int            # Number of bands in A_kn.
    nk: int            # Number of k-points in A_kn (including filtering if any).
    nq: int            # Number of q-points in B_qnu (including filtering if any).
    bstart: int        # First band starts at bstart.
    bstop: int         # Last band (python convention)
    erange: double     # Filtering value (in Ha)
    varpeq: VpqFile

    @classmethod
    def from_vpq(cls, varpeq: VpqFile, spin: int) -> Polaron:
        """
        Build an istance from a VpqFile and the spin index.
        """
        r = varpeq.r
        nstates, nk, nq, nb = r.nstates, r.nk_spin[spin], r.nq_spin[spin], r.nb_spin[spin]
        bstart, bstop = r.brange_spin[spin]
        erange = r.erange_spin[spin]

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

        # Build list of dataframes.
        df_list = []
        for pstate in range(self.nstates):
            n = nstep2cv[pstate]
            dct = {k: scf_hist[pstate, :n, i] * ufact_k(k) for i, k in enumerate(_ALL_ENTRIES)}
            df = pd.DataFrame(dct)
            # Add metadata to the attrs dictionary
            df.attrs["converged"] = bool(cvflag[pstate])
            df.attrs["use_filter"] = bool(abs(self.erange) > 1e-8)
            df.attrs["filter_value"] = self.erange * abu.Ha_eV
            df_list.append(df)

        return df_list

    def get_final_results_df(self, with_params: bool = False) -> pd.DataFrame:
        """
        Return daframe with the last iteration for all polaronic states.
        NB: Energies are in eV.
        """
        row_list = []
        for pstate in range(self.nstates):
            df = self.scf_df_state[pstate]
            row = {"formula": self.structure.reduced_formula,
                   "spgroup": self.structure.get_space_group_info()[1],
                   "polaron": self.varpeq.r.vpq_pkind,
                   "pstate": pstate, "spin": self.spin}
            row.update(df.iloc[-1].to_dict())
            row["converged"] = df.attrs["converged"]
            row["use_filter"] = df.attrs["use_filter"]
            row["filter_value"] = df.attrs["filter_value"]
            if with_params:
                row.update(self.varpeq.params)
            row_list.append(row)

        return pd.DataFrame(row_list)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
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

    def get_title(self, with_gaps: bool = True) -> str:
        """
        Return string with title for matplotlib plots.
        """
        varpeq = self.varpeq
        pre = "" if varpeq.ebands.nsppol == 1 else f"spin={self.spin}"
        if not with_gaps:
            return f"{pre}{varpeq.r.vpq_pkind} polaron"
        else:
            gaps_string = varpeq.ebands.get_gaps_string()
            return f"{pre}{varpeq.r.vpq_pkind} polaron, {gaps_string}"

    def insert_a_inbox(self, fill_value=None) -> tuple:
        """
        Return a_data, ngkpt, shifts where a_data is a
        (nstates, nb, nkx, nky, nkz)) array with A_{pnk} with p the polaron index.
        """
        # Need to know the shape of the k-mesh.
        ngkpt, shifts = self.ngkpt_and_shifts
        k_indices = kpoints_indices(self.kpoints, ngkpt, shifts)
        #print(f"{k_indices=}")
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
        q_indices = kpoints_indices(self.qpoints, ngqpt, shifts)

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

    def get_a2_interpolator_state(self, interp_method: str) -> BzRegularGridInterpolator:
        """
        Build and return an interpolator for |A_nk|^2 for each polaronic state.

        Args:
            interp_method: The method of interpolation. Supported are “linear”, “nearest”,
                “slinear”, “cubic”, “quintic” and “pchip”.
        """
        a_data, ngkpt, shifts = self.insert_a_inbox(fill_value=0)

        return [BzRegularGridInterpolator(self.structure, shifts, np.abs(a_data[pstate])**2, method=interp_method)
                for pstate in range(self.nstates)]

    def get_b2_interpolator_state(self, interp_method: str) -> BzRegularGridInterpolator:
        """
        Build and return an interpolator for |B_qnu|^2 for each polaronic state.

        Args:
            interp_method: The method of interpolation. Supported are “linear”, “nearest”,
                “slinear”, “cubic”, “quintic” and “pchip”.
        """
        b_data, ngqpt, shifts = self.insert_b_inbox(fill_value=0)

        return [BzRegularGridInterpolator(self.structure, shifts, np.abs(b_data[pstate])**2, method=interp_method)
                for pstate in range(self.nstates)]

    def write_a2_bxsf(self, filepath: PathLike, fill_value: float = 0.0) -> None:
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

    def write_b2_bxsf(self, filepath: PathLike, fill_value: float = 0.0) -> None:
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

    @add_fig_kwargs
    def plot_scf_cycle(self, ax_mat=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the SCF cycle.

        Args:
            ax_max: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
        """
        # Build grid of plots.
        nrows, ncols = self.nstates, 3
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        for pstate in range(self.nstates):
            df = self.scf_df_state[pstate]
            niter = len(df)
            xs = np.arange(1, niter + 1)

            for iax, ax in enumerate(ax_mat[pstate]):
                # Create a twin Axes sharing the x-axis

                for ilab, (name, entry) in enumerate(_ALL_ENTRIES.items()):
                    # Convert energies to Hartree. Keep gradient as it is.
                    ys = df[name].to_numpy()

                    if entry.utype == "gradient":
                        # plot only the gradient residual on the 1st panel
                        if iax == 0:
                            energy_like = False
                            ax.plot(xs, ys, label=entry.latex, c='k')
                            ax.set_yscale("log")
                    else:
                        if entry.name == "E_pol":
                        # Solid line for the *variational* quantity, also put it on top
                            ls, zord = '-', 10
                        else:
                        # Dashed lines for non-variational, put them below
                            ls, zord = '--', 0

                        if iax == 1:
                            energy_like =True
                            # Plot values linear scale.
                            ax.plot(xs, ys, label=entry.latex, linestyle=ls, zorder=zord)
                        elif iax == 2:
                            energy_like = True
                            # Plot deltas in logscale.
                            # (remove the last point for pretty-plotting)
                            ax.plot(xs[:-1], np.abs(ys - ys[-1])[:-1], label=entry.latex,
                                    linestyle=ls, zorder=zord)
                            ax.set_yscale("log")

                ax.set_xlim(1, niter)

                if energy_like:
                    ylabel = "Energy (eV)" if iax == 1 else r"$|\Delta E|$ (eV)"
                else:
                    ylabel = r"$|\nabla|$" if iax == 0 else r"$|\Delta |\nabla|$"

                set_grid_legend(ax, fontsize, xlabel="Iteration") #, ylabel=ylabel)
                ax.set_ylabel(ylabel)
                ax.legend()

                if pstate == 0:
                    if iax == 0:
                        ax.set_title("Gradient norm")
                    elif iax == 1:
                        ax.set_title("Energy terms")
                    else:
                        ax.set_title("Log-scale difference")


        fig.suptitle(self.get_title(with_gaps=True))
        fig.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_ank_with_ebands(self, ebands_kpath,
                             ebands_kmesh=None, lpratio: int = 5, with_info = True, with_legend=True,
                             with_ibz_a2dos=True, method="gaussian", step="auto", width="auto",
                             nksmall: int = 20, normalize: bool = False, with_title=True, interp_method="linear",
                             ax_mat=None, ylims=None, scale=50, marker_color="gold", marker_edgecolor="gray",
                             marker_alpha=0.5, fontsize=12, lw_bands=1.0, lw_dos=1.0,
                             filter_value=None, fill_dos=True, **kwargs) -> Figure:
        """
        Plot electron bands with markers whose size is proportional to |A_nk|^2.

        Args:
            ebands_kpath: ElectronBands or netcdf file providing an electronic band structure along a k-path.
            ebands_kmesh: ElectronBands or netcdf file providing an electronic band structure with k-points in the IBZ.
            lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
                The default should be OK in many systems, larger values may be required for accurate derivatives.
            with_ibz_a2dos: True if A2_IBZ(E) should be computed.
            method: Integration scheme for electron DOS.
            step: Energy step (eV) of the linear mesh for DOS computation.
            width: Standard deviation (eV) of the gaussian for DOS computation.
            nksmall: Number of divisions to sample the smallest reciprocal lattice vector when computing electron DOS
            normalize: Rescale electron and A2(E) DOSes to plot them on the same scale.
            with_title: True to add title with chemical formula and gaps.
            interp_method: Interpolation method.
            ax_mat: Matrix |matplotlib-Axes| or None if a new figure should be created.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            scale: Scaling factor for |A_nk|^2.
            marker_color: Color for markers and ADOS.
            marker_edgecolor: Color for marker edges.
            marker_edgecolor: Marker transparency.
            lw_bands: Electronic bands linewidth.
            lw_dos: DOS linewidth.
            filter_value: Energy filter value (eV) wrt the band edge.
            fill_dos: True if the regions defined by ADOS (BZ and IBZ) should be colored.
            fontsize: fontsize for legends and titles
        """
        nrows, ncols = self.nstates, 2
        gridspec_kw = {'width_ratios': [2, 1]}
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=True, squeeze=False, gridspec_kw=gridspec_kw)
        # Get interpolators for |A_nk|^2
        a2_interp_state = self.get_a2_interpolator_state(interp_method)

        df = self.get_final_results_df()

        ebands_kpath = ElectronBands.as_ebands(ebands_kpath)
        ymin, ymax = +np.inf, -np.inf

        pkind = self.varpeq.r.vpq_pkind
        vbm_or_cbm = "vbm" if pkind == "hole" else "cbm"
        bm = self.ebands.get_edge_state(vbm_or_cbm, self.spin).eig
        e0 = self.ebands.fermie

        for pstate in range(self.nstates):
            x, y, s = [], [], []

            a2_max = a2_interp_state[pstate].get_max_abs_data()
            scale *= 1. / a2_max

            for ik, kpoint in enumerate(ebands_kpath.kpoints):
                enes_n = ebands_kpath.eigens[self.spin, ik, self.bstart:self.bstop]
                for e, a2 in zip(enes_n, a2_interp_state[pstate].eval_kpoint(kpoint), strict=True):
                    # Handle filtering
                    allowed = True
                    if filter_value:
                        energy_window = filter_value * 1.1
                        if pkind == "hole" and bm - e > energy_window:
                            allowed = False
                        elif pkind == "electron" and e - bm > energy_window:
                            allowed = False

                    if allowed:
                        x.append(ik); y.append(e); s.append(scale * a2)
                        ymin, ymax = min(ymin, e), max(ymax, e)

            # Plot electron bands with markers.
            ax = ax_mat[pstate, 0]

            points = Marker(x, y, s, color=marker_color, edgecolors=marker_edgecolor,
                            alpha=marker_alpha, label=r'$|A_{n\mathbf{k}}|^2$')

            ebands_kpath.plot(ax=ax, points=points, show=False, linewidth=lw_bands)

            ax.legend(loc="best", shadow=True, fontsize=fontsize)

            # Add energy and convergence status
            if with_info:
                data = (df[df["pstate"] == pstate]).to_dict(orient="list")
                e_pol_ev, converged = float(data["E_pol"][0]), bool(data["converged"][0])
                ax.set_title(f"Formation energy: {e_pol_ev:.3f} eV, {converged=}" , fontsize=8)

            if pstate != self.nstates - 1 or not with_legend:
                set_visible(ax, False, *["legend", "xlabel"])

        vertices_names = [(k.frac_coords, k.name) for k in ebands_kpath.kpoints]

        if ebands_kmesh is None:
            edos_ngkpt = self.structure.calc_ngkpt(nksmall)
            print(f"Computing ebands_kmesh with star-function interpolation and {nksmall=} --> {edos_ngkpt=} ...")
            r = self.ebands.interpolate(lpratio=lpratio, vertices_names=vertices_names, kmesh=edos_ngkpt)
            ebands_kmesh = r.ebands_kmesh
        else:
            ebands_kmesh = ElectronBands.as_ebands(ebands_kmesh)

        # Get electronic DOS from ebands_kmesh.
        # Maybe it's better to set N bandwidth divisions & width factor instead of
        # the step and width arguments themselves?
        bandwidth = ylims[1] - ylims[0] if ylims else 1.2*(ymax - ymin)
        if step == "auto":
            step = bandwidth / 200
        if width == "auto":
            width = 2.0 * step

        edos_kws = dict(method=method, step=step, width=width)
        edos = ebands_kmesh.get_edos(**edos_kws)
        edos_mesh = edos.spin_dos[self.spin].mesh
        e0 = self.ebands.fermie

        ##################
        # Compute A_nk DOS
        ##################
        # NB: A_nk does not necessarily have the symmetry of the lattice so we have to loop over the full BZ.
        # Here we get the mapping BZ --> IBZ needed to obtain the KS eigenvalues e_nk from the IBZ for the DOS.
        kdata = ebands_kmesh.get_bz2ibz_bz_points()

        xmax = -np.inf
        for pstate in range(self.nstates):

            # Compute A^2(E) DOS with A_nk in the full BZ.
            ank_dos = np.zeros(len(edos_mesh))
            for ik_ibz, bz_kpoint in zip(kdata.bz2ibz, kdata.bz_kpoints, strict=True):
                enes_n = ebands_kmesh.eigens[self.spin, ik_ibz, self.bstart:self.bstop]
                a2_n = a2_interp_state[pstate].eval_kpoint(bz_kpoint)
                for band, (e, a2) in enumerate(zip(enes_n, a2_n, strict=True)):
                    ank_dos += a2 * gaussian(edos_mesh, width, center=e-e0)

            ank_dos /= np.prod(kdata.ngkpt)
            ank_dos = Function1D(edos_mesh, ank_dos)
            print(f"For {pstate=}, A^2(E) integrates to:", ank_dos.integral_value, " Ideally, it should be 1.")

            ax = ax_mat[pstate, 1]
            edos_opts = {"color": "black",} if self.spin == 0 else {"color": "red"}
            lines_edos = edos.plot_ax(ax, e0, spin=self.spin, normalize=normalize, exchange_xy=True, label="eDOS(E)", **edos_opts,
                                      linewidth=lw_dos, zorder=3)

            lines_ados = ank_dos.plot_ax(ax, exchange_xy=True, normalize=normalize, label=r"$A^2$(E)", color=marker_color,
                                         linewidth=lw_dos, zorder=2)

            # Computes A2(E) using only k-points in the IBZ. This is just for testing.
            # A2_IBZ(E) should be equal to A2(E) only if A_nk fullfills the lattice symmetries. See notes above.

            if with_ibz_a2dos:
                ank_dos = np.zeros(len(edos_mesh))
                for ik_ibz, ibz_kpoint in enumerate(ebands_kmesh.kpoints):
                    #print("ibz_kpoint:", ibz_kpoint)
                    weight = ibz_kpoint.weight
                    enes_n = ebands_kmesh.eigens[self.spin, ik_ibz, self.bstart:self.bstop]
                    for e, a2 in zip(enes_n, a2_interp_state[pstate].eval_kpoint(ibz_kpoint), strict=True):
                        ank_dos += weight * a2 * gaussian(edos_mesh, width, center=e-e0)

                ank_dos = Function1D(edos_mesh, ank_dos)
                ibz_dos_opts = {"color": "darkred",}
                print(f"For {pstate=}, A2_IBZ(E) integrates to:", ank_dos.integral_value, " Ideally, it should be 1.")
                lines_ados_ibz = ank_dos.plot_ax(ax, exchange_xy=True, normalize=normalize, label=r"$A^2_{IBZ}$(E)", ls="--",
                                                 linewidth=lw_dos, **ibz_dos_opts, zorder=1)

            set_grid_legend(ax, fontsize, xlabel="Arb. unit")
            if pstate != self.nstates - 1 or not with_legend:
                set_visible(ax, False, *["legend", "xlabel"])

            dos_lines = [lines_edos, lines_ados]
            colors = [edos_opts["color"], marker_color]
            span = ymax - ymin

            if with_ibz_a2dos:
                dos_lines.append(lines_ados_ibz)
                colors.append(ibz_dos_opts["color"])

            # determine max x value for auto xlims
            for dos, c in zip(dos_lines, colors):
                for line in dos:
                    x_data, y_data = line.get_xdata(), line.get_ydata()
                    mask = (y_data > ymin-e0-span*0.1) & (y_data < ymax-e0+span*0.1)
                    xmax = max(np.max(x_data[mask]), xmax)

            # fill Ank dos in order
            if fill_dos:
                y_common = np.linspace(ymin-e0-span*0.1, ymax-e0+span*0.1, 100)
                xleft = np.zeros_like(y_common)
                # skip eDOS, fill only ADOS
                for dos, c in zip(dos_lines[1:], colors[1:]):
                    for line in dos:
                        x_data, y_data = line.get_xdata(), line.get_ydata()
                        interp_x = interp1d(y_data, x_data, kind='linear', fill_value='extrapolate')
                        xright = interp_x(y_common)

                        mask = (xright - xleft) > 0
                        y, x0, x1 = y_common[mask], xleft[mask], xright[mask]

                        ax.fill_betweenx(y, x0, x1, alpha=marker_alpha, color=c, linewidth=0)
                        xleft = xright

        # Auto xlims for DOS
        span = xmax
        xmax += 0.1 * span
        for ax in ax_mat[:,1]:
            ax.set_xlim(0, xmax)

        if ylims is None:
            # Automatic ylims.
            span = ymax - ymin
            ymin -= 0.1 * span
            ymax += 0.1 * span
            ylims = [ymin - e0, ymax - e0]

        for ax in ax_mat.ravel():
            set_axlims(ax, ylims, "y")

        # if filtering is used, show the filtering region
        for ax in ax_mat.ravel():
            xmin, xmax = ax.get_xlim()
            xrange = np.linspace(xmin,xmax,100)
            shifted_bm = bm - e0
            if filter_value:
                if pkind == "hole":
                    fill_from, fill_to = shifted_bm - filter_value, shifted_bm
                elif pkind == "electron":
                    fill_from, fill_to = shifted_bm, shifted_bm + filter_value

                ax.axhline(fill_from, c='k', zorder=0, lw=lw_dos)
                ax.axhline(fill_to, c='k', zorder=0, lw=lw_dos)
                ax.fill_between(xrange, ylims[0], fill_from,
                                color='lightgray', linewidth=0, alpha=0.5, zorder=0)
                ax.fill_between(xrange, fill_to, ylims[1],
                                color='lightgray', linewidth=0, alpha=0.5, zorder=0)

        if with_title:
            fig.suptitle(self.get_title(with_gaps=True))

        return fig

    @add_fig_kwargs
    def plot_bqnu_with_ddb(self, ddb, smearing_ev=0.001,
                           with_phdos=True, anaddb_kwargs=None, **kwargs) -> Figure:
        """
        High-level interface to plot phonon energies with markers whose size is proportional to |B_qnu|^2.
        Similar to plot_bqnu_with_phbands but this function receives in input a DdbFile or a
        path to a DDB file and automates the computation of the phonon bands by invoking anaddb.

        Args:
            ddb: DdbFile or path to file.
            with_phdos: True if phonon DOS should be computed and plotter.
            anaddb_kwargs: Optional arguments passed to anaddb.
            kwargs: Optional arguments passed to plot_bqnu_with_phbands.
        """
        ddb = DdbFile.as_ddb(ddb)
        anaddb_kwargs = {} if anaddb_kwargs is None else anaddb_kwargs
        anaddb_kwargs["dos_method"] = f"gaussian:{smearing_ev} eV"

        with ddb.anaget_phbst_and_phdos_files(**anaddb_kwargs) as g:
            phbst_file, phdos_file = g[0], g[1]
            phbands_qpath = phbst_file.phbands
            return self.plot_bqnu_with_phbands(phbands_qpath,
                                               phdos_file=phdos_file if with_phdos else None,
                                               ddb=ddb, **kwargs)

    @add_fig_kwargs
    def plot_bqnu_with_phbands(self, phbands_qpath, with_legend=True,
                               phdos_file=None, ddb=None, width=0.001, normalize: bool = True,
                               verbose=0, anaddb_kwargs=None, with_title=True, interp_method="linear",
                               ax_mat=None, scale=50, marker_color="gold", marker_edgecolor='gray',
                               marker_alpha=0.5, fontsize=12, lw_bands=1.0, lw_dos=1.0,
                               fill_dos=True, **kwargs) -> Figure:
        """
        Plot phonon energies with markers whose size is proportional to |B_qnu|^2.

        Args:
            phbands_qpath: PhononBands or nc file providing a phonon band structure.
            phdos_file:
            ddb: DdbFile or path to file.
            width: Standard deviation (eV) of the gaussian.
            normalize: Rescale the two DOS to plot them on the same scale.
            verbose:
            anaddb_kwargs: Optional arguments passed to anaddb.
            with_title: True to add title with chemical formula and gaps.
            interp_method: Interpolation method.
            ax_mat: List of |matplotlib-Axes| or None if a new figure should be created.
            scale: Scaling factor for |B_qnu|^2.
            marker_color: Color for markers.
            fontsize: fontsize for legends and titles.
        """
        with_phdos = phdos_file is not None and ddb is not None
        nrows, ncols, gridspec_kw = self.nstates, 1, None
        if with_phdos:
            ncols, gridspec_kw = 2, {'width_ratios': [2, 1]}

        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=True, squeeze=False, gridspec_kw=gridspec_kw)

        phbands_qpath = PhononBands.as_phbands(phbands_qpath)

        # Get interpolators for |B_qnu|^2
        b2_interp_state = self.get_b2_interpolator_state(interp_method)

        # TODO: need to fix this hardcoded representation
        units = 'meV'
        units_scale = 1e3 if units == 'meV' else 1

        # Plot phonon bands with markers.
        ymin, ymax = +np.inf, -np.inf
        for pstate in range(self.nstates):
            x, y, s = [], [], []

            b2_max = b2_interp_state[pstate].get_max_abs_data()
            scale *= 1. / b2_max

            for iq, qpoint in enumerate(phbands_qpath.qpoints):
                omegas_nu = phbands_qpath.phfreqs[iq,:]

                for w, b2 in zip(omegas_nu, b2_interp_state[pstate].eval_kpoint(qpoint), strict=True):
                    w *= units_scale
                    x.append(iq); y.append(w); s.append(scale * b2)
                    ymin, ymax = min(ymin, w), max(ymax, w)

            ax = ax_mat[pstate, 0]
            points = Marker(x, y, s, color=marker_color, edgecolors=marker_edgecolor,
                            alpha=marker_alpha, label=r'$|B_{\nu\mathbf{q}}|^2$')
            phbands_qpath.plot(ax=ax, points=points, show=False, linewidth=lw_bands, units=units)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

            if pstate != self.nstates - 1 or not with_legend:
                set_visible(ax, False, *["legend", "xlabel"])

        # determine bandwidth and set ylims. if no negative freqs, set ymin exactly to 0
        if ymin > -1e-6:
            ymin = 0
        bandwidth = ymax - ymin
        ymin -= 0.1 * bandwidth if ymin != 0 else 0
        ymax += 0.1 * bandwidth

        for ax in ax_mat.ravel():
            ax.set_ylim(ymin, ymax)

        if not with_phdos:
            # Return immediately.
            if with_title:
                fig.suptitle(self.get_title(with_gaps=True))
            return fig

        ####################
        # Compute B_qnu DOS
        ####################
        # NB: B_qnu do not necessarily have the symmetry of the lattice so we have to loop over the full BZ.
        # The frequency mesh is in eV, values are in states/eV.
        # (note the units_scale variable before the phbands calculation)
        # Use same q-mesh as phdos
        phdos = phdos_file.phdos
        phdos_ngqpt = np.diagonal(phdos_file.qptrlatt)
        phdos_shifts = [0.0, 0.0, 0.0]
        phdos_nqbz = np.prod(phdos_ngqpt)
        phdos_mesh = phdos.mesh

        with_ibz_b2dos = False

        # Call anaddb to get phonons on the FULL ngqpt mesh.
        # The B_qnu do not necessarily have the symmetry of the lattice so we have to loop over the full BZ.
        anaddb_kwargs = {} if anaddb_kwargs is None else anaddb_kwargs

        bz_qpoints = kmesh_from_mpdivs(phdos_ngqpt, phdos_shifts)
        phbands_bz = ddb.anaget_phmodes_at_qpoints(qpoints=bz_qpoints, ifcflag=1, verbose=verbose, **anaddb_kwargs)
        if len(phbands_bz.qpoints) != np.prod(phdos_ngqpt):
            raise RuntimeError(f"{len(phbands_bz.qpoints)=} != {np.prod(phdos_ngqpt)=}")

        # TODO: Use this new approach so that we can reduce everything to the IBZ:
        #with KmeshFile.from_ngkpt_shifts(structure, phdos_ngqpt, phdos_shifts, kptopt=3, chksymbreak=0) as qmesh:
        #    qibz = qmesh.ibz
        #    qbz2ibz = qmesh.bz2ibz
        #phbands_ibz = ddb.anaget_phmodes_at_qpoints(qpoints=qibz, ifcflag=1, verbose=verbose, **anaddb_kwargs)

        xmax = -np.inf

        for pstate in range(self.nstates):
            # Compute B2(E) by looping over the full BZ.
            bqnu_dos = np.zeros(len(phdos_mesh))
            for iq_bz, qpoint in enumerate(phbands_bz.qpoints):
                #q_weight = 1./phdos_nqbz
                freqs_nu = phbands_bz.phfreqs[iq_bz]
                for w, b2 in zip(freqs_nu, b2_interp_state[pstate].eval_kpoint(qpoint), strict=True):
                    bqnu_dos += b2 * gaussian(phdos_mesh, width, center=w)

            # NB: all the q-weights in PHBST.nc are set to 1.
            bqnu_dos /= phdos_nqbz
            bqnu_dos = Function1D(phdos_mesh, bqnu_dos)

            ax = ax_mat[pstate, 1]
            pdos_opts = {"color": "black"}
            lines_pdos = phdos.plot_dos_idos(ax, exchange_xy=True, units=units, label="phDOS(E)",
                                             normalize=normalize, linewidth=lw_dos, **pdos_opts)
            #phdos.plot_ax(ax, exchange_xy=True, normalize=normalize, label="phDOS(E)", color="black",
            #              linewidth=lw_dos, units=units)
            lines_bdos = bqnu_dos.plot_ax(ax, exchange_xy=True, normalize=normalize,
                                          label=r"$B^2$(E)", color=marker_color, linewidth=lw_dos,
                                          xfactor=units_scale, yfactor=1/units_scale)
            set_grid_legend(ax, fontsize, xlabel="Arb. unit")

            # Get mapping BZ --> IBZ needed to obtain the KS eigenvalues e_nk from the IBZ for the DOS
            # Compute B2(E) using only q-points in the IBZ. This is just for testing.
            # B2_IBZ(E) should be equal to B2(E) only if B_qnu fullfill the lattice symmetries. See notes above.
            ibz_dos_opts = {"color": "darkred"}
            lines_bdos_ibz = None

            dos_lines = [lines_pdos, lines_bdos]
            colors = [pdos_opts["color"], marker_color]
            span = ymax - ymin

            if with_ibz_b2dos:
                dos_lines.append(lines_bdos_ibz)
                colors.append(ibz_dos_opts["color"])

            # determine max x value for auto xlims
            for dos, c in zip(dos_lines, colors):
                for line in dos:
                    x_data, y_data = line.get_xdata(), line.get_ydata()
                    mask = (y_data > ymin) & (y_data < ymax+span*0.1)
                    xmax = max(np.max(x_data[mask]), xmax)

            # fill Bqnu dos in order
            if fill_dos:
                y_common = np.linspace(ymin, ymax+span*0.1, 100)
                xleft = np.zeros_like(y_common)
                # skip phDOS, fill only BDOS
                for dos, c in zip(dos_lines[1:], colors[1:]):
                    for line in dos:
                        x_data, y_data = line.get_xdata(), line.get_ydata()
                        interp_x = interp1d(y_data, x_data, kind='linear', fill_value='extrapolate')
                        xright = interp_x(y_common)

                        mask = (xright - xleft) > 0
                        y, x0, x1 = y_common[mask], xleft[mask], xright[mask]

                        ax.fill_betweenx(y, x0, x1, alpha=marker_alpha, color=c, linewidth=0)
                        xleft = xright

        # Auto xlims for DOS
        span = xmax
        xmax += 0.1 * span
        for ax in ax_mat[:,1]:
            ax.set_xlim(0, xmax)

            if pstate != self.nstates - 1 or not with_legend:
                set_visible(ax, False, *["legend", "xlabel"])

        if with_title:
            fig.suptitle(self.get_title(with_gaps=True))

        return fig


class VpqReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VpqReader
    """

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        # Netcdf Variables

        # int eph_task ;
        # int nkbz ;
        # int nqbz ;
        # int frohl_ntheta ;
        # double vpq_tolgrs ;
        # double e_frohl ;
        # char vpq_pkind(fnlen) ;
        # char vpq_aseed(fnlen) ;
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
        self.vpq_pkind = self.read_string("vpq_pkind")
        #self.vpq_aseed = self.read_string("vpq_aseed")
        self.ngqpt = self.read_value("gstore_ngqpt")
        #self.frohl_ntheta = self.read_value("frohl_ntheta")

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
        self.erange_spin = self.read_value("erange_spin")
        # Total number of k/q points for each spin after filtering (if any)
        #self.glob_spin_nq = self.read_value("gstore_glob_nq_spin")
        #self.glob_nk_spin = self.read_value("gstore_glob_nk_spin")


class VpqRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple VPQ.nc files.

    Usage example:

    .. code-block:: python

        robot = VpqRobot.from_files([
            "out1_VPQ.nc",
            "out2_VPQ.nc",
        ])

        print(robot)
        df = robot.get_final_results_df()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VpqRobot
    """

    EXT = "VPQ"

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append
        df = self.get_final_results_df()
        lines.append(str(df))

        return "\n".join(lines)

    def get_final_results_df(self, spin: int = None, sortby: str = None, with_params: bool = True) -> pd.DataFrame:
        """
        Return dataframe with the last iteration for all polaronic states.
        NB: Energies are in eV.

        Args:
            spin: Spin index, None if all spins should be included.
            sortby: Name to sort by.
            with_params: True if columns with convergence parameters should be added.
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
        if sortby and sortby in df:
            df = df.sort_values(sortby)

        return df

    @add_fig_kwargs
    def plot_erange_conv(self, spin: int = 0, pstate: int = 0, **kwargs) -> list[Figure]:
        """
        Args:
            spin (int, optional): Spin index. Defaults to 0.
            pstate (int, optional): Index of a polaronic state. Defaults to 0.
            **kwargs: Additional arguemtns for plottins functions.

        Returns:
            list: A list of figure objects.
        """

        df = self.get_final_results_df(spin)

        # check if dataframe contains entries with efilter
        df = df[df["use_filter"] & (df["pstate"] == pstate)]
        if df.empty:
            raise RuntimeError("No entries with energy filtering.")

        fig_list = []
        grouped_entries = df.groupby(["formula", "spgroup", "polaron", "ngkpt"])

        # here we iterate over each polaronic group & generate separate figures
        for (formula, spg, pol, ngkpt), group in grouped_entries:

            # check if we have enough filtering values for convergence
            if group["filter_value"].nunique() == 1:
                continue

            group = group.sort_values("filter_value")
            nrows, ncols = 2, 2
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                   sharex=True, sharey=False, squeeze=False)

            for avg_g in [True, False]:
                _df = group[group["avg_g"] == avg_g]
                if _df.empty:
                    continue

                filter_values = _df["filter_value"].to_numpy()
                epol = _df["E_pol"].to_numpy()
                eps = _df["epsilon"].to_numpy()
                frohlich_label = " + LR" if avg_g else ""

                # Convergence
                ax_mat[0,0].plot(filter_values, epol, 's-', label=r'$E_{pol}$' + frohlich_label,
                                 **kwargs)
                ax_mat[1,0].plot(filter_values, eps, 's-', label=r"$\varepsilon$" + frohlich_label,
                                 **kwargs)

                # Relative error
                ax_mat[0,1].plot(filter_values[:-1], np.abs((epol - epol[-1])/epol[-1])[:-1], 's-',
                                   label=r'$E_{pol}$' + frohlich_label, **kwargs)
                ax_mat[1,1].plot(filter_values[:-1], np.abs((eps - eps[-1])/eps[-1])[:-1], 's-',
                                   label=r"$\varepsilon$" + frohlich_label, **kwargs)

            for i, ax_row in enumerate(ax_mat):
                ax_row[1].set_yscale("log")
                ax_row[1].set_ylabel("Relative error (-)")
                ax_row[0].set_ylabel("Energy (eV)")

            for j in range(ncols):
                ax_mat[1, j].set_xlabel("Filter value (eV)")
                ax_mat[0, j].set_title("Binding energy")
                ax_mat[1, j].set_title("Polaron eigenvalue")

            for ax in np.ravel(ax_mat):
                ax.set_xlim(0)
                ax.grid()
                ax.legend()

            fig.suptitle(f"{formula}, space group {spg}, {pol} polaron, k-mesh: {'x'.join(map(str, ngkpt))}")
            fig.tight_layout()

            fig_list.append(fig)

        if not fig_list:
            print("plot_erange_conv: not enough data fro convergence wrt erange")

        return fig_list


    @add_fig_kwargs
    def plot_kconv(self, nfit: int = 3, spin: int = 0, pstate: int = 0, convby: str = "invsc_linsize",
                   add_lr: bool = False, **kwargs) -> list[Figure]:
        """
        Plot the convergence of the results with respect to k-point sampling.

        Args:
            nfit (int, optional): Number of points used in linear extrapolation. Defaults to 3.
            spin (int, optional): Spin index. Defaults to 0.
            pstate (int, optional): Index of a polaronic state. Defaults to 0.
            convby (str, optional): Convergence parameter.
                Possible values are:
                    - `"invsc_linsize"`: Inverse linear size of a supercell (inverse Angstrom).
                    - `"inv_k"`: Inverse linear size of a k-grid (arbitrary units).
                Defaults to `"invsc_linsize"`.
            add_lr (bool, optional): Specifies if LR correction should be added a-posteriori.
                Relevant when LR correction to the matrix elements is not used. Defaults to False.
            **kwargs: Additional arguments for plotting functions.

        Returns:
            list[Figure]: A list of figure objects.
        """

        if convby not in {"invsc_linsize", "inv_k"}:
            raise ValueError(f"Invalid convby value '{convby}'. Choose 'invsc_linsize' or 'inv_k'.")

        df = self.get_final_results_df(spin)
        fig_list = []
        grouped_entries = df.groupby(["formula", "spgroup", "polaron"])

        # Iterate over each polaronic group and generate figures
        for (formula, spg, pol), group in grouped_entries:
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=2, ncols=1, sharex=True, sharey=False, squeeze=False)
            sub_entries = group.groupby(["filter_value", "avg_g"])

            for (filter_value, avg_g), subgroup in sub_entries:
                if subgroup.empty:
                    continue

                # Sort values for consistent plotting
                subgroup = subgroup.sort_values("invsc_linsize")
                params = subgroup[convby].to_numpy()
                e_pol = subgroup["E_pol"].to_numpy()
                eps = subgroup["epsilon"].to_numpy()
                frohlich_label = "+ LR" if avg_g else ""
                filter_label = f"filter {filter_value:.2f} eV," if filter_value > 0 else ""

                # Apply LR correction if needed
                if not avg_g and add_lr:
                    eps_sign = -1 if pol == "hole" else 1
                    e_frohl = subgroup["E_pol"].to_numpy()
                    e_pol += e_frohl
                    eps += e_frohl * eps_sign

                # Compute extrapolation lines
                local_nfit = min(nfit, len(params))
                epol_extr_line = np.poly1d(np.polyfit(params[:local_nfit], e_pol[:local_nfit], 1))
                eps_extr_line = np.poly1d(np.polyfit(params[:local_nfit], eps[:local_nfit], 1))
                xrange = np.linspace(0, np.max(params[:local_nfit]))

                # Plot energy data & extrapolation
                line1, = ax_mat[0, 0].plot(params, e_pol, 'o', **kwargs)
                ax_mat[0, 0].plot(xrange, epol_extr_line(xrange), '--', color=line1.get_color(),
                                  label=rf'{filter_label} $E_{{pol}}$ {frohlich_label}: {epol_extr_line(0):.3f} eV')

                line2, = ax_mat[1, 0].plot(params, eps, 'o', **kwargs)
                ax_mat[1, 0].plot(xrange, eps_extr_line(xrange), '--', color=line2.get_color(),
                                  label=rf'{filter_label} $\varepsilon$ {frohlich_label}: {eps_extr_line(0):.3f} eV')

            # Set axis labels and formatting
            xlabel_map = {
                "invsc_linsize": r"$V_\mathrm{supercell}^{-1/3}$ ($\AA^{-1}$)",
                "inv_k": r"$N_p^{-1/3}$ (-)"
            }
            ax_mat[1, 0].set_xlabel(xlabel_map[convby])

            ax_mat[0, 0].set_title("Binding Energy")
            ax_mat[1, 0].set_title("Polaron Eigenvalue")

            for ax in ax_mat.ravel():
                ax.set_ylabel("Energy (eV)")
                ax.axhline(0, color='k', linestyle='-')
                ax.axvline(0, color='k', linestyle='-')
                ax.grid()
                ax.legend(title="Extrapolation")

            # Set figure title and layout
            fig.suptitle(f"{formula}, space group {spg}, {pol} polaron")
            fig.tight_layout()
            fig_list.append(fig)

        return fig_list

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.plot_scf_cycle(show=False)
        #yield self.plot_kconv()

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            nbv.new_code_cell("robot = abilab.VpqRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

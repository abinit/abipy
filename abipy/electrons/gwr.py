# coding: utf-8
"""Classes to analyze GWR calculations."""
from __future__ import annotations

import numpy as np
import abipy.core.abinit_units as abu

from monty.functools import lazy_property
from monty.string import list_strings, is_string, marquee
#from abipy.core.kpoints import Kpoint, KpointList, Kpath, IrredZone, has_timrev_from_kptopt
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.iotools import ETSF_Reader
#from abipy.tools import duck
from abipy.tools.plotting import (ArrayPlotter, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, Marker,
    set_axlims, set_visible, rotate_ticklabels)

__all__ = [
    "GwrFile",
    "TchimVsSus",
]


class GwrFile(AbinitNcFile, Has_Structure): #, Has_ElectronBands, NotebookWriter):

    @classmethod
    def from_file(cls, filepath: str):
        """Initialize an instance from file."""
        return cls(filepath)

    def __init__(self, filepath):
        """Read data from the netcdf file path."""
        super().__init__(filepath)

        # Keep a reference to the GwrReader.
        self.reader = reader = GwrReader(self.filepath)

    @property
    def structure(self):
        """|Structure| object."""
        return self.reader.structure

    #@property
    #def ebands(self):
    #    """|ElectronBands| with the KS energies."""
    #    return self._ebands

    @lazy_property
    def iw_mesh_ev(self):
        return self.reader.read_value("iw_mesh") * abu.Ha_eV

    @lazy_property
    def params(self) -> dict:
        """:class:`OrderedDict` with parameters that might be subject to convergence studies e.g ecuteps"""
        return {}

    def close(self) -> None:
        """Close the netcdf file."""
        self.reader.close()

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
        #app(self.ebands.to_string(title="Kohn-Sham bands", with_structure=False))
        #app("")

        self.plot_sigma_imag_axis()

        # TODO: Finalize the implementation: add GW metadata.
        #app(marquee("QP direct gaps", mark="="))
        #for kgw in self.gwkpoints:
        #    for spin in range(self.nsppol):
        #        qp_dirgap = self.get_qpgap(spin, kgw)
        #        app("QP_dirgap: %.3f (eV) for K-point: %s, spin: %s" % (qp_dirgap, repr(kgw), spin))
        #        #ks_dirgap =
        #app("")

        # Show QP results
        #strio = StringIO()
        #self.print_qps(precision=3, ignore_imag=verbose == 0, file=strio)
        #strio.seek(0)
        #app("")
        #app(marquee("QP results for each k-point and spin (all in eV)", mark="="))
        #app("".join(strio))
        #app("")

        # TODO: Fix header.
        #if verbose > 1:
        #    app("")
        #    app(marquee("Abinit Header", mark="="))
        #    app(self.hdr.to_string(verbose=verbose))

        return "\n".join(lines)

    @add_fig_kwargs
    def plot_sigma_imag_axis(self, spin=0, fontsize=8, **kwargs):
        """
        Plot ...

        Args:
            include_bands: List of bands to include. None means all.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """

        #"""
        nrows, ncols = 1, 2
        ax_list = None
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()
        re_ax, im_ax = ax_list

        #ik_gw = self.reader.gwkpt2seqindex(kpoint)
        #ik_ibz = self.reader.kpt2fileindex(kpoint)
        #for band in range(self.gwbstart_sk[spin, ik_gw], self.gwbstop_sk[spin, ik_gw]):
        #    ib_gw = band - self.min_gwbstart
        ikcalc = 0

        # nctkarr_t("sigc_iw_diag_kcalc", "dp", "two, ntau, max_nbcalc, nkcalc, nsppol")
        sigc_iw = self.reader.read_value("sigc_iw_diag_kcalc", cmode="c") * abu.Ha_eV

        # nctkarr_t("sigc_it_diag_kcalc", "dp", "two, two, ntau, max_nbcalc, nkcalc, nsppol")
        #sigc_tau = self.reader.read_value("sigc_it_diag_kcalc", cmode="c") * abu.Ha_eV

        #for band in range(0, 7):
        for band in range(0, 2):
            #ib_gw = band - self.min_gwbstart
            ib_gw = band
            fact = 1
            re_ax.plot(self.iw_mesh_ev, sigc_iw[spin, ikcalc, ib_gw].real * fact, label=f"Re band {band}")
            im_ax.plot(self.iw_mesh_ev, sigc_iw[spin, ikcalc, ib_gw].imag * fact, label=f"Im band {band}")

        re_ax.set_ylabel(r"$\Re{\Sigma_c}$ (eV)")
        im_ax.set_ylabel(r"$\Im{\Sigma_c}$ (eV)")
        for ax in ax_list:
            ax.grid(True)
            ax.set_xlabel(r"$i\omega$ (eV)")
            ax.legend(loc="best", fontsize=fontsize, shadow=True)
            set_axlims(ax, [0, 100], "x")

        return fig
        #"""

        nrows, ncols = 1, 2
        ax_list = None
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()
        re_ax, im_ax = ax_list

        # nctkarr_t("e0_kcalc", "dp", "max_nbcalc, nkcalc, nsppol"), &
        # nctkarr_t("sigc_wr_diag_kcalc", "dp", "two, nwr, max_nbcalc, nkcalc, nsppol")
        e0_kcalc = self.reader.read_value("e0_kcalc") * abu.Ha_eV
        sigc_wr = self.reader.read_value("sigc_wr_diag_kcalc", cmode="c") * abu.Ha_eV
        wr_step = self.reader.read_value("wr_step") * abu.Ha_eV
        nwr = self.reader.read_dimvalue("nwr")

        ikcalc = 0
        #ib = 4
        # Exchange at Gamma
        #sigx = np.array([-17.450, -13.026, -13.026, -13.026, -5.665, -5.665, -5.665])
        sigx = np.array([-17.023, -12.6, -12.6, -12.6, -5.665, -5.665, -5.665])
        sigx[:] = 0
        #vxc = np.array([-10.465, -11.330, -11.330, -11.330, -10.030, -10.030, -10.030])
        #sigx -= vxc

        #for ib in [0, 1, 4]:
        for ib in [0, 1]:
            # Mesh is linear and **centered** around e0.
            e0 = e0_kcalc[spin, ikcalc, ib]
            wr_mesh = np.linspace(start=e0 - wr_step * (nwr // 2),
                                  stop=e0 + wr_step * (nwr // 2),  num=nwr)
            # Add sigx manually
            l = re_ax.plot(wr_mesh, sigc_wr[spin, ikcalc, ib].real + sigx[ib], label=f"Re band: {ib}")
            color = l[0].get_color()
            re_ax.axvline(e0, lw=1, color=color, ls="--")
            l = im_ax.plot(wr_mesh, sigc_wr[spin, ikcalc, ib].imag, label=f"Im band: {ib}")
            color = l[0].get_color()
            im_ax.axvline(e0, lw=1, color=color, ls="--")

        re_ax.set_ylabel(r"$\Re{\Sigma}$ (eV)")
        im_ax.set_ylabel(r"$\Im{\Sigma}$ (eV)")
        for ax in ax_list:
            ax.grid(True)
            ax.set_xlabel("$\omega$ (eV)")
            set_axlims(ax, [-60, 60], "x")
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

        if include_bands is not None:
            include_bands = set(include_bands)

        # Build grid of plots.
        nrows, ncols = len(self.sigma_kpoints), 1
        ax_list = None
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        for ik_gw, (kgw, ax) in enumerate(zip(self.sigma_kpoints, ax_list)):
            for spin in range(self.nsppol):
                for band in range(self.gwbstart_sk[spin, ik_gw], self.gwbstop_sk[spin, ik_gw]):
                    if include_bands and band not in include_bands: continue
                    sigw = self.read_sigee_skb(spin, kgw, band)
                    label = r"$A(\omega)$: band: %d, spin: %d" % (spin, band)
                    sigw.plot_ax(ax, what="a", label=label, fontsize=fontsize, **kwargs)

            ax.set_title("K-point: %s" % repr(sigw.kpoint), fontsize=fontsize)

        return fig


class GwrReader(ETSF_Reader):
    r"""
    This object provides method to read data from the GWR file produced ABINIT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GwrReader
    """

    def __init__(self, path: str):
        super().__init__(path)
        # Save important quantities needed to simplify the API.
        #self.ks_bands = ElectronBands.from_file(path)
        #self.nsppol = self.ks_bands.nsppol

        self.structure = self.read_structure()

        """
        # 1) The K-points of the homogeneous mesh.
        self.ibz = self.ks_bands.kpoints

        # 2) The K-points where QPState corrections have been calculated.
        gwred_coords = self.read_value("kptgw")
        self.gwkpoints = KpointList(self.structure.reciprocal_lattice, gwred_coords)
        # Find k-point name
        for kpoint in self.gwkpoints:
            kpoint.set_name(self.structure.findname_in_hsym_stars(kpoint))

        # minbnd[nkptgw, nsppol] gives the minimum band index computed
        # Note conversion between Fortran and python convention.
        self.gwbstart_sk = self.read_value("minbnd") - 1
        self.gwbstop_sk = self.read_value("maxbnd")
        # min and Max band index for GW corrections.
        self.min_gwbstart = np.min(self.gwbstart_sk)
        self.max_gwbstart = np.max(self.gwbstart_sk)

        self.min_gwbstop = np.min(self.gwbstop_sk)
        self.max_gwbstop = np.max(self.gwbstop_sk)
        """


class TchimVsSus:

    def __init__(self, tchim_filepath, sus_filepath):
        from abipy.electrons.scr import SusFile
        self.sus_file = SusFile(sus_filepath)
        #print(self.sus_file)
        self.tchi_reader = ETSF_Reader(tchim_filepath)

    @add_fig_kwargs
    def plot_qpoint_gpairs(self, qpoint, gpairs, fontsize=7, spins=(0, 0), **kwargs):

        gpairs = np.array(gpairs)

        sus_reader = self.sus_file.reader
        _, sus_iq = sus_reader.find_kpoint_fileindex(qpoint)

        # int reduced_coordinates_plane_waves_dielectric_function(number_of_qpoints_dielectric_function,
        # number_of_coefficients_dielectric_function, number_of_reduced_dimensions) ;

        #reduced_coordinates_plane_waves_dielectric_function
        susc_iwmesh = sus_reader.read_value("frequencies_dielectric_function", cmode="c").imag
        # SUSC file uses Gamma-centered G-spere so read at iq = 0
        susc_gvec = sus_reader.read_variable("reduced_coordinates_plane_waves_dielectric_function")[0]

        tchi_qpoints = self.tchi_reader.read_variable("qibz")
        for tchi_iq, tchi_qq in enumerate(tchi_qpoints):
            if np.all(abs(tchi_qq - qpoint) < 1e-6):
                break
        else:
            raise ValueError(f"Cannot find: {qpoint} in TCHIM file")

        tchi_gvec = self.tchi_reader.read_variable("gvecs")[tchi_iq]
        tchi_iwmesh = self.tchi_reader.read_variable("iw_mesh")

        def _find_g(gg, gvec):
            for ig, g_sus in enumerate(gvec):
                if all(gg == g_sus):
                    return ig
            else:
                raise ValueError(f"Cannot find: {gg}")

        sus_inds = [(_find_g(g1, susc_gvec), _find_g(g2, susc_gvec)) for g1, g2 in gpairs]
        chi_inds = [(_find_g(g1, tchi_gvec), _find_g(g2, tchi_gvec)) for g1, g2 in gpairs]

        from abipy.tools.plotting import get_ax_fig_plt, get_axarray_fig_plt

        # Build grid of plots.
        num_plots, ncols, nrows = 2 * len(gpairs), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)

        for i, ((g1, g2), (sus_ig1, sus_ig2), (chi_ig1, chi_ig2)) in enumerate(zip(gpairs, sus_inds, chi_inds)):
            #print(g1, g2)
            ax_re, ax_im = ax_mat[i]

            # number_of_qpoints_dielectric_function, number_of_frequencies_dielectric_function,
            # number_of_spins, number_of_spins, number_of_coefficients_dielectric_function,
            # number_of_coefficients_dielectric_function, complex)

            if sus_reader.path.endswith("SUS.nc"):
               sus_var_name = "polarizability"
            elif sus_reader.path.endswith("SCR.nc"):
                sus_var_name = "inverse_dielectric_function"
            else:
                raise ValueError(f"Cannot detect varname for file: {sus_reader.path}")

            sus_data = sus_reader.read_variable(sus_var_name)[sus_iq, :, 0, 0, sus_ig2, sus_ig1]
            sus_data = sus_data[:, 0] + 1j * sus_data[:, 1]

            # nctkarr_t("mats", "dp", "two, mpw, mpw, ntau, nqibz, nsppol")
            tchi_data = self.tchi_reader.read_variable("mats")[0, tchi_iq, :, chi_ig2, chi_ig1, :]
            tchi_data = tchi_data[:, 0] + 1j * tchi_data[:, 1]

            opts = dict(markersize=4)

            ax_re.plot(susc_iwmesh[1:], sus_data[1:].real, ls="--", marker="o", label="AW Re", **opts)
            ax_re.plot(tchi_iwmesh, tchi_data[0:].real, ls=":", marker="^", label="minimax Re", **opts)
            ax_re.grid(True)

            ax_im.plot(susc_iwmesh[1:], sus_data[1:].imag, ls="--", marker="o", label="AW Im", **opts)
            ax_im.plot(tchi_iwmesh, tchi_data[0:].imag, ls=":", marker="^", label="minimax Im", **opts)
            ax_im.grid(True)

            tex_g1 = r"${\bf{G}}_1$"
            tex_g2 = r"${\bf{G}}_2$"
            tex_qpt = r"${\bf{q}}$"

            ax_re.set_title(f"{tex_g1}: {g1}, {tex_g2}: {g2}, {tex_qpt}: {qpoint}", fontsize=fontsize)
            ax_im.set_title(f"{tex_g1}: {g1}, {tex_g2}: {g2}, {tex_qpt}: {qpoint}", fontsize=fontsize)

            ax_re.legend(loc="best", fontsize=fontsize, shadow=True)
            ax_im.legend(loc="best", fontsize=fontsize, shadow=True)

            #if i == 0:
            if True:
                ax_re.set_ylabel(r'$\Re{\chi^0_{\bf{G}_1 \bf{G}_2}(i\omega)}$')
                ax_im.set_ylabel(r'$\Im{\chi^0_{\bf{G}_1 \bf{G}_2}(i\omega)}$')
                #ax_im.yaxis.set_label_position("right")
                #ax_im.yaxis.tick_right()

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


        with_title = False
        if with_title:
            tex1 = r"$\chi^0(i\omega)$"
            tex2 = r"$\chi^0(i\tau) \rightarrow\chi^0(i\omega)$"
            fig.suptitle(f"Comparison between Adler-Wiser {tex1} and minimax {tex2}\nLeft: real part. Right: imag part",
                         fontsize=10)

        return fig

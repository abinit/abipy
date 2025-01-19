# coding: utf-8
"""
Objects to analyze the TCHIM file produced by the GWR code.
"""
from __future__ import annotations

import numpy as np
#import pandas as pd
import abipy.core.abinit_units as abu

from typing import Any
#from monty.functools import lazy_property
#from monty.string import list_strings, marquee
#from monty.termcolor import cprint
from abipy.core.structure import Structure
#from abipy.core.kpoints import Kpoint, KpointList
#from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.iotools import ETSF_Reader
from abipy.tools import duck
from abipy.tools.typing import Figure, KptSelect
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, Marker,
    set_axlims, set_ax_xylabels, set_visible, rotate_ticklabels, set_grid_legend, hspan_ax_line, Exposer)
from abipy.tools.numtools import data_from_cplx_mode


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

        raise ValueError(f"Cannot find {qpoint=} in TCHIM file")

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
    def plot_matdiff(self,
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

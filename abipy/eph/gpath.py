"""
This module contains objects for analyzing
the PATH.nc file with the e-ph matrix elements along a k/q path
"""
from __future__ import annotations

import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from monty.string import marquee
from monty.functools import lazy_property
#from monty.termcolor import cprint
from abipy.core.structure import Structure
from abipy.core.kpoints import Kpath
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.tools.typing import PathLike
#from abipy.tools.numtools import nparr_to_df
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles, Marker, set_grid_legend, set_axlims)
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.dfpt.phonons import PhononBands
from abipy.dfpt.phtk import NonAnalyticalPh
from abipy.tools.typing import Figure
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader


def k2s(k_vector, fmt=".3f", threshold = 1e-8) -> str:
    k_vector = np.asarray(k_vector)
    k_vector[np.abs(k_vector) < threshold] = 0

    return "[" + ", ".join(f"{x:.3f}" for x in k_vector) + "]"


class GpathFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    This file stores the e-ph matrix elements along a k/q path
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with GpathFile("out_GPATH.nc") as gpath:
            print(gpath)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GpathFile
    """

    @classmethod
    def from_file(cls, filepath: PathLike) -> GpathFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = GpathReader(filepath)

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.r.structure

    @lazy_property
    def ebands_k(self) -> ElectronBands:
        """Electron bands along the k path."""
        return self.r.read_ebands_which_fixed("q")

    @lazy_property
    def ebands_kq(self) -> ElectronBands:
        """Electron bands along the k+q path as a function of q."""
        return self.r.read_ebands_which_fixed("k")

    @lazy_property
    def phbands(self) -> PhononBands:
        """Phonon bands along the q-path (nq_path points)."""
        return self.r.read_phbands()

    def close(self) -> None:
        """Close the file."""
        self.r.close()

    @lazy_property
    def params(self) -> dict:
        """dict with the convergence parameters, e.g. ``nbsum``."""
        #od = OrderedDict([
        #    ("nbsum", self.nbsum),
        #    ("nqibz", self.r.nqibz),
        #])
        ## Add EPH parameters.
        #od.update(self.r.common_eph_params)

        od = {}
        return od

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        app("")
        if self.r.eph_fix_korq == "k":
            app(self.ebands_kq.to_string(with_structure=False, verbose=verbose, title="Electronic Bands (kq)"))
        if self.r.eph_fix_korq == "q":
            app(self.ebands_k.to_string(with_structure=False, verbose=verbose, title="Electronic Bands (k)"))

        #app(f"gstore_cplex: {self.r.cplex}")
        #app(f"gstore_qptopt: {self.r.qptopt}")

        return "\n".join(lines)

    @staticmethod
    def _get_which_g_list(which_g: str) -> list[str]:
        all_choices = ["avg", "raw"]
        if which_g == "all":
            return all_choices

        if which_g not in all_choices:
            raise ValueError(f"Invalid {which=}, should be in {all_choices=}")

        return [which_g]

    def _get_band_range(self, band_range):
        return (self.r.bstart, self.r.bstop) if band_range is None else band_range

    @add_fig_kwargs
    def plot_g_qpath(self, band_range=None, which_g="avg", with_qexp: int=0, scale=1, gmax_mev=250,
                     ph_modes=None, with_phbands=True, with_ebands=False,
                     ax_mat=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the averaged |g(k,q)| in meV units along the q-path

        Args:
            band_range: Band range that will be averaged over (python convention).
            which_g: "avg" to plot the symmetrized |g|, "raw" for unsymmetrized |g|."all" for both.
            with_qexp: Multiply |g(q)| by |q|^{with_qexp}.
            scale: Scaling factor for the marker size used when with_phbands is True.
            gmax_mev: Show results up to gmax in meV.
            ph_modes: List of ph branch indices to show (start from 0). If None all modes are shown.
            with_phbands: False if phonon bands should now be displayed.
            with_ebands: False if electron bands should now be displayed.
            ax_mat: List of |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
        """
        which_g_list = self._get_which_g_list(which_g)
        nrows, ncols = len(which_g_list) + int((np.array([with_ebands, with_phbands]) == True).sum()), self.r.nsppol

        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)
        marker_color = "gold"
        band_range = self._get_band_range(band_range)

        #facts_q, g_label, g_units = self.get_info(which_g, with_qexp)
        facts_q = np.ones(len(self.phbands.qpoints)) if with_qexp == 0 else \
                  np.array([qpt.norm for qpt in self.phbands.qpoints]) ** with_qexp

        q_label = r"$|q|^{%d}$" % with_qexp if with_qexp else ""
        g_units = "(meV)" if with_qexp == 0 else r"(meV $\AA^-{%s}$)" % with_qexp

        for spin in range(self.r.nsppol):
            g_nuq_avg, g_nuq_raw = self.r.get_gnuq_average_spin(spin, band_range)
            ax_cnt = -1

            for which_g in which_g_list:
                # Select ys according to which_g and multiply by facts_q
                g_nuq = dict(avg=g_nuq_avg, raw=g_nuq_raw)[which_g] * facts_q[None,:]

                # Plot g_nu(q)
                ax_cnt += 1
                ax = ax_mat[ax_cnt, spin]
                for nu in range(self.r.natom3):
                    if ph_modes is not None and nu not in ph_modes: continue
                    ax.plot(g_nuq[nu], label=f"{nu=}")
                    self.phbands.decorate_ax(ax, units="meV")
                    g_label = r"$|g^{\text{%s}}_{\mathbf{q}}|$ %s" % (which_g, q_label)
                    set_grid_legend(ax, fontsize, ylabel="%s %s" % (g_label, g_units))

                if gmax_mev is not None and with_qexp == 0:
                    set_axlims(ax, [0, gmax_mev], "y")

            if with_phbands:
                # Plot phonons bands + averaged g(q)  as markers
                ax_cnt += 1
                x, y, s = [], [], []
                for iq, qpoint in enumerate(self.phbands.qpoints):
                    omegas_nu = self.phbands.phfreqs[iq,:]
                    for w, g2 in zip(omegas_nu, g_nuq_avg[:,iq], strict=True):
                        x.append(iq); y.append(w); s.append(scale * g2)

                label = r'$|g^{\text{avg}}_{\mathbf{q}}|$' if with_qexp == 0 else \
                        r'$|g^{\text{avg}}_{\mathbf{q}}| |q|^{%s}$' % with_qexp

                points = Marker(x, y, s, color=marker_color, edgecolors='gray', alpha=0.8, label=label)

                ax = ax_mat[ax_cnt, spin]
                self.phbands.plot(ax=ax, points=points, show=False)
                set_grid_legend(ax, fontsize) #, xlabel=r"Wavevector $\mathbf{q}$")

            if with_ebands:
                # Plot phonons bands + g(q) as markers
                ax_cnt += 1
                ax = ax_mat[ax_cnt, spin]
                self.ebands_kq.plot(ax=ax, spin=spin, band_range=band_range, with_gaps=False, show=False)

        # Add title.
        if (kpt_name := self.structure.findname_in_hsym_stars(self.r.eph_fix_wavec)) is None:
            qpt_name = k2s(self.r.eph_fix_wavec)

        fig.suptitle(f"k = {kpt_name}" + f" m, n = {band_range[0]} - {band_range[1] - 1}")

        return fig

    @add_fig_kwargs
    def plot_g_kpath(self, band_range=None, which_g="avg", scale=1, gmax_mev=250, ph_modes=None,
                    with_ebands=True, ax_mat=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the averaged |g(k,q)| in meV units along the k-path

        Args:
            band_range: Band range that will be averaged over (python convention).
            which_g: "avg" to plot the symmetrized |g|, "raw" for unsymmetrized |g|."all" for both.
            scale: Scaling factor for the marker size used when with_phbands is True.
            gmax_mev: Show results up to gmax in meV.
            ph_modes: List of ph branch indices to show (start from 0). If None all modes are show.
            with_ebands: False if electron bands should now be displayed.
            ax_mat: List of |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
        """
        which_g_list = self._get_which_g_list(which_g)
        nrows, ncols = len(which_g_list) + int((np.array([with_ebands]) == True).sum()), self.r.nsppol

        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)
        marker_color = "gold"
        band_range = self._get_band_range(band_range)

        for spin in range(self.r.nsppol):
            g_nuk_avg, g_nuk_raw = self.r.get_gnuk_average_spin(spin, band_range)
            ax_cnt = -1

            for which_g in which_g_list:
                # Select ys according to which_g
                g_nuk = dict(avg=g_nuk_avg, raw=g_nuk_raw)[which_g]

                # Plot g_nu(q)
                ax_cnt += 1
                ax = ax_mat[ax_cnt, spin]
                for nu in range(self.r.natom3):
                    if ph_modes is not None and nu not in ph_modes: continue
                    ax.plot(g_nuk[nu], label=f"{which_g} {nu=}")

                self.ebands_k.decorate_ax(ax, units="meV")
                set_grid_legend(ax, fontsize, ylabel=r"$|g^{\text{%s}}_{\mathbf{k}}|$ (meV)" % (which_g))
                if gmax_mev is not None:
                    set_axlims(ax, [0, gmax_mev], "y")

            if with_ebands:
                # Plot electron bands
                ax_cnt += 1
                ax = ax_mat[ax_cnt, spin]
                self.ebands_k.plot(ax=ax, spin=spin, band_range=band_range, with_gaps=False, show=False)
                set_grid_legend(ax, fontsize) #, xlabel=r"Wavevector $\mathbf{q}$")

        if (qpt_name := self.structure.findname_in_hsym_stars(self.r.eph_fix_wavec)) is None:
            qpt_name = k2s(self.r.eph_fix_wavec)

        fig.suptitle(f"q = {qpt_name}" + f" m, n = {band_range[0]} - {band_range[1] - 1}")

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        if self.r.eph_fix_korq == "k":
            #yield self.ebands_kq.plot(show=False)
            #yield self.phbands.plot(show=False)
            yield self.plot_g_qpath()

        if self.r.eph_fix_korq == "q":
            #yield self.ebands_k.plot(show=False)
            yield self.plot_g_kpath()

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("gpath = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(gpath)"),
            nbv.new_code_cell("gpath.ebands.plot();"),
            nbv.new_code_cell("gpath.phbands.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class GpathReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GpathReader
    """
    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("nsppol")
        self.natom = self.read_dimvalue("natom")
        self.natom3 = self.read_dimvalue("natom3")
        self.nb_in_g = self.read_dimvalue("nb_in_g")
        self.nk_path = self.read_dimvalue("nk_path")
        self.nq_path = self.read_dimvalue("nq_path")

        self.structure = self.read_structure()

        # eigens are in Ha, phfreq are in eV for historical reason
        self.phfreqs_ha = self.read_value("phfreqs")
        self.all_eigens_k = self.read_value("all_eigens_k")
        self.all_eigens_kq = self.read_value("all_eigens_kq")

        # Read important variables.
        self.eph_fix_korq = self.read_string("eph_fix_korq")
        if self.eph_fix_korq not in {"k", "q"}:
            raise ValueError(f"Invalid value for {self.eph_fix_korq=}")
        self.eph_fix_wavec = self.read_value("eph_fix_wavevec")
        self.dbdb_add_lr = self.read_value("dvdb_add_lr")
        #self.used_ftinterp = self.read_value("used_ftinterp")
        #self.completed = self.read_value("gstore_completed")

        # Note conversion Fortran --> C for the bstart index.
        nband = self.read_dimvalue("nband")
        self.bstart = self.read_value("bstart") - 1
        self.bstop = self.read_value("bstop")
        self.band_range = [self.bstart, self.bstop]

    def read_ebands_which_fixed(self, which_fixed: str):
        """
        Overrides method of superclass as we cannot rely of the etsf-io file format,
        and we have to build the ebands manually.
        """
        nspinor = self.read_dimvalue("nspinor")
        nspden = self.read_dimvalue("nspden")
        nelect = self.read_value("nelect")
        fermie = self.read_value("fermie") * abu.Ha_eV

        structure = self.read_structure()
        kpath_frac_coords = self.read_value("kpoints")
        qpath_frac_coords = self.read_value("qpoints")

        # eigens are in Ha
        if which_fixed == "k":
            frac_coords = kpath_frac_coords + qpath_frac_coords
            all_eigens = self.read_value("all_eigens_kq") * abu.Ha_eV
        elif which_fixed == "q":
            frac_coords = kpath_frac_coords
            all_eigens = self.read_value("all_eigens_k") * abu.Ha_eV
        else:
            raise ValueError(f"Invalid value {which_fixed=}")

        occfacts = np.zeros_like(all_eigens)
        path = Kpath(structure.lattice.reciprocal_lattice, frac_coords)

        #print(f"Before ElectronBands {len(path)=}, {all_eigens.shape=}")
        #print(path)

        return ElectronBands(structure, path, all_eigens, fermie, occfacts, nelect, nspinor, nspden)

    def read_phbands(self) -> PhononBands:
        """
        Read the phonon band structure along the q-path.
        """
        amu_list = self.read_value("atomic_mass_units")
        atomic_numbers = self.read_value("atomic_numbers")
        amu = {at: a for at, a in zip(atomic_numbers, amu_list)}

        # phfreqs are in eV for historical reason
        phfreqs = self.read_value("phfreqs")
        # Complex array with the Cartesian displacements in **Angstrom**
        phdispl_cart = self.read_value("phdispl_cart", cmode="c") * abu.Bohr_Ang

        structure = self.read_structure()
        qpath_frac_coords = self.read_value("qpoints")

        path_qq = Kpath(structure.lattice.reciprocal_lattice, qpath_frac_coords, weights=None, names=None, ksampling=None)

        non_anal_ph = NonAnalyticalPh.from_ncreader(self) if "non_analytical_directions" in self.rootgrp.variables else None

        return PhononBands(structure=structure, qpoints=path_qq, phfreqs=phfreqs, phdispl_cart=phdispl_cart, amu=amu,
                           non_anal_ph=non_anal_ph,
                           # TODO ?
                           #epsinf=epsinf,
                           #zcart=zcart,
                           )

    def get_gnuq_average_spin(self, spin: int, band_range: list|tuple|None, eps_mev: float=0.01) -> tuple:
        """
        Average e-matrix elements over phonon modes, and k- k+q electrons when the matrix elements
        have been computed along a q-path.

        Args:
            spin: spin index
            band_range: Band range that will be averaged over (python convention).
            eps_mev: Tolerance in meV used to detect degeneracies for phonons and electrons.

        Return:
            tuple with two numpy array
        """
        # Consistency check
        if self.nk_path != 1:
            raise ValueError(f"{self.nk_path=} != 1. In this case, one cannot ask for q-dependent g(k,q)!")

        # Tolerences to detect degeneracies.
        eps_ha = eps_mev / abu.Ha_meV
        eps_ev = eps_ha * abu.Ha_eV

        # Number of m, n bands in g_mn, the first band starts at bstart.
        nb_in_g = self.nb_in_g
        bstart, bstop = self.bstart, self.bstop
        nsppol, natom3 = self.nsppol, self.natom3

        # double all_eigens_k(nsppol, nk_path, nband) ;
        # double all_eigens_kq(nsppol, nq_path, nband) ;
        all_eigens_k, all_eigens_kq = self.all_eigens_k, self.all_eigens_kq  # eV units
        phfreqs_ha = self.phfreqs_ha                                        # Ha units

        # Now read the e-ph matrix elements. On disk we have
        #                                                  n-index, m-index
        # double gkq2_nu(nsppol, nk_path, nq_path, natom3, nb_in_g, nb_in_g) ;
		#  gkq2_nu:_FillValue = -1. ;
        #
        # in Ha^2 with nk_path == 1
        #                                      m-index, n-index
        # In memory we want: (nq_path, natom3, nb_in_g, nb_in_g)

        absg = np.sqrt(self.read_variable("gkq2_nu")[spin, 0][:].transpose(0, 1, 3, 2).copy()) * abu.Ha_meV
        absg_raw = absg.copy()

        # Average over degenerate phonon modes.
        absg_avg = np.zeros_like(absg)
        for iq in range(self.nq_path):
            for nu in range(natom3):
                # Sum the squared values of absg over the degenerate phonon mu indices.
                mask_nu = np.abs(phfreqs_ha[iq, :] - phfreqs_ha[iq, nu]) < eps_ha
                g2_mn = np.sum(absg[iq, mask_nu, :, :]**2, axis=0)
                # Compute the symmetrized value and divide by the number of degenerate ph-modes for this iq.
                absg_avg[iq, nu, :, :] = np.sqrt(g2_mn / np.sum(mask_nu))

        # MG FIXME: Note the difference with a similar function in gkq here I use absg and not absgk
        # Average over degenerate k electrons taking bstart into account.
        absg = absg_avg.copy()
        g2_nu = np.zeros((natom3), dtype=float)
        for iq in range(self.nq_path):
            for m_kq in range(nb_in_g):
                for n_k in range(nb_in_g):
                    w_1 = all_eigens_k[spin, 0, n_k + bstart]
                    g2_nu[:], nn = 0.0, 0
                    for bsum_k in range(nb_in_g):
                        w_2 = all_eigens_k[spin, 0, bsum_k + bstart]
                        if abs(w_2 - w_1) >= eps_ev: continue
                        nn += 1
                        g2_nu += absg[iq,:,m_kq,bsum_k] ** 2
                    absg_avg[iq,:,m_kq,n_k] = np.sqrt(g2_nu / nn)

        # Average over degenerate k+q electrons taking bstart into account.
        absg = absg_avg.copy()
        for iq in range(self.nq_path):
          for n_k in range(nb_in_g):
              for m_kq in range(nb_in_g):
                  w_1 = all_eigens_kq[spin, iq, m_kq + bstart]
                  g2_nu[:], nn = 0.0, 0
                  for bsum_kq in range(nb_in_g):
                      w_2 = all_eigens_kq[spin, iq, bsum_kq + bstart]
                      if abs(w_2 - w_1) >= eps_ev: continue
                      nn += 1
                      g2_nu += absg[iq,:,bsum_kq,n_k] ** 2
                  absg_avg[iq,:,m_kq,n_k] = np.sqrt(g2_nu / nn)

        # Transpose the data: (nq_path, natom3, nb_in_g, nb_in_g) -> (natom3, nq_path, nb_in_g, nb_in_g)
        absg_avg, absg_raw = absg_avg.transpose(1, 0, 2, 3).copy(), absg_raw.transpose(1, 0, 2, 3).copy()

        # Slice the last two band dimensions if band_range is given in input.
        nb = nb_in_g
        if band_range is not None:
            nb = band_range[1] - band_range[0]
            b0, b1 = band_range[0] - bstart, band_range[1] - bstart
            absg_avg, absg_raw = absg_avg[..., b0:b1, b0:b1], absg_raw[..., b0:b1, b0:b1]

        # Average over bands 1/n_b**2 sum_{mn}
        return np.sum(absg_avg, axis=(-2, -1)) / nb**2, np.sum(absg_raw, axis=(-2, -1)) / nb**2

    def get_gnuk_average_spin(self, spin: int, band_range: list|tuple|None, eps_mev: float=0.01) -> tuple:
        """
        Average g elements over phonon modes, and k- k+q electrons when the matrix elements
        have been computed along a k-path.

        Args:
            spin: spin index
            band_range: Band range that will be averaged over (python convention).
            eps_mev: Tolerance in meV used to detect degeneracies for phonons and electrons.

        Return:
            tuple with two numpy array
        """
        # Consistency check
        if self.nq_path != 1:
            raise ValueError(f"{self.nq_path=} != 1. In this case, one cannot ask for l-dependent g(k,q)!")

        # Tolerences to detect degeneracies.
        eps_ha = eps_mev / abu.Ha_meV
        eps_ev = eps_ha * abu.Ha_eV

        # Number of m, n bands in g_mn, the first band starts at bstart.
        nb_in_g = self.nb_in_g
        bstart, bstop = self.bstart, self.bstop
        nsppol, natom3 = self.nsppol, self.natom3

        # double all_eigens_k(nsppol, nk_path, nband) ;
        # double all_eigens_kq(nsppol, nq_path, nband) ;
        all_eigens_k, all_eigens_kq = self.all_eigens_k, self.all_eigens_kq  # eV units
        phfreqs_ha = self.phfreqs_ha                                         # Ha units

        # Now read the e-ph matrix elements. On disk we have
        #                                                  n-index, m-index
        # double gkq2_nu(nsppol, nk_path, nq_path, natom3, nb_in_g, nb_in_g) ;
		#  gkq2_nu:_FillValue = -1. ;
        #
        # in Ha^2 with nq_path == 1
        #                                      m-index, n-index
        # In memory we want: (nk_path, natom3, nb_in_g, nb_in_g)

        absg = np.sqrt(self.read_variable("gkq2_nu")[spin, :, 0, :, :,:][:].transpose(0, 1, 3, 2).copy()) * abu.Ha_meV
        absg_raw = absg.copy()

        # Average over degenerate phonon modes for this q
        iq = 0
        absg_avg = np.zeros_like(absg)
        for ik in range(self.nk_path):
            for nu in range(natom3):
               # Sum the squared values of absg over the degenerate phonon mu indices.
               mask_nu = np.abs(phfreqs_ha[iq, :] - phfreqs_ha[iq, nu]) < eps_ha
               g2_mn = np.sum(absg[ik, mask_nu, :, :]**2, axis=0)
               # Compute the symmetrized value and divide by the number of degenerate ph-modes for this iq.
               absg_avg[ik, nu, :, :] = np.sqrt(g2_mn / np.sum(mask_nu))

        # MG FIXME: Note the difference with a similar function in gkq here I use absg and not absgk
        # Average over degenerate k electrons taking bstart into account.
        absg = absg_avg.copy()
        g2_nu = np.zeros((natom3), dtype=float)
        for ik in range(self.nk_path):
            for m_kq in range(nb_in_g):
                for n_k in range(nb_in_g):
                    w_1 = all_eigens_k[spin, ik, n_k + bstart]
                    g2_nu[:], nn = 0.0, 0
                    for bsum_k in range(nb_in_g):
                        w_2 = all_eigens_k[spin, ik, bsum_k + bstart]
                        if abs(w_2 - w_1) >= eps_ev: continue
                        nn += 1
                        g2_nu += absg[ik,:,m_kq,bsum_k] ** 2
                    absg_avg[ik,:,m_kq,n_k] = np.sqrt(g2_nu / nn)

        # Average over degenerate k+q electrons taking bstart into account.
        absg = absg_avg.copy()
        for n_k in range(nb_in_g):
            for m_kq in range(nb_in_g):
                w_1 = all_eigens_kq[spin, 0, m_kq + bstart]
                g2_nu[:], nn = 0.0, 0
                for bsum_kq in range(nb_in_g):
                    w_2 = all_eigens_kq[spin, 0, bsum_kq + bstart]
                    if abs(w_2 - w_1) >= eps_ev: continue
                    nn += 1
                    g2_nu += absg[ik,:,bsum_kq,n_k] ** 2
                absg_avg[ik,:,m_kq,n_k] = np.sqrt(g2_nu / nn)

        # Transpose the data: (nk_path, natom3, nb_in_g, nb_in_g) -> (natom3, nk_path, nb_in_g, nb_in_g)
        absg_avg, absg_raw = absg_avg.transpose(1, 0, 2, 3).copy(), absg_raw.transpose(1, 0, 2, 3).copy()

        # Slice the last two band dimensions if band_range is given in input.
        nb = nb_in_g
        if band_range is not None:
            nb = band_range[1] - band_range[0]
            b0, b1 = band_range[0] - bstart, band_range[1] - bstart
            absg_avg, absg_raw = absg_avg[..., b0:b1, b0:b1], absg_raw[..., b0:b1, b0:b1]

        # Average over bands 1/n_b**2 sum_{mn}
        return np.sum(absg_avg, axis=(-2, -1)) / nb**2, np.sum(absg_raw, axis=(-2, -1)) / nb**2


class GpathRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple GPATH.nc files.

    Usage example:

    .. code-block:: python

        robot = GpathRobot.from_files([
            "t04o_GPATH.nc",
            "t05o_GPATH.nc",
            ])

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GpathRobot
    """
    EXT = "GPATH"

    @add_fig_kwargs
    def plot_g_qpath(self, which_g="avg", gmax_mev=250, ph_modes=None,
                    colormap="jet", **kwargs) -> Figure:
        """
        Compare the g-matrix along a q-path.

        Args
            which_g: "avg" to plot the symmetrized |g|, "raw" for unsymmetrized |g|."all" for both.
            gmax_mev: Show results up to gmax in me
            ph_modes: List of ph branch indices to show (start from 0). If None all modes are show.
            colormap: Color map. Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
        """
        nsppol, nq_path, natom3, eph_fix_wavec, eph_fix_korq = self.getattrs_alleq(
            "nsppol", "nq_path", "natom3", "eph_fix_wavec", "eph_fix_korq"
        )
        xs = np.arange(nq_path)

        nrows, ncols = 1, nsppol
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)
        cmap = plt.get_cmap(colormap)

        # TODO: Compute common band range.
        band_range = None
        ref_ifile= 0
        #q_label = r"$|q|^{%d}$" % with_qexp if with_qexp else ""
        #g_units = "(meV)" if with_qexp == 0 else r"(meV $\AA^-{%s}$)" % with_qexp

        for spin in range(nsppol):
            ax_cnt = 0
            ax = ax_mat[ax_cnt, spin]

            for ifile, gpath in enumerate(self.abifiles):
                g_nuq_avg, g_nuq_raw = gpath.r.get_gnuq_average_spin(spin, band_range)
                # Select ys according to which_g and multiply by facts_q
                g_nuq = dict(avg=g_nuq_avg, raw=g_nuq_raw)[which_g] # * facts_q[None,:]

                for nu in range(natom3):
                    if ph_modes is not None and nu not in ph_modes: continue
                    color = cmap(nu / natom3)
                    if ifile == ref_ifile:
                        ax.scatter(xs, g_nuq[nu], color=color, label=f"{nu=}", marker="o")
                        gpath.phbands.decorate_ax(ax, units="meV")
                        #g_label = r"$|g^{\text{%s}}_{\mathbf{q}}|$ %s" % (which_g, q_label)
                        #set_grid_legend(ax, fontsize, ylabel="%s %s" % (g_label, g_units))
                    else:
                        ax.plot(g_nuq[nu], color=color, label=f"{nu=}")

            #if gmax_mev is not None and with_qexp == 0:
            if gmax_mev is not None:
                set_axlims(ax, [0, gmax_mev], "y")

        return fig

    #@add_fig_kwargs
    #def plot_g_kpath(self, **kwargs) --> Figure

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        #for fig in self.get_ebands_plotter().yield_figs(): yield fig

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.GstoreRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        # Mixins
        #nb.cells.extend(self.get_baserobot_code_cells())
        #nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

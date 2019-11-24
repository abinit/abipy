# coding: utf-8
"""
Objects to compute electronic effective masses via finite differences starting from an |ElectronBands| object.
"""

import numpy as np
import pandas as pd
import pymatgen.core.units as units

from collections import OrderedDict
from monty.termcolor import cprint
from abipy.core.mixins import Has_Structure, Has_ElectronBands
from abipy.tools.derivatives import finite_diff
from abipy.tools.printing import print_dataframe
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_visible


class EffMassAnalyzer(Has_Structure, Has_ElectronBands):
    """
    This objects provides a high-level API to compute electronic effective masses
    via finite differences starting from an |ElectronBands| object.

    Usage example:

    .. code-block:: python

        emana = EffmassAnalyzer.from_file("out_GSR.nc")
        print(emana)
        emana.select_vbm()
        #emana.select_kpoint_band(kpoint=[0, 0, 0], band=3)
        emans.plot_emass()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EffMassAnalyzer
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf file with an ebands object."""
        from abipy.abilab import abiopen
        with abiopen(filepath) as ncfile:
            return cls(ncfile.ebands, copy=False)

    def __init__(self, ebands, copy=True):
        """Initialize the object from an ebands object with k-points along a path."""
        if not ebands.kpoints.is_path:
            raise ValueError("EffmassAnalyzer requires k-points along a path. Got:\n %s" % repr(ebands.kpoints))

        # Copy ebands before changing fermie because we don't want side effects!
        self._ebands = ebands.deepcopy() if copy else ebands
        self._ebands.set_fermie_to_vbm()
        self.segments = []

    def __repr__(self):
        """Invoked by repr"""
        return repr(self.ebands)

    def __str__(self):
        """Invoked by str"""
        return self.ebands.to_string()

    def to_string(self, verbose):
        """
        Human-readable string with useful info such as band gaps, position of HOMO, LOMO...

        Args:
            verbose: Verbosity level.
        """
        return self.ebands.to_string(with_structure=True, with_kpoints=True, verbose=verbose)

    def select_kpoint_band(self, kpoint, band, spin=0, etol_ev=1e-3):
        """
        Construct line segments based on the k-point ``kpoint`` and the band index ``band``.
        This is the most flexible interface and allows one to include extra bands
        by increasing the value of `etol_ev`.

        Args:
            kpoint: Fractional coordinates or |Kpoint| object.
            band: Band index
            spin: Spin index.
            etol_ev: Include all bands at this k-point whose energy differ from the one at
                (kpoint, band) less that ``etol_ev`` in eV.
        """
        # Find all k-indices associated to the input kpoint.
        ik_indices = self.kpoints.get_all_kindices(kpoint)
        ik = ik_indices[0]

        # Find band indices of "degenerate" states within tolerance etol_ev.
        e0 = self.ebands.eigens[spin, ik, band]
        band_inds_k = [be[0] for be in enumerate(self.ebands.eigens[spin, ik]) if abs(be[1] - e0) <= etol_ev]
        band_inds_k = len(ik_indices) * [band_inds_k]
        self._build_segments(spin, ik_indices, band_inds_k)

    def select_cbm(self, spin=0, etol_ev=1e-3):
        """Select conduction band minimum."""
        ik_indices, band_inds_k = self._select(["cbm"], spin, etol_ev)
        self._build_segments(spin, ik_indices, band_inds_k)

    def select_vbm(self, spin=0, etol_ev=1e-3):
        """Select valence band maximum."""
        ik_indices, band_inds_k = self._select(["vbm"], spin, etol_ev)
        self._build_segments(spin, ik_indices, band_inds_k)

    def select_band_edges(self, spin=0, etol_ev=1e-3):
        """Select conduction band minimum and valence band maximum."""
        ik_indices, band_inds_k = self._select(["cbm", "vbm"], spin, etol_ev)
        self._build_segments(spin, ik_indices, band_inds_k)

    def _select(self, what_list, spin, etol_ev):
        ik_indices, band_inds_k = [], []

        if "vbm" in what_list:
            homo = self.ebands.homos[spin]
            homo_iks = self.kpoints.get_all_kindices(homo.kpoint)
            ik = homo_iks[0]
            e0 = self.ebands.eigens[spin, ik, homo.band]
            homo_bands = [be[0] for be in enumerate(self.ebands.eigens[spin, ik]) if abs(be[1] - e0) <= etol_ev]
            ik_indices.extend(homo_iks)
            for i in range(len(homo_iks)):
                band_inds_k.append(homo_bands)

        if "cbm" in what_list:
            lumo = self.ebands.lumos[spin]
            lumo_iks = self.kpoints.get_all_kindices(lumo.kpoint)
            ik = lumo_iks[0]
            e0 = self.ebands.eigens[spin, ik, lumo.band]
            lumo_bands = [be[0] for be in enumerate(self.ebands.eigens[spin, ik]) if abs(be[1] - e0) <= etol_ev]
            ik_indices.extend(lumo_iks)
            for i in range(len(lumo_iks)):
                band_inds_k.append(lumo_bands)

        return ik_indices, band_inds_k

    def _build_segments(self, spin, ik_indices, band_inds_k):
        """
        Build list of segments.

        Args:
            spin: Spin index.
            ik_indices: List of k-point indices [nk][kids]
            band_inds_k: [nk][kids]
        """
        #print("in build_segments with:\n\tik_indices:", ik_indices, "\n\tband_inds_k:", band_inds_k)
        self.spin = spin

        dims = len(ik_indices), len(band_inds_k)
        if dims[0] != dims[1]:
            raise RuntimeError("len(ik_indices) %s != len(band_inds_k) %s" % (dims[0], dims[1]))

        self.segments = []
        for ik, bids in zip(ik_indices, band_inds_k):
            for iline, line in enumerate(self.kpoints.lines):
                if line[-1] >= ik >= line[0]: break
            else:
                raise ValueError("Cannot find k-index `%s` in lines: `%s`" % (ik, self.kpoints.lines))
            #print("line[-1]", line[-1], "ik", ik, "line[0]", line[0])

            self.segments.append(Segment(ik, spin, line, bids, self.ebands))

    @property
    def ebands(self):
        """|ElectronBands| object."""
        return self._ebands

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    def _consistency_check(self):
        if not self.segments:
            methods = ", ".join(["select_cbm", "select_vbm", "select_band_edges", "select_kpoint_band"])
            raise RuntimeError("You must call one among: `%s`\nto build segments before analyzing data." % methods)

    def summarize(self, acc=4):
        self._consistency_check()
        for segment in self.segments:
            df = segment.get_dataframe_with_accuracies(acc_list=[acc])
            title = "k: %s, spin: %s, nbands: %d" % (repr(segment.k0), segment.spin, segment.nb)
            #label="direction: %s" % segment.kdir.tos(m="fracart") if ib == 0 else None)
            print_dataframe(df, title=title)
            #print("")

    #def print_segments(self):
    #    self._consistency_check()
    #    for segment in self.segments

    @add_fig_kwargs
    def plot_emass(self, acc=4, fontsize=6, colormap="viridis", **kwargs):
        """
        Plot electronic dispersion and quadratic curve based on the
        effective masses computed along each segment.

        Args:
            acc:
            fontsize: legend and title fontsize.
            colormap: matplotlib colormap
        """
        self._consistency_check()

        # Build grid of plots for this spin.
        num_plots, ncols, nrows = len(self.segments), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        for iseg, (segment, ax) in enumerate(zip(self.segments, ax_list)):
            irow, icol = divmod(iseg, ncols)
            segment.plot_emass(ax=ax, acc=acc, fontsize=fontsize, colormap=colormap, show=False)
            if iseg != 0: set_visible(ax, False, "ylabel")
            if irow != nrows - 1: set_visible(ax, False, "xticklabels")

        # don't show the last ax if numeb is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        return fig

    @add_fig_kwargs
    def plot_all_segments(self, ax=None, colormap="viridis", fontsize=8, **kwargs):
        """

        Args:
            colormap: matplotlib colormap
            fontsize: legend and title fontsize.

        Return: |matplotlib-Figure|
        """
        self._consistency_check()

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        markers = ["o", "s", "x", "D", "+", "v", ">", "<"] * 8

        ax.grid(True)
        ax.set_ylabel('Energy (eV)')
        pad = 0
        for iseg, segment in enumerate(self.segments):
            color = cmap(float(iseg / len(self.segments)))
            for ib in range(segment.nb):
                ax.plot(segment.kpoint_indices + pad, segment.energies_bk[ib],
                        linestyle=":", marker=markers[ib], markersize=2, color=color,
                        label="direction: %s" % segment.kdir.tos(m="fracart") if ib == 0 else None)
            pad += 10

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        #title = "k: %s, spin: %s, nband: %d" % (repr(self.efm_kpoint), self.spin, segment.nb)
        #ax.set_title(title, fontsize=fontsize)

        return fig


class Segment:
    """
        .. attribute:: Z

            Atomic number
    """
    def __init__(self, ik, spin, line, band_inds, ebands):
        """
        Args:
            kpoint_indices (list(int)): the kpoint indices of the segment
        """
        self.ik = ik
        self.spin = spin
        self.kpos = line.index(ik)
        self.kpos_type = "central"
        if self.kpos == 0: self.kpos_type = "right"
        if self.kpos == len(line) - 1: self.kpos_type = "left"

        self.kdir = ebands.kpoints.versors[line[0]]
        self.dk = ebands.kpoints.ds[line[0]]
        if not np.allclose(self.dk, ebands.kpoints.ds[line[:-1]]):
            raise ValueError("For finite difference derivatives, the path must be homogeneous!\n" +
                             str(ebands.kpoints.ds[line[:-1]]))

        self.kpoint_indices = np.asarray(line)
        self.band_inds = band_inds
        # The reference point and |k-k0|^2
        self.k0 = ebands.kpoints[ik]
        self.kmk0_2 = np.array([(k - self.k0).norm ** 2 for k in (ebands.kpoints[l] for l in line)])

        self.energies_bk = []
        for band in band_inds:
            self.energies_bk.append(ebands.eigens[spin, line, band])

        self.energies_bk = np.reshape(self.energies_bk, (len(band_inds), len(line)))
        self.nb = len(self.energies_bk)
        self.ebands = ebands

    def __repr__(self):
        return "k0: %s, kdir: %s, dk: %.3f (Ang-1)" % (repr(self.k0), repr(self.kdir), self.dk)

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app("k-point: %s, nband: %s, spin: %d" % (self.k0.to_string(verbose=verbose), self.nb, self.spin))
        return "\n".join(lines)

    def __str__(self):
        return self.to_string()

    def get_fd_emass_d2(self, enes_kline, acc):
        d2 = finite_diff(enes_kline, self.dk, order=2, acc=acc, index=self.kpos)
        emass = 1. / (d2.value * (units.eV_to_Ha / units.bohr_to_ang ** 2))
        return emass, d2

    def get_dataframe_with_accuracies(self, acc_list=(2, 4, 6, 8)):
        """
        Build and return a |pandas-Dataframe| with effective masses computed with different accuracies (npts)
        """
        rows = []
        for acc in acc_list:
            emass_dict = OrderedDict()
            for ib, enes_kline in enumerate(self.energies_bk):
                emass, d2 = self.get_fd_emass_d2(enes_kline, acc)
                emass_dict["m%d" % ib] = emass

            od = OrderedDict([
                ("acc", acc), ("npts", d2.npts),
            ])
            od.update(emass_dict)
            rows.append(od)

        return pd.DataFrame(rows, columns=list(rows[0].keys()))

    @add_fig_kwargs
    def plot_emass(self, acc=4, ax=None, fontsize=8, colormap="viridis", **kwargs):
        """

        Args:
            acc:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: legend and title fontsize.
            colormap: matplotlib colormap

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)
        cmap = plt.get_cmap(colormap)

        for ib, enes_kline in enumerate(self.energies_bk):
            # Plot KS-DFT points.
            xs = range(len(enes_kline))
            ax.scatter(xs, enes_kline, marker="o", color=cmap(float(ib) / self.nb))

            # Compute effective masses
            try:
                emass, d2 = self.get_fd_emass_d2(enes_kline, acc)
            except Exception as exc:
                cprint("Eception for segment: %s" % str(self), "red")
                raise exc

            ys = ((self.kmk0_2 * units.bohr_to_ang ** 2) / (2 * emass)) * units.Ha_to_eV + self.energies_bk[ib, self.kpos]
            label = r"$m^*$ = %.3f, %d-pts %s finite-diff" % (emass, d2.npts, d2.mode)
            ax.plot(xs, ys, linestyle="--", color=cmap(float(ib) / self.nb), label=label)

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        title = r"${\bf k}_0$: %s, direction: %s" % (repr(self.k0), self.kdir.tos(m="fracart"))
        ax.set_title(title, fontsize=fontsize)
        ax.set_ylabel('Energy (eV)')

        return fig


#class MultipleEffMassAnalyzer:
#
#    @classmethod
#    def from_files(cls, filepaths):
#        """Initialize the object from a list of files providing an |ElectronBands| object."""
#        return cls([ElectronBands.from_file(f, copy=False) for f in list_strings(filepaths)])
#
#    def __init__(self, eband_list)
#        self.eband_list = eband_list

# coding: utf-8
"""
"""

import numpy as np
import pandas as pd

from collections import OrderedDict
import pymatgen.core.units as units
from abipy.core.mixins import Has_Structure, Has_ElectronBands #, NotebookWriter
from abipy.tools.derivatives import finite_diff
from abipy.tools.printing import print_dataframe
#from abipy.tools import duck
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt #, set_axlims,


#class EffMassAnalyzer(Has_Structure, Has_ElectronBands, NotebookWriter):
class EffMassAnalyzer(Has_Structure, Has_ElectronBands):
    """
    Usage example:

    .. code-block:: python

        emana = EffmassAnalyzer.from_file("out_GSR.nc")
        print(emana)
        emana.set_kpoint_band(kpoint=[0, 0, 0], band=3)
        emans.plot_emass()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EffmassAnalyzer
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf file with an ebands object."""
        from abipy.abilab import abiopen
        with abiopen(filepath) as ncfile:
            return cls(ncfile.ebands)

    def __init__(self, ebands):
        """Initialize the object from an ebands object with k-points along a path."""
        self._ebands = ebands
        if not ebands.kpoints.is_path:
            raise ValueError("EffmassAnalyzer requires k-points along a path. Got:\n %s" % repr(ebands.kpoints))

        self.efm_kpoint = None

    def set_kpoint_band(self, kpoint, band, degtol_ev=1e-3):
        all_iks = self.kpoints.get_all_kindexes(kpoint)
        ik = all_iks[0]
        self.efm_kpoint = self.kpoints[ik]
        #self.efm_band = band
        #self.efm_degtol_ev = degtol_ev

        # Find indices of degenerate states.
        band_indices_spin = []
        for spin in range(self.nsppol):
            e0 = self.ebands.eigens[spin, ik, band]
            inds = [enum[0] for enum in enumerate(self.ebands.eigens[spin, ik]) if abs(enum[1] - e0) <= degtol_ev]
            band_indices_spin.append(np.array(inds))

        self.segments_spin = self.nsppol * [[]]
        for ik in self.kpoints.get_all_kindexes(kpoint):
            for iline, line in enumerate(self.kpoints.lines):
                if line[-1] >= ik >= line[0]: break
            else:
                raise ValueError("Cannot find k-index `%s` in lines: `%s`" % (ik, self.kpoints.lines))

            for spin in range(self.nsppol):
                segment = Segment(ik, spin, line, band_indices_spin[spin], self.ebands)
                self.segments_spin[spin].append(segment)

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        app(self.kpoints.to_string(verbose=verbose, title="List of k-points"))

        return "\n".join(lines)

    @property
    def ebands(self):
        """|ElectronBands| object."""
        return self._ebands

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    def _consistency_check(self):
        if not hasattr(self, "segments_spin"):
            raise RuntimeError("You must call set_kpoint_band to select the k-point and the bands")

    def summarize(self, acc=4):
        self._consistency_check()
        for spin in range(self.nsppol):
            for segment in self.segments_spin[spin]:
                df = segment.get_dataframe_with_accuracies(acc_list=[acc])
                title = "k: %s, spin: %s, nbands: %d" % (repr(self.efm_kpoint), spin, segment.nb)
                #label="direction: %s" % segment.kdir.tos(m="fracart") if ib == 0 else None)
                print_dataframe(df, title=title)
                #print("")

    @add_fig_kwargs
    def plot_emass(self, acc=4, fontsize=8, colormap="viridis", **kwargs):
        """
        Plot electronic dispersion and quadratic expression based on effective mass computed along each segment.

        Args:
            acc:
            fontsize: legend and title fontsize.
            colormap: matplotlib colormap
        """
        self._consistency_check()

        # Build grid of plots.
        nrows, ncols = max(len(segs) for segs in self.segments_spin), self.nsppol
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=True, squeeze=False)

        for spin in range(self.nsppol):
            for iseg, segment in enumerate(self.segments_spin[spin]):
                segment.plot_emass(ax=ax_mat[iseg, spin], acc=acc, fontsize=fontsize, colormap=colormap, show=False)

        return fig

    @add_fig_kwargs
    def plot_all_segments(self, colormap="viridis", fontsize=8, **kwargs):
        """

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: matplotlib colormap
            fontsize: legend and title fontsize.

        Return: |matplotlib-Figure|
        """
        self._consistency_check()

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=1, ncols=self.nsppol,
                                               sharex=False, sharey=True, squeeze=False)
        cmap = plt.get_cmap(colormap)
        markers = ["o", "s", "x", "D", "+", "v", ">", "<"] * 8

        for spin, ax in enumerate(ax_mat.ravel()):
            ax.grid(True)
            ax.set_ylabel('Energy (eV)')
            pad = 0
            for iseg, segment in enumerate(self.segments_spin[spin]):
                color = cmap(float(iseg / len(self.segments_spin[spin])))
                for ib in range(segment.nb):
                    ax.plot(segment.kpoint_indices + pad, segment.energies_bk[ib],
                            linestyle=":", marker=markers[ib], markersize=4, color=color,
                            label="direction: %s" % segment.kdir.tos(m="fracart") if ib == 0 else None)
                pad += 10

            ax.legend(loc="best", fontsize=fontsize, shadow=True)
            title = "k: %s, spin: %s, nband: %d" % (repr(self.efm_kpoint), spin, segment.nb)
            ax.set_title(title, fontsize=fontsize)

        return fig

    #def yield_figs(self, **kwargs):  # pragma: no cover
    #    """
    #    This function *generates* a predefined list of matplotlib figures with minimal input from the user.
    #    """
    #    #for fig in self.yield_structure_figs(**kwargs): yield fig
    #    for fig in self.yield_ebands_figs(**kwargs): yield fig

    #def write_notebook(self, nbpath=None):
    #    """
    #    Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
    #    working directory is created. Return path to the notebook.
    #    """
    #    nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

    #    nb.cells.extend([
    #        nbv.new_code_cell("eskw = abilab.abiopen('%s')" % self.filepath),
    #        nbv.new_code_cell("print(eskw)"),
    #    ])

    #    return self._write_nb_nbpath(nb, nbpath)


class Segment:

    def __init__(self, ik, spin, line, band_indices, ebands):
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
        self.band_indices = band_indices
        # The reference point and |k - k0|^2
        self.k0 = ebands.kpoints[ik]
        self.kmk0_2 = np.array([(k - self.k0).norm ** 2 for k in (ebands.kpoints[l] for l in line)])

        self.energies_bk = []
        for band in band_indices:
            self.energies_bk.append(ebands.eigens[spin, line, band])

        self.energies_bk = np.reshape(self.energies_bk, (len(band_indices), len(line)))
        self.nb = len(self.energies_bk)
        self.ebands = ebands

    #def __repr__(self):
    #    return "emass_left: %s, emass_right: %s" % (self.em_left, self.em_right)

    def to_string(self, verbose=0):
        """String representation."""
        lines = ["foo"]; app = lines.append
        #app("For spin: %s, band: %s, k-point: %s, eig: %.3f [eV], accuracy: %s" % (
        #    self.spin, self.band, repr(self.kpoint), self.eig, self.acc))
        #app("K-point: %s, eigenvalue: %s (eV)" % (repr(self.kpoint), self.eig))
        return "\n".join(lines)

    def __str__(self):
        return self.to_string()

    def get_fd_emass_d2(self, enes_kline, acc):
        d2 = finite_diff(enes_kline, self.dk, order=2, acc=acc, index=self.kpos)
        emass = 1. / (d2.value * (units.eV_to_Ha / units.bohr_to_ang ** 2))
        return emass, d2

    def get_dataframe_with_accuracies(self, acc_list=(2, 4, 6, 8)):
        """
        Build and return a pandas dataframe with effective masses computed with different accuracies (npts)
        """
        rows = []
        for acc in acc_list:
            emass_dict = OrderedDict()
            for ib, enes_kline in enumerate(self.energies_bk):
                emass, d2 = self.get_fd_emass_d2(enes_kline, acc)
                emass_dict["m%d" % ib] = emass

            od = OrderedDict([
                ("acc", acc),
                ("npts", d2.npts),
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
            emass, d2 = self.get_fd_emass_d2(enes_kline, acc)
            ys = ((self.kmk0_2 * units.bohr_to_ang ** 2) / (2 * emass)) * units.Ha_to_eV + self.energies_bk[ib, self.kpos]
            label = r"$m^*$ = %.3f, %d-pts %s finite-diff" % (emass, d2.npts, d2.mode)
            ax.plot(xs, ys, linestyle="--", color=cmap(float(ib) / self.nb), label=label)

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        title = "$m^*$  at k: %s, along direction: %s" % (repr(self.k0), self.kdir.tos(m="fracart"))
        ax.set_title(title, fontsize=fontsize)
        ax.set_ylabel('Energy (eV)')
        #ax1.set_xlabel('Energy (eV)')

        return fig

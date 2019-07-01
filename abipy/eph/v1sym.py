# coding: utf-8
"""
Object to analyze the results stored in the V1SYM.nc file (mainly for debugging purposes)
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict
from monty.string import marquee
from monty.functools import lazy_property
from abipy.tools.plotting import add_fig_kwargs, get_axarray_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.core.kpoints import KpointList, Kpoint
from abipy.iotools import ETSF_Reader
from abipy.tools import duck


class V1symFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super(V1symFile, self).__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        # Read dimensions.
        self.nfft = r.read_dimvalue("nfft")
        self.nspden = r.read_dimvalue("nspden")
        self.natom3 = len(self.structure) * 3
        self.symv1scf = r.read_value("symv1scf")
        # Read FFT mesh.
        #self.ngfft = r.read_value("ngfft")

    @lazy_property
    def structure(self):
        """|Structure| object."""
        return self.reader.read_structure()

    @lazy_property
    def pertsy_qpt(self):
        """
        Determine the symmetrical perturbations. Meaning of pertsy:

           0 for non-target perturbations
           1 for basis perturbations
          -1 for perturbations that can be found from basis perturbations
        """
        # Fortran array: nctkarr_t("pertsy_qpt", "int", "three, mpert, nqpt")))
        return self.reader.read_value("pertsy_qpt")

    def close(self):
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app("symv1scf: %s" % self.symv1scf)

        return "\n".join(lines)

    @lazy_property
    def qpoints(self):
        return KpointList(self.structure.reciprocal_lattice, frac_coords=self.reader.read_value("qpts"))

    def _find_iqpt_qpoint(self, qpoint):
        if duck.is_intlike(qpoint):
            iq = qpoint
            qpoint = self.qpoints[iq]
        else:
            qpoint = Kpoint.as_kpoint(qpoint, self.structure.reciprocal_lattice)
            iq = self.qpoints.index(qpoint)

        return iq, qpoint

    def read_v1_at_iq(self, key, iq, reshape_nfft_nspden=False):
        # Fortran array ("two, nfft, nspden, natom3, nqpt")
        v1 = self.reader.read_variable(key)[iq]
        v1 = v1[..., 0] + 1j * v1[..., 1]  
        # reshape (nspden, nfft) dims because we are not interested in the spin dependence.
        if reshape_nfft_nspden: v1 = np.reshape(v1, (self.natom3, self.nspden * self.nfft))
        return v1

    @add_fig_kwargs
    def plot_diff_at_qpoint(self, qpoint=0, fontsize=8, **kwargs):
        """
        Args:
            qpoint:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        iq, qpoint = self._find_iqpt_qpoint(qpoint)

        # complex arrays with shape: (natom3, nspden * nfft)
        origin_v1 = self.read_v1_at_iq("origin_v1scf", iq, reshape_nfft_nspden=True)
        symm_v1 = self.read_v1_at_iq("recons_v1scf", iq, reshape_nfft_nspden=True)

        num_plots, ncols, nrows = self.natom3, 3, self.natom3 // 3
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)

        for nu, ax in enumerate(ax_list.ravel()):
            idir = nu % 3
            ipert = (nu - idir) // 3

            # l1_rerr(f1, f2) = \int |f1 - f2| dr / (\int |f2| dr
            abs_diff = np.abs(origin_v1[nu] - symm_v1[nu])
            l1_rerr = np.sum(abs_diff) / np.sum(np.abs(origin_v1[nu]))

            stats = OrderedDict([
                ("max", abs_diff.max()),
                ("min", abs_diff.min()),
                ("mean", abs_diff.mean()),
                ("std", abs_diff.std()),
                ("L1_rerr", l1_rerr),
            ])

            xs = np.arange(len(abs_diff))
            ax.hist(abs_diff, facecolor='g', alpha=0.75)
            ax.grid(True)
            ax.set_title("idir: %d, iat: %d, pertsy: %d" % (idir, ipert, self.pertsy_qpt[iq, ipert, idir]), 
                         fontsize=fontsize)

            ax.axvline(stats["mean"], color='k', linestyle='dashed', linewidth=1)
            _, max_ = ax.get_ylim()
            ax.text(0.7, 0.7, "\n".join("%s = %.1E" % item for item in stats.items()), 
                    fontsize=fontsize, horizontalalignment='center', verticalalignment='center', 
                    transform=ax.transAxes)

        fig.suptitle("qpoint: %s" % repr(qpoint))
        return fig

    @add_fig_kwargs
    def plot_pots_at_qpoint(self, qpoint=0, fontsize=8, **kwargs):
        """
        Args:
            qpoint:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        iq, qpoint = self._find_iqpt_qpoint(qpoint)

        # complex arrays with shape: (natom3, nspden * nfft)
        origin_v1 = self.read_v1_at_iq("origin_v1scf", iq, reshape_nfft_nspden=True)
        symm_v1 = self.read_v1_at_iq("recons_v1scf", iq, reshape_nfft_nspden=True)

        num_plots, ncols, nrows = self.natom3, 3, self.natom3 // 3
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)

        natom = len(self.structure)
        xs = np.arange(self.nspden * self.nfft)
        for nu, ax in enumerate(ax_list.ravel()):
            idir = nu % 3
            ipert = (nu - idir) // 3

            # l1_rerr(f1, f2) = \int |f1 - f2| dr / (\int |f2| dr
            abs_diff = np.abs(origin_v1[nu] - symm_v1[nu])
            l1_rerr = np.sum(abs_diff) / np.sum(np.abs(origin_v1[nu]))

            stats = OrderedDict([
                ("max", abs_diff.max()),
                ("min", abs_diff.min()),
                ("mean", abs_diff.mean()),
                ("std", abs_diff.std()),
                ("L1_rerr", l1_rerr),
            ])

            ax.grid(True)
            ax.set_title("idir: %d, iat: %d, pertsy: %d" % (idir, ipert, self.pertsy_qpt[iq, ipert, idir]), 
                         fontsize=fontsize)
            # Plot absolute error
            #ax.plot(xs, abs_diff, linestyle="-", color="red", alpha=1.0, label="Abs diff" if nu == 0 else None)

            # Plot absolute values
            #ax.plot(xs, np.abs(origin_v1[nu]), linestyle="--", color="red", alpha=0.4, label="Origin" if nu == 0 else None)
            #ax.plot(xs, -np.abs(symm_v1[nu]), linestyle="--", color="blue", alpha=0.4, label="-Symm" if nu == 0 else None)

            # Plot real and imag
            #ax.plot(xs, origin_v1[nu].real, linestyle="--", color="red", alpha=0.4, label="Re Origin" if nu == 0 else None)
            #ax.plot(xs, -symm_v1[nu].real, linestyle="--", color="blue", alpha=0.4, label="Re Symm" if nu == 0 else None)

            data = np.angle(origin_v1[nu], deg=True) - np.angle(symm_v1[nu], deg=True)
            #data = data[abs_diff > stats["mean"]]
            data = data[np.abs(origin_v1[nu]) > 1e-5]
            ax.plot(np.arange(len(data)), data, 
                    linestyle="--", color="red", alpha=0.4, label="diff angle degrees" if nu == 0 else None)

            #ax.plot(xs, origin_v1[nu].real, linestyle="--", color="red", alpha=0.4, label="Re Origin" if nu == 0 else None)
            #ax.plot(xs, -symm_v1[nu].real, linestyle="--", color="blue", alpha=0.4, label="Re Symm" if nu == 0 else None)

            #ax.plot(xs, origin_v1[nu].real - symm_v1[nu].real, linestyle="--", color="red", alpha=0.4, 
            #        label="Re Origin" if nu == 0 else None)

            #ax.plot(xs, origin_v1[nu].imag, linestyle=":", color="red", alpha=0.4, label="Imag Origin" if nu == 0 else None)
            #ax.plot(xs, -symm_v1[nu].imag, linestyle=":", color="blue", alpha=0.4, label="Imag Symm" if nu == 0 else None)

            #ax.plot(xs, origin_v1[nu].imag - symm_v1[nu].imag, linestyle="--", color="blue", alpha=0.4, 
            #        label="Re Origin" if nu == 0 else None)

            if nu == 0: 
                ax.set_ylabel(r"Abs diff")
                ax.legend(loc="best", fontsize=fontsize, shadow=True)
            if ipert == natom -1: 
                ax.set_xlabel(r"FFT index")

            #ax.axvline(stats["mean"], color='k', linestyle='dashed', linewidth=1)
            _, max_ = ax.get_ylim()
            ax.text(0.7, 0.7, "\n".join("%s = %.1E" % item for item in stats.items()), 
                    fontsize=fontsize, horizontalalignment='center', verticalalignment='center', 
                    transform=ax.transAxes)

            #ax2 = ax.twinx()
            #rerr = 100 * abs_diff / np.abs(origin_v1[nu]) 
            #ax2.plot(xs, rerr, linestyle="--", color="blue", alpha=0.4, 
            #          label=r"|V_{\mathrm{origin}}|" if nu == 0 else None)

        fig.suptitle("qpoint: %s" % repr(qpoint))
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        maxnq = 3
        for iq, qpoint in enumerate(self.qpoints):
            if iq > maxnq:
                print("Only the first %d q-points are show..." % maxnq)
                break
            #yield self.plot_diff_at_qpoint(qpoint=iq, **kwargs, show=False)
            yield self.plot_pots_at_qpoint(qpoint=iq, **kwargs, show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
        ])

        for iq, qpoint in enumerate(self.qpoints):
            nb.cells.append(nbv.new_code_cell("ncfile.plot_diff_at_qpoint(qpoint=%d);" % iq))
            #nb.cells.append(nbv.new_code_cell("ncfile.plot_diff_at_qpoint(qpoint=%d);" % iq))

        return self._write_nb_nbpath(nb, nbpath)

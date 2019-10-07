# coding: utf-8
"""
Object to plot DFPT potentials in the phonon mode representation.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict
from monty.string import marquee
from monty.functools import lazy_property
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.core.kpoints import KpointList, Kpoint
from abipy.tools import duck
from abipy.iotools import Visualizer, xsf, ETSF_Reader #, cube


class V1qnuFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super(V1qnuFile, self).__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        # Read dimensions.
        self.nfft = r.read_dimvalue("nfft")
        self.nspden = r.read_dimvalue("nspden")
        self.natom3 = len(self.structure) * 3
        # Read FFT grid.
        self.ngfft = r.read_value("ngfft")

    @lazy_property
    def structure(self):
        """|Structure| object."""
        return self.reader.read_structure()

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
        for iq, qpt in enumerate(self.qpoints):
            app("[%d] %s" % (iq, repr(qpt)))

        return "\n".join(lines)

    @lazy_property
    def qpoints(self):
        return KpointList(self.structure.reciprocal_lattice, frac_coords=self.reader.read_value("qlist"))

    def _find_iqpt_qpoint(self, qpoint):
        if duck.is_intlike(qpoint):
            iq = qpoint
            qpoint = self.qpoints[iq]
        else:
            qpoint = Kpoint.as_kpoint(qpoint, self.structure.reciprocal_lattice)
            iq = self.qpoints.index(qpoint)

        return iq, qpoint

    def visualize_qpoint_nu(self, qpoint, nu, spin=0, appname="vesta"):
        iq, qpoint = self._find_iqpt_qpoint(qpoint)

        # Fortran array nctkarr_t("v1_qnu", "dp", "two, nfft, nspden, natom3, nqlist")])
        v1_qnu = self.reader.read_variable("v1_qnu")[iq, nu, spin]
        v1_qnu = v1_qnu[:, 0] + 1j * v1_qnu[:, 1]  
        #wqnu = self.reader.read_variable["phfreqs"][nu]
        #v1_qnu /= np.sqrt(2 * wqnu)
        datar = np.reshape(np.abs(v1_qnu), self.ngfft)

        visu = Visualizer.from_name(appname)
        ext = "xsf"
        if ext not in visu.supported_extensions():
            raise ValueError("Visualizer %s does not support XSF files" % visu)
        from abipy.core.globals import abinb_mkstemp
        _, filename = abinb_mkstemp(suffix="." + ext, text=True)

        with open(filename, mode="wt") as fh:
            if ext == "xsf":
                xsf.xsf_write_structure(fh, self.structure)
                xsf.xsf_write_data(fh, self.structure, datar, add_replicas=True)
            else:
                raise NotImplementedError("extension %s is not supported." % ext)
            
        return visu(filename)

    @add_fig_kwargs
    def plot_v1qnu_vs_lr(self, ax=None, fontsize=8, **kwargs):
        """
        Plot the difference between the ab-initio v1_qnu and the potential obtained with Verdi's model

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Label and title fontsize.

        Return: |matplotlib-Figure|
        """
        # Fortran array nctkarr_t("v1_qnu", "dp", "two, nfft, nspden, natom3, nqlist")])
        v1_qnu = self.reader.read_value("v1_qnu", cmode="c")
        v1lr_qnu = self.reader.read_value("v1lr_qnu", cmode="c")

        stats = OrderedDict([
            ("min", []),
            ("max", []),
            ("mean", []),
            ("std", []),
        ])

        qnorms = []
        for iq, qpt in enumerate(self.qpoints):
            qnorms.append(qpt.norm)
            abs_diff = np.abs(v1_qnu[iq] - v1lr_qnu[iq])
            for key in stats.keys():
                stats[key].append(getattr(abs_diff, key)())

        # Sort values by |q|.
        qindex = np.arange(len(qnorms))
        items = sorted(list(zip(qnorms, qindex)),  key=lambda t: t[0])
        qnorm = np.array([t[0] for t in items])
        qindex = np.array([t[1] for t in items])
        #print(qnorm, "\n", qindex)
        for key in stats.keys():
            stats[key] = np.array(stats[key])[qindex]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for key, values in stats.items():
            ax.plot(values, label=key)

        ax.grid(True)
        ax.set_xlabel(r"$|\bf{q}|$ 1/Ang")
        #ax.set_ylabel(r"$|g_{\bf q}|$ (meV)")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        #title = "band_kq: %s, band_k: %s, kpoint: %s" % (band_kq, band_k, repr(kpoint))
        #ax.set_title(title, fontsize=fontsize)
        
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot_v1qnu_vs_lr(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("# ncfile.visualize_qpoint_nu(qpoint, nu, spin=0, appname='vesta'")
        ])

        return self._write_nb_nbpath(nb, nbpath)

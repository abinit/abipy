# coding: utf-8
"""
Object to plot the lattice representation of the DFPT potentials (phonon perturbations).
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
        #self.has_dielt_zeff = bool(r.read_value("has_dielt_zeff"))
        #self.add_lr_part = bool(r.read_value("add_lr_part"))

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
        #NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("pertsy_qpt", "int", "three, mpert, nqpt")))
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
        #app("has_dielt_zeff: %s" % self.has_dielt_zeff)
        #app("add_lr_part: %s" % self.add_lr_part)
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

        # Fortran array "origin_v1scf", "dp", ("two, nfft, nspden, natom3, nqpt")
        origin_v1 = self.reader.read_variable("origin_v1scf")[iq]
        origin_v1 = origin_v1[..., 0] + 1j * origin_v1[..., 1]  
        # reshape nspden * nfft because we are not interested in the spin dependence.
        origin_v1 = np.reshape(origin_v1, (self.natom3, self.nspden * self.nfft))

        symm_v1 = self.reader.read_variable("recons_v1scf")[iq]
        symm_v1 = symm_v1[..., 0] + 1j * symm_v1[..., 1]  
        symm_v1 = np.reshape(symm_v1, (self.natom3, self.nspden * self.nfft))

        num_plots, ncols, nrows = self.natom3, 3, self.natom3 // 3
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)

        for nu, ax in enumerate(ax_list.ravel()):
            idir = nu % 3
            ipert = (nu - idir) // 3

            abs_diff = np.abs(origin_v1[nu] - symm_v1[nu]).ravel()
            #rel_diff = abs_diff / (np.abs(origin_v1[nu]) + 1e-10)

            stats = OrderedDict([
                ("max", abs_diff.max()),
                ("min", abs_diff.min()),
                ("mean", abs_diff.mean()),
                ("std", abs_diff.std()),
            ])
            #stats_title = ", ".join("%s = %.1E" % item for item in stats.items())

            xs = np.arange(len(abs_diff))
            ax.hist(abs_diff, facecolor='g', alpha=0.75)
            ax.grid(True)
            ax.set_title("idir: %d, iat: %d, pertsy: %d" % (idir, ipert, self.pertsy_qpt[iq, ipert, idir]), 
                         fontsize=fontsize)

            ax.axvline(stats["mean"], color='k', linestyle='dashed', linewidth=1)
            _, max_ = ax.get_ylim()
            ax.text(0.7, 0.7, 
                    "\n".join("%s = %.1E" % item for item in stats.items()), 
                    fontsize=fontsize, horizontalalignment='center', verticalalignment='center', 
                    transform=ax.transAxes)

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
            yield self.plot_diff_at_qpoint(qpoint=iq, **kwargs, show=False)

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

        return self._write_nb_nbpath(nb, nbpath)

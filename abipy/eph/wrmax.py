# coding: utf-8
"""
Object to plot the lattice representation of the DFPT potentials (phonon perturbations).
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from monty.string import marquee
from monty.functools import lazy_property
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, ax_share)
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.iotools import ETSF_Reader


class WRmaxFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super(WRmaxFile, self).__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        self.nrpt = r.read_dimvalue("nrpt")

        self.has_dielt_zeff = bool(r.read_value("has_dielt_zeff"))
        self.add_lr_part = bool(r.read_value("add_lr_part"))
        self.symv1 = bool(r.read_value("symv1"))

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
        return self.to_string(with_doc=False)

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app("ngqpt: %s" % str(self.ngqpt))
        app("has_dielt_zeff: %s" % self.has_dielt_zeff)
        app("add_lr_part: %s" % self.add_lr_part)
        app("symv1: %s" % self.symv1)

        return "\n".join(lines)

    @lazy_property
    def ngqpt(self):
        """Division of coarse q-mesh (ab-initio points)."""
        return self.reader.read_value("ngqpt")

    @lazy_property
    def rmod(self):
        """|R| in Bohr"""
        return self.reader.read_value("rmod")

    @add_fig_kwargs
    def plot(self, fontsize=12, **kwargs):
        """
        Args:
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        # Build grid of plots.
        natom = len(self.structure)
        ncols, nrows = (2, natom //2) if natom % 2 == 0 else (1, natom)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # nctkarr_t("maxw", "dp", "nrpt, natom3")
        maxw = np.reshape(self.reader.read_value("maxw"), (natom, 3, self.nrpt))
        for iatom, ax in enumerate(ax_list.ravel()):
            site = self.structure[iatom]
            title = "{} [{:.4f} {:.4f} {:.4f}]".format(site.specie.symbol, *site.frac_coords)
            ax.set_title(title, fontsize=fontsize)
            ax.semilogy(self.rmod, maxw[iatom, 0], marker="o", ls=":", lw=0, label="$L_x$" if iatom == 0 else None)
            ax.semilogy(self.rmod, maxw[iatom, 1], marker="o", ls=":", lw=0, label="$L_y$" if iatom == 0 else None)
            ax.semilogy(self.rmod, maxw[iatom, 2], marker="o", ls=":", lw=0, label="$L_z$" if iatom == 0 else None)
            ax.grid(True)
            if iatom == 0: 
                ax.set_ylabel(r"$Max_{\bf{r}} \| W(\bf{r}, \bf{R}) \|$")
                ax.legend(loc="best", fontsize=fontsize, shadow=True)
            if iatom == len(ax_list) - 1: ax.set_xlabel(r"$\|\bf{R}\|$ (Bohr)")

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        #for fig in self.yield_structure_figs(**kwargs): yield fig
        yield self.plot(**kwargs)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("wrmax = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(wrmax)"),
            nbv.new_code_cell("wrmax.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

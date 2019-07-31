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
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader


class WRmaxFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super(WRmaxFile, self).__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        self.nrpt = r.read_dimvalue("nrpt")
        self.has_zeff = bool(r.read_value("has_zeff"))
        self.has_dielt = bool(r.read_value("has_dielt"))
        self.dvdb_add_lr = r.read_value("dvdb_add_lr")
        self.symv1 = r.read_value("symv1")
        #self.symv1 = -1
        self.qdamp = r.read_value("qdamp")

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
        app("has_dielt: %s" % self.has_dielt)
        app("has_zeff: %s" % self.has_zeff)
        app("dvdb_add_lr: %s" % self.dvdb_add_lr)
        app("symv1: %s" % self.symv1)

        return "\n".join(lines)

    @lazy_property
    def ngqpt(self):
        """Division of coarse q-mesh (ab-initio points)."""
        return self.reader.read_value("ngqpt")

    @lazy_property
    def rmod(self):
        """|R| in Bohr."""
        return self.reader.read_value("rmod")

    @add_fig_kwargs
    def plot(self, scale="semilogy", ax=None, fontsize=8, **kwargs):
        """
        Plot the decay of max_{r,idir,ipert} |W(R,r,idir,ipert)|

        Args:
            scale: "semilogy", "loglog" or "plot".
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # nctkarr_t("maxw", "dp", "nrpt, natom3")
        maxw = self.reader.read_value("maxw")
        data = np.max(maxw, axis=0)
        f = {"plot": ax.plot, "semilogy": ax.semilogy, "loglog": ax.loglog}[scale]
        f(self.rmod, data, marker="o", ls=":", lw=0, **kwargs)
        ax.grid(True)
        ax.set_ylabel(r"$Max_{({\bf{r}}, idir, ipert)} \| W({\bf{r}}, {\bf{R}}, idir, ipert) \|$")
        ax.set_xlabel(r"$\|{\bf{R}}\|$ (Bohr)")

        if kwargs.pop("with_title", True):
            ax.set_title("dvdb_add_lr %d, qdamp: %s, symv1: %d" % (self.dvdb_add_lr, self.qdamp, self.symv1), 
                         fontsize=fontsize)
        return fig

    @add_fig_kwargs
    def plot_perts(self, scale="semilogy", fontsize=8, **kwargs):
        """
        Plot the decay of max_r |W(R,r,idir,ipert)| for the individual atomic perturbations

        Args:
            scale: "semilogy", "loglog" or "plot".
            fontsize: fontsize for legends and titles.

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
            f = {"plot": ax.plot, "semilogy": ax.semilogy, "loglog": ax.loglog}[scale]
            f(self.rmod, maxw[iatom, 0], marker="o", ls=":", lw=0, label="$L_x$" if iatom == 0 else None)
            f(self.rmod, maxw[iatom, 1], marker="o", ls=":", lw=0, label="$L_y$" if iatom == 0 else None)
            f(self.rmod, maxw[iatom, 2], marker="o", ls=":", lw=0, label="$L_z$" if iatom == 0 else None)
            ax.grid(True)
            if iatom == 0: 
                ax.set_ylabel(r"$Max_{{\bf{r}}} \| W({\bf{r}}, {\bf{R}}) \|$")
                ax.legend(loc="best", fontsize=fontsize, shadow=True)
            if iatom == len(ax_list) - 1: ax.set_xlabel(r"$\|{\bf{R}}\|$ (Bohr)")

        fig.suptitle("dvdb_add_lr %d, qdamp: %s, symv1: %d" % (self.dvdb_add_lr, self.qdamp, self.symv1), 
                     fontsize=fontsize)
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function generates a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)
        yield self.plot_perts(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.plot();"),
            nbv.new_code_cell("ncfile.plot_perts();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class WrmaxRobot(Robot):
    """
    This robot analyzes the results contained in multiple WRmax.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: WrmaxRobot
    """
    EXT = "Wrmax"

    @add_fig_kwargs
    def plot(self, ax=None, **kwargs):
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for label, abifile in self.items():
            abifile.plot(ax=ax, label=label, with_title=False, show=False, **kwargs)
        return fig

    def yield_figs(self, **kwargs): # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.plot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to `nbpath`. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.WrmaxRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

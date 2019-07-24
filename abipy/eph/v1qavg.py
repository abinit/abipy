# coding: utf-8
"""

"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict
from monty.string import marquee
from monty.functools import lazy_property
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, ax_share)
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.core.kpoints import Kpath
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader


class V1qAvgFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super(V1qAvgFile, self).__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        #self.has_zeff = bool(r.read_value("has_zeff"))
        #self.has_dielt = bool(r.read_value("has_dielt"))
        #self.dvdb_add_lr = r.read_value("dvdb_add_lr")
        #self.symv1 = bool(r.read_value("symv1"))

    @lazy_property
    def structure(self):
        """|Structure| object."""
        return self.reader.read_structure()

    @lazy_property
    def qpoints(self):
        """List of Q-points."""
        frac_coords = self.reader.read_value('qpoints')
        return Kpath(self.structure.reciprocal_lattice, frac_coords, ksampling=None)

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
        app(self.qpoints.to_string(verbose=verbose, title="Q-path"))
        #app("has_dielt: %s" % self.has_dielt)
        #app("has_zeff: %s" % self.has_zeff)
        #app("dvdb_add_lr: %s" % self.dvdb_add_lr)
        #app("symv1: %s" % self.symv1)

        return "\n".join(lines)

    def _make_ticks_and_labels(self):
        # Find the k-point names in the pymatgen database.

        od = OrderedDict()
        # If the first or the last k-point are not recognized in findname_in_hsym_stars
        # matplotlib won't show the full band structure along the k-path
        # because the labels are not defined. So we have to make sure that
        # the labels for the extrema of the path are always defined.
        od[0] = " "

        for idx, kpoint in enumerate(self.qpoints):
            name = kpoint.name if kpoint.name is not None else self.structure.findname_in_hsym_stars(kpoint)
            if name:
                od[idx] = name
                if kpoint.name is None: kpoint.set_name(name)

        last = len(self.qpoints) - 1
        if last not in od: od[last] = " "

        return list(od.keys()), list(od.values())

    @add_fig_kwargs
    def plot(self, fontsize=8, ispden=0, **kwargs):
        """
        Plot

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        # nctkarr_t("v1scf_avg", "dp", "two, nspden, three, natom, nqpt")
        # nctkarr_t("v1lr_avg",  "dp", "two, nspden, three, natom, nqpt")
        v1scf_avg = self.reader.read_value("v1scf_avg")
        v1lr_avg = self.reader.read_value("v1lr_avg")
        #v1scf_abs_avg = self.reader.read_value("v1scf_abs_avg")

        natom = len(self.structure)
        nrows, ncols = natom, 3
        num_plots = ncols * nrows

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)

        xs = np.arange(len(self.qpoints))
        ticks, labels = self._make_ticks_and_labels()

        for ip, ax in enumerate(ax_mat.ravel()):
            idir = ip % 3
            iat = (ip - idir) // 3

            #ax.plot(xs, v1scf_abs_avg[:, iat, idir, ispden, 0], label="Re v1scf_abs_avg")
            #ax.plot(xs, v1scf_abs_avg[:, iat, idir, ispden, 1], label="Im v1scf_abs_avg")

            ax.plot(xs, v1scf_avg[:, iat, idir, ispden, 0], label=r"$\Re \langle v1q \rangle$")
            ax.plot(xs, v1scf_avg[:, iat, idir, ispden, 1], label=r"$\Im \langle v1q \rangle$")
            ax.plot(xs, v1lr_avg[:, iat, idir, ispden, 0], label=r"$\Re \langle v1lr \rangle$")
            ax.plot(xs, v1lr_avg[:, iat, idir, ispden, 1], label=r"$\Im \langle v1lr \rangle$")

            ax.grid(True)
            if iat == natom - 1: ax.set_xlabel("Wave vector")
            #if idir == 0: ax.set_ylabel(r"$\|{\bf{R}}\|$ (Bohr)")

            if ticks:
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.pop("qlabel_size", "large"))

            if ip == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

            site = self.structure[iat]
            ax.set_title("idir: %d, iat: %d, %s" % (idir, iat, site.specie.symbol), fontsize=fontsize)

        #if kwargs.pop("with_title", True):
        #    ax.set_title("dvdb_add_lr %d, alpha_gmin: %s, symv1: %d" % (self.dvdb_add_lr, self.alpha_gmin, self.symv1), 
        #                 fontsize=fontsize)

        #fig.suptitle("qpoint: %s" % repr(self.qpoint), fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function generates a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)

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
        ])

        return self._write_nb_nbpath(nb, nbpath)


class V1qAvgRobot(Robot):
    """
    This robot analyzes the results contained in multiple V1qAvgFile.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: V1qAvgRobot
    """
    EXT = "V1QAVG"

    #@add_fig_kwargs
    #def plot(self, ax=None, **kwargs):
    #    ax, fig, plt = get_ax_fig_plt(ax=ax)
    #    for label, abifile in self.items():
    #        abifile.plot(ax=ax, label=label, with_title=False, show=False, **kwargs)
    #    return fig

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
            nbv.new_code_cell("robot = abilab.V1qAvgRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

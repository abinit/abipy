# coding: utf-8
"""

"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict
from monty.string import list_strings, marquee
from monty.functools import lazy_property
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, linestyles
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.core.kpoints import Kpath
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader


def _get_style(reim, what):
    symb, linestyle, linewidth, alpha = {
        "v1scf_avg": (r"v1_{\bf q}", "-", 1, 2.0),
        "v1lr_avg": (r"v1_{\mathrm{lr}}", "--", 1, 2.0),
        "v1scfmlr_avg": (r"(v1_{\bf q} - v1_{\mathrm{lr}})", "-.", 1, 2.0),
    }[what]

    return dict(
        #marker="x", markersize=1, 
        label=r"$\langle \%s %s \rangle$" % ({0: "Re", 1: "Im"}[reim], symb), 
        color={0: "blue", 1: "red"}[reim],
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=alpha,
    )


class V1qAvgFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super(V1qAvgFile, self).__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        #self.has_zeff = bool(r.read_value("has_zeff"))
        #self.has_dielt = bool(r.read_value("has_dielt"))
        #self.has_quadrupoles = bool(r.read_value("has_quadrupoles"))
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
        return self.to_string()

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

    def make_ticks_and_labels(self):
        # Find the k-point names in the pymatgen database.

        od = OrderedDict()
        # If the first or the last k-point are not recognized in findname_in_hsym_stars
        # matplotlib won't show the full band structure along the k-path
        # because the labels are not defined. So we have to make sure that
        # the labels for the extrema of the path are always defined.
        od[0] = " "

        for idx, qpoint in enumerate(self.qpoints):
            name = qpoint.name if qpoint.name is not None else self.structure.findname_in_hsym_stars(qpoint)
            if name:
                od[idx] = name
                if qpoint.name is None: qpoint.set_name(name)

        last = len(self.qpoints) - 1
        if last not in od: od[last] = " "

        return list(od.keys()), list(od.values())

    @add_fig_kwargs
    def plot(self, what_list="all", ispden=0, fontsize=8, sharey=False, **kwargs):
        """
        Plot

        Args:
            what_list:
            ispden: Spin density component to plot.
            fontsize: fontsize for legends and titles
            sharey: True to share y-axes.

        Return: |matplotlib-Figure|
        """
        #all_varnames = ["v1scf_avg", "v1lr_avg", "v1scfmlr_avg", "v1scf_abs_avg"]
        what_list = list_strings(what_list) if what_list != "all" else ["v1scf_avg", "v1lr_avg"]
        data = {}
        for vname in what_list:
            # nctkarr_t("v1scf_avg", "dp", "two, nspden, three, natom, nqpt")
            data[vname] = self.reader.read_value(vname)

        natom = len(self.structure)
        nrows, ncols = natom, 3
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)

        xs = np.arange(len(self.qpoints))
        ticks, labels = self.make_ticks_and_labels()

        for ip, ax in enumerate(ax_mat.ravel()):
            idir = ip % 3
            iat = (ip - idir) // 3
            for reim in (0, 1):
                for vname in what_list:
                    ys = data[vname][:, iat, idir, ispden, reim]
                    ax.plot(xs, ys, **_get_style(reim, vname))

            ax.grid(True)
            if iat == natom - 1: ax.set_xlabel("Q Wave Vector")
            if idir == 0: ax.set_ylabel(r"Hartree/Bohr")

            if ticks:
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.get("qlabel_size", "large"))
                ax.set_xlim(ticks[0], ticks[-1])

            if ip == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

            site = self.structure[iat]
            s = "%s [%.3f, %.3f, %.3f]" % (site.specie.symbol, site.frac_coords[0], site.frac_coords[1], site.frac_coords[2])
            ax.set_title("idir: %d, iat: %d, %s" % (idir, iat, s), fontsize=fontsize)

        #if kwargs.pop("with_title", True):
        #    ax.set_title("dvdb_add_lr %d, alpha_gmin: %s, symv1: %d" % (self.dvdb_add_lr, self.alpha_gmin, self.symv1), 
        #                 fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function generates a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(what_list="v1scfmlr_avg", title=r"$v1_{\bf q} - v1_{\bf q}^{\mathrm{LR}}$", show=False)
        yield self.plot(title=r"$v1_{\bf q}\,vs\,v1_{\bfq}^{\mathrm{LR}}$", show=False)
        #import os
        #yield self.plot(title=os.path.basename(self.filepath), show=False)

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

    @lazy_property
    def qpoints(self):
        """List of Q-points."""
        if len(self) == 1: return self.abifiles[0].qpoints

        if (any(len(ncfile.qpoints) != len(self.abifiles[0].qpoints) for ncfile in self.abifiles)):
            raise RuntimeError("Assuming ncfiles with same number of q-points.")

        for abifile in self.abifiles[1:]:
            if np.any(np.abs(abifile.qpoints.frac_coords - self.abifiles[0].qpoints.frac_coords) > 1e-3):
                raise RuntimeError("Found different q-points!")

        return self.abifiles[0].qpoints

    @add_fig_kwargs
    def plot(self, ispden=0, sharey=False, fontsize=8, **kwargs):
        """
        Plot

        Args:
            ispden: Spin density component to plot.
            sharey: True to share y-axes.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        # Caveat: No check is done on the consistency among structures
        ref_file = self.abifiles[0]
        structure = ref_file.structure

        natom = len(structure)
        nrows, ncols = natom, 3
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)

        xs = np.arange(len(self.qpoints))
        ticks, labels = ref_file.make_ticks_and_labels()

        vname = "v1scf_avg"
        data_file = {abilabel: abifile.reader.read_value(vname) for abilabel, abifile in self.items()}
        style = dict(marker=".", markersize=2)

        for ip, ax in enumerate(ax_mat.ravel()):
            idir = ip % 3
            iat = (ip - idir) // 3

            for ifile, (abilabel, abifile) in enumerate(self.items()):
                v1scf_avg = data_file[abilabel]
                for reim in (0, 1):
                    style = _get_style(reim, vname)
                    label = "%s %s" % (style.pop("label"), abilabel)
                    style["linestyle"] = {0: "-", 1: "--", 2: "-.", 3: ":"}[ifile]
                    ax.plot(xs, v1scf_avg[:, iat, idir, ispden, reim], label=label if ip == 0 else None, **style)

            if ip == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

            ax.grid(True)
            if iat == natom - 1: ax.set_xlabel("Q Wave Vector")
            if idir == 0: ax.set_ylabel(r"Hartree/Bohr")

            site = ref_file.structure[iat]
            s = "%s [%.3f, %.3f, %.3f]" % (site.specie.symbol, site.frac_coords[0], site.frac_coords[1], site.frac_coords[2])
            ax.set_title("idir: %d, iat: %d, %s" % (idir, iat, s), fontsize=fontsize)
                                                                                                                 
            if ticks:
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.get("qlabel_size", "large"))
                ax.set_xlim(ticks[0], ticks[-1])
                
        return fig

    def yield_figs(self, **kwargs): # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
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

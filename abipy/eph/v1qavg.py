# coding: utf-8
"""
Tools to analyze the V1QAVG file produced by the E-PH code (eph_task +15 or -15)
"""
import numpy as np

from collections import OrderedDict
from monty.string import list_strings, marquee
from monty.functools import lazy_property
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.core.kpoints import Kpath
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader


def _get_style(reim, what, marker=None, markersize=None, alpha=2.0):
    lw = 1
    symbol, linestyle, linewidth = {
        "v1scf_avg": (r"v1_{\bf q}", "-", lw),
        "v1scf_abs_avg": (r"|v1_{\bf q}|", "-", lw),
        "v1lr_avg": (r"v1_{\mathrm{lr}}", "--", lw),
        "v1lr_abs_avg": (r"|v1_{\mathrm{lr}}|", "--", lw),
        "v1scfmlr_avg": (r"(v1_{\bf q} - v1_{\mathrm{lr}})", "-.", lw),
        "v1scfmlr_abs_avg": (r"(|v1_{\bf q} - v1_{\mathrm{lr}|})", "-.", lw),
        "v1scf_gsmall": (r"v1_{\bf q}(G)", "-", lw),
        "v1lr_gsmall": (r"v1lr_{\bf q}(G)", "-.", lw),
        "v1scfmlr_gsmall": (r"(v1_{\bf q}(G) - v1_{\mathrm{lr}}(G))", "-.", lw),
    }[what]

    return dict(
        marker=marker,
        markersize=markersize,
        label=r"$\langle \%s %s \rangle$" % ({0: "Re", 1: "Im"}[reim], symbol),
        color={0: "blue", 1: "red"}[reim],
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=alpha,
    )


class V1qAvgFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    The V1QAVG.nc file contains the average over the unit cell of the periodic part the DFPT scattering potential.
    This file is produced by the E-PH code by setting eph_task to +15 or -15.
    If eph_task is +15, the input DVDB contains a q-mesh and the potentials are interpolated on a list of q-points
    (usually a q-path) specified by the user. In this case the V1QAVG.nc file also contains an extra array
    with Max_r |W(R, r)|, useful to study the decay of the scattering potentials in R-space.
    If eph_task is -15, the netcdf file contains the average for the q-points found in the DVDB file.
    This option is usually used to visualize the ab-initio potentials and compare then with the model for the LR part.
    """

    def __init__(self, filepath):
        super().__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        # Read medadata
        self.has_zeff = bool(r.read_value("has_zeff"))
        self.has_dielt = bool(r.read_value("has_dielt"))
        self.has_quadrupoles = bool(r.read_value("has_quadrupoles"))
        self.has_efield = bool(r.read_value("has_efield", default=False))
        self.dvdb_add_lr = r.read_value("dvdb_add_lr")
        self.symv1scf = r.read_value("symv1scf")
        self.interpolated = bool(r.read_value("interpolated"))
        self.qdamp = r.read_value("qdamp")

    @lazy_property
    def structure(self):
        """|Structure| object."""
        return self.reader.read_structure()

    @lazy_property
    def qpoints(self):
        """List of Q-points."""
        frac_coords = self.reader.read_value('qpoints')
        return Kpath(self.structure.reciprocal_lattice, frac_coords, ksampling=None)

    @lazy_property
    def has_maxw(self):
        """True if ncfile contains Max_r |W(R, r)|"""
        return "maxw" in self.reader.rootgrp.variables

    def close(self):
        self.reader.close()

    @lazy_property
    def params(self):
        """Dict with parameters that might be subject to convergence studies."""
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
        app("")
        app("has_dielt: %s, has_zeff: %s, has_quadrupoles: %s, has_efield: %s" % (
            self.has_dielt, self.has_zeff, self.has_quadrupoles, self.has_efield))
        app("dvdb_add_lr: %s, symv1scf: %s, interpolated: %s, qdamp: %s" % (
            self.dvdb_add_lr, self.symv1scf, self.interpolated, self.qdamp))

        return "\n".join(lines)

    def make_ticks_and_labels(self):
        """Find the k-point names in the pymatgen database."""

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

    #def xsf_write(self, idir=0, ipert=0, ispden=0):
    #    natom = len(self.structure)
    #    ip = idir + 3 * ipert
    #    r = self.reader
    #    ngfft = r.read_value("ngfft")

    #    from abipy.tools.numtools import transpose_last3dims
    #    from abipy.iotools import xsf

    #    def write(datar, filename):
    #        with open(filename, mode="wt") as fh:
    #            xsf.xsf_write_structure(fh, self.structure)
    #            xsf.xsf_write_data(fh, self.structure, datar, add_replicas=True, cplx_mode="abs")

    #    # nctkarr_t("v1r_interpolated", "dp", "two, nfft, nspden, natom3")
    #    var = r.read_variable("v1r_interpolated")
    #    datar = var[ip, ispden, :, 0] + 1j * var[ip, ispden, :, 1]
    #    # Because we are reading Fortran data (z,y,x) and xsf_write_data expects (..., x, y, z)
    #    #print(datar.shape)
    #    datar = np.reshape(datar, np.flip(ngfft))
    #    datar = transpose_last3dims(datar)
    #    write(datar, "v1r_interp_idir%d_ipert%d.xsf" % (idir, ipert))

    #    var = r.read_variable("v1r_lrmodel")
    #    datar = var[ip, ispden, :, 0] + 1j * var[ip, ispden, :, 1]
    #    datar = np.reshape(datar, np.flip(ngfft))
    #    datar = transpose_last3dims(datar)
    #    write(datar, "v1r_lr_idir%d_ipert%d.xsf" % (idir, ipert))

    @add_fig_kwargs
    def plot(self, what_list="all", ispden=0, fontsize=6, sharey=False, **kwargs):
        """
        Plot

        Args:
            what_list:
            ispden: Spin density component to plot.
            fontsize: fontsize for legends and titles
            sharey: True to share y-axes.

        Return: |matplotlib-Figure|
        """
        #all_varnames = ["v1scf_avg", "v1lr_avg", "v1scfmlr_avg", "v1scfmlr_abs_avg", "v1scf_abs_avg"]
        what_list = list_strings(what_list) if what_list != "all" else ["v1scf_avg", "v1lr_avg"]
        data = {}
        for vname in what_list:
            # Fortran array: nctkarr_t("v1scf_avg", "dp", "two, nspden, three, natom, nqpt")
            data[vname] = self.reader.read_value(vname)

        # Build [natom, 3] grid plot.
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
            if idir == 0: ax.set_ylabel(r"(Hartree/Bohr)")

            if ticks:
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.get("qlabel_size", "large"))
                ax.set_xlim(ticks[0], ticks[-1])

            if ip == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

            site = self.structure[iat]
            s = "%s [%.3f, %.3f, %.3f]" % (site.specie.symbol, site.frac_coords[0], site.frac_coords[1], site.frac_coords[2])
            ax.set_title("idir: %d, iat: %d, %s" % (idir, iat, s), fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_gvec(self, gvec, ispden=0, fontsize=6, sharey=False, **kwargs):
        """
        Plot

        Args:
            gvec: G vector in reduced coordinates.
            ispden: Spin density component to plot.
            fontsize: fontsize for legends and titles
            sharey: True to share y-axes.

        Return: |matplotlib-Figure|
        """
        gsmall = self.reader.read_value("gsmall")
        for ig, g in enumerate(gsmall):
            if np.all(g == gvec): break
        else:
            raise RuntimeError("Cannot find gvec %s in gsmall array" % str(gvec))

        # nctkarr_t("v1scf_gsmall", "dp", "two, ngsmall, nspden, three, natom, nqpt"), &
        # nctkarr_t("v1lr_gsmall", "dp", "two, ngsmall, nspden, three, natom, nqpt"), &
        what_list = ["v1scf_gsmall", "v1lr_gsmall"]
        data = {}
        for vname in what_list:
            ncvar = self.reader.read_variable(vname)
            data[vname] = ncvar[..., ig, :]

        # Build [natom, 3] grid plot.
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

                # plot difference.
                ys = (data["v1scf_gsmall"][:, iat, idir, ispden, reim]
                     - data["v1lr_gsmall"][:, iat, idir, ispden, reim]) # * 10
                ax.plot(xs, ys, **_get_style(reim, "v1scfmlr_gsmall"))

            ax.grid(True)
            if iat == natom - 1: ax.set_xlabel("Q Wave Vector")
            if idir == 0: ax.set_ylabel(r"(Hartree/Bohr)")

            if ticks:
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.get("qlabel_size", "large"))
                ax.set_xlim(ticks[0], ticks[-1])

            if ip == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

            site = self.structure[iat]
            s = "%s [%.3f, %.3f, %.3f]" % (site.specie.symbol, site.frac_coords[0], site.frac_coords[1], site.frac_coords[2])
            ax.set_title("idir: %d, iat: %d, %s" % (idir, iat, s), fontsize=fontsize)

        fig.suptitle("G = %s" % str(gvec))

        return fig

    @add_fig_kwargs
    def plot_maxw(self, scale="semilogy", ax=None, fontsize=8, **kwargs):
        """
        Plot the decay of max_{r,idir,ipert} |W(R,r,idir,ipert)|

        Args:
            scale: "semilogy", "loglog" or "plot".
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        if not self.has_maxw: return None
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Fortran array: nctkarr_t("maxw", "dp", "nrpt, natom3")
        rmod = self.reader.read_value("rmod")
        maxw = self.reader.read_value("maxw")
        data = np.max(maxw, axis=0)
        f = {"plot": ax.plot, "semilogy": ax.semilogy, "loglog": ax.loglog}[scale]
        f(rmod, data, marker="o", ls=":", lw=0, **kwargs)

        ax.grid(True)
        ax.set_ylabel(r"$Max_{({\bf{r}}, idir, ipert)} \| W({\bf{r}}, {\bf{R}}, idir, ipert) \|$")
        ax.set_xlabel(r"$\|{\bf{R}}\|$ (Bohr)")

        #if kwargs.pop("with_title", True):
        #    ax.set_title("dvdb_add_lr %d, qdamp: %s, symv1scf: %d" % (self.dvdb_add_lr, self.qdamp, self.symv1scf),
        #                 fontsize=fontsize)
        return fig

    @add_fig_kwargs
    def plot_maxw_perts(self, scale="semilogy", sharey=False, fontsize=8, **kwargs):
        """
        Plot the decay of max_r |W(R,r,idir,ipert)| for the individual atomic perturbations.

        Args:
            scale: "semilogy", "loglog" or "plot".
            sharey: True is y-axes should be shared.
            fontsize: fontsize for legends and titles.

        Return: |matplotlib-Figure|
        """
        if not self.has_maxw: return None
        # Build grid of plots.
        natom = len(self.structure)
        ncols, nrows = (2, natom // 2) if natom % 2 == 0 else (1, natom)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()

        # Fortran array: nctkarr_t("maxw", "dp", "nrpt, natom3")
        nrpt = self.reader.read_dimvalue("nrpt")
        rmod = self.reader.read_value("rmod")
        maxw = np.reshape(self.reader.read_value("maxw"), (natom, 3, nrpt))

        for iatom, ax in enumerate(ax_list.ravel()):
            site = self.structure[iatom]
            title = "{} [{:.4f} {:.4f} {:.4f}]".format(site.specie.symbol, *site.frac_coords)
            ax.set_title(title, fontsize=fontsize)
            f = {"plot": ax.plot, "semilogy": ax.semilogy, "loglog": ax.loglog}[scale]
            f(rmod, maxw[iatom, 0], marker="o", ls=":", lw=0, label="$L_x$" if iatom == 0 else None)
            f(rmod, maxw[iatom, 1], marker="o", ls=":", lw=0, label="$L_y$" if iatom == 0 else None)
            f(rmod, maxw[iatom, 2], marker="o", ls=":", lw=0, label="$L_z$" if iatom == 0 else None)
            ax.grid(True)
            if iatom == 0:
                ax.set_ylabel(r"$Max_{{\bf{r}}} \| W({\bf{r}}, {\bf{R}}) \|$")
                ax.legend(loc="best", fontsize=fontsize, shadow=True)
            if iatom == len(ax_list) - 1: ax.set_xlabel(r"$\|{\bf{R}}\|$ (Bohr)")

        #fig.suptitle("dvdb_add_lr %d, qdamp: %s, symv1scf: %d" % (self.dvdb_add_lr, self.qdamp, self.symv1scf),
        #             fontsize=fontsize)
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function generates a predefined list of matplotlib figures with minimal input from the user.
        """
        title = r"$\langle v1_{\bf q} \rangle \,vs\, \langle v1_{\bfq}^{\mathrm{LR}} \rangle$"
        yield self.plot(title=title, show=False)
        yield self.plot(what_list="v1scfmlr_avg", title=r"$v1_{\bf q} - v1_{\bf q}^{\mathrm{LR}}$", show=False)
        #yield self.plot(what_list=["v1scf_abs_avg", "v1lr_abs_avg"], title=r"ABS", show=False)
        #yield self.plot(what_list="v1scfmlr_abs_avg", title=r"$|v1_{\bf q} - v1_{\bf q}^{\mathrm{LR}|}$", show=False)

        for gvec in [[0, 0, 0], [1, 0, 0], [0, 1, 1], [1, 1, 1], [2, 2, 2]]:
            yield self.plot_gvec(gvec, show=False)

        if self.has_maxw:
            if kwargs.get("verbose", 0) > 0:
                yield self.plot_maxw(show=False)
                yield self.plot_maxw_perts(show=False)
            else:
                print("Use verbose to print the decay of W(r,R)")

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.plot();"),
        ])

        if self.has_maxw:
            nbv.new_code_cell("ncfile.plot_maxw();"),
            nbv.new_code_cell("ncfile.plot_maxw_perts();"),

        return self._write_nb_nbpath(nb, nbpath)


class V1qAvgRobot(Robot):
    """
    This robot analyzes the results contained in multiple V1qAvgFile.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: V1qAvgRobot
    """
    EXT = "V1QAVG"

    # Absolute tolerance used to compare q-points.
    atol = 1e-2

    @lazy_property
    def qpoints(self):
        """List of q-points."""
        if len(self) == 1: return self.abifiles[0].qpoints

        if (any(len(ncfile.qpoints) != len(self.abifiles[0].qpoints) for ncfile in self.abifiles)):
            raise RuntimeError("Assuming ncfiles with same number of q-points.\nFound %s" % (
                str([len(ncfile.qpoints) for ncfile in self.abifiles])))

        for abifile in self.abifiles[1:]:
            if np.any(np.abs(abifile.qpoints.frac_coords - self.abifiles[0].qpoints.frac_coords) > self.atol):
                for q1, q2 in zip(self.abifiles[0].qpoints, abifile.qpoints):
                    print("q1:", q1, ", q2:", q2)
                raise RuntimeError("Found different q-points with tolerance: %s!" % self.atol)

        return self.abifiles[0].qpoints

    @add_fig_kwargs
    def plot(self, ispden=0, vname="v1scf_avg", sharey=False, fontsize=8, **kwargs):
        """
        Plot

        Args:
            ispden: Spin density component to plot.
            sharey: True to share y-axes.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        # Caveat: No check is done on the consistency among structures.
        ref_file = self.abifiles[0]
        structure = ref_file.structure

        natom = len(structure)
        nrows, ncols = natom, 3
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=sharey, squeeze=False)

        xs = np.arange(len(self.qpoints))
        ticks, labels = ref_file.make_ticks_and_labels()

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
            if idir == 0: ax.set_ylabel(r"(Hartree/Bohr)")

            site = ref_file.structure[iat]
            s = "%s [%.3f, %.3f, %.3f]" % (site.specie.symbol, site.frac_coords[0], site.frac_coords[1], site.frac_coords[2])
            ax.set_title("idir: %d, iat: %d, %s" % (idir, iat, s), fontsize=fontsize)

            if ticks:
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.get("qlabel_size", "large"))
                ax.set_xlim(ticks[0], ticks[-1])

        return fig

    @add_fig_kwargs
    def plot_maxw(self, ax=None, **kwargs):
        """
        Plot

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        if any(not abifile.has_maxw for abifile in self.abifiles): return None
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for label, abifile in self.items():
            abifile.plot_maxw(ax=ax, label=label, with_title=False, show=False, **kwargs)

        return fig

    def yield_figs(self, **kwargs): # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        for vname in ("v1scf_avg", "v1scf_abs_avg"):
            yield self.plot(vname=vname, title=vname, show=False)

        if all(abifile.has_maxw for abifile in self.abifiles):
            yield self.plot_maxw(self, show=False, **kwargs)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to `nbpath`. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.V1qAvgRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


#class V1qAvgPlotter(object):
#
#    def __init__(self, v1q_dfpt_path, v1q_frohl_path, v1q_frohl_quad_path, v1q_froh_quad_efield_path):
#        self.v1q_dfpt = V1qAvgFile.from_file(v1q_dfpt_path)
#        #v1q_frohl = V1qAvgFile.from_file(v1q_frohl_path) if v1q_frohl_path is not None else None
#        #v1q_frohl_quad_path =
#        #v1q_frohl_quad_efield_path =
#
#    #def close(self):
#    #    for k in ("v1q_frohl", "v1q_frohl_quad_path", "v1q_frohl_quad_efield_path"):
#    #        ncfile = getattr(self, k)
#    #        if ncfile is not None: ncfile.close()
#
#    @add_fig_kwargs
#    def plot(self, iatom, idir, ax=None, ispden=0, fontsize=6, **kwargs):
#
#        ax, fig, plt = get_ax_fig_plt(ax=ax)
#
#        return fig

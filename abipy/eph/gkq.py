"""
Interface to the GKQ.nc file storing the e-ph matrix elements 
in the atomic representation (idir, ipert) for a single q-point.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import abipy.core.abinit_units as abu

from collections import OrderedDict
from monty.string import marquee
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.kpoints import Kpoint
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    get_ax3d_fig_plt, rotate_ticklabels, set_visible, plot_unit_cell, set_ax_xylabels)
from abipy.tools import duck
from abipy.abio.robots import Robot
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands
from abipy.eph.common import glr_frohlich, EPH_WTOL

class GkqFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(GkqFile, self).__init__(filepath)
        self.reader = GkqReader(filepath)
        #self.alpha_gmin = self.reader.read_value("alpha_gmin")

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        app("qpoint: %s" % str(self.qpoint))
        app("Macroscopic dielectric tensor in Cartesian coordinates")
        app(str(self.epsinf_cart))
        app("")
        app("Born effective charges in Cartesian coordinates:")
        for i, (site, bec) in enumerate(zip(self.structure, self.becs_cart)):
            app("[%d]: %s" % (i, repr(site)))
            app(str(bec))
            app("")

        #if verbose > 1:
        #    app("")
        #    app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    def close(self):
        self.reader.close()

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @lazy_property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def uses_interpolated_dvdb(self):
        """True if the matrix elements have been computed with an interpolated potential."""
        return int(self.reader.read_value("interpolated")) == 1

    @lazy_property
    def params(self):
        """Dict with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @lazy_property
    def qpoint(self):
        """Q-point object."""
        return Kpoint(self.reader.read_value('qpoint'), self.structure.reciprocal_lattice)

    @lazy_property
    def phfreqs(self):
        """(3 * natom) array with phonon frequencies in Ha."""
        return self.reader.read_value("phfreqs")

    @lazy_property
    def phdispl_cart(self):
        """(natom3_nu, natom3) complex array with phonon displacement in cartesian coordinates."""
        return self.reader.read_value("phdispl_cart", cmode="c")

    @lazy_property
    def phdispl_red(self):
        """(natom3_nu, natom3) complex array with phonon displacement in reduced coordinates."""
        return self.reader.read_value("phdispl_red", cmode="c")

    @lazy_property
    def becs_cart(self):
        """(natom, 3, 3) array with the Born effective charges in Cartesian coordinates."""
        return self.reader.read_value("becs_cart").transpose(0, 2, 1).copy()

    @lazy_property
    def epsinf_cart(self):
        """(3, 3) array with electronic macroscopic dielectric tensor in Cartesian coordinates."""
        return self.reader.read_value("emacro_cart").T.copy()

    def read_all_gkq(self, mode="atom"):

        if mode not in ("atom", "phonon"):
            raise ValueError("Invalid mode: %s" % mode)

        # Read e-ph matrix element in the atomic representation (idir, ipert)
        # Fortran array on disk has shape:
        # nctkarr_t('gkq', "dp", &
        # 'complex, max_number_of_states, max_number_of_states, number_of_phonon_modes, number_of_kpoints, number_of_spins') &
        gkq_atm = self.reader.read_value("gkq", cmode="c")
        if mode == "atom": return gkq_atm

        # Convert from atomic to phonon representation.
        # May use np.einsum for better efficiency but oh well!
        nband = gkq_atm.shape[-1]
        nb2 = nband ** 2
        assert nband == gkq_atm.shape[-2] and nband == self.ebands.nband
        natom = len(self.structure)
        natom3 = natom * 3
        phfreqs, phdispl_red = self.phfreqs, self.phdispl_red
        gkq_nu = np.empty_like(gkq_atm)
        cwork = np.empty((natom3, nb2), dtype=np.complex)
        for spin in range(self.ebands.nsppol):
            for ik in range(self.ebands.nkpt):
                g = np.reshape(gkq_atm[spin, ik], (-1, nb2))
                for nu in range(natom3):
                    if phfreqs[nu] > EPH_WTOL:
                        cwork[nu] = np.dot(phdispl_red[nu], g) / np.sqrt(2.0 * phfreqs[nu]) 
                    else:
                        cwork[nu] = 0.0
                gkq_nu[spin, ik] = np.reshape(cwork, (natom3, nband, nband))

        return gkq_nu

    @add_fig_kwargs
    def plot(self, mode="phonon", with_glr=True, ax=None, fontsize=8, **kwargs):
        """
        Plot the gkq matrix elements for a given q-point.

            mode: "phonon" to plot eph matrix elements in the phonon representation, 
                "atom" for atomic representation.
            with_glr: True to plot the long-range component estimated from Verdi's model.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Label and title fontsize.

        Return: |matplotlib-Figure|
        """
        gkq = np.abs(self.read_all_gkq(mode=mode))
        if mode == "phonon": gkq *= abu.Ha_meV

        stats = OrderedDict([
            ("min", gkq.min()),
            ("max", gkq.max()),
            ("mean", gkq.mean()),
            ("std", gkq.std()),
        ])

        eigens_kq = self.reader.read_value("eigenvalues_kq") * abu.Ha_eV

        c = np.empty_like(gkq)
        for spin in range(self.ebands.nsppol):
            for ik in range(self.ebands.nkpt):
                for ib_kq in range(self.ebands.mband):
                    for ib_k in range(self.ebands.mband):
                        c[spin, ik, :, ib_k, ib_kq] = abs(eigens_kq[spin, ik, ib_kq] - self.ebands.eigens[spin, ik, ib_k]) 
                       

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        colormap = "viridis"
        cmap = plt.get_cmap(colormap)
        #from matplotlib.colors import Normalize
        #norm = Normalize(vmin=abs_ediff.min(), vmax=abs_ediff.max())

        xs = np.arange(len(gkq.ravel()))
        sc = ax.scatter(xs, gkq.ravel(), alpha=0.9, s=30, c=c.ravel(), cmap=cmap)
                        #facecolors='none', edgecolors='orange')
        plt.colorbar(sc)

        if with_glr and mode == "phonon":
            # Add horizontal bar with matrix elements computed from model (only G = 0, no wavefunctions).
            gkq_lr = glr_frohlich(self.qpoint, self.becs_cart, self.epsinf_cart, 
                                  self.phdispl_cart, self.phfreqs, self.structure)
            for hval in np.abs(gkq_lr) * abu.Ha_meV:
                ax.axhline(hval, color='k', linestyle='dashed', linewidth=1)

        ax.grid(True)
        ax.set_xlabel("Matrix element index")
        ylabel = r"$|g^{atm}_{\bf q}|$" if mode == "atom" else r"$|g_{\bf q}|$ (meV)"
        ax.set_ylabel(ylabel)
        #ax.legend(loc="best", fontsize=fontsize, shadow=True)

        _, max_ = ax.get_ylim()
        ax.text(0.7, 0.7, 
                "\n".join("%s = %.1E" % item for item in stats.items()), 
                fontsize=fontsize, horizontalalignment='center', verticalalignment='center', 
                transform=ax.transAxes)

        return fig

    @add_fig_kwargs
    def plot_diff_with_other(self, other, mode="phonon", labels=None, fontsize=8, **kwargs):
        """
        Compare the gkq matrix elements.

            other: other GkqFile instance.
            mode: "phonon" to plot eph matrix elements in the phonon representation, 
                "atom" for atomic representation.
            labels: Labels associated to self and other
            fontsize: Label and title fontsize.

        Return: |matplotlib-Figure|
        """
        if self.qpoint != other.qpoint:
            raise ValueError("Found different q-points: %s and %s" % (self.qpoint, other.qpoint))

        labels = ["this", "other"] if labels is None else labels

        this_gkq = np.abs(self.read_all_gkq(mode=mode))
        other_gkq = np.abs(other.read_all_gkq(mode=mode))
        if mode == "phonon": 
            this_gkq *= abu.Ha_meV
            other_gkq *= abu.Ha_meV

        absdiff_gkq = np.abs(this_gkq - other_gkq)

        stats = OrderedDict([
            ("min", absdiff_gkq.min()),
            ("max", absdiff_gkq.max()),
            ("mean", absdiff_gkq.mean()),
            ("std", absdiff_gkq.std()),
        ])

        num_plots, ncols, nrows = 2, 2, 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Downsample datasets. Show only points with error > mean.
        threshold = stats["mean"]
        this_gkq = this_gkq[absdiff_gkq > threshold]
        xs = np.arange(len(this_gkq.ravel()))

        ax = ax_list[0]
        ax.scatter(xs, this_gkq.ravel(), alpha=0.9, s=30, label=labels[0], 
                   facecolors='none', edgecolors='orange')

        other_gkq = other_gkq[absdiff_gkq > threshold]
        ax.scatter(xs, other_gkq.ravel(), alpha=0.3, s=10, marker="x", 
                   facecolors="g", edgecolors="none", label=labels[1])

        ax.grid(True)
        ax.set_xlabel("Matrix element index")
        ylabel = r"$\Delta |g^{atm}_{\bf q}|$" if mode == "atom" else r"$\Delta |g_{\bf q}|$ (meV)"
        ax.set_ylabel(ylabel)
        ax.set_title("qpoint: %s" % repr(self.qpoint), fontsize=fontsize)
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        ax = ax_list[1]
        ax.hist(absdiff_gkq.ravel(), facecolor='g', alpha=0.75)
        ax.grid(True)
        ax.set_xlabel("Absolute Error" if mode == "atom" else "Absolute Error (meV)")
        ax.set_ylabel("Count")

        ax.axvline(stats["mean"], color='k', linestyle='dashed', linewidth=1)
        _, max_ = ax.get_ylim()
        ax.text(0.7, 0.7, 
                "\n".join("%s = %.1E" % item for item in stats.items()), 
                fontsize=fontsize, horizontalalignment='center', verticalalignment='center', 
                transform=ax.transAxes)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        for fig in self.yield_structure_figs(**kwargs): yield fig
        for fig in self.yield_ebands_figs(**kwargs): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("gkq = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(gkq)"),
            nbv.new_code_cell("gkq.ebands.plot();"),
            nbv.new_code_cell("gkq.epsinf_cart;"),
            nbv.new_code_cell("gkq.becs_cart;"),
            nbv.new_code_cell("""
              #other = abilab.abiopen('other_GKW.nc')
              #gkq.plot_diff_with_other(other);
            """)
        ])

        return self._write_nb_nbpath(nb, nbpath)


class GkqReader(ElectronsReader):
     """
     This object reads the results stored in the GKQ file produced by ABINIT.
     It provides helper function to access the most important quantities.
 
     .. rubric:: Inheritance Diagram
     .. inheritance-diagram:: GkqReader
     """


class GkqRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple GSR.nc_ files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GkqRobot
    """
    EXT = "GKQ"

    @lazy_property
    def kpoints(self):
        # Consistency check: kmesh should be the same in each file.
        ref_kpoints = self.abifiles[0].ebands.kpoints
        for i, abifile in enumerate(self.abifiles):
            if i == 0: continue
            if abifile.kpoints != ref_kpoints:
                for k1, k2 in zip(ref_kpoints, abifile.kpoints):
                    print("k1:", k1, "--- k2:", k2)
                raise ValueError("Found different list of kpoints in %s" % str(abifile.filepath))
        return ref_kpoints

    #@lazy_property
    #def qpoints(self):
    #    return [abifile.qpoint for abifile in self.abifiles]

    @add_fig_kwargs
    def plot_gkq2_qpath(self, band_kq, band_k, kpoint=0, with_glr=True, nu_list=None, # spherical_average=False, 
                        ax=None, fontsize=12, **kwargs):
        r"""
        Plot the magnitude of the electron-phonon matrix elements <k+q, band_kq| Delta_{q\nu} | k, band_k>
        for a given set of (band_kq, band, k) as a function of the q-point.

        Args:
            band_ks: Band index of the k+q states (starts at 0)
            band_k: Band index of the k state (starts at 0)
            kpoint: |Kpoint| object or index.
            with_glr: True to plot the long-range component estimated from Verdi's model.
            nu_list: List of phonons modes to be selected (starts at 0). None to select all modes.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Label and title fontsize.

        Return: |matplotlib-Figure|
        """
        if duck.is_intlike(kpoint):
            ik = kpoint
            kpoint = self.kpoints[ik]
        else:
            kpoint = Kpoint.as_kpoint(kpoint, self.abifiles[0].structure.reciprocal_lattice)
            ik = self.kpoints.index(kpoint)

        # Assume abifiles are already ordered according to q-path.
        xs = list(range(len(self.abifiles)))
        natom3 = len(self.abifiles[0].structure) * 3
        nsppol = self.abifiles[0].nsppol
        nqpt = len(self.abifiles)
        gkq_snuq = np.empty((nsppol, natom3, nqpt), dtype=np.complex)
        if with_glr: gkq_lr = np.empty((nsppol, natom3, nqpt), dtype=np.complex)

        xticks, xlabels = [], []
        for iq, abifile in enumerate(self.abifiles):
            qpoint = abifile.qpoint
            
            name = qpoint.name if qpoint.name is not None else abifile.structure.findname_in_hsym_stars(qpoint)
            if qpoint.name is not None:
                xticks.append(iq)
                xlabels.append(name)

            phfreqs, phdispl_red = abifile.phfreqs, abifile.phdispl_red
            ncvar = abifile.reader.read_variable("gkq")
            for spin in range(nsppol):
                gkq_atm = ncvar[spin, ik, :, band_k, band_kq] 
                gkq_atm = gkq_atm[:, 0] + 1j * gkq_atm[:, 1]
                #gkq_snuq[spin, :, iq] = np.abs(gkq_atm)

                # Transform the gkk matrix elements from (atom, red_direction) basis to phonon-mode basis.
                gkq_snuq[spin, :, iq] = 0.0
                for nu in range(natom3):
                    if phfreqs[nu] < EPH_WTOL: continue
                    #fact = one if not spherical_average else np.sqrt(4 * np.pi) * qpoint.norm
                    gkq_snuq[spin, nu, iq] = np.dot(phdispl_red[nu], gkq_atm) / np.sqrt(2.0 * phfreqs[nu]) 

            if with_glr:
                # Compute long range part with (simplified) generalized Frohlich model.
                #fact = one if not spherical_average else np.sqrt(4 * np.pi) * qpoint.norm
                gkq_lr[spin, :, iq] = glr_frohlich(qpoint, abifile.becs_cart, abifile.epsinf_cart, 
                                                   abifile.phdispl_cart, phfreqs, abifile.structure)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        nu_list = list(range(natom3)) if nu_list is None else list(nu_list)
        for spin in range(nsppol):
            for nu in nu_list:
                ys = np.abs(gkq_snuq[spin, nu]) * abu.Ha_meV
                label = "nu: %s" % nu if nsppol == 1 else "nu: %s, spin: %s" % (nu, spin)
                ax.plot(xs, ys, linestyle="--",  label=label)
                if with_glr:
                    ys = np.abs(gkq_lr[spin, nu]) * abu.Ha_meV
                    label = r"$g_{\bf q}^{\text{lr}} nu: %s" % (
                            nu if nsppol == 1 else "nu: %s, spin: %s" % (nu, spin))
                    ax.plot(xs, ys, linestyle="", marker="o") #, label=label)

        ax.grid(True)
        ax.set_xlabel("Wave Vector")
        ax.set_ylabel(r"$|g_{\bf q}|$ (meV)")
        if xticks:
            ax.set_xticks(xticks, minor=False)
            ax.set_xticklabels(xlabels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        title = "band_kq: %s, band_k: %s, kpoint: %s" % (band_kq, band_k, repr(kpoint))
        ax.set_title(title, fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs): # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        for fig in self.get_ebands_plotter().yield_figs(): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.GkqRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

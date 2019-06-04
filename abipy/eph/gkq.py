"""
Interface to the GKQ.nc_ file storing the e-ph matrix elements for a single q-point.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

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


EPH_WTOL = 1e-6


class GkqFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands): #, NotebookWriter):

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(GkqFile, self).__init__(filepath)
        self.reader = GkqReader(filepath)

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
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @lazy_property
    def qpoint(self):
        """Q-point object."""
        return Kpoint(self.reader.read_value('qpoint'), self.structure.reciprocal_lattice)

    @lazy_property
    def phfreqs(self):
        """array with 3 * natom phonon frequencies in Ha."""
        return self.reader.read_value("phfreqs")

    @lazy_property
    def phdispl_cart(self):
        """(natom3_nu, natom3) complex array with ph displacement in cartesian coordinates."""
        return self.reader.read_value("phdispl_cart", cmode="c")

    @lazy_property
    def phdispl_red(self):
        """(natom3_nu, natom3) complex array with ph displacement in reduced coordinates."""
        return self.reader.read_value("phdispl_red", cmode="c")

    @lazy_property
    def becs_cart(self):
        """(natom, 3, 3) array with the Born effective charges in Cartesian coordinates."""
        return self.reader.read_value("becs_cart").T.copy()

    @lazy_property
    def epsinf_cart(self):
        """(3, 3) array with macroscopic dielectric tensor in Cartesian coordinates."""
        return self.reader.read_value("emacro_cart").T.copy()

    #def get_all_gkkp(self):
    #    """
    #    Get gkq matrix elements from the file
    #    """
    #    #nctkarr_t("gkq_representation", "char", "character_string_length"), &
    #    #nctkarr_t('gkq', "dp", &
    #    # 'complex, max_number_of_states, max_number_of_states, number_of_phonon_modes, number_of_kpoints, number_of_spins') &
    #    return self.reader.read_value("gkq", cmode="c")

    @add_fig_kwargs
    def plot_scatter_with_other(self, other, ax=None, fontsize=12, **kwargs):
        """
        Compare gkq_atm matrix elements for a given q-point.

            other: other GkqFile instance.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Return: |matplotlib-Figure|
        """
        if self.qpoint != other.qpoint:
            raise ValueError("Found different q-points: %s and %s" % (self.qpoint, other.qpoint))

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        this_gkq_atm = np.abs(self.reader.read_value("gkq", cmode="c"))
        other_gkq_atm = np.abs(other.reader.read_value("gkq", cmode="c"))

        adiff_gkq_atm = np.abs(this_gkq_atm - other_gkq_atm)
        reldiff_gkq_atm = np.abs(this_gkq_atm - other_gkq_atm) / (0.5 * np.abs(this_gkq_atm + other_gkq_atm) + 1e-20)

        print("max:", adiff_gkq_atm.max())
        print("min:", adiff_gkq_atm.min())
        print("mean:", adiff_gkq_atm.mean())
        print("std:", adiff_gkq_atm.std())

        xs = np.arange(len(this_gkq_atm.ravel()))
        #ax.scatter(xs, this_gkq_atm.ravel(), label="this")
        #ax.scatter(xs, other_gkq_atm.ravel(), label="other")
        ax.scatter(xs, adiff_gkq_atm.ravel(), label="abs_diff")
        ax.scatter(xs, reldiff_gkq_atm.ravel(), label="rel_diff")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig


class GkqReader(ElectronsReader):
     """
     This object reads the results stored in the GKQ file produced by ABINIT.
     It provides helper function to access the most important quantities.
 
     .. rubric:: Inheritance Diagram
     .. inheritance-diagram:: GsrReader
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
    def plot_gkq2_qpath(self, band_kq, band_k, kpoint=0, nu_list=None, with_glr=True, 
                        ax=None, fontsize=12, **kwargs):
        """
        Plot the magnitude of the electron-phonon matrix elements along a path
        """
        if duck.is_intlike(kpoint):
            ik = kpoint
            kpoint = self.kpoints[ik]
        else:
            kpoint = Kpoint.as_kpoint(kpoint, self.abifiles[0].structure.reciprocal_lattice)
            ik = self.kpoints.index(kpoint)

        # Assume abifiles are already ordered according to q-path
        xs = list(range(len(self.abifiles)))
        natom3 = len(self.abifiles[0].structure) * 3
        nsppol = self.abifiles[0].nsppol
        nqpt = len(self.abifiles)
        gkq_snuq = np.empty((nsppol, natom3, nqpt), dtype=np.complex)
        if with_glr: gkq_lr = np.empty((nsppol, natom3, nqpt), dtype=np.complex)


        xticks, xlabels = [], []
        for iq, abifile in enumerate(self.abifiles):
            qpoint = abifile.qpoint
            #qnorm = qpoint.norm
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
                    gkq_snuq[spin, nu, iq] = np.dot(phdispl_red[nu], gkq_atm) / np.sqrt(2.0 * phfreqs[nu]) 

            if with_glr:
                # Compute g long range with (simplified) generalized Frohlich model.
                gkq_lr[spin, :, iq] = glr_frohlich(qpoint, abifile.becs_cart, abifile.epsinf_cart, 
                                                   phdispl_cart, phfreqs, abifile.structure)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        
        nu_list = list(range(natom3)) if nu_list is None else list(nu_list)
        for spin in range(nsppol):
            for nu in nu_list:

                ys = np.abs(gkq_snuq[spin, nu]) * abu.Ha_meV
                label = "nu: %s" % nu if nsppol == 1 else "nu: %s, spin: %s" % (nu, spin)
                ax.plot(xs, ys, ls="-o", label=label)

                if with_glr:
                    ys = np.abs(gkq_lr[spin, nu]) * abu.Ha_meV
                    label = "glr nu: %s" % nu if nsppol == 1 else "nu: %s, spin: %s" % (nu, spin)
                    ax.plot(xs, ys, ls="-o", label=label)

        ax.grid(True)
        ax.set_xlabel("Wave Vector")
        ax.set_ylabel("gkq (meV)")
        if xticks:
            ax.set_xticks(xticks, minor=False)
            ax.set_xticklabels(xlabels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        title = "band_kq: %s, band_k: %s, kpoint: %s" % (band_kq, band_k, repr(kpoint) )
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
            #nbv.new_code_cell("robot = abilab.GsrRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)


def glr_frohlich(qpoint, zeff_cart, eps_inf, phdispl_cart, phfreqs, structure, tol_qnorm=1e-6):
    """
    Compute the long-range part of the e-ph matrix element with the simplified Frohlich model 
    i.e. we include only G = 0 and the <k+q,b1|e^{i(q+G).r}|b2,k> coefficient is replaced by delta_{b1,b2}

    Args:
        qpoint:
        zeff_cart:
        eps_inf:
        phdispl_cart:
        phfreqs:
        structure:

    Return:
        (natom3) complex array with gkq_LR.
    """
    natom = len(structure)
    natom3 = natom * 3
    qeq = np.dot(qpoint.cart_coords, np.matmul(eps_inf, qpoint.cart_coords))
    #if qpoint.is_gamma

    glr = np.zeros(natom3)
    for nu in range(3 if qpt.norm < tol_qnorm else 0, natom3):
        if phfreqs[nu] < EPH_WTOL: continue
        num = 0.0
        for iat in range(natom):
            num += np.dot(qpoint.cart_coords, np.matmul(zeff_cart[iat], phdispl_cart[nu, iat]))
        glr[nu] = num / (qeq * np.sqrt(two * phfreqs[nu]))

    return glr * 4j * np.pi / structure.volume

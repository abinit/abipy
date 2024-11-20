"""
Interface to the GKQ.nc file storing the e-ph matrix elements
in the atomic representation (idir, ipert) for a single q-point.
This file is produced by the eph code with eph_task -4.

To analyze the e-ph scattering potentials, use v1qavg and eph_task 15 or -15
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from collections import OrderedDict
from monty.string import marquee
from monty.functools import lazy_property
from abipy.core.structure import Structure
from abipy.core.kpoints import Kpoint
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.tools import duck
from abipy.tools.typing import Figure, PathLike
from abipy.abio.robots import Robot
from abipy.electrons.ebands import ElectronBands, ElectronsReader, RobotWithEbands
from abipy.eph.common import glr_frohlich, EPH_WTOL


class GkqFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):

    @classmethod
    def from_file(cls, filepath: PathLike):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = GkqReader(filepath)

    def __str__(self) -> str:
        """String representation."""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        app("qpoint in g(k,q): %s" % str(self.qpoint))
        app("uses_interpolated_dvdb: %s" % str(self.uses_interpolated_dvdb))
        app("phonon frequencies in Ha %s:" % str(self.phfreqs_ha))
        app("Macroscopic dielectric tensor in Cartesian coordinates:")
        app(str(self.epsinf_cart))
        app("")
        app("Born effective charges in Cartesian coordinates:")
        for i, (site, bec) in enumerate(zip(self.structure, self.becs_cart)):
            app("[%d]: %s" % (i, repr(site)))
            app(str(bec))
            app("")

        app(r"Fulfillment of charge neutrality, F_{ij} = \sum_{atom} Z^*_{ij,atom}:")
        f = np.sum(self.becs_cart, axis=0)
        app(str(f) + "\n")

        return "\n".join(lines)

    def close(self) -> None:
        self.r.close()

    @lazy_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @lazy_property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def uses_interpolated_dvdb(self) -> bool:
        """True if the matrix elements have been computed with an interpolated potential."""
        return int(self.r.read_value("interpolated")) == 1

    @lazy_property
    def params(self) -> dict:
        """Dict with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @lazy_property
    def qpoint(self) -> Kpoint:
        """Q-point object."""
        return Kpoint(self.r.read_value('qpoint'), self.structure.reciprocal_lattice)

    @lazy_property
    def phfreqs_ha(self) -> np.ndarray:
        """(3 * natom) array with phonon frequencies in Ha."""
        return self.r.read_value("phfreqs")

    @lazy_property
    def phdispl_cart_bohr(self) -> np.ndarray:
        """(natom3_nu, natom3) complex array with the phonon displacement in cartesian coordinates in Bohr."""
        return self.r.read_value("phdispl_cart", cmode="c")

    @lazy_property
    def phdispl_red(self) -> np.ndarray:
        """(natom3_nu, natom3) complex array with the phonon displacement in reduced coordinates."""
        return self.r.read_value("phdispl_red", cmode="c")

    @lazy_property
    def becs_cart(self) -> np.ndarray:
        """(natom, 3, 3) array with the Born effective charges in Cartesian coordinates."""
        return self.r.read_value("becs_cart").transpose(0, 2, 1).copy()

    @lazy_property
    def epsinf_cart(self) -> np.ndarray:
        """(3, 3) array with electronic macroscopic dielectric tensor in Cartesian coordinates."""
        return self.r.read_value("emacro_cart").T.copy()

    @lazy_property
    def eigens_kq(self) -> np.ndarray:
        """(spin, nkpt, mband) array with eigenvalues on the k+q grid in eV."""
        return self.r.read_value("eigenvalues_kq") * abu.Ha_eV

    @staticmethod
    def _check_mode(mode: str) -> None:
        """Check whether mode is allowed."""
        allowed_modes = ("atom", "phonon")
        if mode not in allowed_modes:
            raise ValueError(f"Invalid {mode=}, it should be in {allowed_modes=}")

    def read_gkq_kpoint(self, kpoint, mode: str = "phonon") -> np.ndarray:
        """
        Read e-ph matrix stored on disk for all spins and the given k-point

        Args:
            kpoint:
            mode: "phonon" for e-ph matrix elements in phonon representation,
                  "atom" for e-ph matrix elements in the atomic representation (idir, iatom).

        Return: complex array with shape: (nsppol, 3*natom, mband, mband)
                                                            m_kq, n_k    <-- band indices.
        """
        self._check_mode(mode)

    def read_all_gkq(self, mode: str = "phonon") -> np.ndarray:
        """
        Read all e-ph matrix stored on disk.

        Args:
            mode: "phonon" for e-ph matrix elements in phonon representation,
                  "atom" for e-ph matrix elements in the atomic representation (idir, iatom).

        Return: complex array with shape: (nsppol, nkpt, 3*natom, mband, mband)
                                                                   m_kq, n_k    <-- band indices.
        """
        self._check_mode(mode)

        # Read the e-ph matrix element in the atomic representation (idir, ipert). Fortran array on disk has shape:
        # nctkarr_t('gkq', "dp",
        #           'complex, max_number_of_states, max_number_of_states, number_of_phonon_modes, number_of_kpoints, number_of_spins')
        # The first band index in Fortran refers to m_kq, the second one to n_k.
        # hence we have to transpose the (nb_kq, nb_k) submatrix written by ABINIT.
        gkq_atm = self.r.read_value("gkq", cmode="c").transpose(0, 1, 2, 4, 3).copy()
        if mode == "atom":
            return gkq_atm

        # Convert from atomic to phonon representation. May use np.einsum for better efficiency but oh well!
        nband = gkq_atm.shape[-1]
        assert nband == gkq_atm.shape[-2] and nband == self.ebands.nband
        nb2, natom3 = nband ** 2, 3 * len(self.structure)
        phfreqs_ha, phdispl_red = self.phfreqs_ha, self.phdispl_red

        gkq_nu, cwork = np.empty_like(gkq_atm), np.empty((natom3, nb2), dtype=complex)
        for spin in range(self.ebands.nsppol):
            for ik in range(self.ebands.nkpt):
                gc = np.reshape(gkq_atm[spin, ik], (-1, nb2))
                for nu in range(natom3):
                    cwork[nu] = np.dot(phdispl_red[nu], gc) / np.sqrt(2.0 * phfreqs_ha[nu]) if (phfreqs_ha[nu] > EPH_WTOL) else 0.0
                gkq_nu[spin, ik] = np.reshape(cwork, (natom3, nband, nband))

        return gkq_nu

    def get_absg_kpoint(self, kpoint, eps_mev: float=0.01) -> tuple[np.ndarray, np.ndarray, int, Kpoint]:
        """
        Args:
            kpoint: |Kpoint| object or list/tuple with reduced coordinates or integer with the index
            eps_mev: Tolerance in mev used to detect degeneracies
        """
        if duck.is_intlike(kpoint):
            ik = kpoint
            kpoint = self.kpoints[ik]
        else:
            kpoint = Kpoint.as_kpoint(kpoint, self.structure.reciprocal_lattice)
            ik = self.kpoints.index(kpoint)

        eps_ha = eps_mev / abu.Ha_meV
        eps_ev = eps_ha * abu.Ha_eV

        nsppol = self.ebands.nsppol
        natom3 = len(self.structure) * 3
        nb = self.ebands.nband

        phfreqs_ha = self.phfreqs_ha
        eigens_k = self.ebands.eigens
        eigens_kq = self.eigens_kq

        # (nsppol, nkpt, 3*natom, mband, mband) real array.
        absg = np.abs(self.read_all_gkq(mode="phonon")) * abu.Ha_meV
        absgk = absg[:,ik].copy()
        absg_unsym = absg[:,ik].copy()
        absg_sym = np.zeros_like(absgk)

        # Average over phonons.
        for spin in range(nsppol):
            g2_mn = np.zeros((nb, nb), dtype=float)
            for nu in range(natom3):
                w_1 = phfreqs_ha[nu]
                g2_mn[:], nn = 0.0, 0
                for mu in range(natom3):
                    w_2 = phfreqs_ha[mu]
                    if abs(w_1 - w_2) >= eps_ha: continue
                    nn += 1
                    g2_mn += absgk[spin,mu,:,:] ** 2
                absg_sym[spin,nu,:,:] = np.sqrt(g2_mn / nn)

        # Average over k electrons.
        absg = absg_sym.copy()
        g2_nu = np.zeros((natom3), dtype=float)
        for spin in range(nsppol):
            for jbnd in range(nb):
                for ibnd in range(nb):
                    w_1 = eigens_k[spin, ik, ibnd]
                    g2_nu[:], nn = 0.0, 0
                    for pbnd in range(nb):
                        w_2 = eigens_k[spin, ik, pbnd]
                        if abs(w_2 - w_1) >= eps_ev: continue
                        nn += 1
                        # MG FIXME: Why absgk and not absg here as done below for k+q?
                        g2_nu += absgk[spin,:,jbnd,pbnd] ** 2
                    absg_sym[spin,:,jbnd,ibnd] = np.sqrt(g2_nu / nn)

        # Average over k+q electrons.
        absgk = absg_sym.copy()
        for spin in range(nsppol):
            for ibnd in range(nb):
                for jbnd in range(nb):
                    w_1 = eigens_kq[spin, ik, jbnd]
                    g2_nu[:], nn = 0.0, 0
                    for pbnd in range(nb):
                        w_2 = eigens_kq[spin, ik, pbnd]
                        if abs(w_2 - w_1) >= eps_ev: continue
                        nn += 1
                        g2_nu += absgk[spin,:,pbnd,ibnd] ** 2
                    absg_sym[spin,:,jbnd,ibnd] = np.sqrt(g2_nu / nn)

        return absg_sym, absg_unsym, ik, kpoint

    def get_qe_dataframe(self, kpoint) -> pd.DataFrame:
        """
        Build and return a dataframe with |g(k,q)|^2 for the given k-point and all bands.

        Args:
            kpoint:
        """
        absg, absg_unsym, ik, kpoint = self.get_absg_kpoint(kpoint)

        # Now insert absg array in a pandas dataframe.
        # Flatten the array, get the indices and combine indices and values into a DataFrame
        shape, ndim = absg.shape, absg.ndim
        indices = np.indices(shape).reshape(ndim, -1).T
        df = pd.DataFrame(indices, columns=["spin", "imode", "m_kq", "n_k"])
        df["|g|[meV]"] = absg.flatten()
        df["ik"] = ik

        # Add columns with phonon frequencies and electron energies in meV at k and k+q.
        imodes = df["imode"].to_numpy()
        df["omega(q)[meV]"] = (self.phfreqs_ha * abu.Ha_meV)[imodes]
        spin_inds, mkq_inds, nk_inds = df["spin"].to_numpy(), df["m_kq"].to_numpy(), df["n_k"].to_numpy()
        #print(self.ebands.eigens[spin_inds,ik,nk_inds].shape)
        df["e_nk[eV]"] = self.ebands.eigens[spin_inds, ik, nk_inds]
        df["e_mkq[eV]"] = self.eigens_kq[spin_inds, ik, mkq_inds]

        # Reorder the columns and drop the index
        new_order = ["n_k", "m_kq", "spin", "imode", "e_nk[eV]", "e_mkq[eV]", "omega(q)[meV]", "|g|[meV]"]
        return df[new_order].reset_index(drop=True)

    @add_fig_kwargs
    def plot(self, mode="phonon", with_glr=True, fontsize=8, colormap="viridis", sharey=True, **kwargs) -> Figure:
        """
        Plot the gkq matrix elements for a given q-point.

        Args:
            mode: "phonon" to plot eph matrix elements in the phonon representation,
                  "atom" for atomic representation.
            with_glr: True to plot the long-range component estimated from Verdi's model.
            fontsize: Label and title fontsize.
            colormap: matplotlib colormap
            sharey: True if yaxis should be shared among axes.

        Return: |matplotlib-Figure|
        """
        gkq = np.abs(self.read_all_gkq(mode=mode))
        if mode == "phonon": gkq *= abu.Ha_meV

        # Compute e_{k+q} - e_k for all possible (b, b')
        ediffs = np.empty_like(gkq)
        for spin in range(self.ebands.nsppol):
            for ik in range(self.ebands.nkpt):
                for ib_kq in range(self.ebands.mband):
                    for ib_k in range(self.ebands.mband):
                        ed = np.abs(self.eigens_kq[spin, ik, ib_kq] - self.ebands.eigens[spin, ik, ib_k])
                        ediffs[spin, ik, :, ib_k, ib_kq] = ed

        if with_glr and mode == "phonon":
            # Add horizontal bar with matrix elements computed from Verdi's model (only G = 0, \delta_nm in bands).
            dcart_bohr = self.phdispl_cart_bohr
            #dcart_bohr = self.r.read_value("phdispl_cart_qvers", cmode="c").real
            gkq_lr = glr_frohlich(self.qpoint, self.becs_cart, self.epsinf_cart,
                                  dcart_bohr, self.phfreqs_ha, self.structure)
            # self.phdispl_cart_bohr, self.phfreqs_ha, self.structure)
            gkq2_lr = np.abs(gkq_lr) * abu.Ha_meV

        natom = len(self.structure)
        num_plots, ncols, nrows = 3 * natom, 3, natom
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()
        cmap = plt.get_cmap(colormap)

        for nu, ax in enumerate(ax_list):
            idir = nu % 3
            iat = (nu - idir) // 3
            data, c = gkq[:, :, nu, :, :].ravel(), ediffs[:,:,nu,:,:].ravel()
            # Filter items according to ediff
            index = c <= 1.2 * self.phfreqs_ha.max() * abu.Ha_eV
            data, c = data[index], c[index]
            sc = ax.scatter(np.arange(len(data)), data, alpha=0.9, s=30, c=c, cmap=cmap)
                            #facecolors='none', edgecolors='orange')
            plt.colorbar(sc, ax=ax)

            ax.grid(True)
            if iat == natom - 1:
                ax.set_xlabel("Matrix element index")
            if idir == 0:
                ylabel = r"$|g^{atm}_{\bf q}|$" if mode == "atom" else r"$|g_{\bf q}|$ (meV)"
                ax.set_ylabel(ylabel)

            ax.set_title(r"$\nu$: %d, $\omega_{{\bf q}\nu}$ = %.2E (meV)" %
                         (nu, self.phfreqs_ha[nu] * abu.Ha_meV), fontsize=fontsize)

            if with_glr:
                ax.axhline(gkq2_lr[nu], color='k', linestyle='dashed', linewidth=2)

        fig.suptitle("qpoint: %s" % repr(self.qpoint), fontsize=fontsize)
        return fig

    @add_fig_kwargs
    def plot_diff_with_other(self, other, mode="phonon", ax_list=None, labels=None, fontsize=8, **kwargs) -> Figure:
        """
        Produce scatter plot and histogram to compare the gkq matrix elements stored in two files.

            other: other GkqFile instance.
            mode: "phonon" to plot eph matrix elements in the phonon representation,
                  "atom" for atomic representation.
            ax_list: List with 2 matplotlib axis. None if new ax_list should be created.
            labels: Labels associated to self and other
            fontsize: Label and title fontsize.

        Return: |matplotlib-Figure|
        """
        if self.qpoint != other.qpoint:
            raise ValueError("Found different q-points: %s and %s" % (self.qpoint, other.qpoint))

        if labels is None:
            labels = ["this (interpolated: %s)" % self.uses_interpolated_dvdb,
                      "other (interpolated: %s)" % other.uses_interpolated_dvdb]

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
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Downsample datasets. Show only points with error > threshold.
        ntot = absdiff_gkq.size
        threshold = stats["mean"] + stats["std"]
        data = this_gkq[absdiff_gkq > threshold].ravel()
        nshown = len(data)
        xs = np.arange(len(data))

        ax = ax_list[0]
        ax.scatter(xs, data, alpha=0.9, s=30, label=labels[0], facecolors='none', edgecolors='orange')

        data = other_gkq[absdiff_gkq > threshold].ravel()
        ax.scatter(xs, data, alpha=0.3, s=10, marker="x", label=labels[1], facecolors="g", edgecolors="none")

        ax.grid(True)
        ax.set_xlabel("Matrix element index")
        ylabel = r"$|g^{atm}_{\bf q}|$" if mode == "atom" else r"$|g_{\bf q}|$ (meV)"
        ax.set_ylabel(ylabel)
        ax.set_title(r"qpt: %s, $\Delta$ > %.1E (%.1f %%)" % (
                     repr(self.qpoint), threshold, 100 * nshown / ntot), fontsize=fontsize)
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        ax = ax_list[1]
        ax.hist(absdiff_gkq.ravel(), facecolor='g', alpha=0.75)
        ax.grid(True)
        ax.set_xlabel("Absolute Error" if mode == "atom" else "Absolute Error (meV)")
        ax.set_ylabel("Count")

        ax.axvline(stats["mean"], color='k', linestyle='dashed', linewidth=1)
        _, max_ = ax.get_ylim()
        ax.text(0.7, 0.7,  "\n".join("%s = %.1E" % item for item in stats.items()),
                fontsize=fontsize, horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.plot()

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
              #with abilab.abiopen('other_GKQ.nc') as other:
              #     gkq.plot_diff_with_other(other);
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
    This robot analyzes the results contained in multiple GKQ.nc files.

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

    def _check_qpoints_equal(self):
        """Raises ValueError if different `qpoint` in files."""
        ref_qpoint = self.abifiles[0].qpoint
        for i, abifile in enumerate(self.abifiles):
            if i == 0: continue
            if abifile.qpoint != ref_qpoint:
                raise ValueError("Found different qpoint in %s" % str(abifile.filepath))

    #@add_fig_kwargs
    #def plot_gkq2_qpath(self, band_kq, band_k, kpoint=0, with_glr=False, qdamp=None, nu_list=None, # spherical_average=False,
    #                    ax=None, fontsize=8, eph_wtol=EPH_WTOL, **kwargs):
    #    ncols, nrows = 2, len(self) - 1
    #    num_plots = ncols * nrows
    #    ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,

    @add_fig_kwargs
    def plot_gkq2_qpath(self, band_kq, band_k, kpoint=0, with_glr=False, qdamp=None, nu_list=None, # spherical_average=False,
                        ax=None, fontsize=8, eph_wtol=EPH_WTOL, kq_labels=False, **kwargs) -> Figure:
        r"""
        Plot the magnitude of the electron-phonon matrix elements <k+q, band_kq| Delta_{q\nu} V |k, band_k>
        for a given set of (band_kq, band, k) as a function of the q-point.

        Args:
            band_ks: Band index of the k+q states (starts at 0)
            band_k: Band index of the k state (starts at 0)
            kpoint: |Kpoint| object or index.
            with_glr: True to plot the long-range component estimated from Verdi's model.
            qdamp:
            nu_list: List of phonons modes to be selected (starts at 0). None to select all modes.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Label and title fontsize.
            kq_labels: If True, use the label associated to k+q instead of q.

        Return: |matplotlib-Figure|
        """
        if duck.is_intlike(kpoint):
            ik = kpoint
            kpoint = self.kpoints[ik]
        else:
            kpoint = Kpoint.as_kpoint(kpoint, self.abifiles[0].structure.reciprocal_lattice)
            ik = self.kpoints.index(kpoint)

        # Assume abifiles in the robot are already ordered according to q-path.
        xs = list(range(len(self.abifiles)))
        natom3 = len(self.abifiles[0].structure) * 3
        nsppol = self.abifiles[0].nsppol
        nqpt = len(self.abifiles)
        gkq_snuq = np.empty((nsppol, natom3, nqpt), dtype=complex)
        if with_glr: gkq_lr = np.empty((nsppol, natom3, nqpt), dtype=complex)

        # TODO: Should take into account possible degeneracies in k and k+q and phonon modes.
        xticks, xlabels = [], []
        for iq, abifile in enumerate(self.abifiles):
            qpoint = abifile.qpoint
            #d3q_fact = one if not spherical_average else np.sqrt(4 * np.pi) * qpoint.norm

            if kq_labels:
                name = abifile.structure.findname_in_hsym_stars(kpoint + qpoint)
            else:
                name = qpoint.name if qpoint.name is not None else abifile.structure.findname_in_hsym_stars(qpoint)

            if name is not None:
                xticks.append(iq)
                xlabels.append(name)

            phfreqs_ha, phdispl_red = abifile.phfreqs_ha, abifile.phdispl_red
            ncvar = abifile.r.read_variable("gkq")
            for spin in range(nsppol):
                gkq_atm = ncvar[spin, ik, :, band_k, band_kq]
                gkq_atm = gkq_atm[:, 0] + 1j * gkq_atm[:, 1]

                # Transform the gkq matrix elements from (atom, red_direction) basis to phonon-mode basis.
                gkq_snuq[spin, :, iq] = 0.0
                for nu in range(natom3):
                    if phfreqs_ha[nu] < eph_wtol: continue
                    gkq_snuq[spin, nu, iq] = np.dot(phdispl_red[nu], gkq_atm) / np.sqrt(2.0 * phfreqs_ha[nu])

            if with_glr:
                # Compute long range part with (simplified) generalized Frohlich model.
                gkq_lr[spin, :, iq] = glr_frohlich(qpoint, abifile.becs_cart, abifile.epsinf_cart,
                                                   abifile.phdispl_cart_bohr, phfreqs_ha, abifile.structure, qdamp=qdamp)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        nu_list = list(range(natom3)) if nu_list is None else list(nu_list)
        for spin in range(nsppol):
            for nu in nu_list:
                ys = np.abs(gkq_snuq[spin, nu]) * abu.Ha_meV
                pre_label = kwargs.pop("pre_label",r"$g_{\bf q}$")
                if nsppol == 1: label = r"%s $\nu$: %s" % (pre_label, nu)
                if nsppol == 2: label = r"%s $\nu$: %s, spin: %s" % (pre_label, nu, spin)
                ax.plot(xs, ys, linestyle="--", label=label)
                if with_glr:
                    # Plot model with G = 0 and delta_nn'
                    ys = np.abs(gkq_lr[spin, nu]) * abu.Ha_meV
                    label = r"$g_{\bf q}^{\mathrm{lr0}}$ $\nu$: %s" % nu
                    ax.plot(xs, ys, linestyle="", marker="o", label=label)

        ax.grid(True)
        ax.set_xlabel("Wave Vector")
        ax.set_ylabel(r"$|g_{\bf q}|$ (meV)")
        if xticks:
            ax.set_xticks(xticks, minor=False)
            ax.set_xticklabels(xlabels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))

        ax.legend(loc="best", fontsize=fontsize, shadow=True)
        title = r"$band_{{\bf k} + {\bf q}: %s, band_{\bf{k}}: %s, kpoint: %s" % (band_kq, band_k, repr(kpoint))
        ax.set_title(title, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_gkq2_diff(self, iref=0, **kwargs) -> Figure:
        """
        Wraps gkq.plot_diff_with_other
        Produce scatter and histogram plot to compare the gkq matrix elements stored in all the files
        contained in the robot. Assume all files have the same q-point. Compare the `iref` file with others.
        kwargs are passed to `plot_diff_with_other`.
        """
        if len(self) <= 1: return None
        self._check_qpoints_equal()

        ncols, nrows = 2, len(self) - 1
        num_plots = ncols * nrows
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        ref_gkq, ref_label = self.abifiles[iref], self.labels[iref]
        cnt = -1
        for ifile, (other_label, other_gkq) in enumerate(zip(self.labels, self.abifiles)):
            if ifile == iref: continue
            cnt += 1
            labels = [ref_label, other_label]
            ref_gkq.plot_diff_with_other(other_gkq, ax_list=ax_mat[cnt], labels=labels, show=False, **kwargs)

        return fig

    def yield_figs(self, **kwargs): # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        for fig in self.get_ebands_plotter().yield_figs(): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to `nbpath`. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.GkqRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("# robot.plot_gkq2_diff();"),
            nbv.new_code_cell("# robot.plot_gkq2_qpath(band_kq=0, band_k=0, kpoint=0, with_glr=True, qdamp=None);")
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

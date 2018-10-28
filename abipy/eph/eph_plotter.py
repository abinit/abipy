# coding: utf-8
"""
Objects to plot electronic, vibrational and e-ph properties.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from monty.string import marquee, list_strings
from monty.termcolor import cprint
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, ax_share)
from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.phonons import PhbstFile, PhdosFile
from abipy.eph.sigeph import SigEPhFile


class EphPlotter(object):
    """
    This object provides methods to plot electron and phonons for a single system.
    An EphPlotter has:

        - An |ElectronBands| on a k-path
        - An |ElectronBands| on a k-mesh (optional)
        - A |PhbstFile| with phonons along a q-path.
        - A |PhdosFile| with different kinds of phonon DOSes.

    EphPlotter uses these objects/files and other inputs/files provided by
    the user to generate matplotlib plots related to e-ph interaction.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EphPlotter
    """
    @classmethod
    def from_ddb(cls, ddb, ebands_kpath, ebands_kmesh=None, **kwargs):
        """
        Build the object from the ddb file, invoke anaddb to get phonon properties.
        This entry point is needed to have phonon plots with LO-TO splitting
        as AbiPy will generate an anaddb input with the different q --> 0 directions
        required in phbands.plot to plot the LO-TO splitting correctly.

        Args:
            ddb: |DdbFile| or filepath.
            ebands_kpath: |ElectronBands| with energies on a k-path or filepath.
            ebands_kpath: (optional) |ElectronBands| with energies on a k-mesh or filepath.
            kwargs: Passed to anaget_phbst_and_phdos_files
        """
        ddb = DdbFile.as_ddb(ddb)
        phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(**kwargs)
        return cls(ebands_kpath, phbst_file, phdos_file, ebands_kmesh=ebands_kmesh)

    def __init__(self, ebands_kpath, phbst_file, phdos_file, ebands_kmesh=None):
        self.eb_kpath = ElectronBands.as_ebands(ebands_kpath)
        self.eb_kmesh = ElectronBabds.as_ebands(ebands_kmesh) if ebands_kmesh is not None else None

        self.phbst_file = phbst_file
        if duck.is_string(self.phbst_file):
            self.phbst_file = PhbstFile(self.phbst_file)
        self.phb_qpath = self.phbst_file.phbands

        self.phdos_file = phdos_file
        if duck.is_string(self.phdos_file):
            self.phdos_file = PhdosFile(phdos_file)

    @add_fig_kwargs
    def plot(self, eb_ylims=None, **kwargs):
        """
        Plot electrons with possible (phonon-mediated) scattering channels for the CBM and VBM.
        Also plot phonon band structure and phonon PJDOS.

        Args:
            eb_ylims: Set the data limits for the y-axis of the electron band. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If None, limits are selected automatically.

        Return: |matplotlib-Figure|
        """
        # Build grid. Share y-axis for Phbands and Phdos
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax0 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
        ax1 = plt.subplot2grid((3, 3), (2, 0), colspan=2, rowspan=1)
        ax2 = plt.subplot2grid((3, 3), (2, 2), colspan=1, rowspan=1)
        ax1.get_shared_y_axes().join(ax1, ax2)

        # Plot electrons with possible e-ph scattering channels.
        self.eb_kpath.plot(ax=ax0, with_gaps=True, ylims=eb_ylims, max_phfreq=self.phb_qpath.maxfreq, show=False)

        # Plot phonon bands
        self.phb_qpath.plot(ax=ax1, show=False)
        #ax1.yaxis.set_visible(False)
        #set_visible(ax1, False, "ylabel")

        # Plot phonon PJDOS
        self.phdos_file.plot_pjdos_type(ax=ax2, fontsize=8, exchange_xy=True, show=False)
        set_visible(ax2, False, "ylabel")
        ax2.tick_params("y", left=False, labelleft=False)
        ax2.tick_params("y", right=True, labelright=True)

        # Adjust y-limits for phonons
        ylims = self.phb_qpath.minfreq, self.phb_qpath.maxfreq + 0.1 * abs(self.phb_qpath.maxfreq)
        for ax in (ax1, ax2):
            set_axlims(ax, ylims, "y")

        return fig

    @add_fig_kwargs
    def plot_phonons_occ(self, temps=(100, 200, 300, 400), **kwargs):
        """
        Plot phonon band structure with markers proportional to the occupation
        of each phonon mode for different temperatures.

        Args:
            temps: List of temperatures in Kelvin.

        Return: |matplotlib-Figure|
        """
        temps = np.array(temps)
        ntemp = len(temps)

        # Build plot grid.
        num_plots, ncols, nrows = ntemp, 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        for ax, temp in zip(ax_list, temps.ravel()):
            self.phb_qpath.plot(ax=ax, units="eV", temp=temp, fontsize=8, show=False)

        return fig

    @add_fig_kwargs
    def plot_linewidths_sigeph(self, sigeph, eb_ylims=None, **kwargs):
        """
        Plot e-bands + e-DOS + Im(Sigma_{eph}) + phonons + gkq^2

        Args:
            sigeph: |SigephFile| or string with path to file.
            eb_ylims: Set the data limits for the y-axis of the electron band. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If None, limits are selected automatically.
        """
        closeit = False
        if duck.is_string(sigeph):
            sigeph = SigEPhFile.from_file(sigeph)
            closeit = True

        # Build grid. share y-axis for Phbands and Phdos
        import matplotlib.pyplot as plt
        fig = plt.figure()

        # Electrons
        ax0 = plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=1)
        ax1 = plt.subplot2grid((2, 4), (0, 2), colspan=1, rowspan=1)
        ax2 = plt.subplot2grid((2, 4), (0, 3), colspan=1, rowspan=1)
        # Share y-axis
        ax_share("y", ax0, ax1, ax2)

        # Phonons
        ax3 = plt.subplot2grid((2, 4), (1, 0), colspan=2, rowspan=1)
        ax4 = plt.subplot2grid((2, 4), (1, 2), colspan=1, rowspan=1)
        ax5 = plt.subplot2grid((2, 4), (1, 3), colspan=1, rowspan=1)
        # Share y-axis
        ax_share("y", ax3, ax4, ax5)

        e0 = "fermie"

        # Plot electrons with possible e-ph scattering channels.
        self.eb_kpath.plot(ax=ax0, e0=e0, with_gaps=True, ylims=eb_ylims, max_phfreq=self.phb_qpath.maxfreq, show=False)
        sigeph.plot_lws_vs_e0(ax=ax1, e0=e0, exchange_xy=True, show=False)
        sigeph.edos.plot(ax=ax2, e0=e0, exchange_xy=True, show=False)

        # Plot phonon bands
        self.phb_qpath.plot(ax=ax3, show=False)
        sigeph.plot_a2fw_skb_sum(ax=ax4, what="gkq2", exchange_xy=True, fontsize=8, show=False)
        # Plot phonon PJDOS
        self.phdos_file.plot_pjdos_type(ax=ax5, fontsize=8, exchange_xy=True, show=False)
        #set_visible(ax4, False, "ylabel")
        #ax4.tick_params("y", left=False, labelleft=False)
        #ax4.tick_params("y", right=True, labelright=True)

        if closeit: sigeph.close()

        return fig

    #def close(self):
    #    self.phbst_file.close()
    #    self.phdos_file.close()


#class EphMultiPlotter(object):

# coding: utf-8
"""Classes for the analysis of electronic structures."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
#import pymatgen.core.units as units

from collections import OrderedDict
#from monty.string import is_string
#from monty.functools import lazy_property
from pymatgen.core.periodic_table import Element
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands
from abipy.electrons.ebands import ElectronsReader
from abipy.tools import gaussian

import logging
logger = logging.getLogger(__name__)


class FatBandsFile(AbinitNcFile, Has_Structure, Has_ElectronBands):
    """
    File with  ...

    Usage example:

    .. code-block:: python

        with FatBandsFile("foo_GSR.nc") as fb:
            fb.plot()
    """
    l2color = {0: "red", 1: "blue", 2: "green", 3: "yellow"}

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(FatBandsFile, self).__init__(filepath)
        self.reader = r = ElectronsReader(filepath)

        # Initialize the electron bands from file
        self._ebands = ebands = r.read_ebands()
        self.natom = len(self.structure)
        #self.nsppol, self.mband, self.nkpt = ebands.nsppol, ebands.mband, self.ebands.nkpt
        self.nkpt = self.ebands.nkpt # FIXME

        # Read metadata so that we know how to handle the content of the file.
        self.ecut = r.read_value("kinetic_energy_cutoff")
        self.usepaw = r.read_value("usepaw")
        self.prtdos = r.read_value("prtdos")
        self.prtdosm = r.read_value("prtdosm")
        self.pawprtdos = r.read_value("pawprtdos")
        self.pawfatbnd = r.read_value("pawfatbnd")
        self.natsph = r.read_dimvalue("natsph")
        self.iatsph = r.read_value("iatsph") - 1 # F --> C
        self.ndosfraction = r.read_dimvalue("ndosfraction")
        #self.ratsph = r.read_value("ratsph")
        #self.natsph_extra = r.read_dimvalue("natsph_extra")
        #self.ratsph_extra = r.read_value("ratsph_extra")
        #self.xredsph_extra = r.read_value("xredsph_extra")
        self.mbesslang = r.read_dimvalue("mbesslang")

        # pawtab_l_size(ntypat): Maximum value of l+1 leading to non zero Gaunt coeffs: l_size=2*l_max+1
        self.pawtab_l_size = r.read_value("pawtab_l_size")
        typat = r.read_value("atom_species") - 1 # F --> C
        self.lmax_atom = np.empty(self.natom, dtype=np.int)
        for iat in range(self.natom):
            self.lmax_atom[iat] = (self.pawtab_l_size[typat[iat]] - 1) // 2
        self.lsize = self.lmax_atom.max() + 1

        self.symbols = sorted(self.structure.symbol_set, key=lambda s: Element[s].Z)
        self.symbol2indices, self.lmax_symbol = OrderedDict(), OrderedDict()
        for symbol in self.symbols:
            self.symbol2indices[symbol] = np.array(self.structure.indices_from_symbol(symbol))
            self.lmax_symbol[symbol] = self.lmax_atom[self.symbol2indices[symbol][0]]

        self.symbol2color = {}
        if len(self.symbols) < 5:
            for i, symb in enumerate(self.symbols):
                self.symbol2color[symb] = self.l2color[i]
        else:
            # color will now be an RGBA tuple
            import matplotlib.pyplot as plt
            cm = plt.get_cmap('jet')
            nsymb = len(self.symbols)
            for i, symb in enumerate(self.symbols):
                self.symbol2color[symb] = cm(i/nsymb)

        #lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
        #dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang)
        #integer,intent(in) :: m_dos_flag,mbesslang,ndosfraction,pawfatbnd
        #dos_fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang*m_dos_flag)
        #             = m-resolved projected dos inside PAW sphere.

        #dtset
        #pawfatbnd = keyword for fatbands
        #mbesslang = maximum angular momentum for Bessel function expansion
        #m_dos_flag = option for the m-contributions to the partial DOS
        #ndosfraction = natsph*mbesslang

        # LDA+U params.
        #self.usepawu = r.read_value("usepawu")
        #self.lpawu = r.read_value("lpawu")
        #self.upawu = r.read_value("upawu")
        #self.jpawu = r.read_value("jpawu")

        # Array dimensiones as natom. Set to true if iatom has been calculated
        #self.has_atom = np.zeros(self.natom, dtype=bool)
        #self.has_atom[self.iatsph] = True

        # Read dos_fraction_m and reshape
        # [ndosfraction * mbesslang, nsppol, mband, nkpt]
        dos_fractions = self.reader.read_value("dos_fractions_m")

        if len(self.iatsph) == self.natom and np.all(self.iatsph == np.arange(self.natom)):
            self.walm_sbk = np.reshape(dos_fractions,
                                      (self.natom, self.mbesslang**2, self.nsppol, self.mband, self.nkpt))
        else:
            raise NotImplementedError("Should transfer data.")

    @property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self._ebands

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    def close(self):
        return self.reader.close()

    def __str__(self):
        lines = []
        app = lines.append
        app("lmax_symbol = %s" % str(self.lmax_symbol))
        return "\n".join(lines)

    def wl_atom(self, iatom, spin=None, band=None):
        if spin is None and band is None:
            wl = np.zeros((self.lsize, self.nsppol, self.mband, self.nkpt))
            for l in range(self.lmax_atom[iatom]+1):
                padl = l**2
                for m in range(2*l + 1):
                    wl[l] += self.walm_sbk[iatom, padl + m]
        else:
            assert spin is not None and band is not None
            wl = np.zeros((self.lsize, self.nkpt))
            for l in range(self.lmax_atom[iatom]+1):
                padl = l**2
                for m in range(2*l + 1):
                    wl[l] += self.walm_sbk[iatom, padl + m, spin, band, :]

        return wl

    def wl_symbol(self, symbol, spin=None, band=None):

        if spin is None and band is None:
            wl = np.zeros((self.lsize, self.nsppol, self.mband, self.nkpt))
            for iat in self.symbol2indices[symbol]:
                for l in range(self.lmax_atom[iat]+1):
                    padl = l**2
                    for m in range(2*l + 1):
                        wl[l] += self.walm_sbk[iat, padl + m]
        else:
            assert spin is not None and band is not None
            wl = np.zeros((self.lsize, self.nkpt))
            for iat in self.symbol2indices[symbol]:
                for l in range(self.lmax_atom[iat]+1):
                    padl = l**2
                    for m in range(2*l + 1):
                        wl[l, :] += self.walm_sbk[iat, padl + m, spin, band, :]

        return wl

    @add_fig_kwargs
    def plot_fatbands(self, view="all", e0="fermie", fact=3.0, alpha=0.9, **kwargs):
        """
        Plot the electronic fatbands.

        Args:
            view: "all", "types", "ineq"
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            fact:  float used to scale the stripe size.
            alpha

        Returns:
            `matplotlib` figure
        """
        # Build grid of plots.
        num_plots = self.natom
        ncols, nrows = 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)

        for iat, ax in enumerate(ax_list):
            # Plot the energies.
            ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
            ebands.decorate_ax(ax, title=str(self.structure[iat]))

            # Add width around each band.
            for spin in ebands.spins:
                for band in range(ebands.mband):
                    wlk = self.wl_atom(iat, spin=spin, band=band) * fact
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup.copy()
                    for l in range(self.lmax_atom[iat]+1):
                        w = wlk[l,:] / 2
                        y1, y2 = yup + w, ydown - w
                        ax.fill_between(x, yup, y1, alpha=alpha, facecolor=self.l2color[l])
                        ax.fill_between(x, ydown, y2, alpha=alpha, facecolor=self.l2color[l])
                        yup, ydown = y1, y2

        return fig

    @add_fig_kwargs
    def plot_fatbands_lview(self, view="all", e0="fermie", fact=3.0, alpha=0.9, **kwargs):
        """
        """
        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)

        import matplotlib.pyplot as plt
        ncols = self.lsize
        fig, ax_list = plt.subplots(nrows=1, ncols=ncols, sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        for l, ax in enumerate(ax_list):
            ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
            ebands.decorate_ax(ax, title="l=%s" % l)

            for spin in range(self.nsppol):
                for band in range(self.mband):
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup.copy()
                    for symbol in self.symbols:
                        wlk = self.wl_symbol(symbol, spin=spin, band=band)
                        w = wlk[l] / 2
                        y1, y2 = yup + w, ydown - w
                        # Add width around each band.
                        ax.fill_between(x, yup, y1, alpha=alpha, facecolor=self.symbol2color[symbol])
                        ax.fill_between(x, ydown, y2, alpha=alpha, facecolor=self.symbol2color[symbol])
                        yup, ydown = y1, y2

        return fig

    def get_edos_pjdos(self, method="gaussian", step=0.1, width=0.2):
        edos = self.ebands.get_edos(method=method, step=step, width=width)
        mesh = edos.spin_dos[0].mesh

        # Compute l-decomposed PJ DOS for each type of atom.
        if method == "gaussian":
            pjdos_symbls = OrderedDict()
            for symbol in self.symbols:
                lmax = self.lmax_symbol[symbol]
                wlsbk = self.wl_symbol(symbol)
                lso = np.zeros((self.lsize, self.nsppol, len(mesh)))
                for spin in range(self.nsppol):
                    for k, kpoint in enumerate(ebands.kpoints):
                        weight = kpoint.weight
                        for band in range(ebands.nband_sk[spin, k]):
                            e = ebands.eigens[spin,k,band]
                            for l in range(lmax+1):
                                lso[l, spin] += wlsbk[l, spin, band, k] * weight * gaussian(mesh, width, center=e)
                pjdos_symbls[symbol] = lso

        else:
            raise ValueError("Method %s is not supported" % method)

        return edos, pjdos_symbls

    @add_fig_kwargs
    def plot_pjdos(self, e0="fermie", what="l", method="gaussian", step=0.1, width=0.2, **kwargs):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns:
            :class:`ElectronDos` object.
        """
        edos, dos_symbls = self.get_edos_pjdos(self, method=method, step=step, width=width)
        mesh = edos.spin_dos[0].mesh
        # Define the zero of energy.
        e0 = ebands.get_e0(e0)
        mesh -= e0

        # Plot data.
        import matplotlib.pyplot as plt
        ncols = self.lsize
        fig, ax_list = plt.subplots(nrows=1, ncols=ncols, sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        for symbol in self.symbols:
            lso = dos_symbls[symbol]
            for spin in self.ebands.spins:
                spin_sign = +1 if spin == 0 else -1
                for l in range(self.lmax_symbol[symbol]+1):
                    tot_line = ax_list[l].plot(mesh, spin_sign * edos.spin_dos[spin].values, color="k", label="Total")
                    line_l, = ax_list[l].plot(mesh, spin_sign * lso[l, spin], color=self.symbol2color[symbol], label=symbol)

        for l, ax in enumerate(ax_list):
            ax.grid(True)
            ax.set_title("L = %d" % l)
            ax.set_xlabel("Energy [eV]")
            if l == 0:
                ax.legend(loc="best")
                ax.set_ylabel('DOS [states/eV]')

        return fig

    @add_fig_kwargs
    def plot_fatbands_with_pjdos(self, e0="fermie", what="l", method="gaussian", step=0.1, width=0.2, **kwargs):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns:
            :class:`ElectronDos` object.
        """

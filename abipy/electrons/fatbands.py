# coding: utf-8
"""Classes for the analysis of fatbands and PJDOS."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict, defaultdict
from tabulate import tabulate
from pymatgen.core.periodic_table import Element
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons.ebands import ElectronsReader
from abipy.tools import gaussian


class FatBandsFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    Provides methods to analyze the data stored in the FATBANDS.nc file.

    Usage example:

    .. code-block:: python

        with FatBandsFile("foo_FATBANDS.nc") as fb:
            fb.plot_fatbands_lview()

    Alternatively, one can use::

        with abiopen("foo_FATBANDS.nc") as fb:
            fb.plot_fatbands_lview()
    """
    # Mapping L --> color used in plots.
    l2color = {0: "red", 1: "blue", 2: "green", 3: "yellow"}

    # Markers used for up/down bands.
    marker_spin = {0: "^", 1: "v"}

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
        # TODO: nspden, nspinor,

        # Read metadata so that we know how to handle the content of the file.
        self.ecut = r.read_value("kinetic_energy_cutoff")
        self.usepaw = r.read_value("usepaw")
        self.prtdos = r.read_value("prtdos")
        self.prtdosm = r.read_value("prtdosm")
        self.pawprtdos = r.read_value("pawprtdos")
        #self.pawfatbnd = r.read_value("pawfatbnd")
        self.natsph = r.read_dimvalue("natsph")
        self.iatsph = r.read_value("iatsph") - 1 # F --> C
        self.ndosfraction = r.read_dimvalue("ndosfraction")
        self.mbesslang = r.read_dimvalue("mbesslang")
        self.ratsph_type = r.read_value("ratsph")
        # natsph_extra is present only if we have an extra sphere.
        self.natsph_extra = r.read_dimvalue("natsph_extra", default=0)
        if self.natsph_extra:
            self.ratsph_extra = r.read_value("ratsph_extra")
            self.xredsph_extra = r.read_value("xredsph_extra")
        if self.natsph_extra != 0:
            raise NotImplementedError("natsph_extra is not implemented, "
              "but it's just a matter of using natom + natsph_extra")

        # This is a tricky part. Note the following:
        # If usepaw == 0, lmax_type gives the max l included in the non-local part of Vnl
        #   The wavefunction can have l-components > lmax_type
        # If usepaw == 1, lmax_type represents the max l included in the PAW basis set.
        #   The AE wavefunction cannot have more ls than l-max if pawprtdos == 2 and
        #   the cancellation between PS-onsite and FFT part is perfect
        # TODO: Decide how to change lmax_type at run-time: API or global self.set_lmax?
        self.lmax_type = r.read_value("lmax_type")
        if not self.usepaw or self.pawprtdos != 0:
            self.lmax_type[:] = self.mbesslang - 1

        self.typat = r.read_value("atom_species") - 1 # F --> C
        self.lmax_atom = np.empty(self.natom, dtype=np.int)
        for iat in range(self.natom):
            self.lmax_atom[iat] = self.lmax_type[self.typat[iat]]
        # lsize is used to dimension arrays that depend on L.
        self.lsize = self.lmax_type.max() + 1

        # Sort the chemical symbols and use OrderedDict because we are gonna use these dicts for looping.
        self.symbols = sorted(self.structure.symbol_set, key=lambda s: Element[s].Z)
        self.symbol2indices, self.lmax_symbol = OrderedDict(), OrderedDict()
        for symbol in self.symbols:
            self.symbol2indices[symbol] = np.array(self.structure.indices_from_symbol(symbol))
            self.lmax_symbol[symbol] = self.lmax_atom[self.symbol2indices[symbol][0]]
        self.ntypat = len(self.symbols)

        # Mapping chemical symbol --> color used in plots.
        self.symbol2color = {}
        if len(self.symbols) < 5:
            for i, symb in enumerate(self.symbols):
                self.symbol2color[symb] = self.l2color[i]
        else:
            # Use colormap. Color will now be an RGBA tuple
            import matplotlib.pyplot as plt
            cm = plt.get_cmap('jet')
            nsymb = len(self.symbols)
            for i, symb in enumerate(self.symbols):
                self.symbol2color[symb] = cm(i/nsymb)

        # Array dimensiones as natom. Set to true if iatom has been calculated
        self.has_atom = np.zeros(self.natom, dtype=bool)
        self.has_atom[self.iatsph] = True

        #lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
        #dos_fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang*m_dos_flag)
        #             = m-resolved projected dos inside PAW sphere.
        #mbesslang = maximum angular momentum for Bessel function expansion

        # Read dos_fraction_m from file and build walm_sbk array of shape [natom, lmax**2, nsppol, mband, nkpt].
        # Note that Abinit allows the users to select a subset of atoms with iatsph. Moreover the order
        # of the atoms could differ from the one in the structure even when natom == natsph (unlikely but possible).
        # To keep it simple, the code always operate on an array dimensioned with the total number of atoms
        # Entries that are not computed are set to zero and a warning is issued.
        wshape = (self.natom, self.mbesslang**2, self.nsppol, self.mband, self.nkpt)

        if self.natsph == self.natom and np.all(self.iatsph == np.arange(self.natom)):
            self.walm_sbk = np.reshape(r.read_value("dos_fractions_m"), wshape)

        else:
            # Need to tranfer data. Note np.zeros.
            self.walm_sbk = np.zeros(wshape)
            if self.natsph == self.natom and np.any(self.iatsph != np.arange(self.natom)):
                print("Will rearrange filedata since iatsp != [1, 2, ...])")
                filedata = np.reshape(r.read_value("dos_fractions_m"), wshape)
                for i, iatom in enumerate(self.iatsph):
                    self.walm_sbk[iatom] = filedata[i]

            else:
                print("natsph < natom. Will set to zero the PJDOS contributions for the atoms that are not included.")
                assert self.natsph < self.natom
                filedata = np.reshape(r.read_value("dos_fractions_m"),
                                     (self.natsph, self.mbesslang**2, self.nsppol, self.mband, self.nkpt))
                for i, iatom in enumerate(self.iatsph):
                    self.walm_sbk[iatom] = filedata[i]

        # In principle, this should never happen (unless there's a bug in Abinit or a
        # very bad cancellation between the FFT and the PS-PAW term (pawprtden=0).
        num_neg = np.sum(self.walm_sbk < 0)
        if num_neg:
            print("WARNING: There are %d (%.1f%%) negative entries in LDOS weights" % (
                  num_neg, 100 * num_neg / self.walm_sbk.size))

    @property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self._ebands

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    def close(self):
        """Called at the end of the `with` context manager."""
        return self.reader.close()

    def __str__(self):
        """String representation"""
        lines = []
        app = lines.append
        app("usepaw=%d, prtdos=%d, pawprtdos=%d, prtdosm=%d, mbesslang=%d" % (
            self.usepaw, self.prtdos, self.pawprtdos, self.prtdosm, self.mbesslang))
        app("nsppol=%d, nkpt=%d, mband=%d" % (self.nsppol, self.nkpt, self.mband))
        app("")

        table = [["Idx", "Symbol", "Reduced_Coords", "Lmax", "Ratsph [Bohr]", "Has_Atom"]]
        for iatom, site in enumerate(self.structure):
            table.append([
                iatom,
                site.specie.symbol,
                "%.5f %.5f %.5f" % tuple(site.frac_coords),
                self.lmax_atom[iatom],
                self.ratsph_type[self.typat[iatom]],
                bool(self.has_atom[iatom]),
            ])

        app(tabulate(table, headers="firstrow"))

        return "\n".join(lines)

    def get_wl_atom(self, iatom, spin=None, band=None):
        """
        Return the l-dependent DOS weights for atom index `iatom`. The weights are summed over m.
        If `spin` and `band` are not specified, the method returns the weights
        for all spin and bands else the contribution for (spin, band)
        """
        if spin is None and band is None:
            wl = np.zeros((self.lsize, self.nsppol, self.mband, self.nkpt))
            for l in range(self.lmax_atom[iatom]+1):
                for m in range(2*l + 1):
                    wl[l] += self.walm_sbk[iatom, l**2 + m]
        else:
            assert spin is not None and band is not None
            wl = np.zeros((self.lsize, self.nkpt))
            for l in range(self.lmax_atom[iatom]+1):
                for m in range(2*l + 1):
                    wl[l] += self.walm_sbk[iatom, l**2 + m, spin, band, :]

        return wl

    def get_wl_symbol(self, symbol, spin=None, band=None):
        """
        Return the l-dependent DOS weights for a given type specified in terms of the
        chemical symbol `symbol`. The weights are summed over m and over all atoms of the same type.
        If `spin` and `band` are not specified, the method returns the weights
        for all spins and bands else the contribution for (spin, band).
        """
        if spin is None and band is None:
            wl = np.zeros((self.lsize, self.nsppol, self.mband, self.nkpt))
            for iat in self.symbol2indices[symbol]:
                for l in range(self.lmax_atom[iat]+1):
                    for m in range(2*l + 1):
                        wl[l] += self.walm_sbk[iat, l**2 + m]
        else:
            assert spin is not None and band is not None
            wl = np.zeros((self.lsize, self.nkpt))
            for iat in self.symbol2indices[symbol]:
                for l in range(self.lmax_atom[iat]+1):
                    for m in range(2*l + 1):
                        wl[l, :] += self.walm_sbk[iat, l**2 + m, spin, band, :]

        return wl

    def get_w_symbol(self, symbol, spin=None, band=None):
        """
        Return the DOS weights for a given type specified in terms of the
        chemical symbol `symbol`. The weights are summed over m and lmax[symbol] and over all atoms of the same type.
        If `spin` and `band` are not specified, the method returns the weights
        for all spins and bands else the contribution for (spin, band).
        """
        if spin is None and band is None:
            wl = self.get_wl_symbol(symbol)
            w = np.zeros((self.nsppol, self.mband, self.nkpt))
            for l in range(self.lmax_symbol[symbol]+1):
                w += wl[l]

        else:
            assert spin is not None and band is not None
            wl = self.get_wl_symbol(symbol, spin=spin, band=spin)
            w = np.zeros((self.nkpt))
            for l in range(self.lmax_symbol[symbol]+1):
                w += wl[l]

        return w

    def get_spilling(self, spin=None, band=None):
        """
        """
        if spin is None and band is None:
            sp = np.zeros((self.nsppol, self.mband, self.nkpt))
            for iatom in range(self.natom):
                for l in range(self.lmax_atom[iatom]+1):
                    for m in range(2*l + 1):
                        sp += self.walm_sbk[iatom, l**2 + m]
        else:
            assert spin is not None and band is not None
            sp = np.zeros((self.nkpt))
            for iatom in range(self.natom):
                for l in range(self.lmax_atom[iatom]+1):
                    for m in range(2*l + 1):
                        sp += self.walm_sbk[iatom, l**2 + m, spin, band, :]

        return 1.0 - sp

    @add_fig_kwargs
    def plot_fatbands_siteview(self, e0="fermie", view="inequivalent", fact=2.0, alpha=0.6, **kwargs):
        """
        Plot fatbands for each atom in the unit cell. By default, only the "inequivalent" atoms are shown.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            view: "inequivalent", "all"
            fact:  float used to scale the stripe size.
            alpha:

        Returns:
            `matplotlib` figure
        """
        # Define num_plots and ax2atom depending on view.
        # ax2natom[1:num_plots] --> iatom index in structure.
        # TODO
        #if view == "inequivalent" and self.nspden == 2 and self.nsppol == 1:
        #    print("The system contains magnetic symmetries but the spglib API used does not handle them.")
        #    print("Setting view to `all`")
        #    view = "all"

        if view == "all" or self.natom == 1:
            num_plots, ax2iatom = self.natom, np.arange(self.natom)

        elif view == "inequivalent":
            print("Calling spglib to find inequivalent sites. Note: magnetic symmetries are not taken into account")
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            spgan = SpacegroupAnalyzer(self.structure)
            spgdata = spgan.get_symmetry_dataset()
            equivalent_atoms = spgdata["equivalent_atoms"]
            ax2iatom = []
            eqmap = defaultdict(list)
            for pos, eqpos in enumerate(equivalent_atoms):
                if pos == eqpos: ax2iatom.append(pos)
                eqmap[eqpos].append(pos)
            ax2iatom = np.array(ax2iatom)
            num_plots = len(ax2iatom)
            print("Found %d inequivalent position(s)." % num_plots)
            for i, irr_pos in enumerate(sorted(eqmap.keys())):
                print("Irred_Site: %s" % str(self.structure[irr_pos]))
                for eqind in eqmap[irr_pos]:
                    if eqind == irr_pos: continue
                    print("\tSymEq: %s" % str(self.structure[eqind]))
        else:
            raise ValueError("Wrong value for view: %s" % str(view))

        # Build grid of plots.
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

        for i, ax in enumerate(ax_list):
            iat = ax2iatom[i]
            # Plot the energies.
            # ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
            for spin in range(self.nsppol):
                ebands.plot_ax(ax, e0, spin=spin,
                              color="grey", linewidth=0.2, markersize=2.0, marker=self.marker_spin[spin])

            ebands.decorate_ax(ax, title=str(self.structure[iat]))

            # Add width around each band.
            for spin in ebands.spins:
                for band in range(ebands.mband):
                    wlk = self.get_wl_atom(iat, spin=spin, band=band) * fact
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup.copy()
                    for l in range(self.lmax_atom[iat]+1):
                        w = wlk[l,:] / 2
                        y1, y2 = yup + w, ydown - w
                        ax.fill_between(x, yup, y1, alpha=alpha, facecolor=self.l2color[l])
                        ax.fill_between(x, ydown, y2, alpha=alpha, facecolor=self.l2color[l],
                                        label="L=%d" % l if (i, spin, band) == (0, 0, 0) else None)
                                        # Note: could miss a label in the other plots if lmax is not large enough!
                        yup, ydown = y1, y2

        ax_list[0].legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_fatbands_lview(self, e0="fermie", fact=2.0, alpha=0.6, ax_list=None, lmax=None, **kwargs):
        """
        Plot the electronic fatbands grouped by l.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            fact:  float used to scale the stripe size.
            alpha:
            ax_list:

        Returns:
            `matplotlib` figure
        """
        mylsize = self.lsize if lmax is None else lmax + 1
        # Get ax_list and fig.
        import matplotlib.pyplot as plt
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=1, ncols=mylsize, sharex=True, sharey=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            if len(ax_list) != mylsize:
                raise ValueError("len(ax_list) != self.lsize")
            fig = plt.gcf()

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)

        for l, ax in enumerate(ax_list):
            # ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
            for spin in range(self.nsppol):
                ebands.plot_ax(ax, e0, spin=spin,
                              color="grey", linewidth=0.2, markersize=2.0, marker=self.marker_spin[spin])

            ebands.decorate_ax(ax, title="L = %s" % l)
            if l != 0:
                ax.set_ylabel("")

            for spin in range(self.nsppol):
                for band in range(self.mband):
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup.copy()
                    for symbol in self.symbols:
                        wlk = self.get_wl_symbol(symbol, spin=spin, band=band) * fact
                        w = wlk[l] / 2
                        y1, y2 = yup + w, ydown - w
                        # Add width around each band.
                        ax.fill_between(x, yup, y1, alpha=alpha, facecolor=self.symbol2color[symbol])
                        ax.fill_between(x, ydown, y2, alpha=alpha, facecolor=self.symbol2color[symbol],
                                        label=symbol if (l, spin, band) == (0, 0, 0) else None)
                        yup, ydown = y1, y2

        ax_list[0].legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_fatbands_mview(self, iatom, e0="fermie", fact=6.0, alpha=0.6, lmax=None, **kwargs):
        """
        Plot the electronic fatbands grouped by l.

        Args:
            iatom: Index of the atom in the structure.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            fact:  float used to scale the stripe size.
            alpha:

        Returns:
            `matplotlib` figure
        """
        mylmax = self.lmax_atom[iatom] if lmax is None else lmax

        # Build plot grid.
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure()
        nrows, ncols = 2 * (mylmax+1), mylmax + 1
        gspec = GridSpec(nrows=nrows, ncols=ncols)
        gspec.update(wspace=0.1, hspace=0.1)

        ax_lim = {}
        for im in range(nrows):
            for l in range(ncols):
                k = (l, im)
                ax00 = None if l == 0 else ax_lim[(0, 0)]
                ax = plt.subplot(gspec[im, l], sharex=ax00, sharey=ax00)
                if im < 2*l + 1:
                    #ax.set_title("l=%d, m=%d" % (l, im - l))
                    ax_lim[k] = ax
                else:
                    ax.axis("off")

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)

        for lim, ax in ax_lim.items():
            l, im = lim[0], lim[1]
            # ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
            for spin in range(self.nsppol):
                ebands.plot_ax(ax, e0, spin=spin,
                              color="grey", linewidth=0.2, markersize=2.0, marker=self.marker_spin[spin])

            if im == 2 * l:
               ebands.decorate_ax(ax)
            #if l > 0:
            #    ax.set_ylabel("")

            for spin in range(self.nsppol):
                for band in range(self.mband):
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup.copy()

                    w = self.walm_sbk[iatom,  l**2 + im, spin, band, :] * fact / 2
                    y1, y2 = yup + w, ydown - w
                    # Add width around each band.
                    ax.fill_between(x, yup, y1, alpha=alpha, facecolor=self.l2color[l])
                    ax.fill_between(x, ydown, y2, alpha=alpha, facecolor=self.l2color[l])
                    yup, ydown = y1, y2

        #fig.tight_layout()
        return fig

    @add_fig_kwargs
    def plot_fatbands_typeview(self, e0="fermie", fact=2.0, alpha=0.6, ax_list=None, **kwargs):
        """
        Plot the electronic fatbands

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            fact:  float used to scale the stripe size.
            alpha:
            ax_list:

        Returns:
            `matplotlib` figure
        """
        # Get ax_list and fig.
        import matplotlib.pyplot as plt
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=1, ncols=self.ntypat, sharex=True, sharey=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            if len(ax_list) != self.ntypat:
                raise ValueError("len(ax_list) != self.ntypat")
            fig = plt.gcf()

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)

        for i, (symbol, ax) in enumerate(zip(self.symbols, ax_list)):
            # ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
            for spin in range(self.nsppol):
                ebands.plot_ax(ax, e0, spin=spin,
                              color="grey", linewidth=0.2, markersize=2.0, marker=self.marker_spin[spin])

            ebands.decorate_ax(ax, title=symbol)
            if i != 0:
                ax.set_ylabel("")

            w_sbk = self.get_w_symbol(symbol) * fact
            for spin in range(self.nsppol):
                for band in range(self.mband):
                    y = ebands.eigens[spin, :, band] - e0
                    w = w_sbk[spin, band, :] / 2
                    # Add width around each band.
                    ax.fill_between(x, y, y + w, alpha=alpha, facecolor=self.symbol2color[symbol])
                    ax.fill_between(x, y, y - w, alpha=alpha, facecolor=self.symbol2color[symbol])

        return fig

    @add_fig_kwargs
    def plot_spilling(self, e0="fermie", fact=5.0, alpha=0.6, ax=None, **kwargs):
        """
        Plot the electronic fatbands

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            fact:  float used to scale the stripe size.
            alpha:

        Returns:
            `matplotlib` figure
        """
        # Get ax_list and fig.
        import matplotlib.pyplot as plt
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)

        spill_sbk = self.get_spilling() * fact

        ebands.decorate_ax(ax)
        for spin in range(self.nsppol):
            ebands.plot_ax(ax, e0, spin=spin,
                           color="grey", linewidth=0.2, markersize=2.0, marker=self.marker_spin[spin])

            for band in range(self.mband):
                y = ebands.eigens[spin, :, band] - e0
                w = spill_sbk[spin, band, :] / 2

                # Handle negative spilling values.
                wispos = w >= 0.0
                wisneg = np.logical_not(wispos)
                num_neg = np.sum(wisneg)

                # Add width around each band.
                ax.fill_between(x, y, y + w, where=wispos, alpha=alpha, facecolor="blue")
                ax.fill_between(x, y, y - w, where=wispos, alpha=alpha, facecolor="blue")

                # Show regions with negative spilling in red.
                if num_neg:
                    print("For spin:", spin, "band:", band,
                          "There are %d (%.1f%%) k-points with negative spilling. Min: %.2E" % (
                           num_neg, 100 * num_neg / self.nkpt, w.min()))

                    absw = np.abs(w)
                    ax.fill_between(x, y, y + absw, where=wisneg, alpha=alpha, facecolor="red")
                    ax.fill_between(x, y, y - absw, where=wisneg, alpha=alpha, facecolor="red")

        return fig

    def nelect_in_spheres(self, energy=None, method="gaussian", step=0.1, width=0.2):
        """
        Print the number of electrons inside each atom-centered sphere.
        Note that this is a very crude estimate of the charge density distribution.

        Args:
            energy: PJDOS is integrated up to this energy [eV]. If None, the Fermi level is used.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
        """
        raise NotImplementedError("")
        if energy is None: energy = self.ebands.fermie
        #edos, pjdos_symbls, cumdos_ls = self._get_edos_pjdos_cumdos(method=method, step=step, width=width)

        for iatm, site in enumerate(self.structure):
            # ?? site_dos: [natom, spin, nomega]
            if iatm == 0:
                # Find the index of the first point in the mesh whose value is >= value. -1 if not found
                stop_spin = {}
                for spin in self.spins:
                    stop = site_dos[spin].find_mesh_index(energy)
                    if stop == -1:
                        raise ValueError("For spin %d: cannot find index in mesh such that mesh[i] >= energy." % spin)
                    stop_spin[spin] = stop

            nel_spin = {}
            for spin in self.spins:
                nel_spin[spin] = site_dos[spin].integral(stop=stop_spin[spin])
            print("iatom", iatm, "site", site, nel_spin)

    def _get_edos_pjdos_cumdos(self, method="gaussian", step=0.1, width=0.2):
        ebands = self.ebands
        edos = ebands.get_edos(method=method, step=step, width=width)
        mesh = edos.spin_dos[0].mesh

        # Compute l-decomposed PJ DOS for each type of atom.
        if method == "gaussian":
            pjdos_symbls = OrderedDict()
            for symbol in self.symbols:
                lmax = self.lmax_symbol[symbol]
                wlsbk = self.get_wl_symbol(symbol)
                lso = np.zeros((self.lsize, self.nsppol, len(mesh)))
                for spin in range(self.nsppol):
                    for k, kpoint in enumerate(ebands.kpoints):
                        weight = kpoint.weight
                        for band in range(ebands.nband_sk[spin, k]):
                            e = ebands.eigens[spin, k, band]
                            for l in range(lmax+1):
                                lso[l, spin] += wlsbk[l, spin, band, k] * weight * gaussian(mesh, width, center=e)
                pjdos_symbls[symbol] = lso

        else:
            raise ValueError("Method %s is not supported" % method)

        cumdos_ls = None
        if True:
            # Compute cumdos_ls datastructure for stacked DOS.
            # cumdos_ls maps (l, spin) onto a numpy array [nsymbols, nfreqs] where
            # [isymb, :] contains the cumulative sum of the PJDOS up to symbol isymb.
            from itertools import product
            dls = defaultdict(dict)
            for symbol, lso in pjdos_symbls.items():
                for l, spin in product(range(self.lmax_symbol[symbol]+1), range(self.nsppol)):
                    dls[(l, spin)][symbol] = lso[l, spin]

            cumdos_ls = {}
            nsymb = len(self.symbols)
            for (l, spin), dvals in dls.items():
                arr = np.zeros((nsymb, len(mesh)))
                for isymb, symbol in enumerate(self.symbols):
                    arr[isymb] = dvals[symbol]
                cumdos_ls[(l, spin)] = np.cumsum(arr, axis=0)

        return edos, pjdos_symbls, cumdos_ls

    @add_fig_kwargs
    def plot_pjdos_lview(self, e0="fermie", method="gaussian", step=0.1, width=0.2, stacked=True, alpha=0.6,
                         ax_list=None, exchange_xy=False, **kwargs):
        """
        Plot the PJ-DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            stacked:
            alpha
            ax_list:
            exchange_xy: True if the dos should be plotted on the x axis instead of y.

        Returns:
            `matplotlib` figure
        """
        edos, pjdos_symbls, cumdos_ls = self._get_edos_pjdos_cumdos(method=method, step=step, width=width)

        # Get energy mesh from total DOS and define the zero of energy
        # Note that the mesh is not not spin-dependent.
        mesh = edos.spin_dos[0].mesh
        e0 = self.ebands.get_e0(e0)
        mesh -= e0

        # Plot data.
        import matplotlib.pyplot as plt
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=1, ncols=self.lsize, sharex=True, sharey=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            if len(ax_list) != self.lsize:
                raise ValueError("len(ax_list) != self.lsize")
            fig = plt.gcf()

        if not stacked:
            # Plot PJDOS as lines.
            for isymb, symbol in enumerate(self.symbols):
                for spin in range(self.nsppol):
                    spin_sign = +1 if spin == 0 else -1
                    # Loop over the columns of the grid.
                    for l in range(self.lmax_symbol[symbol]+1):
                        x, y = mesh, spin_sign * edos.spin_dos[spin].values
                        if exchange_xy: x, y = y, x
                        ax_list[l].plot(x, y, color="k", label="Tot" if (l, spin, isymb) == (0, 0, 0) else None)
                        x, y = mesh, spin_sign * pjdos_symbls[symbol][l, spin]
                        if exchange_xy: x, y = y, x
                        ax_list[l].plot(x, y, color=self.symbol2color[symbol],
                                        label=symbol if (l, spin, isymb) == (0, 0, 0) else None)
        else:
            # Plot stacked PJDOS
            # Loop over the columns of the grid.
            for l in range(self.lsize):
                for spin in self.ebands.spins:
                    spin_sign = +1 if spin == 0 else -1

                    # Plot total DOS.
                    x, y = mesh, spin_sign * edos.spin_dos[spin].values
                    if exchange_xy: x, y = y, x
                    ax_list[l].plot(x, y, color="k", label="Tot" if (l, spin) == (0, 0) else None)

                    # Plot cumulative PJ-DOS(l, spin)
                    cumdos = cumdos_ls[(l, spin)] * spin_sign
                    for isymb, symbol in enumerate(self.symbols):
                        yup = cumdos[isymb]
                        ydown = cumdos[isymb-1] if isymb != 0 else np.zeros(len(mesh))
                        fill = ax_list[l].fill_between if not exchange_xy else ax_list[l].fill_betweenx
                        fill(mesh, yup, ydown, alpha=alpha, facecolor=self.symbol2color[symbol],
                             label="%s (stacked)" % symbol if (l, spin) == (0, 0) else None)

        # Decorate axis.
        for l, ax in enumerate(ax_list):
            ax.set_title("L = %d" % l)
            ax.grid(True)

            if not exchange_xy:
                # Display yticklabels on the first plot and last plot only.
                # and display the legend only on the first plot.
                ax.set_xlabel("Energy [eV]")
                if l == 0:
                    ax.legend(loc="best")
                    ax.set_ylabel('DOS [states/eV]')
                elif l == self.lsize - 1:
                    ax.yaxis.set_ticks_position("right")
                    ax.yaxis.set_label_position("right")
                else:
                    # Plots in the middle: don't show labels.
                    # Trick: Don't change the labels but set their fontsize to 0 otherwise
                    # also the other axes are affecred (likely due to sharey=True).
                    #ax.set_yticklabels([])
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(0)
            else:
                if l == 0:
                    ax.legend(loc="best")
                    ax.set_xlabel('DOS [states/eV]')
                elif l == self.lsize - 1:
                    ax.yaxis.set_ticks_position("right")
                    ax.yaxis.set_label_position("right")
                else:
                    # See trick used above.
                    #ax.set_yticklabels([])
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(0)

        return fig

    @add_fig_kwargs
    def plot_fatbands_with_pjdos(self, e0="fermie", fact=2.0, alpha=0.6,
                                 pjdosfile=None, edos_kwargs=None, stacked=True, **kwargs):
        """
        Compute the fatbands and the PJDOS.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            fact:  float used to scale the stripe size.
            alpha:
            pjdosfile: FATBANDS file used to compute the PJDOS. If None, the PJDOS is taken from self.
            edos_kwargs
            stacked:

        Returns:
            `matplotlib` figure
        """
        closeit = False
        if pjdosfile is not None:
            if not isinstance(pjdosfile, FatBandsFile):
                # String --> open the file here and close it before returning.
                pjdosfile = FatBandsFile(pjdosfile)
                closeit = True
        else:
            # Compute PJDOS from self.
            pjdosfile = self

        # Build plot grid.
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure()
        gspec = GridSpec(nrows=1, ncols=self.lsize)
        fatbands_axlist, pjdos_axlist = [], []
        for l in range(self.lsize):
            subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gspec[l], width_ratios=[2, 1], wspace=0.05)
            # Get axes and align bands and DOS.
            prev_ax = None if l == 0 else fatbands_axlist[0]
            ax1 = plt.subplot(subgrid[0], sharex=prev_ax, sharey=prev_ax)
            ax2 = plt.subplot(subgrid[1], sharey=ax1)
            fatbands_axlist.append(ax1)
            pjdos_axlist.append(ax2)

        # Plot bands on fatbands_axlist
        self.plot_fatbands_lview(e0=e0, fact=fact, alpha=alpha, ax_list=fatbands_axlist, show=False)

        # Plot PJDOS on pjdos_axlist.
        if edos_kwargs is None: edos_kwargs = {}
        pjdosfile.plot_pjdos_lview(e0=e0, ax_list=pjdos_axlist, exchange_xy=True,
                                   stacked=stacked, show=False, **edos_kwargs)

        if closeit: pjdosfile.close()
        return fig

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        workind directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("fbfile = abilab.abiopen('%s')\nprint(fbfile)" % self.filepath),
            nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands_lview()"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands_mview(iatom=0)"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands_siteview()"),
            nbv.new_code_cell("fig = fbfile.plot_pjdos_lview()"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands_with_pjdos(pjdosfile=None)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

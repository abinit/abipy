# coding: utf-8
"""Classes for the analysis of electronic structures."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict, defaultdict
from pymatgen.core.periodic_table import Element
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons.ebands import ElectronsReader
from abipy.tools import gaussian


class FatBandsFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    File with  ...

    Usage example:

    .. code-block:: python

        with FatBandsFile("foo_FATBANDS.nc") as fb:
            fb.plot_fatbands_lview()
    """
    # Mapping L --> color used in plots.
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

        # Sort the chemical symbols and use OrderedDict because we are gonna use these dicts for looping.
        self.symbols = sorted(self.structure.symbol_set, key=lambda s: Element[s].Z)
        self.symbol2indices, self.lmax_symbol = OrderedDict(), OrderedDict()
        for symbol in self.symbols:
            self.symbol2indices[symbol] = np.array(self.structure.indices_from_symbol(symbol))
            self.lmax_symbol[symbol] = self.lmax_atom[self.symbol2indices[symbol][0]]

        # Mapping chemical symbol --> color used in plots.
        self.symbol2color = {}
        if len(self.symbols) < 5:
            for i, symb in enumerate(self.symbols):
                self.symbol2color[symb] = self.l2color[i]
        else:
            # Ues colormap. Color will now be an RGBA tuple
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

        # Read dos_fraction_m from file and build walm_sbk array of shape [natom, lmax**2, nsppol, mband, nkpt].
        # Note that Abinit allows the users to select a subset of atoms with iatsph. Moreover the order
        # of the atoms could differ from the one in the structure.
        # To keep it simple, the code always operate on an array dimensioned with the total number of atoms
        # Entries that are not computed are set to zero and a warning is issued.
        wshape = (self.natom, self.mbesslang**2, self.nsppol, self.mband, self.nkpt))

        if self.natsph == self.natom and np.all(self.iatsph == np.arange(self.natom)):
            self.walm_sbk = np.reshape(r.read_value("dos_fractions_m"), wshape)

        else:
            # Need to tranfer data. Note np.zeros.
            self.walm_sbk = np.zeros(wshape)
            if self.natsph == self.natom and np.any(self.iatsph != np.arange(self.natom)):
                print("Will rearrange filedata since iatsp != [1, 2, ...])")
                filedata = np.reshape(r.read_value("dos_fractions_m"), wshape)
                for i, iatom in enumerate(self.iatsph):
                    self.walm_sbl[iatom] = filedata[i]

            else
                print("natsph < natom. Will set to zero the PJDOS contributions for the atoms that are not included.")
                assert self.natsph < self.natom
                filedata = np.reshape(r.read_value("dos_fractions_m"),
                                     (self.natsph, self.mbesslang**2, self.nsppol, self.mband, self.nkpt))
                for i, iatom in enumerate(self.iatsph):
                    self.walm_sbl[iatom] = filedata[i]

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
        app("lmax_symbol = %s" % str(self.lmax_symbol))
        return "\n".join(lines)

    def wl_atom(self, iatom, spin=None, band=None):
        """
        Return the l-dependent DOS weights for atom index `iatom`. The weights are summed over m.
        If `spin` and `band` are not specified, the method returns the weights
        for all spin and bands else the contribution for (spin, band)
        """
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
        """
        Return the l-dependent DOS weights for a given type specified in terms of the
        chemical symbol `symbol`. The weights are summed over m and over all atoms of the same type.
        If `spin` and `band` are not specified, the method returns the weights
        for all spin and bands else the contribution for (spin, band).
        """
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
    def plot_fatbands(self, e0="fermie", fact=3.0, alpha=0.6, **kwargs):
        """
        Plot the electronic fatbands for all atoms in the unit cell.

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
    def plot_fatbands_lview(self, e0="fermie", fact=3.0, alpha=0.6, ax_list=None, **kwargs):
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
            fig, ax_list = plt.subplots(nrows=1, ncols=self.lsize, sharex=True, sharey=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            if len(ax_list) != self.lsize:
                raise ValueError("len(ax_list) != self.lsize")
            fig = plt.gcf()

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)

        for l, ax in enumerate(ax_list):
            ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
            ebands.decorate_ax(ax, title="L = %s" % l)
            if l != 0:
                ax.set_ylabel("")

            for spin in range(self.nsppol):
                for band in range(self.mband):
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup.copy()
                    for symbol in self.symbols:
                        wlk = self.wl_symbol(symbol, spin=spin, band=band) * fact
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
    def plot_fatbands_m(self, iatom, e0="fermie", fact=5.0, alpha=0.6, **kwargs):
        """
        """
        lmax = self.lmax_atom[iatom]
        # Build plot grid.
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure()
        nrows, ncols = 2 * (lmax+1), lmax + 1
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
            ebands.plot_ax(ax, e0, color="grey", linewidth=0.2, marker="o", markersize=2.0)
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

    def get_edos_pjdos(self, method="gaussian", step=0.1, width=0.2):
        ebands = self.ebands
        edos = ebands.get_edos(method=method, step=step, width=width)
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
    def plot_pjdos(self, e0="fermie", method="gaussian", step=0.1, width=0.2, stacked=True, alpha=0.6,
                   ax_list=None, exchange_xy=False, **kwargs):
        """
        Plot the PJ-DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            ax_list:
            exchange_xy: True if the dos should be plotted on the x axis insted of y.

        Returns:
            :class:`ElectronDos` object.
        """
        edos, dos_symbls = self.get_edos_pjdos(method=method, step=step, width=width)
        mesh = edos.spin_dos[0].mesh
        # Define the zero of energy.
        e0 = self.ebands.get_e0(e0)
        mesh -= e0

        # Plot data.
        import matplotlib.pyplot as plt
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=1, ncols=self.lsize, sharex=True, sharey=False, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            if len(ax_list) != self.lsize:
                raise ValueError("len(ax_list) != self.lsize")
            fig = plt.gcf()

        if not stacked:
            for isymb, symbol in enumerate(self.symbols):
                for spin in self.ebands.spins:
                    spin_sign = +1 if spin == 0 else -1
                    # Loop over the columns of the grid.
                    for l in range(self.lmax_symbol[symbol]+1):
                        x, y = mesh, spin_sign * edos.spin_dos[spin].values
                        if exchange_xy: x, y = y, x
                        ax_list[l].plot(x, y, color="k", label="Tot" if (l, spin, isymb) == (0, 0, 0) else None)
                        x, y = mesh, spin_sign * dos_symbls[symbol][l, spin]
                        if exchange_xy: x, y = y, x
                        ax_list[l].plot(x, y, color=self.symbol2color[symbol],
                                        label=symbol if (l, spin, isymb) == (0, 0, 0) else None)
        else:
            # Compute datastructure for stacked DOS.
            from itertools import product
            dls = defaultdict(dict)
            for symbol, lso in dos_symbls.items():
                for l, spin in product(range(self.lmax_symbol[symbol]+1), range(self.nsppol)):
                    dls[(l, spin)][symbol] = lso[l, spin]
            cumdos_ls = {}
            nsymb = len(self.symbols)
            for (l, spin), dvals in dls.items():
                arr = np.zeros((nsymb, len(mesh)))
                for isymb, symbol in enumerate(self.symbols):
                    arr[isymb] = dvals[symbol]
                cumdos_ls[(l, spin)] = np.cumsum(arr, axis=0)

            # Loop over the columns of the grid.
            for l in range(self.lsize):
                for spin in self.ebands.spins:
                    spin_sign = +1 if spin == 0 else -1

                    x, y = mesh, spin_sign * edos.spin_dos[spin].values
                    if exchange_xy: x, y = y, x
                    ax_list[l].plot(x, y, color="k", label="Tot" if (l, spin) == (0, 0) else None)

                    cumdos = cumdos_ls[(l, spin)] * spin_sign
                    for isymb, symbol in enumerate(self.symbols):
                        yup = cumdos[isymb]
                        ydown = cumdos[isymb-1] if isymb != 0 else np.zeros(len(mesh))
                        fill = ax_list[l].fill_between if not exchange_xy else ax_list[l].fill_betweenx
                        fill(mesh, yup, ydown, alpha=alpha, facecolor=self.symbol2color[symbol],
                             label="%s (stacked)" % symbol if (l, spin) == (0, 0) else None)

        # Decorate axis.
        for l, ax in enumerate(ax_list):
            ax.grid(True)
            ax.set_title("L = %d" % l)

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
                    ax.set_yticklabels([])

            else:
                if l == 0:
                    ax.legend(loc="best")
                    ax.set_xlabel('DOS [states/eV]')
                elif l == self.lsize - 1:
                    ax.yaxis.set_ticks_position("right")
                    ax.yaxis.set_label_position("right")
                else:
                    ax.set_yticklabels([])

        return fig

    @add_fig_kwargs
    def plot_fatbands_with_pjdos(self, e0="fermie", fact=3.0, alpha=0.6,
                                 ncfile=None, edos_kwargs=None, stacked=True, **kwargs):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            exchange_xy: True to exchange x-y axes.

        Returns:
        """
        closeit = False
        if ncfile is not None:
            if not isinstance(ncfile, FatBandsFile):
                ncfile = FatBandsFile(ncfile)
                closeit = True
        else:
            ncfile = self

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

        self.plot_fatbands_lview(e0=e0, fact=fact, alpha=alpha, ax_list=fatbands_axlist, show=False)

        if edos_kwargs is None: edos_kwargs = {}
        ncfile.plot_pjdos(e0=e0, ax_list=pjdos_axlist, exchange_xy=True, stacked=stacked, show=False, **edos_kwargs)
        if closeit: ncfile.close()

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file is created.
        Return path to the notebook.
        """
        import io, tempfile
        if nbpath is None:
            _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("fbfile = abilab.abiopen('%s')\nprint(fbfile)" % self.filepath),
            nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands()"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands_lview()"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands_m(iatom=0)"),
            nbv.new_code_cell("fig = fbfile.plot_pjdos()"),
            nbv.new_code_cell("fig = fbfile.plot_fatbands_with_pjdos(ncfile=None)"),
        ])

        with io.open(nbpath, 'wt', encoding="utf8") as fh:
            nbformat.write(nb, fh)
        return nbpath

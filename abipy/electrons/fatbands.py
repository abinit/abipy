# coding: utf-8
"""Classes for the analysis of electronic fatbands and projected DOSes."""
from __future__ import print_function, division, unicode_literals, absolute_import

import traceback
import numpy as np

from collections import OrderedDict, defaultdict
from tabulate import tabulate
from monty.termcolor import cprint
from monty.functools import lazy_property
from monty.string import marquee, list_strings
from pymatgen.core.periodic_table import Element
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons.ebands import ElectronsReader
from abipy.tools import gaussian
from abipy.tools.plotting import set_axlims, get_axarray_fig_plt, add_fig_kwargs, get_ax_fig_plt


def gaussians_dos(dos, mesh, width, values, energies, weights):
    assert len(dos) == len(mesh) and len(values) == len(energies) == len(weights)
    for vw, e, w in zip(values * weights, energies, weights):
        dos += vw * gaussian(mesh, width, center=e)
    return dos


class FatBandsFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    Provides methods to analyze the data stored in the FATBANDS.nc_ file.

    Usage example:

    .. code-block:: python

        with abiopen("out_FATBANDS.nc") as fb:
            print(fb.structure)
            fb.plot_fatbands_lview()

    Alternatively, one can use the `abiopen.py` script to open the file in an ipython notebook with::

        abiopen.py out_FATBANDS.nc -nb

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: FatBandsFile
    """
    # These class attributes can be redefined in order to customize the plots.

    # Mapping L --> color used in plots.
    l2color = {0: "red", 1: "blue", 2: "green", 3: "yellow", 4: "orange"}

    # Mapping L --> title used in subplots that depend on L.
    l2tex = {0: "$l=s$", 1: "$l=p$", 2: "$l=d$", 3: "$l=f$", 4: "$l=g$"}

    # Markers used for up/down bands (collinear spin)
    marker_spin = {0: "^", 1: "v"}

    # \U starts an eight-character Unicode escape. raw strings do not work in python2.7
    # and we need a latex symbol to avoid errors in matplotlib --> replace myuparrow --> uparrow

    # Mapping spin --> title used in subplots that depend on (collinear) spin.
    spin2tex = {k: v.replace("myuparrow", "uparrow") for k, v in
            {0: r"$\sigma=\myuparrow$", 1: r"$\sigma=\downarrow$"}.items()}

    # Mappings used for non-collinear spins.
    spinors2tex = {k: v.replace("myuparrow", "uparrow") for k, v in
        {"up-up": r"$\myuparrow,\myuparrow$", "up-down": r"$\myuparrow,\downarrow$",
                   "down-up": r"$\downarrow,\myuparrow$", "down-down": r"$\downarrow,\downarrow$",
                   "sigma_x": r"$\sigma_{x}$", "sigma_y": r"$\sigma_{y}$", "sigma_z": r"$\sigma_{z}$",
                  }.items()}

    spinors2color = {"up-up": "black", "up-down": "brown",
                     "down-up": "violet", "down-down": "yellow",
                     "sigma_x": "green", "sigma_y": "blue", "sigma_z": "red",
                   }

    # The alpha blending value, between 0 (transparent) and 1 (opaque)
    alpha = 0.6

    # Options passed to ebands.plot_ax. See also eb_plotax_kwargs
    # Either you change these values or subclass `FatBandsFile` and redefine eb_plotax_kwargs.
    marker_size = 3.0
    linewidth = 0.1
    linecolor = "grey"
    #klabel_size = None

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(FatBandsFile, self).__init__(filepath)
        self.reader = r = ElectronsReader(filepath)

        # Initialize the electron bands from file
        self._ebands = r.read_ebands()
        self.natom = len(self.structure)

        # Read metadata so that we know how to handle the content of the file.
        self.prtdos = r.read_value("prtdos")
        self.prtdosm = r.read_value("prtdosm")
        self.pawprtdos = r.read_value("pawprtdos")
        self.usepaw = r.read_value("usepaw")
        self.natsph = r.read_dimvalue("natsph")
        self.iatsph = r.read_value("iatsph") - 1 # F --> C
        self.ndosfraction = r.read_dimvalue("ndosfraction")
        # mbesslang: 1 + maximum angular momentum for Bessel function expansion.
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
        #   The wavefunction can have l-components > lmax_type, especially if vloc = vlmax.
        # If usepaw == 1, lmax_type represents the max l included in the PAW basis set.
        #   The AE wavefunction cannot have more ls than l-max if pawprtdos == 2 and
        #   the cancellation between PS-onsite and the smooth part is exact.
        # TODO: Decide how to change lmax_type at run-time: API or global self.set_lmax?
        self.lmax_type = r.read_value("lmax_type")
        if self.usepaw == 0 or (self.usepaw == 1 and self.pawprtdos != 2):
            self.lmax_type[:] = self.mbesslang - 1

        self.typat = r.read_value("atom_species") - 1 # F --> C
        self.lmax_atom = np.empty(self.natom, dtype=np.int)
        for iat in range(self.natom):
            self.lmax_atom[iat] = self.lmax_type[self.typat[iat]]
        # lsize is used to dimension arrays that depend on L.
        self.lsize = self.lmax_type.max() + 1

        # Sort the chemical symbols and use OrderedDict because we are gonna use these dicts for looping.
        # Note that we don't have arrays dimensioned with ntypat in the nc file so we can define
        # our own ordering for symbols.
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
            cmap = plt.get_cmap('jet')
            nsymb = len(self.symbols)
            for i, symb in enumerate(self.symbols):
                self.symbol2color[symb] = cmap(i/nsymb)

        # Array dimensioned with natom. Set to true if iatom has been calculated
        if self.prtdos == 3:
            self.has_atom = np.zeros(self.natom, dtype=bool)
            self.has_atom[self.iatsph] = True

    @lazy_property
    def wal_sbk(self):
        """
        |numpy-array| of shape [natom, mbesslang, nsppol, mband, nkpt]
        with the L-contributions. Present only if prtdos == 3.
        """
        return self._read_wal_sbk()

    @lazy_property
    def walm_sbk(self):
        """
        |numpy-array| of shape [natom, mbesslang**2, nsppol, mband, nkpt]
        with the LM-contribution. Present only if prtdos == 3 and prtdosm != 0
        """
        return self._read_walm_sbk(key="dos_fractions_m")

    def _read_wal_sbk(self, key="dos_fractions"):
        # Read dos_fraction_m from file and build wal_sbk array of shape
        # [natom, lmax, nsppol, mband, nkpt].
        #
        # In abinit the **Fortran** array has shape
        #   dos_fractions(nkpt,mband,nsppol,ndosfraction)
        #
        # Note that Abinit allows the users to select a subset of atoms with iatsph. Moreover the order
        # of the atoms could differ from the one in the structure even when natom == natsph (unlikely but possible).
        # To keep it simple, the code always operate on an array dimensioned with the total number of atoms
        # Entries that are not computed are set to zero and a warning is issued.
        if self.prtdos != 3:
            raise RuntimeError("The file does not contain L-DOS since prtdos=%i" % self.prtdos)

        wshape = (self.natom, self.mbesslang, self.nsppol, self.mband, self.nkpt)

        if self.natsph == self.natom and np.all(self.iatsph == np.arange(self.natom)):
            # All atoms have been calculated and the order if ok.
            wal_sbk = np.reshape(self.reader.read_value(key), wshape)

        else:
            # Need to tranfer data. Note np.zeros.
            wal_sbk = np.zeros(wshape)
            if self.natsph == self.natom and np.any(self.iatsph != np.arange(self.natom)):
                print("Will rearrange filedata since iatsp != [1, 2, ...])")
                filedata = np.reshape(self.reader.read_value(key), wshape)
                for i, iatom in enumerate(self.iatsph):
                    wal_sbk[iatom] = filedata[i]
            else:
                print("natsph < natom. Will set to zero the PJDOS contributions for the atoms that are not included.")
                assert self.natsph < self.natom
                filedata = np.reshape(self.reader.read_value(key),
                                     (self.natsph, self.mbesslang, self.nsppol, self.mband, self.nkpt))
                for i, iatom in enumerate(self.iatsph):
                    wal_sbk[iatom] = filedata[i]

        # In principle, this should never happen (unless there's a bug in Abinit or a
        # very bad cancellation between the FFT and the PS-PAW term (pawprtden=0).
        num_neg = np.sum(wal_sbk < 0)
        if num_neg:
            print("WARNING: There are %d (%.1f%%) negative entries in LDOS weights" % (
                  num_neg, 100 * num_neg / wal_sbk.size))

        return wal_sbk

    def _read_walm_sbk(self, key="dos_fraction_m"):
        # Read dos_fraction_m from file and build walm_sbk array of shape
        # [natom, lmax**2, nsppol, mband, nkpt].
        #
        # In abinit the **Fortran** array has shape
        #   dos_fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang*m_dos_flag)
        #
        # Note that Abinit allows the users to select a subset of atoms with iatsph. Moreover the order
        # of the atoms could differ from the one in the structure even when natom == natsph (unlikely but possible).
        # To keep it simple, the code always operate on an array dimensioned with the total number of atoms
        # Entries that are not computed are set to zero and a warning is issued.
        if self.prtdos != 3:
            raise RuntimeError("The file does not contain L-DOS since prtdos=%i" % self.prtdos)
        if self.prtdosm == 0:
            raise RuntimeError("The file does not contain LM-DOS since prtdosm=%i" % self.prtdosm)

        wshape = (self.natom, self.mbesslang**2, self.nsppol, self.mband, self.nkpt)

        if self.natsph == self.natom and np.all(self.iatsph == np.arange(self.natom)):
            # All atoms have been calculated and the order if ok.
            walm_sbk = np.reshape(self.reader.read_value(key), wshape)

        else:
            # Need to tranfer data. Note np.zeros.
            walm_sbk = np.zeros(wshape)
            if self.natsph == self.natom and np.any(self.iatsph != np.arange(self.natom)):
                print("Will rearrange filedata since iatsp != [1, 2, ...])")
                filedata = np.reshape(self.reader.read_value(key), wshape)
                for i, iatom in enumerate(self.iatsph):
                    walm_sbk[iatom] = filedata[i]
            else:
                print("natsph < natom. Will set to zero the PJDOS contributions for the atoms that are not included.")
                assert self.natsph < self.natom
                filedata = np.reshape(r.read_value(key),
                                     (self.natsph, self.mbesslang**2, self.nsppol, self.mband, self.nkpt))
                for i, iatom in enumerate(self.iatsph):
                    walm_sbk[iatom] = filedata[i]

        # In principle, this should never happen (unless there's a bug in Abinit or a
        # very bad cancellation between the FFT and the PS-PAW term (pawprtden=0).
        num_neg = np.sum(walm_sbk < 0)
        if num_neg:
            print("WARNING: There are %d (%.1f%%) negative entries in LDOS weights" % (
                  num_neg, 100 * num_neg / walm_sbk.size))

        return walm_sbk

    @property
    def ebands(self):
        """|ElectronBands| object."""
        return self._ebands

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    def close(self):
        """Called at the end of the ``with`` context manager."""
        return self.reader.close()

    def __str__(self):
        """String representation"""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=True, title="Electronic Bands"))
        app("")
        app(marquee("Fatbands Info", mark="="))
        app("prtdos=%d, prtdosm=%d, mbesslang=%d, pawprtdos=%d, usepaw=%d" % (
            self.prtdos, self.prtdosm, self.mbesslang, self.pawprtdos, self.usepaw))
        app("nsppol=%d, nkpt=%d, mband=%d" % (self.nsppol, self.nkpt, self.mband))
        app("")

        if self.prtdos == 3:
            # Print table with info on atoms.
            table = [["Idx", "Symbol", "Reduced_Coords", "Lmax", "Ratsph [Bohr]", "Has_Atom"]]
            for iatom, site in enumerate(self.structure):
                table.append([
                    iatom,
                    site.specie.symbol,
                    "%.5f %.5f %.5f" % tuple(site.frac_coords),
                    self.lmax_atom[iatom],
                    self.ratsph_type[self.typat[iatom]],
                    "Yes" if self.has_atom[iatom] else "No",
                ])
            app(tabulate(table, headers="firstrow"))

        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    def get_wl_atom(self, iatom, spin=None, band=None):
        """
        Return the l-dependent DOS weights for atom index ``iatom``. The weights are summed over m.
        If ``spin`` and ``band`` are not specified, the method returns the weights
        for all spin and bands else the contribution for (spin, band)
        """
        if spin is None and band is None:
            return self.wal_sbk[iatom]
        else:
            assert spin is not None and band is not None
            return self.wal_sbk[iatom, :, spin, band, :]

    def get_wl_symbol(self, symbol, spin=None, band=None):
        """
        Return the l-dependent DOS weights for a given type specified in terms of the
        chemical symbol ``symbol``. The weights are summed over m and over all atoms of the same type.
        If ``spin`` and ``band`` are not specified, the method returns the weights
        for all spins and bands else the contribution for (spin, band).
        """
        if spin is None and band is None:
            wl = np.zeros((self.lsize, self.nsppol, self.mband, self.nkpt))
            for iat in self.symbol2indices[symbol]:
                for l in range(self.lmax_atom[iat]+1):
                    wl[l] += self.wal_sbk[iat, l]
        else:
            assert spin is not None and band is not None
            wl = np.zeros((self.lsize, self.nkpt))
            for iat in self.symbol2indices[symbol]:
                for l in range(self.lmax_atom[iat]+1):
                    wl[l, :] += self.wal_sbk[iat, l, spin, band, :]

        return wl

    def get_w_symbol(self, symbol, spin=None, band=None):
        """
        Return the DOS weights for a given type specified in terms of the
        chemical symbol ``symbol``. The weights are summed over m and lmax[symbol] and
        over all atoms of the same type.
        If ``spin`` and ``band`` are not specified, the method returns the weights
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
        Return the spilling parameter
        If ``spin`` and ``band`` are not specified, the method returns the spilling for all states
        as a [nsppol, mband, nkpt] numpy array else the spilling for (spin, band) with shape [nkpt].
        """
        if spin is None and band is None:
            sp = np.zeros((self.nsppol, self.mband, self.nkpt))
            for iatom in range(self.natom):
                for l in range(self.lmax_atom[iatom]+1):
                    sp += self.wal_sbk[iatom, l]
        else:
            assert spin is not None and band is not None
            sp = np.zeros((self.nkpt))
            for iatom in range(self.natom):
                for l in range(self.lmax_atom[iatom]+1):
                    sp += self.wal_sbk[iatom, l, spin, band, :]

        return 1.0 - sp

    def eb_plotax_kwargs(self, spin):
        """
        Dictionary with the options passed to ``ebands.plot_ax``
        when plotting a band line with spin index ``spin``.
        Subclasses can redefine the implementation to customize the plots.
        """
        return dict(
            color=self.linecolor,
            linewidth=self.linewidth,
            markersize=self.marker_size,
            marker=self.marker_spin[spin],
            #klabel_size=self.klabel_size,
        )

    @add_fig_kwargs
    def plot_fatbands_siteview(self, e0="fermie", view="inequivalent", fact=1.0, fontsize=12,
                               ylims=None, blist=None, **kwargs):
        """
        Plot fatbands for each atom in the unit cell. By default, only the **inequivalent** atoms are shown.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            view: "inequivalent", "all"
            fact:  float used to scale the stripe size.
            fontsize: fontsize for titles and legend
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            blist: List of band indices for the fatband plot. If None, all bands are included

        Returns: |matplotlib-Figure|
        """
        # Define num_plots and ax2atom depending on view.
        # ax2natom[1:num_plots] --> iatom index in structure.
        # TODO: ebands.used_magnetic_symmetries?
        if view == "inequivalent" and (self.nspden == 2 and self.nsppol == 1) or (self.nspinor == 2 and self.nspden != 4):
            cprint("The system with magnetic symmetries but the spglib API used by pymatgen does not support them.", "yellow")
            cprint("    nsppol: %s, nspden: %s, nspinor: %s" % (self.nsppol, self.nspden, self.nspinor), "yellow")
            cprint("Setting view to `all`", "yellow")
            view = "all"

        # TODO: spin
        if view == "all" or self.natom == 1:
            num_plots, ax2iatom = self.natom, np.arange(self.natom)

        elif view == "inequivalent":
            print("Calling spglib to find inequivalent sites.")
            print("Note that `symafm` magnetic symmetries (if any) are not taken into account.")
            ea = self.structure.spget_equivalent_atoms(printout=True)
            num_plots = len(ea.irred_pos)
            ax2iatom = ea.irred_pos
        else:
            raise ValueError("Wrong value for view: %s" % str(view))

        # Build plot grid.
        ncols, nrows = 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_mat[-1, -1].axis("off")

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)
        mybands = range(ebands.mband) if blist is None else blist

        for iax, iatom in enumerate(ax2iatom):
            ax = ax_mat.flat[iax]
            # Plot the energies.
            for spin in range(self.nsppol):
                ebands.plot_ax(ax, e0, spin=spin, **self.eb_plotax_kwargs(spin))

            site = self.structure[iatom]
            ebands.decorate_ax(ax, title=str(site))

            # Add width around each band.
            for spin in ebands.spins:
                for band in mybands:
                    wlk = self.get_wl_atom(iatom, spin=spin, band=band) * (fact / 2)
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup
                    for l in range(self.lmax_atom[iatom]+1):
                        w = wlk[l,:]
                        y1, y2 = yup + w, ydown - w
                        ax.fill_between(x, yup, y1, alpha=self.alpha, facecolor=self.l2color[l])
                        ax.fill_between(x, ydown, y2, alpha=self.alpha, facecolor=self.l2color[l],
                                        label=self.l2tex[l] if (spin, band) == (0, 0) else None)
                                        # Note: could miss a label in the other plots if lmax is not large enough!
                        yup, ydown = y1, y2

            set_axlims(ax, ylims, "y")

        ax_mat[0, 0].legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_fatbands_lview(self, e0="fermie", fact=1.0, ax_mat=None, lmax=None,
                            ylims=None, blist=None, fontsize=12, **kwargs):
        """
        Plot the electronic fatbands grouped by L.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            fact:  float used to scale the stripe size.
            ax_mat: Matrix of axes, if None a new figure is produced.
            lmax: Maximum L included in plot. None means full set available on file.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            blist: List of band indices for the fatband plot. If None, all bands are included
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        mylsize = self.lsize if lmax is None else lmax + 1
        # Build or get grid with (nsppol, mylsize) axis.
        nrows, ncols = self.nsppol, mylsize
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)
        ax_mat = np.reshape(ax_mat, (nrows, ncols))

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)
        mybands = range(ebands.mband) if blist is None else blist

        for spin in range(self.nsppol):
            for l in range(mylsize):
                ax = ax_mat[spin, l]
                ebands.plot_ax(ax, e0, spin=spin, **self.eb_plotax_kwargs(spin))
                title = "%s, %s" % (self.l2tex[l], self.spin2tex[spin]) if self.nsppol == 2 else "%s" % self.l2tex[l]
                ebands.decorate_ax(ax, title=title)

                if l != 0:
                    ax.set_ylabel("")
                    # Only the first column show labels.
                    # Trick: Don't change the labels but set their fontsize to 0 otherwise
                    # also the other axes are affected (likely due to sharey=True).
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(0)

                for ib, band in enumerate(mybands):
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup
                    for symbol in self.symbols:
                        wlk = self.get_wl_symbol(symbol, spin=spin, band=band) * (fact / 2)
                        w = wlk[l]
                        y1, y2 = yup + w, ydown - w
                        # Add width around each band. Only the [0,0] plot has the legend.
                        ax.fill_between(x, yup, y1, alpha=self.alpha, facecolor=self.symbol2color[symbol])
                        ax.fill_between(x, ydown, y2, alpha=self.alpha, facecolor=self.symbol2color[symbol],
                                        label=symbol if (l, spin, ib) == (0, 0, 0) else None)
                        yup, ydown = y1, y2

                set_axlims(ax, ylims, "y")

        ax_mat[0, 0].legend(loc="best", fontsize=fontsize, shadow=True)
        return fig

    @add_fig_kwargs
    def plot_fatbands_mview(self, iatom, e0="fermie", fact=1.0, lmax=None,
                            ylims=None, blist=None, **kwargs):
        """
        Plot the electronic fatbands grouped by LM.

        Args:
            iatom: Index of the atom in the structure.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            fact:  float used to scale the stripe size.
            lmax: Maximum L included in plot. None means full set available on file.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            blist: List of band indices for the fatband plot. If None, all bands are included

        Returns: |matplotlib-Figure|
        """
        # TODO
        #if self.nsppol == 2:
        #    raise NotImplementedError("To be tested with spin")

        if not (self.prtdos == 3 and self.prtdosm != 0):
            cprint("Fatbands plots with LM-character require `prtdos = 3 and prtdosm != 0`", "red")
            return None

        mylmax = self.lmax_atom[iatom] if lmax is None else lmax

        # Build plot grid.
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure()
        nrows, ncols = 2 * (mylmax+1), mylmax + 1
        gspec = GridSpec(nrows=nrows, ncols=ncols, wspace=0.1, hspace=0.1)

        # Build plot grid (L along the column, each L has 2L+1 subplots).
        # ax_lim[(l, im)] gives the axis.
        ax_lim = {}
        for im in range(nrows):
            for l in range(ncols):
                k = (l, im)
                ax00 = None if l == 0 else ax_lim[(0, 0)]
                ax = plt.subplot(gspec[im, l], sharex=ax00, sharey=ax00)
                if im < 2*l + 1:
                    #ax.set_title("l=%d, m=%d" % (l, im - l))
                    ax_lim[k] = ax
                    ax.grid(True)
                else:
                    ax.axis("off")

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)
        mybands = range(ebands.mband) if blist is None else blist

        for lim, ax in ax_lim.items():
            l, im = lim[0], lim[1]
            for spin in range(self.nsppol):
                ebands.plot_ax(ax, e0, spin=spin, **self.eb_plotax_kwargs(spin))

            if im == 2 * l:
               ebands.decorate_ax(ax)
            #if l > 0:
            #    ax.set_ylabel("")

            for spin in range(self.nsppol):
                for band in mybands:
                    #print("band:", band)
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup

                    w = self.walm_sbk[iatom,  l**2 + im, spin, band, :] * (fact / 2)
                    y1, y2 = yup + w, ydown - w
                    # Add width around each band.
                    ax.fill_between(x, yup, y1, alpha=self.alpha, facecolor=self.l2color[l])
                    ax.fill_between(x, ydown, y2, alpha=self.alpha, facecolor=self.l2color[l])
                    yup, ydown = y1, y2

            set_axlims(ax, ylims, "y")

        return fig

    @add_fig_kwargs
    def plot_fatbands_typeview(self, e0="fermie", fact=1.0, lmax=None, ax_mat=None, ylims=None,
                              blist=None, fontsize=8, **kwargs):
        """
        Plot the electronic fatbands grouped by atomic type.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            fact:  float used to scale the stripe size.
            lmax: Maximum L included in plot. None means full set available on file.
            ax_mat: Matrix of axis. None if a new figure should be created.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            blist: List of band indices for the fatband plot. If None, all bands are included
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        mylsize = self.lsize if lmax is None else lmax + 1

        # Get ax_mat and fig.
        nrows, ncols = self.nsppol, self.ntypat
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)
        ax_mat = np.reshape(ax_mat, (nrows, ncols))

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)
        mybands = range(ebands.mband) if blist is None else blist

        for itype, symbol in enumerate(self.symbols):
            wl_sbk = self.get_wl_symbol(symbol) * (fact / 2)
            for spin in range(self.nsppol):
                ax = ax_mat[spin, itype]
                ebands.plot_ax(ax, e0, spin=spin, **self.eb_plotax_kwargs(spin))

                title = ("type=%s, %s" % (symbol, self.spin2tex[spin]) if self.nsppol == 2
                         else "type=%s" % symbol)
                ebands.decorate_ax(ax, title=title)
                if itype != 0:
                    ax.set_ylabel("")

                # Plot fatbands for given (symbol, spin) and all angular momenta.
                for band in mybands:
                    yup = ebands.eigens[spin, :, band] - e0
                    ydown = yup
                    for l in range(min(self.lmax_symbol[symbol] + 1, mylsize)):
                        # Add width around each band.
                        w = wl_sbk[l, spin, band]
                        y1, y2 = yup + w, ydown - w
                        ax.fill_between(x, yup, y1, alpha=self.alpha, facecolor=self.l2color[l])
                        ax.fill_between(x, ydown, y2, alpha=self.alpha, facecolor=self.l2color[l],
                                        label=self.l2tex[l] if (itype, spin, band) == (0, 0, 0) else None)
                                        # Note: could miss a label in the other plots if lmax is not large enough!
                        yup, ydown = y1, y2

                set_axlims(ax, ylims, "y")

        ax_mat[0, 0].legend(loc="best", fontsize=fontsize, shadow=True)
        return fig

    @add_fig_kwargs
    def plot_spilling(self, e0="fermie", fact=1.0, ax_list=None, ylims=None, blist=None, **kwargs):
        """
        Plot the electronic fatbands

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            fact:  float used to scale the stripe size.
            ax_list: List of matplotlib axes for plot. If None, new figure is produced
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            blist: List of band indices for the fatband plot. If None, all bands are included

        Returns: |matplotlib-Figure|
        """
        nrows, ncols = 1, self.nsppol
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        ebands = self.ebands
        e0 = ebands.get_e0(e0)
        x = np.arange(self.nkpt)
        mybands = range(ebands.mband) if blist is None else blist
        spill_sbk = self.get_spilling() * (fact / 2)

        for spin in range(self.nsppol):
            ax = ax_list[spin]
            ebands.plot_ax(ax, e0, spin=spin, **self.eb_plotax_kwargs(spin))
            ebands.decorate_ax(ax)

            for band in mybands:
                y = ebands.eigens[spin, :, band] - e0
                w = spill_sbk[spin, band, :]

                # Handle negative spilling values.
                wispos = w >= 0.0
                wisneg = np.logical_not(wispos)
                num_neg = np.sum(wisneg)

                # Add width around each band.
                ax.fill_between(x, y, y + w, where=wispos, alpha=self.alpha, facecolor="blue")
                ax.fill_between(x, y, y - w, where=wispos, alpha=self.alpha, facecolor="blue")

                # Show regions with negative spilling in red.
                if num_neg:
                    print("For spin:", spin, "band:", band,
                          "There are %d (%.1f%%) k-points with negative spilling. Min: %.2E" % (
                           num_neg, 100 * num_neg / self.nkpt, w.min()))

                    absw = np.abs(w)
                    ax.fill_between(x, y, y + absw, where=wisneg, alpha=self.alpha, facecolor="red")
                    ax.fill_between(x, y, y - absw, where=wisneg, alpha=self.alpha, facecolor="red")

            set_axlims(ax, ylims, "y")

        return fig

    # TODO: THIS CODE IS STILL UNDER DEVELOPMENT
    #@add_fig_kwargs
    #def plot_fatbands_spinor(self, terms=("sigma_z",), e0="fermie", fact=1,
    #                         ylims=None, blist=None, ax_list=None, **kwargs):
    #    """
    #    Plot spinor-resolved electronic fatbands. Require prtdos = 5 and nspinor = 2.

    #    Args:
    #        terms: List of strings defining the quantities to plot. Possible values in:

    #            {"up-up", "up-down", "down-up", "down-down", "sigma_x", "sigma_y", "sigma_z"}

    #        e0: Option used to define the zero of energy in the band structure plot. Possible values:
    #            - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy.
    #            -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
    #            -  None: Don't shift energies, equivalent to e0=0
    #        fact:  float used to scale the marker size.
    #        ax_list: List of matplotlib axes for plot. If None, new figure is produced
    #        ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
    #               or scalar e.g. `left`. If left (right) is None, default values are used
    #        blist: List of band indices for the fatband plot. If None, all bands are included

    #    Returns: |matplotlib-Figure|
    #    """
    #    if not (self.prtdos == 5 and self.nspinor == 2):
    #        print("This method assumes prtdos=5 and nspinor=2 but file has prtdos=%s and nspinor=%s"
    #            % (self.prtdos, self.nspinor))
    #        return None

    #    # Build ax_list
    #    import matplotlib.pyplot as plt
    #    if ax_list is None:
    #        fig, ax_list = plt.subplots(nrows=len(terms), ncols=1, sharex=True, sharey=True, squeeze=False)
    #        ax_list = ax_list.ravel()
    #    else:
    #        ax_list = np.reshape(ax_list, (1, len(terms))).ravel()
    #        fig = plt.gcf()

    #    ebands = self.ebands
    #    e0 = ebands.get_e0(e0)
    #    xvals = np.arange(self.nkpt)
    #    mybands = range(ebands.mband) if blist is None else blist
    #    spin0 = 0

    #    # Dos fractions contains: [ndosfraction, nsppol, mband, nkpt].
    #    # where the first dimension stores: (up,up  up,down  down,up  down,down  sigma_x sigma_y sigma_z)
    #    w_sbk = self.reader.read_value("dos_fractions")
    #    term2idx = {"up-up":0, "up-down":1, "down-up":2, "down-down":3, "sigma_x":4, "sigma_y":5, "sigma_z":6}

    #    # TODO: To be tested.
    #    terms = list_strings(terms)
    #    for i, (term, ax) in enumerate(zip(terms, ax_list)):
    #        ebands.plot_ax(ax, e0, spin=spin0, color=self.linecolor, linewidth=self.linewidth)
    #        ebands.decorate_ax(ax, title=term)
    #        if i != 0:
    #            ax.set_ylabel("")
    #            # Only the first column show labels.
    #            # Trick: Don't change the labels but set their fontsize to 0 otherwise
    #            # also the other axes are affecred (likely due to sharey=True).
    #            for tick in ax.yaxis.get_major_ticks():
    #                tick.label.set_fontsize(0)

    #        idx = term2idx[term]
    #        color = self.spinors2color[term]
    #        # Rescale weights in [0, 1]
    #        data = w_sbk[idx, spin0, :, :].copy()
    #        data -= data.min()
    #        dmax = data.max()
    #        if dmax != 0.0: data /= dmax
    #        #print(data)

    #        for ib, band in enumerate(mybands):
    #            yvals = ebands.eigens[spin0, :, band] - e0
    #            ax.scatter(xvals, yvals, c=data[band], s=20, alpha=self.alpha,
    #                       #edgecolor=plt.get_cmap('jet')(data[band]),
    #                       marker="o", label=self.spinors2tex[term] if ib == 0 else None)

    #            #ws = w_sbk[idx, spin0, band, :] * fact
    #            #pos_inds, neg_inds = np.where(ws >= 0), np.where(ws < 0)
    #            #ax.scatter(xvals[pos_inds], yvals[pos_inds], c=color, s=ws[pos_inds],
    #            #           marker="^", label=term if ib == 0 else None)
    #            #ax.scatter(xvals[neg_inds], yvals[neg_inds], c=color, s=np.abs(ws[neg_inds]),
    #            #           marker="v", label=term if ib == 0 else None)

    #        set_axlims(ax, ylims, "y")

    #    ax_list[0].legend(loc="best", fontsize=fontsize, shadow=True)
    #    return fig

    # TODO: THIS CODE IS STILL UNDER DEVELOPMENT
    #def nelect_in_spheres(self, start_energy=None, stop_energy=None,
    #                     method="gaussian", step=0.1, width=0.2):
    #    """
    #    Print the number of electrons inside each atom-centered sphere.
    #    Note that this is a very crude estimate of the charge density distribution.

    #    Args:
    #        start_energy: PJDOS is integrated from this energy (eV). If None, the lower bound in used.
    #        stop_energy: PJDOS is integrated up to this energy (eV). If None, the Fermi level is used.
    #        method: String defining the method for the computation of the DOS.
    #        step: Energy step (eV) of the linear mesh.
    #        width: Standard deviation (eV) of the gaussian.
    #    """
    #    intg = self.get_dos_integrator(method, step, width)
    #    raise NotImplementedError("")
    #    site_edos = intg.site_edos
    #    if stop_energy is None: stop_energy = self.ebands.fermie

    #    # find_mesh_indes returns the first point in the mesh whose value is >= value. -1 if not found
    #    edos = site_edos[0]
    #    start_spin, stop_spin = {}, {}
    #    for spin in self.spins:
    #        start = 0
    #        if start_energy is not None:
    #            start = edos[spin].find_mesh_index(start_energy)
    #        if start == -1:
    #            raise ValueError("For spin %d: cannot find index in mesh such that mesh[i] >= start." % spin)
    #        if start > 0: start -= 1
    #        start_spin[spin] = start

    #        stop = edos[spin].find_mesh_index(stop_energy)
    #        if stop == -1:
    #            raise ValueError("For spin %d: cannot find index in mesh such that mesh[i] >= energy." % spin)
    #        stop_spin[spin] = stop

    #    for iatm, site in enumerate(self.structure):
    #        edos = site_edos[iatm]
    #        nel_spin = {}
    #        for spin in self.spins:
    #            nel_spin[spin] = edos[spin].integral(start=start_spin[spin], stop=stop_spin[spin])
    #        print("iatom", iatm, "site", site, nel_spin)

    def get_dos_integrator(self, method, step, width):
        """
        FatBandsFile can use differerent integrators that are cached in self._cached_dos_integrators
        """
        if not hasattr(self, "_cached_dos_integrators"): self._cached_dos_integrators = {}
        key = (method, step, width)
        intg = self._cached_dos_integrators.get(key, None)
        if intg is not None: return intg
        # Build integrator, cache it and return it.
        intg = _DosIntegrator(self, method, step, width)
        self._cached_dos_integrators[key] = intg
        return intg

    @add_fig_kwargs
    def plot_pjdos_lview(self, e0="fermie", lmax=None, method="gaussian", step=0.1, width=0.2,
                         stacked=True, combined_spins=True, ax_mat=None, exchange_xy=False,
                         with_info=True, with_spin_sign=True, xlims=None, ylims=None, fontsize=8, **kwargs):
        """
        Plot the PJ-DOS on a linear mesh.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            lmax: Maximum L included in plot. None means full set available on file.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            combined_spins: Define how up/down DOS components should be plotted when nsppol==2.
                If True, up/down DOSes are plotted on the same figure (positive values for up,
                negative values for down component)
                If False, up/down components are plotted on different axes.
            ax_mat:
            exchange_xy: True if the dos should be plotted on the x axis instead of y.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ylims: Same meaning as ``xlims`` but for the y-axis
            fontsize: Legend and label fontsize

        Returns: |matplotlib-Figure|
        """
        try:
            intg = self.get_dos_integrator(method, step, width)
        except Exception:
            msg = traceback.format_exc()
            msg += ("Error while trying to compute the DOS.\n"
                    "Verify that the k-points form a homogenous sampling of the BZ.\n"
                    "Returning None\n")
            cprint(msg, "red")
            return None

        # Get energy mesh from total DOS and define the zero of energy
        # Note that the mesh is not not spin-dependent.
        e0 = self.ebands.get_e0(e0)
        mesh = intg.mesh.copy()
        mesh -= e0
        edos, symbols_lso = intg.edos, intg.symbols_lso

        # Get grid of axes.
        mylsize = self.lsize if lmax is None else lmax + 1
        nrows = self.nsppol if not combined_spins else 1
        ncols = mylsize
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)
        ax_mat = np.reshape(ax_mat, (nrows, ncols))

        # The code below expectes a matrix of axes of shape[nsppol, self.lsize]
        # If spins are plotted on the same graph (combined_spins), I build a new matrix so that
        # ax_mat[spin=0] is ax_mat[spin=1] and aliased_axis is set to True
        aliased_axis = False
        if self.nsppol == 2 and combined_spins:
            aliased_axis = True
            ax_mat = np.array([ax_mat.ravel(), ax_mat.ravel()])

        spin_sign = +1
        if not stacked:
            # Plot PJDOS as lines.
            for isymb, symbol in enumerate(self.symbols):
                for spin in range(self.nsppol):
                    if with_spin_sign: spin_sign = +1 if spin == 0 else -1
                    # Loop over the columns of the grid.
                    for l in range(min(self.lmax_symbol[symbol] + 1, mylsize)):
                        ax = ax_mat[spin, l]

                        # Plot total DOS.
                        x, y = mesh, spin_sign * edos.spin_dos[spin].values
                        if exchange_xy: x, y = y, x
                        label = "Tot" if (l, spin, isymb) == (0, 0, 0) else None
                        ax.plot(x, y, color="k", label=label if with_info else None)

                        # Plot PJ-DOS(l, spin)
                        x, y = mesh, spin_sign * symbols_lso[symbol][l, spin]
                        if exchange_xy: x, y = y, x
                        label = symbol if (l, spin, isymb) == (0, 0, 0) else None
                        ax.plot(x, y, color=self.symbol2color[symbol], label=label if with_info else None)

                    set_axlims(ax, xlims, "x")
                    set_axlims(ax, ylims, "y")

        else:
            # Plot stacked PJDOS
            # Loop over the columns of the grid.
            ls_stackdos = intg.ls_stackdos
            spin_sign = +1
            zerodos = np.zeros(len(mesh))
            for l in range(mylsize):
                for spin in self.ebands.spins:
                    if with_spin_sign: spin_sign = +1 if spin == 0 else -1
                    ax = ax_mat[spin, l]

                    # Plot total DOS.
                    x, y = mesh, spin_sign * edos.spin_dos[spin].values
                    if exchange_xy: x, y = y, x
                    label = "Tot" if (l, spin) == (0, 0) else None
                    ax.plot(x, y, color="k", label=label if with_info else None)

                    # Plot cumulative PJ-DOS(l, spin)
                    stack = ls_stackdos[(l, spin)] * spin_sign
                    for isymb, symbol in enumerate(self.symbols):
                        yup = stack[isymb]
                        ydown = stack[isymb-1] if isymb != 0 else zerodos
                        label ="%s (stacked)" % symbol if (l, spin) == (0, 0) else None
                        fill = ax.fill_between if not exchange_xy else ax.fill_betweenx
                        fill(mesh, yup, ydown, alpha=self.alpha, facecolor=self.symbol2color[symbol],
                             label=label if with_info else None)

                    set_axlims(ax, xlims, "x")
                    set_axlims(ax, ylims, "y")

        # Decorate axis.
        for spin in range(self.nsppol):
            if aliased_axis and spin == 1: break # Don't repeat yourself!

            for l in range(mylsize):
                ax = ax_mat[spin, l]
                if with_info:
                    if combined_spins:
                        title = self.l2tex[l]
                    else:
                        title = "%s, %s" % (self.l2tex[l], self.spin2tex[spin]) if self.nsppol == 2 else \
                                 self.l2tex[l]
                    ax.set_title(title, fontsize=fontsize)
                ax.grid(True)

                # Display yticklabels on the first plot and last plot only.
                # and display the legend only on the first plot.
                ax.set_xlabel("Energy (eV)")
                if l == 0:
                    if with_info:
                        ax.legend(loc="best", fontsize=fontsize, shadow=True)
                        if exchange_xy:
                            ax.set_xlabel('DOS (states/eV)')
                        else:
                            ax.set_ylabel('DOS (states/eV)')
                elif l == mylsize - 1:
                    ax.yaxis.set_ticks_position("right")
                    ax.yaxis.set_label_position("right")
                else:
                    # Plots in the middle: don't show labels.
                    # Trick: Don't change the labels but set their fontsize to 0 otherwise
                    # also the other axes are affecred (likely due to sharey=True).
                    # ax.set_yticklabels([])
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(0)

        return fig

    @add_fig_kwargs
    def plot_pjdos_typeview(self, e0="fermie", lmax=None, method="gaussian", step=0.1, width=0.2,
                            stacked=True, combined_spins=True, ax_mat=None, exchange_xy=False,
                            with_info=True, with_spin_sign=True, xlims=None, ylims=None, fontsize=8, **kwargs):
        """
        Plot the PJ-DOS on a linear mesh.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            lmax: Maximum L included in plot. None means full set available on file.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            combined_spins: Define how up/down DOS components should be plotted when nsppol==2.
                If True, up/down DOSes are plotted on the same figure (positive values for up,
                negative values for down component)
                If False, up/down components are plotted on different axes.
            ax_mat:
            exchange_xy: True if the dos should be plotted on the x axis instead of y.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ylims: Same meaning as ``xlims`` but for the y-axis
            fontsize: Legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        mylsize = self.lsize if lmax is None else lmax + 1

        try:
            intg = self.get_dos_integrator(method, step, width)
        except Exception:
            msg = traceback.format_exc()
            msg += ("Error while trying to compute the DOS.\n"
                    "Verify that the k-points form a homogenous sampling of the BZ.\n"
                    "Returning None\n")
            cprint(msg, "red")
            return None

        # Get energy mesh from total DOS and define the zero of energy
        # Note that the mesh is not not spin-dependent.
        e0 = self.ebands.get_e0(e0)
        mesh = intg.mesh.copy()
        mesh -= e0
        edos, symbols_lso = intg.edos, intg.symbols_lso

        # Get grid of axes.
        nrows = self.nsppol if not combined_spins else 1
        ncols = self.ntypat
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)
        ax_mat = np.reshape(ax_mat, (nrows, ncols))

        # The code below expectes a matrix of axes of shape[nsppol, self.ntypat]
        # If spins are plotted on the same graph (combined_spins), I build a new matrix so that
        # ax_mat[spin=0] is ax_mat[spin=1] and aliased_axis is set to True
        aliased_axis = False
        if self.nsppol == 2 and combined_spins:
            aliased_axis = True
            ax_mat = np.array([ax_mat.ravel(), ax_mat.ravel()])

        spin_sign = +1
        if not stacked:
            for spin in range(self.nsppol):
                if with_spin_sign: spin_sign = +1 if spin == 0 else -1
                # Loop over the columns of the grid.
                for isymb, symbol in enumerate(self.symbols):
                    ax = ax_mat[spin, isymb]

                    # Plot total DOS.
                    x, y = mesh, spin_sign * edos.spin_dos[spin].values
                    if exchange_xy: x, y = y, x
                    label = "Tot" if (spin, isymb) == (0, 0) else None
                    ax.plot(x, y, color="k", label=label if with_info else None)

                    for l in range(min(self.lmax_symbol[symbol] + 1, mylsize)):
                        # Plot PJ-DOS(l, spin)
                        x, y = mesh, spin_sign * symbols_lso[symbol][l, spin]
                        if exchange_xy: x, y = y, x
                        label = self.l2tex[l] if (spin, isymb) == (0, 0) else None
                        ax.plot(x, y, color=self.l2color[l], label=label if with_info else None)

                    set_axlims(ax, xlims, "x")
                    set_axlims(ax, ylims, "y")

        else:
            # Plot stacked PJDOS
            # Loop over the columns of the grid.
            #ls_stackdos = intg.ls_stackdos
            spin_sign = +1
            zerodos = np.zeros(len(mesh))
            for spin in range(self.nsppol):
                if with_spin_sign: spin_sign = +1 if spin == 0 else -1
                for isymb, symbol in enumerate(self.symbols):
                    ax = ax_mat[spin, isymb]

                    # Plot total DOS.
                    x, y = mesh, spin_sign * edos.spin_dos[spin].values
                    if exchange_xy: x, y = y, x
                    label = "Tot" if (spin, isymb) == (0, 0) else None
                    ax.plot(x, y, color="k", label=label if with_info else None)

                    # Plot cumulative PJ-DOS(l, spin)
                    stack = intg.get_lstack_symbol(symbol, spin) * spin_sign
                    for l in range(min(self.lmax_symbol[symbol] + 1, mylsize)):
                        yup = stack[l]
                        ydown = stack[l-1] if l != 0 else zerodos
                        label ="%s (stacked)" % self.l2tex[l] if (isymb, spin) == (0, 0) else None
                        fill = ax.fill_between if not exchange_xy else ax.fill_betweenx
                        fill(mesh, yup, ydown, alpha=self.alpha, facecolor=self.l2color[l],
                             label=label if with_info else None)

                    set_axlims(ax, xlims, "x")
                    set_axlims(ax, ylims, "y")

        # Decorate axis.
        for spin in range(self.nsppol):
            if aliased_axis and spin == 1: break # Don't repeat yourself!

            for itype, symbol in enumerate(self.symbols):
                ax = ax_mat[spin, itype]
                if with_info:
                    if combined_spins:
                        title = "Type: %s" % symbol
                    else:
                        title = "%s, %s" % (symbol, self.spin2tex[spin]) if self.nsppol == 2 else symbol
                    ax.set_title(title, fontsize=fontsize)
                ax.grid(True)

                # Display yticklabels on the first plot and last plot only.
                # and display the legend only on the first plot.
                ax.set_xlabel("Energy (eV)")
                if itype == 0:
                    if with_info:
                        ax.legend(loc="best", fontsize=fontsize, shadow=True)
                        if exchange_xy:
                            ax.set_xlabel('DOS (states/eV)')
                        else:
                            ax.set_ylabel('DOS (states/eV)')
                elif itype == self.ntypat - 1:
                    ax.yaxis.set_ticks_position("right")
                    ax.yaxis.set_label_position("right")
                else:
                    # Plots in the middle: don't show labels.
                    # Trick: Don't change the labels but set their fontsize to 0 otherwise
                    # also the other axes are affected (likely due to sharey=True).
                    # ax.set_yticklabels([])
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(0)

        return fig

    @add_fig_kwargs
    def plot_fatbands_with_pjdos(self, e0="fermie", fact=1.0, lmax=None, blist=None, view="type",
                                 pjdosfile=None, edos_kwargs=None, stacked=True, width_ratios=(2, 1),
                                 fontsize=8, ylims=None, **kwargs):
        """
        Compute the fatbands and the PJDOS on the same figure, a.k.a the Sistine Chapel.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy.
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            fact: float used to scale the stripe size.
            lmax: Maximum L included in plot. None means full set available on file.
            blist: List of band indices for the fatband plot. If None, all bands are included
            pjdosfile: FATBANDS file used to compute the PJDOS. If None, the PJDOS is taken from self.
            edos_kwargs:
            stacked: True if DOS partial contributions should be stacked on top of each other.
            width_ratios: Defines the ratio between the band structure plot and the dos plot.
            fontsize: Legend and label fontsize.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used

        Returns: |matplotlib-Figure|
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

        if not pjdosfile.ebands.kpoints.is_ibz:
            cprint("DOS requires k-points in the IBZ but got pjdosfile: %s" % repr(pjdosfile), "yellow")
            cprint("Returning None", "yellow")
            return None

        if edos_kwargs is None: edos_kwargs = {}

        # Build plot grid.
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure()
        # Define number of columns depending on view
        mylsize = self.lsize if lmax is None else lmax + 1
        #ncols = dict(type=self.ntypat, lview=self.lsize)[view]
        ncols = dict(type=self.ntypat, lview=mylsize)[view]
        gspec = GridSpec(nrows=self.nsppol, ncols=ncols)
        fatbands_axmat = np.empty((self.nsppol, ncols), dtype=object)
        pjdos_axmat = fatbands_axmat.copy()

        for spin in range(self.nsppol):
            for icol in range(ncols):
                subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gspec[spin, icol],
                                                  width_ratios=width_ratios, wspace=0.05)
                # Similar plots share the x-axis. All plots share the y-axis.
                prev_fatbax = None if (icol == 0 and spin == 0) else fatbands_axmat[0, 0]
                ax1 = plt.subplot(subgrid[0], sharex=prev_fatbax, sharey=prev_fatbax)
                prev_pjdosax = None if (icol == 0 and spin == 0) else pjdos_axmat[0, 0]
                ax2 = plt.subplot(subgrid[1], sharex=prev_pjdosax, sharey=ax1)
                fatbands_axmat[spin, icol] = ax1
                pjdos_axmat[spin, icol] = ax2

        # Plot bands on fatbands_axmat and PJDOS on pjdos_axmat.
        if view == "lview":
            self.plot_fatbands_lview(e0=e0, fact=fact, lmax=lmax, blist=blist, ax_mat=fatbands_axmat, ylims=ylims, show=False)
            pjdosfile.plot_pjdos_lview(e0=e0, lmax=lmax, ax_mat=pjdos_axmat, exchange_xy=True,
                                       stacked=stacked, combined_spins=False, fontsize=fontsize,
                                       with_info=False, with_spin_sign=False, show=False, ylims=ylims,
                                       **edos_kwargs)

        elif view == "type":
            self.plot_fatbands_typeview(e0=e0, fact=fact, lmax=lmax, blist=blist, ax_mat=fatbands_axmat, ylims=ylims, show=False)
            pjdosfile.plot_pjdos_typeview(e0=e0, lmax=lmax, ax_mat=pjdos_axmat, exchange_xy=True,
                                          stacked=stacked, combined_spins=False, fontsize=fontsize,
                                          with_info=False, with_spin_sign=False, show=False, ylims=ylims,
                                          **edos_kwargs)
        else:
            raise ValueError("Don't know how to handle view=%s" % str(view))

        # Remove labels from DOS plots.
        for ax in pjdos_axmat.ravel():
            ax.set_xlabel("")
            ax.set_ylabel("")
            for xtick, ytick in zip(ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()):
                xtick.label.set_fontsize(0)
                ytick.label.set_fontsize(0)

        if closeit: pjdosfile.close()
        return fig

    @add_fig_kwargs
    def plot_pawdos_terms(self, lmax=None, method="gaussian", step=0.1, width=0.2, xlims=None, *kwargs):
        """
        Plot ...

        Args:
            lmax: Maximum L included in plot. None means full set available on file.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used

        Returns: |matplotlib-Figure|
        """
        if self.usepaw != 1:
            print("This is not a PAW calculation!")
            return None
        if self.prtdos != 3:
            print("prtdos == 3 is required!")
            return None
        # TODO: More tests.
        #if self.pawprtdos != 0:

        mylsize = self.lsize if lmax is None else lmax + 1

        # Onsite contributions.
        # fracts_paw1,(dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction))
        #wshape = (self.natom, self.mbesslang, self.nsppol, self.mband, self.nkpt)
        wal_sbk = self.wal_sbk
        paw1_wal_sbk = self._read_wal_sbk(key="dos_fractions_paw1")
        pawt1_wal_sbk = self._read_wal_sbk(key="dos_fractions_pawt1")

        ebands = self.ebands
        kpoints = ebands.kpoints
        nband_sk = ebands.nband_sk
        eigens = ebands.eigens

        # Compute the linear mesh for DOS.
        epad = 1.0
        e_min = ebands.enemin() - epad
        e_max = ebands.enemax() + epad

        nw = int(1 + (e_max - e_min) / step)
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)
        totdos_al = np.zeros((self.natom, self.lsize, self.nsppol, nw))
        paw1dos_al = np.zeros((self.natom, self.lsize, self.nsppol, nw))
        pawt1dos_al = np.zeros((self.natom, self.lsize, self.nsppol, nw))

        if method == "gaussian":
            for spin in range(self.nsppol):
                for k, kpoint in enumerate(kpoints):
                    weight = kpoint.weight
                    for band in range(nband_sk[spin, k]):
                        gs = gaussian(mesh, width, center=eigens[spin,k,band])
                        for iatom in range(self.natom):
                            if not self.has_atom[iatom]: continue
                            for l in range(min(self.lmax_atom[iatom] + 1, mylsize)):
                                totdos_al[iatom, l, spin] += weight * gs * wal_sbk[iatom, l, spin, band, k]
                                paw1dos_al[iatom, l, spin] += weight * gs * paw1_wal_sbk[iatom, l, spin, band, k]
                                pawt1dos_al[iatom, l, spin] += weight * gs * pawt1_wal_sbk[iatom, l, spin, band, k]

        else:
            raise ValueError("Method %s is not supported" % method)

        # TOT = PW + AE - PS
        pwdos_al = totdos_al - paw1dos_al + pawt1dos_al

        # Build plot grid.
        nrows, ncols = np.count_nonzero(self.has_atom), self.lsize
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)
        ax_mat = np.reshape(ax_mat, (nrows, ncols))

        irow = -1
        for iatom in range(self.natom):
            if not self.has_atom[iatom]: continue
            irow += 1
            #for l in range(min(self.lmax_atom[iatom] + 1, mylsize)):
            for l in range(min(self.lsize, mylsize)):
                ax = ax_mat[irow, l]
                if l >= self.lmax_atom[iatom]+1:
                    # don't show this plots and cycle
                    ax.axis("off")
                    continue
                ax.grid(True)
                if l != 0:
                    ax.set_ylabel("")
                    # Only the first column show labels.
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(0)

                for spin in range(self.nsppol):
                    spin_sign = +1 if spin == 0 else -1
                    ax.plot(mesh, totdos_al[iatom, l, spin] * spin_sign, color="k",
                            label="Total" if (irow, l, spin) == (0, 0, 0) else None)
                    ax.plot(mesh, pwdos_al[iatom, l, spin] * spin_sign, color="r",
                            label="PW part" if (irow, l, spin) == (0, 0, 0) else None)
                    ax.plot(mesh, paw1dos_al[iatom, l, spin] * spin_sign, color="b",
                            label="AE-onsite" if (irow, l, spin) == (0, 0, 0) else None)
                    ax.plot(mesh, pawt1dos_al[iatom, l, spin] * spin_sign, color="g",
                            label="PS-onsite" if (irow, l, spin) == (0, 0, 0) else None)

        for ax in ax_mat[-1, :]:
            ax.set_xlabel('Energy (eV)')
            set_axlims(ax, xlims, "x")

    # TODO: THIS CODE IS STILL UNDER DEVELOPMENT
    #@add_fig_kwargs
    #def plot_pjdos_spinor(self, terms=("sigma_x", "sigma_y", "sigma_z"),
    #                      e0="fermie", method="gaussian", step=0.1, width=0.2,
    #                      blist=None,
    #                      #stacked=True,
    #                      ax=None, exchange_xy=False, xlims=None,
    #                      #with_info=True,
    #                      **kwargs):
    #    """
    #    Plot the PJ-DOS on a linear mesh.

    #    Args:
    #        terms:
    #        e0: Option used to define the zero of energy in the band structure plot. Possible values:
    #            - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy.
    #            -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
    #            -  None: Don't shift energies, equivalent to e0=0
    #        method: String defining the method for the computation of the DOS.
    #        step: Energy step (eV) of the linear mesh.
    #        width: Standard deviation (eV) of the gaussian.
    #        stacked: True if DOS partial contributions should be stacked on top of each other.
    #        ax: matplotlib axis, if None a new figure is generated.
    #        exchange_xy: True if the dos should be plotted on the x axis instead of y.
    #        xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
    #               or scalar e.g. `left`. If left (right) is None, default values are used

    #    Returns:
    #        `matplotlib` figure
    #    """
    #    # TODO: To be tested.
    #    if self.prtdos != 5 and self.nspinor != 2:
    #        print("This method assumes prtdos=5 and nspinor=2 but file has prtdos=%s and nspinor=%s"
    #            % (self.prtdos, self.nspinor))
    #        return None

    #    try:
    #        intg = self.get_dos_integrator(method, step, width)
    #    except Exception:
    #        msg = traceback.format_exc()
    #        msg += ("Error while trying to compute the DOS.\n"
    #                "Verify that the k-points form a homogenous sampling of the BZ.\n"
    #                "Returning None\n")
    #        print(msg)
    #        return None

    #    # Get energy mesh from total DOS and define the zero of energy
    #    # Note that the mesh is not not spin-dependent.
    #    e0 = self.ebands.get_e0(e0)
    #    mesh = intg.mesh.copy()
    #    mesh -= e0
    #    ebands = self.ebands
    #    mybands = range(ebands.mband) if blist is None else blist

    #    # Read data from file.
    #    # Dos fractions contains: [ndosfraction, nsppol, mband, nkpt].
    #    # where the first dimension stores: (up,up  up,down  down,up  down,down  sigma_x sigma_y sigma_z)
    #    w_sbk = self.reader.read_value("dos_fractions")
    #    term2idx = {"up-up":0, "up-down":1, "down-up":2, "down-down":3, "sigma_x":4, "sigma_y":5, "sigma_z":6}

    #    spin0 = 0
    #    terms = list_strings(terms)
    #    dos_terms = {term: np.zeros(len(mesh)) for term in terms}
    #    weights = np.array([k.weight for k in ebands.kpoints])
    #    for i, term in enumerate(terms):
    #        idx = term2idx[term]
    #        for band in mybands:
    #            gaussians_dos(dos_terms[term], mesh, width,
    #                          w_sbk[idx, spin0, band, :], ebands.eigens[spin0, :, band], weights)

    #    # Plot doses.
    #    ax, fig, plt = get_ax_fig_plt(ax=ax)
    #    for term, yvals in dos_terms.items():
    #        ax.plot(mesh, yvals, color=self.spinors2color[term], label=self.spinors2tex[term])

    #    ax.grid(True)
    #    ax.set_xlabel("Energy (eV)")
    #    set_axlims(ax, xlims, "x")
    #    ax.legend(loc="best", fontsize=fontsize, shadow=True)

    #    return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        #for fig in self.yield_ebands_figs(): yield fig
        if self.ebands.kpoints.is_path:
            yield self.ebands.kpoints.plot(show=False)
            yield self.plot_fatbands_lview(show=False)
            yield self.plot_fatbands_typeview(show=False)
        else:
            yield self.plot_pjdos_lview(show=False)
            yield self.plot_pjdos_typeview(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("fbnc = abilab.abiopen('%s')\nprint(fbnc)" % self.filepath),
            nbv.new_code_cell("fbnc.structure"),
            nbv.new_code_cell("#fbnc.ebands.kpoints.plot();"),
            nbv.new_code_cell("xlims = (None, None)\nylims = (None, None)"),
            nbv.new_code_cell("#fbnc.ebands.plot(ylims=ylims);"),
            nbv.new_code_cell("fbnc.ebands.boxplot();"),
        ])

        if self.prtdos == 3:
            nb.cells.extend([
                nbv.new_markdown_cell("## Fatbands plots with L-character \n(require `prtdos=3`)"),
                nbv.new_code_cell("fbnc.plot_fatbands_typeview(ylims=ylims);"),
                nbv.new_code_cell("fbnc.plot_fatbands_lview(ylims=ylims);"),
                nbv.new_code_cell("fbnc.plot_fatbands_siteview(ylims=ylims);"),
            ])

        if self.prtdos == 3 and self.prtdosm != 0:
            nb.cells.extend([
                nbv.new_markdown_cell("## Fatbands plots with LM-character \n(require `prtdos = 3 and prtdosm != 0`)"),
                nbv.new_code_cell("fbnc.plot_fatbands_mview(iatom=0, ylims=ylims);"),
            ])

        if self.prtdos == 5 and self.nspinor == 2:
            nb.cells.extend([
                nbv.new_markdown_cell("## Fatbands plots with Spin-character \n(require `prtdos = 5 and nspinor == 2`)"),
                nbv.new_code_cell("fbnc.plot_fatbands_spinor(ylims=ylims);"),
            ])

        if self.prtdos == 3 and self.ebands.kpoints.is_mpmesh:
            nb.cells.extend([
                nbv.new_markdown_cell("## L-DOSes plots \n(require `prtdos = 3` and BZ sampling)"),
                nbv.new_code_cell("fbnc.plot_pjdos_lview(xlims=xlims);"),
                nbv.new_code_cell("fbnc.plot_pjdos_typeview(xlims=xlims);"),
            ])

        if self.prtdos == 3 and self.ebands.kpoints.is_path:
            nb.cells.extend([
                nbv.new_markdown_cell("## L-DOSes with fatbands\n"
                                     "(require `prtdos=3`, `fbnc` must contain a k-path, "
                                     "`pjdosfile` is a `FATBANDS.nc` file with a BZ sampling)"),
                nbv.new_code_cell("fbnc.plot_fatbands_with_pjdos(pjdosfile=None, ylims=ylims, view='type');"),
            ])

        if self.usepaw == 1 and self.prtdos != 3:
            nb.cells.extend([
                nbv.new_markdown_cell("## PAW L-DOS decomposed into smooth PW part, AE and PS terms"),
                nbv.new_code_cell("fbnc.plot_pawdos_terms();"),
            ])

        return self._write_nb_nbpath(nb, nbpath)


class _DosIntegrator(object):
    """
    This object is responsible for the integration of the DOS/PJDOS.
    It's an internal object that should not be instantiated directly outside of this module.
    PJDOSes are computed lazily and stored in the integrator so that we can reuse the results
    if needed.
    """
    def __init__(self, fbfile, method, step, width):
        """
        """
        self.fbfile, self.method, self.step, self.width = fbfile, method, step, width

        # Compute Total DOS from ebands and define energy mesh.
        self.edos = fbfile.ebands.get_edos(method=method, step=step, width=width)
        self.mesh = self.edos.spin_dos[0].mesh

    #@lazy_property
    #def site_edos(self):
    #    """Array [natom, nsppol, lmax**2]"""

    @lazy_property
    def symbols_lso(self):
        """
        """
        fbfile, ebands = self.fbfile, self.fbfile.ebands

        # Compute l-decomposed PJDOS for each type of atom.
        symbols_lso = OrderedDict()
        if self.method == "gaussian":

            for symbol in fbfile.symbols:
                lmax = fbfile.lmax_symbol[symbol]
                wlsbk = fbfile.get_wl_symbol(symbol)
                lso = np.zeros((fbfile.lsize, fbfile.nsppol, len(self.mesh)))
                for spin in range(fbfile.nsppol):
                    for k, kpoint in enumerate(ebands.kpoints):
                        weight = kpoint.weight
                        for band in range(ebands.nband_sk[spin, k]):
                            e = ebands.eigens[spin, k, band]
                            for l in range(lmax + 1):
                                lso[l, spin] += wlsbk[l, spin, band, k] *\
                                                weight * gaussian(self.mesh, self.width, center=e)
                symbols_lso[symbol] = lso

        else:
            raise ValueError("Method %s is not supported" % self.method)

        return symbols_lso

    @lazy_property
    def ls_stackdos(self):
        """
        Compute ``ls_stackdos`` datastructure for stacked DOS.
        ls_stackdos maps (l, spin) onto a numpy array [nsymbols, nfreqs] where
        [isymb, :] contains the cumulative sum of the PJDOS(l,s) up to symbol isymb.
        """
        fbfile = self.fbfile
        from itertools import product
        dls = defaultdict(dict)
        for symbol, lso in self.symbols_lso.items():
            for l, spin in product(range(fbfile.lmax_symbol[symbol]+1), range(fbfile.nsppol)):
                dls[(l, spin)][symbol] = lso[l, spin]

        ls_stackdos = {}
        nsymb = len(fbfile.symbols)
        for (l, spin), dvals in dls.items():
            arr = np.zeros((nsymb, len(self.mesh)))
            for isymb, symbol in enumerate(fbfile.symbols):
                arr[isymb] = dvals[symbol]
            ls_stackdos[(l, spin)] = arr.cumsum(axis=0)

        return ls_stackdos

    def get_lstack_symbol(self, symbol, spin):
        """
        Return |numpy-array| with the cumulative sum over l for a given
        atom type (specified by the chemical symbol ``symbol``) and spin.
        """
        lso = self.symbols_lso[symbol][:, spin]
        return lso.cumsum(axis=0)

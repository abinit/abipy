"""Phonon density of states."""
from __future__ import print_function, division

import os
import collections
import numpy as np

from abipy.core.func1d import Function1D
from abipy.iotools import ETSF_Reader, AbinitNcFile

__all__ = [
    "PhononDOS",
    "PHDOS_Reader",
    "PHDOS_File",
]


class PhononDOS(object):
    """This object stores the phonon density of states."""

    def __init__(self, mesh, values):
        """
        Args:
            mesh:
                array-like object with the points of the mesh.
            mesh:
                array-like object with the DOS values.

        .. note:
            mesh is given in eV, values are in states/eV.
        """
        self.dos = Function1D(mesh, values)
        self.idos = self.dos.integral()

    def plot_ax(self, ax, what="d", exchange_xy=False, *args, **kwargs):
        """
        Helper function to plot the data on the axis ax.

        Args:
            ax:
                matplotlib axis
            what:
                string selecting the quantity to plot:
                "d" for DOS, "i" for IDOS. chars can be concatenated
                hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy:
                True to exchange exis
            args, kwargs:
                Options passes to matplotlib.

        Return value is a list of lines that were added.
        """
        opts = [c.lower() for c in what]

        cases = {"d": self.dos,
                 "i": self.idos}

        lines = list()
        for c in opts:
            f = cases[c]
            ls = f.plot_ax(ax, exchange_xy=exchange_xy, *args, **kwargs)
            lines.extend(ls)

        return lines

    def plot(self, *args, **kwargs):
        """
        Plot DOS and IDOS.

        Args:
            args:
                Positional arguments passed to :mod:`matplotlib`.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure.
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        gspec = GridSpec(2, 1, height_ratios=[1,2])
        ax1 = plt.subplot(gspec[0])
        ax2 = plt.subplot(gspec[1])

        for ax in (ax1, ax2):
            ax.grid(True)

        ax2.set_xlabel('Energy [eV]')

        ax1.set_ylabel("IDOS")
        ax2.set_ylabel("DOS")

        if title:
            ax1.set_title(title)

        self.plot_ax(ax1, what="i", *args, **kwargs)
        self.plot_ax(ax2, what="d", *args, **kwargs)

        if show:
            plt.show()

        fig = plt.gcf()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


class PHDOS_Reader(ETSF_Reader):
    """
    This object reads data from the PHDOS.nc file produced by anaddb.

    .. note:
            Frequencies are in eV, DOSes are in states/eV.
    """
    def _lazy_get(self, varname):
        """Helper function used to create lazy properties."""
        hiddename = "__" + varname
        try:
            return getattr(self, hiddename)
        except AttributeError:
            setattr(self, hiddename, self.read_value(varname))
            return self._lazy_get(varname)

    @property
    def wmesh(self):
        """The frequency mesh in eV."""
        return self._lazy_get("wmesh")

    @property
    def pjdos_type(self):
        """DOS projected over atom types e.g. pjdos_type(ntypat,nomega)."""
        return self._lazy_get("pjdos_type")

    @property
    def pjdos_rc_type(self):
        """DOS projected over atom types and reduced directions e.g. pjdos_type(3,ntypat,nomega)."""
        return self._lazy_get("pjdos__rc_type")

    @property
    def pjdos(self):
        """DOS projected over atoms and reduced directions pjdos(natom,3,nomega)."""
        return self._lazy_get("pjdos")

    @property
    def structure(self):
        """The crystalline structure."""
        if not hasattr(self, "__structure"):
            self.__structure = self.read_structure()
        return self.__structure

    def read_phdos(self):
        """Return the PhononDOS."""
        return PhononDOS(self.wmesh, self.read_value("phdos"))

    def read_pjdos_type(self, symbol):
        """
        The contribution to the DOS due to the atoms of given chemical symbol.
        pjdos_type(ntypat,nomega)
        """
        type_idx = self.typeidx_from_symbol(symbol)
        return PhononDOS(self.wmesh, self.pjdos_type[type_idx])

    # def read_pjdos(self, atom_idx=None):
    #     """
    #     projected DOS (over atoms)
    #     """
    #     return self.read_value("phonon_frequencies")

    # def read_pjdos_rc_type(self, symbol=None):
    #     """
    #     phdos(3,ntypat,nomega)
    #     phonon DOS contribution arising from a particular atom-type
    #     decomposed along the three reduced directions.
    #     """
    #     return self.read_value("phonon_frequencies")


class PHDOS_File(AbinitNcFile):
    """
    Container object storing the different DOSes stored in the
    PHDOS.nc file produced by anaddb. Provides helper function
    to visualize/extract data.
    """
    def __init__(self, filepath):
        # Open the file, read data and create objects.
        super(PHDOS_File, self).__init__(filepath)

        pjdos_type_dict = collections.OrderedDict()
        with PHDOS_Reader(filepath) as r:
            self.structure = r.structure
            self.wmesh = r.wmesh
            self.phdos = r.read_phdos()

            for symbol in r.chemical_symbols:
                #print(symbol, ncdata.typeidx_from_symbol(symbol))
                pjdos_type_dict[symbol] = r.read_pjdos_type(symbol)

        self.pjdos_type_dict = pjdos_type_dict

    def get_structure(self):
        """Returns the `Structure` object."""
        return self.structure

    def plot_pjdos_type(self, colormap="jet", **kwargs):
        """
        Stacked Plot of the  projected DOS (projection is for atom types)

        Args:
            colormap
                Have a look at the colormaps here and decide which one you'd like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            matplotlib figure.
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt
        fig = plt.figure()

        ax = fig.add_subplot(1,1,1)
        ax.grid(True)

        xlim = kwargs.pop("xlim", None)
        if xlim is not None:
            ax.set_xlim(xlim)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None:
            ax.set_ylim(ylim)

        ax.set_xlabel('Frequency [eV]')
        ax.set_ylabel('PJDOS [states/eV]')

        if title:
            ax.set_title(title)

        # Type projected DOSes.
        num_plots = len(self.pjdos_type_dict)
        cumulative = np.zeros(len(self.wmesh))
        for i, (symbol, pjdos) in enumerate(self.pjdos_type_dict.items()):
            f = pjdos.dos
            x, y = f.mesh, f.values
            color = plt.get_cmap(colormap)(float(i)/(num_plots-1))
            ax.plot(x, cumulative + y, lw=2, label=symbol, color=color)
            ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=0.7)
            cumulative += y

        # Total PHDOS
        f = self.phdos.dos
        x, y = f.mesh, f.values
        ax.plot(x, y, lw=2, label="Total PHDOS", color='black')

        ax.legend(loc="best")

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

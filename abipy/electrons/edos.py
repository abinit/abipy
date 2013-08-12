"""Electronic density of states."""
from __future__ import division, print_function

import numpy as np
import collections

from abipy.core.func1d import Function1D

__all__ = [
    "ElectronDOS",
    "ElectronDOSPlotter",
]


class ElectronDOS(object):
    """This object stores the electronic density of states."""

    def __init__(self, mesh, spin_dos):
        """
        Args:
            mesh:
                array-like object with the mesh points.
            spin_dos:
                array-like object with the DOS value for the different spins.
                Example::

                     spin_dos[nw] if spin-unpolarized.
                     spin_dos[nsppol, nw] if spin-polarized case.

        .. note:
            mesh is given in eV, spin_dos is in states/eV.
        """
        spin_dos = np.atleast_2d(spin_dos)
        assert spin_dos.ndim == 2
        self.nsppol = len(spin_dos)

        # Save DOS and IDOS for each spin.
        sumv = np.zeros(len(mesh))
        self.spin_dos, self.spin_idos = [], []
        for values in spin_dos:
            sumv += values
            f = Function1D(mesh, values)
            self.spin_dos.append(f)
            self.spin_idos.append(f.integral())

        # Total DOS and IDOS.
        if self.nsppol == 1: sumv = 2 * sumv
        self.tot_dos = Function1D(mesh, sumv)
        self.tot_idos = self.tot_dos.integral()

    def dos_idos(self, spin=None):
        """
        Returns DOS and IDOS for given spin. Total DOS and IDOS if spin is None.
        """
        if spin is None:
            return self.tot_dos, self.tot_idos
        else:
            return self.spin_dos[spin], self.spin_idos[spin]

    def find_mu(self, nelect, spin=None, num=500, atol=1.e-5):
        """Finds the chemical potential given the number of electrons."""
        idos = self.tot_idos if spin is None else self.spin_idos[spin]

        # Cannot use bisection because DOS might be negative due to smearing.
        # This one is safer albeit slower.
        for i, (ene, intg) in enumerate(idos):
            if intg > nelect:
                break
        else:
            raise ValueError("Cannot find I(e) such that I(e) > nelect")

        # Now use spline to get a more accurate mu (useful if mesh is coarse)
        e0, e1 = idos.mesh[i-1], idos.mesh[i]
        idos_spline = idos.spline
        for mu in np.linspace(e0, e1, num=num):
            if abs(idos_spline(mu) - nelect) < atol:
                return mu
        else:
            raise RuntimeError("Cannot find mu, try to increase num and/or atol")

    def plot_ax(self, ax, spin=None, what="d", exchange_xy=False, *args, **kwargs):
        """
        Helper function to plot the data on the axis ax.

        Args:
            ax:
                matplotlib axis.
            spin:
                selects the spin component, None for total DOS, IDOS.
            what:
                string selecting what will be plotted:
                "d" for DOS, "i" for IDOS. chars can be concatenated
                hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy:
                True to exchange axis.
            args, kwargs:
                Options passes to matplotlib.

        Return value is a list of lines that were added.
        """
        #print("spin",spin)
        dosf, idosf = self.dos_idos(spin=spin)

        opts = [c.lower() for c in what]

        lines = []
        for c in opts:
            if c == "d": f = dosf
            if c == "i": f = idosf
            ls = f.plot_ax(ax, exchange_xy=exchange_xy, *args, **kwargs)
            lines.extend(ls)
        return lines

    def plot(self, spin=None, *args, **kwargs):
        """
        Plot DOS and IDOS.

        Args:
            spin:
                Selects the spin component, None if total DOS is wanted.
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
            matplotlib figure.
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

        ax1.set_ylabel("TOT IDOS" if spin is None else "IDOS (spin %s)" % spin)
        ax2.set_ylabel("TOT DOS" if spin is None else "DOS (spin %s)" % spin)

        if title:
            ax1.set_title(title)

        self.plot_ax(ax1, spin=spin, what="i", *args, **kwargs)
        self.plot_ax(ax2, spin=spin, what="d", *args, **kwargs)

        if show:
            plt.show()

        fig = plt.gcf()
        if savefig is not None:
            fig.savefig(savefig)

        return fig


class ElectronDOSPlotter(object):
    """
    Class for plotting DOSes.
    """
    def __init__(self):
        self._doses = collections.OrderedDict()

    def add_dos(self, label, dos):
        """
        Adds a DOS for plotting.

        Args:
            label:
                label for the MDF. Must be unique.
            dos:
                MacroscopicDielectricFunction object.
        """
        if label in self._doses:
            raise ValueError("label %s is already in %s" % (label, self._doses.keys()))

        self._doses[label] = dos

    def add_dos_dict(self, dos_dict, key_sort_func=None):
        """
        Add a dictionary of DOSes, with an optional sorting function for the keys.

        Args:
            dos_dict:
                dict of {label: dos}
            key_sort_func:
                function used to sort the dos_dict keys.
        """
        if key_sort_func:
            keys = sorted(dos_dict.keys(), key=key_sort_func)
        else:
            keys = dos_dict.keys()

        for label in keys:
            self.add_dos(label, dos_dict[label])

    def plot(self, *args, **kwargs):
        """
        Get a matplotlib plot showing the DOSes.



        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        xlim            x-axis limits. None (default) for automatic determination.
        ylim            y-axis limits.  None (default) for automatic determination.
        ==============  ==============================================================
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)
        import matplotlib.pyplot as plt
        fig = plt.figure()

        ax = fig.add_subplot(1,1,1)
        ax.grid(True)

        xlim = kwargs.pop("xlim", None)
        if xlim is not None: ax.set_xlim(xlim)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None: ax.set_ylim(ylim)

        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('DOS [states/eV]')

        if title is not None:
            ax.set_title(title)

        lines, legends = [], []
        for (label, dos) in self._doses.items():
            l = dos.plot_ax(ax, *args, **kwargs)[0]

            lines.append(l)
            legends.append("DOS: %s" % label)

        # Set legends.
        ax.legend(lines, legends, 'upper right', shadow=True)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

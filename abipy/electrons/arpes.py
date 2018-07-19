# coding: utf-8
"""
Arpese Plotter (still under development)
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

#from collections import OrderedDict
from scipy.interpolate import UnivariateSpline
#from monty.string import marquee # is_string, list_strings
#from monty.functools import lazy_property
from monty.collections import dict2namedtuple
#from monty.termcolor import cprint
from abipy.core.mixins import Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons import ElectronBands
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_ax3d_fig_plt, get_axarray_fig_plt #set_axlims,


class ArpesPlotter(Has_Structure, Has_ElectronBands, NotebookWriter):
    """

    Usage example:

    .. code-block:: python

        with abilab.abiopen("foo_ABIWAN.nc") as abiwan:
            print(abiwan)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ArpesPlotter
    """
    @classmethod
    def model_from_ebands(cls, ebands, tmesh=(0, 300, 600), poorman_polaron=False):
        ebands = ElectronBands.as_ebands(ebands)

        ntemp = len(tmesh)
        nwr =  1000
        wr_step = 0.01

        aw = np.empty((ebands.nsppol, ebands.nkpt, ebands.mband, ntemp, nwr))
        aw_meshes = np.empty((ebands.nsppol, ebands.nkpt, ebands.mband, nwr))
        #aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
        #aw_meshes: [max_nbcalc, nkcalc, nsppol] array with energy mesh in eV
        from abipy.tools.numtools import lorentzian
        from scipy.integrate import cumtrapz
        for spin in ebands.spins:
            for ik, kpt in enumerate(ebands.kpoints):
                for band in range(ebands.nband_sk[spin, ik]):
                    e0 = ebands.eigens[spin, ik, band]
                    emin = e0 - wr_step * (nwr // 2)
                    emax = e0 + wr_step * (nwr // 2)
                    emesh = np.linspace(emin, emax, num=nwr)
                    aw_meshes[spin, ik, band] = emesh
                    # Naive model: lorentzian centered on KS energy with T-dep broadening
                    for itemp, temp in enumerate(tmesh):
                        width = 0.2 + (temp / 300) * 0.2
                        avals = lorentzian(emesh, width=width, center=e0, height=None)
                        if poorman_polaron:
                            if band in (1, 2, 3) and kpt.norm < 0.3:
                                avals += 1.1 * lorentzian(emesh, width=0.1 * width, center=e0 - 0.4, height=None)
                            avals /= cumtrapz(avals, x=emesh)[-1]
                        aw[spin, ik, band, itemp] = avals

        return cls(ebands, aw, aw_meshes, tmesh)

    def __init__(self, ebands, aw, aw_meshes, tmesh):
        """
        Args:
            ebands: |ElectronBands| object
            aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
            aw_meshes: [max_nbcalc, nkcalc, nsppol] array with energy mesh in eV
            tmesh: Temperature mesh in Kelvin.

        .. note::

            - Each (spin, k-point, band) can have different frequency meshes
              (usually centered on the initial KS energy).
              The plotter will spline the data.

            - The order of the k-points in ebands and aw should be the same.

            - The treatment of bands if complicated by the fact that we can have
              different nband(k) whose indices are not necessarily aligned.
              Consider, for example. what happes if we use symsigma or arbitrary bdgw.
              Use MaskedArray or metadata with start, stop?
        """
        self._ebands = ebands
        self.aw = aw
        self.aw_meshes = aw_meshes
        self.tmesh = tmesh
        self.ntemp = len(tmesh)
        #assert

        # Options passed to UnivariateSpline
        self.ext, self.k, self.s = "zeros", 3, 0

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @property
    def ebands(self):
        """|ElectronBands| object."""
        return self._ebands

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        lines = []; app = lines.append
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))

        #if verbose > 1:
        return "\n".join(lines)

    def with_points_along_path(self, frac_bounds=None, knames=None, dist_tol=1e-12):
        """
        Args:
            frac_bounds: [M, 3] array  with the vertexes of the k-path in reduced coordinates.
                If None, the k-path is automatically selected from the structure.
            knames: List of strings with the k-point labels defining the k-path. It has precedence over frac_bounds.
            dist_tol: A point is considered to be on the path if its distance from the line
                is less than dist_tol.
        """
        r = self.ebands.with_points_along_path(frac_bounds=frac_coords, knames=knames, dist_tol=dist_tol)
        # Transfer data using r.ik_new2prev table.
        return self.__class__(r.ebands,
                              aw=self.aw[:, :, :, r.ik_new2prev, :].copy(),
                              aw_meshes=self.aw_meshes[:, r.ik_new2prev, :].copy(),
                              tmesh=self.tmesh)

    def interpolate(self):
        new_ebands = self.ebands.interpolate(lpratio=5, vertices_names=None, line_density=20,
                            kmesh=None, is_shift=None, filter_params=None, verbose=0)

        # Build interpolator.
        #from abipy.core.skw import SkwInterpolator
        #my_kcoords = [k.frac_coords for k in self.kpoints]
        #cell = (self.structure.lattice.matrix, self.structure.frac_coords,
        #        self.structure.atomic_numbers)

        #skw = SkwInterpolator(lpratio, my_kcoords, self.eigens, self.fermie, self.nelect,
        #                      cell, fm_symrel, self.has_timrev,
        #                      filter_params=filter_params, verbose=verbose)

        # Interpolate energies.
        #eigens_kpath = skw.interp_kpts(kfrac_coords).eigens

        return self.__class__(new_ebands, new_aw, aw_meshes, self.tmesh)

    def get_emesh_eminmax(self, estep):
        """Compute linear mesh covering entire energy range."""
        emin = self.ebands.enemin()
        emin -= 0.1 * abs(emin)
        emax = self.ebands.enemax()
        emax += 0.1 * abs(emax)
        return np.arange(emin, emax, estep), emin, emax

    def get_data_nmtuple(self, itemp, estep, spins=None):
        nkpt = self.ebands.nkpt
        spins = range(self.ebands.nsppol) if spins is None else spins

        emesh, emin, emax = self.get_emesh_eminmax(estep)
        nene = len(emesh)
        #print("nkpt", nkpt, "nene", nene)
        data = np.zeros((nkpt, nene))

        # aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
        for spin in spins:
            for ik in range(nkpt):
                for band in range(self.ebands.nband_sk[spin, ik]):
                    w = self.aw_meshes[spin, ik, band]
                    aw = self.aw[spin, ik, band, itemp]
                    data[ik] += UnivariateSpline(w, aw, k=self.k, s=self.s, ext=self.ext)(emesh)

        return dict2namedtuple(data=data, emesh=emesh, emin=emin, emax=emax, spins=spins, nkpt=nkpt)

    def get_atw(self, wmesh, spin, ikpt, band_inds, temp_inds):
        ntemp, nene = len(temp_inds), len(wmesh)
        atw = np.zeros((ntemp, nene))
        for band in range(self.ebands.nband_sk[spin, ikpt]):
            if band_inds is not None and band not in band_inds: continue
            w = self.aw_meshes[spin, ikpt, band]
            for it, itemp in enumerate(temp_inds):
                aw = self.aw[spin, ikpt, band, itemp]
                atw[it] += UnivariateSpline(w, aw, k=self.k, s=self.s, ext=self.ext)(wmesh)

        return atw

    @add_fig_kwargs
    def plot_ekmap_temps(self, temp_inds=None, spins=None, estep=0.02, with_colorbar=True,
                        ylims=None, fontsize=8, **kwargs):
        """
        Plot (k, e) color maps for different temperatures.

        Args:
            fontsize (int): fontsize for titles and legend

        Return: |matplotlib-Figure|
        """
        temp_inds = range(self.ntemp) if temp_inds is None else temp_inds
        # Build plot grid.
        num_plots, ncols, nrows = len(temp_inds), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        # Don't show the last ax if numeb is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for itemp, ax in zip(temp_inds, ax_list):
            self.plot_ekmap_itemp(itemp=itemp, spins=spins, estep=estep, ax=ax, ylims=ylims,
                    with_colorbar=with_colorbar, show=False, **kwargs)
            ax.set_title("T = %.1f K" % self.tmesh[itemp], fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_ekmap_itemp(self, itemp=0, spins=None, estep=0.02, ax=None, ylims=None, with_colorbar=True, **kwargs):
        """
        Plot (k, e) color map for given temperature.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                - Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                - None: Don't shift energies, equivalent to ``e0 = 0``.
            spins: Selects the spin to be plotted, None if all spins.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            with_colorbar: True to add color bar.
            kwargs: options passed to ``ax.imshow``.

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        a = self.get_data_nmtuple(itemp, estep, spins=spins)

        cmap = "jet"
        img = ax.imshow(a.data.T, origin="lower", extent=[0, a.nkpt, a.emin, a.emax],
                  cmap=cmap,
                  interpolation="bilinear",
                  #interpolation="spline36",
                  #interpolation="bicubic",
                  #vmin=0, vmax=np.abs(data).max()
                  )
        self.ebands.plot(ax=ax, e0=0, show=False, color="w", lw=1, ls="--")

        #ax.set_zlabel(r"$A(\omega)$")
        #self.ebands.plot(ax=ax, e0=0, color="r", lw=1)
        #ax.imshow(data, cmap=None, norm=None, aspect=None, interpolation=None, alpha=None, vmin=None, vmax=None,
        #       origin=None, extent=None, shape=None, filternorm=1, filterrad=4.0, imlim=None, resample=None,
        #       url=None, hold=None, data=None, **kwargs)

        if with_colorbar:
            # Make a color bar
            #plt.colorbar(img, cmap=cmap)
            # https://stackoverflow.com/questions/13310594/positioning-the-colorbar
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            #cax = divider.new_vertical(size="5%", pad=0.1, pack_start=True)
            #cax = divider.new_horizontal(size="5%", pad=0.1, pack_start=True)

            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            # https://matplotlib.org/2.0.2/mpl_toolkits/axes_grid/users/overview.html#axesdivider
            cax = divider.append_axes("right", size="5%", pad=0.05)

            #fig.add_axes(cax)
            #fig.colorbar(img, cax=cax, ax=ax, orientation="horizontal")
            fig.colorbar(img, cax=cax, ax=ax)

        return fig

    @add_fig_kwargs
    def plot_ak_vs_temp(self, temp_inds=None, spins=None, band_inds=None, kpt_inds=None,
                        apad=1.0, estep=0.02, colormap="jet", fontsize=8, **kwargs):
        """

        Args:
            temp_inds:
            spins:
            band_inds:
            kpt_inds:
            estep:
            colormap
            fontsize (int): fontsize for titles and legend

        Return: |matplotlib-Figure|
        """
        temp_inds = range(self.ntemp) if temp_inds is None else temp_inds
        ntemp = len(temp_inds)
        spins = range(self.ebands.nsppol) if spins is None else spins
        kpt_inds = range(self.ebands.nkpt) if kpt_inds is None else kpt_inds
        nkpt = len(kpt_inds)

        xs, emin, emax = self.get_emesh_eminmax(estep)
        nene = len(xs)

        num_plots, ncols, nrows = nkpt, 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        # Build plot grid.
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = np.array(ax_list).ravel()
        cmap = plt.get_cmap(colormap)

        for isp, spin in enumerate(spins):
            spin_sign = +1 if spin == 0 else -1
            for ik, (ikpt, ax) in enumerate(zip(kpt_inds, ax_list)):
                ax.grid(True)
                atw = self.get_atw(xs, spin, ikpt, band_inds, temp_inds)
                for it, itemp in enumerate(temp_inds):
                    ys = spin_sign * atw[it] + (it * apad)
                    ax.plot(xs, ys, lw=2, alpha=0.8, color=cmap(float(it) / ntemp),
                            label = "T = %.1f K" % self.tmesh[itemp] if (ik, isp) == (0, 0) else None)

                if spin == 0:
                    kpt = self.ebands.kpoints[ikpt]
                    ax.set_title("k:%s" % (repr(kpt)), fontsize=fontsize)

                if (ik, isp) == (0, 0):
                    ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    #@add_fig_kwargs
    #def plot_ak(self, temp_inds=None, spins=None, band_inds=None, kpt_inds=None,
    #           apad=1.0, estep=0.02, colormap="jet", fontsize=8, **kwargs):
    #    """

    #    Args:
    #        temp_inds:
    #        spins:
    #        band_inds:
    #        kpt_inds:
    #        estep:
    #        colormap
    #        fontsize (int): fontsize for titles and legend

    #    Return: |matplotlib-Figure|
    #    """
    #    temp_inds = range(self.ntemp) if temp_inds is None else temp_inds
    #    ntemp = len(temp_inds)
    #    spins = range(self.ebands.nsppol) if spins is None else spins
    #    kpt_inds = range(self.ebands.nkpt) if kpt_inds is None else kpt_inds
    #    nkpt = len(kpt_inds)

    #    xs, emin, emax = self.get_emesh_eminmax(estep)
    #    nene = len(xs)

    #    num_plots, ncols, nrows = nkpt, 1, 1
    #    if num_plots > 1:
    #        ncols = 2
    #        nrows = (num_plots // ncols) + (num_plots % ncols)

    #    # Build plot grid.
    #    ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
    #                                            sharex=True, sharey=True, squeeze=False)
    #    ax_list = np.array(ax_list).ravel()
    #    cmap = plt.get_cmap(colormap)

    #    for isp, spin in enumerate(spins):
    #        spin_sign = +1 if spin == 0 else -1
    #        for ik, (ikpt, ax) in enumerate(zip(kpt_inds, ax_list)):
    #            ax.grid(True)
    #            atw = self.get_atw(xs, spin, ikpt, band_inds, temp_inds)
    #            for it, itemp in enumerate(temp_inds):
    #                ys = spin_sign * atw[it] + (it * apad)
    #                ax.plot(xs, ys, lw=2, alpha=0.8, color=cmap(float(it) / ntemp),
    #                        label = "T = %.1f K" % self.tmesh[itemp] if (ik, isp) == (0, 0) else None)

    #            if spin == 0:
    #                kpt = self.ebands.kpoints[ikpt]
    #                ax.set_title("k:%s" % (repr(kpt)), fontsize=fontsize)

    #            if (ik, isp) == (0, 0):
    #                ax.legend(loc="best", fontsize=fontsize, shadow=True)

    #    return fig

    @add_fig_kwargs
    def plot_3dlines(self, itemp=0, estep=0.02, spins=None, band_inds=None, ax=None, **kwargs):
        ax, fig, plt = get_ax3d_fig_plt(ax=ax)

        xs, emin, emax = self.get_emesh_eminmax(estep)
        nene = len(xs)
        nkpt = self.ebands.nkpt
        cmap = plt.get_cmap("jet")
        # TODO: Reshift everything?
        if band_inds is not None:
            band_inds = set(band_inds)

        # aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
        spins = range(self.ebands.nsppol) if spins is None else spins
        for spin in spins:
            for ik in range(nkpt):
                ys = np.ones(nene) * ik
                zs = np.zeros(nene)
                for band in range(self.ebands.nband_sk[spin, ik]):
                    if band_inds is not None and band not in band_inds: continue
                    w = self.aw_meshes[spin, ik, band]
                    aw = self.aw[spin, ik, band, itemp]
                    zs += UnivariateSpline(w, aw, k=self.k, s=self.s, ext=self.ext)(xs)

                ax.plot(ys, xs, zs, color="k", lw=1, alpha=0.8) #cmap(float(ik) / nkpt))

                # Code to convert data in 3D polygons
                # See https://stackoverflow.com/questions/33641551/vertically-fill-3d-matplotlib-plot
                #h = 0.0
                #from mpl_toolkits.mplot3d.art3d import Poly3DCollection
                #xs = xs.copy()
                #ys = xs.copy()
                #zs = zs.copy()
                #v = []
                #for k in range(0, len(xs) - 1):
                #    x = [xs[k], xs[k+1], xs[k+1], xs[k]]
                #    y = [ys[k], ys[k+1], ys[k+1], ys[k]]
                #    z = [zs[k], zs[k+1],       h,     h]
                #    v.append(list(zip(x, y, z)))
                #poly3dCollection = Poly3DCollection(v)
                #ax.add_collection3d(poly3dCollection)

        ax.set_zlabel(r"$A(\omega)$")
        self.ebands.plot(ax=ax, e0=0, color="r", lw=1)

        return fig

    @add_fig_kwargs
    def plot_surface(self, itemp=0, estep=0.02, spins=None, ax=None, **kwargs):
        ax, fig, plt = get_ax3d_fig_plt(ax=ax)

        xs, emin, emax = self.get_emesh_eminmax(estep)
        nene = len(xs)
        nkpt = self.ebands.nkpt
        cmap = plt.get_cmap("jet")

        # aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
        spins = range(self.ebands.nsppol) if spins is None else spins
        spin = 0
        zs = np.zeros((nkpt, nene))
        ys = np.arange(nkpt)
        for ik in range(nkpt):
            for band in range(self.ebands.nband_sk[spin, ik]):
                w = self.aw_meshes[spin, ik, band]
                aw = self.aw[spin, ik, band, itemp]
                zs[ik] += UnivariateSpline(w, aw, k=self.k, s=self.s, ext=self.ext)(xs)

        # Plot the surface.
        xs, ys = np.meshgrid(xs, ys)

        #surf = ax.plot_surface(ys, xs, zs,
        #                       rstride=1, cstride=1,
        #                       linewidth=0, antialiased=True,
        #                       #cmap=cmap,
        #                       )
        #cb = fig.colorbar(surf, shrink=0.5)

        ax.plot_wireframe(ys, xs, zs, rstride=1, cstride=5)

        ax.set_zlabel(r"$A(\omega)$")
        self.ebands.plot(ax=ax, e0=0, color="r", lw=1)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield None
        # TODO
        #yield self.combiplot(show=False)
        #yield self.gridplot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        # Use pickle files for data persistence.
        tmpfile = self.pickle_dump()

        # TODO
        nb.cells.extend([
            nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("plotter = abilab.ArpesPlotter.pickle_load('%s')" % tmpfile),
            nbv.new_code_cell("print(plotter)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


if __name__ == "__main__":
    import sys
    plotter = ArpesPlotter.model_from_ebands(sys.argv[1]) #, aw, aw_meshes, tmesh)
    print(plotter.to_string(verbose=2))
    #plotter.plot_ekmap_itemp(itemp=0, estep=0.05, with_colorbar=True)
    plotter.plot_ekmap_temps(estep=0.05, with_colorbar=True)
    #plotter.plot_3dlines(itemp=0, estep=0.05, band_inds=[1, 2, 3])
    #plotter.plot_ak()
    #plotter.plot_ak_vs_temp()
    #plotter.plot_surface(istep=0, estep=0.05)

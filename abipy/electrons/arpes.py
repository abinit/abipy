# coding: utf-8
""""""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

#from collections import OrderedDict
from scipy.interpolate import UnivariateSpline
#from monty.string import marquee # is_string, list_strings,
#from monty.functools import lazy_property
from monty.collections import dict2namedtuple
#from monty.termcolor import cprint
#import pymatgen.core.units as units
from abipy.core.mixins import Has_Structure, Has_ElectronBands #, NotebookWriter
#from abipy.core.kpoints import KpointList, is_diagonal, find_points_along_path
from abipy.electrons import ElectronBands
from abipy.tools.plotting import set_axlims, add_fig_kwargs, get_ax_fig_plt, get_ax3d_fig_plt


class ArpesPlotter(Has_Structure, Has_ElectronBands):

    @classmethod
    def model_from_ebands(cls, ebands):
        ebands = ElectronBands.as_ebands(ebands)

        tmesh = [0]
        ntemp = len(tmesh)
        nwr =  1000
        wr_step = 0.01

        aw = np.empty((ebands.nsppol, ebands.nkpt, ebands.mband, ntemp, nwr))
        aw_meshes = np.empty((ebands.nsppol, ebands.nkpt, ebands.mband, nwr))
        #aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
        #aw_meshes: [max_nbcalc, nkcalc, nsppol] array with energy mesh in eV
        from abipy.tools.numtools import gaussian, lorentzian
        from scipy.integrate import cumtrapz
        for spin in ebands.spins:
            for ik, kpt in enumerate(ebands.kpoints):
                for band in range(ebands.nband_sk[spin, ik]):
                    e0 = ebands.eigens[spin, ik, band]
                    emin = e0 - wr_step * (nwr // 2)
                    emax = e0 + wr_step * (nwr // 2)
                    emesh = np.linspace(emin, emax, num=nwr)
                    aw_meshes[spin, ik, band] = emesh
                    avals = lorentzian(emesh, width=0.2, center=e0, height=None)
                    #if band in (1, 2, 3) and kpt.norm < 0.3:
                    #    avals += 1.1 * lorentzian(emesh, width=0.1, center=e0 - 0.5, height=None)
                    #avals /= cumtrapz(avals, x=emesh)[-1]
                    itemp = 0
                    aw[spin, ik, band, itemp] = avals

        return cls(ebands, aw, aw_meshes, tmesh)

    def __init__(self, ebands, aw, aw_meshes, tmesh):
        """
        Args:
            ebands: |ElectronBands| object
            aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
            aw_meshes: [max_nbcalc, nkcalc, nsppol] array with energy mesh in eV

        .. note::

            Each (spin, k-point, band) can have different frequency meshes
            (usually centered on the initial KS energy).
            The plotter will spline the data.

            The order of the k-points in ebands and aw should be the same.

            The treatment of bands if complicated by the fact that we can have
            different nband(k) whose indices are not necessarily aligned.
            Consider, for example. what happes if we use symsigma or arbitrary gw.
            Use MaskedArray?
        """
        self._ebands = ebands
        self.aw = aw
        self.aw_meshes = aw_meshes
        self.tmesh = tmesh
        #assert

    @property
    def structure(self):
        return self.ebands.structure

    @property
    def ebands(self):
        return self._ebands

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        lines = []
        app = lines.append
        return "\n".join(lines)

    def with_points_along_path(self, frac_bounds=None, knames=None, dist_tol=1e-12):
        p = self.ebands.with_points_along_path()
        # Transfer data using p.ik_new2prev table.
        p.ik_new2prev
        return self.__class__(p.new_ebands, aw, aw_meshes, self.tmesh)

    def interpolate(self):
        new_ebands = self.ebands.interpolate(self, lpratio=5, vertices_names=None, line_density=20,
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

    def get_data_nmtuple(self, itemp, estep, spins=None, k=3, s=0, ext="extrapolate"):
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
                    data[ik] += UnivariateSpline(w, aw, k=k, s=s, ext=ext)(emesh)

        return dict2namedtuple(data=data, emesh=emesh, emin=emin, emax=emax, spins=spins, nkpt=nkpt)

    @add_fig_kwargs
    def plot_itemp(self, itemp=0, spins=None, estep=0.02, ax=None, ylims=None, **kwargs):
        """
        Plot electronic DOS

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                - Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                - None: Don't shift energies, equivalent to ``e0 = 0``.
            spins: Selects the spin to be plotted, None if all spins.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            kwargs: options passed to ``ax.imshow``.

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        a = self.get_data_nmtuple(itemp, estep, spins=spins)

        ax.imshow(a.data.T, origin="lower", extent=[0, a.nkpt, a.emin, a.emax],
                  cmap="jet",
                  interpolation="bilinear",
                  #interpolation="spline36",
                  #interpolation="bicubic",
                  #vmin=0, vmax=np.abs(data).max()
                  )
        self.ebands.decorate_ax(ax)
        #ax.set_zlabel(r"$A(\omega)$")
        #self.ebands.plot(ax=ax, e0=0, color="r", lw=1)

        #ax.imshow(data, cmap=None, norm=None, aspect=None, interpolation=None, alpha=None, vmin=None, vmax=None,
        #       origin=None, extent=None, shape=None, filternorm=1, filterrad=4.0, imlim=None, resample=None,
        #       url=None, hold=None, data=None, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_3dlines(self, itemp=0, estep=0.02, spins=None, band_list=None,
                     k=3, s=0, ext="extrapolate", ax=None, **kwargs):
        ax, fig, plt = get_ax3d_fig_plt(ax=ax)

        xs, emin, emax = self.get_emesh_eminmax(estep)
        nene = len(xs)
        nkpt = self.ebands.nkpt
        cmap = plt.get_cmap("jet")
        # TODO: Reshift everything?
        if band_list is not None:
            band_list = set(band_list)

        # aw: [nwr, ntemp, max_nbcalc, nkcalc, nsppol] array
        spins = range(self.ebands.nsppol) if spins is None else spins
        for spin in spins:
            for ik in range(nkpt):
                ys = np.ones(nene) * ik
                zs = np.zeros(nene)
                for band in range(self.ebands.nband_sk[spin, ik]):
                    if band_list is not None and band not in band_list: continue
                    w = self.aw_meshes[spin, ik, band]
                    aw = self.aw[spin, ik, band, itemp]
                    zs += UnivariateSpline(w, aw, k=k, s=s, ext=ext)(xs)

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
    def plot_surface(self, itemp=0, estep=0.02, spins=None, k=3, s=0, ext="extrapolate", ax=None, **kwargs):
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
                zs[ik] += UnivariateSpline(w, aw, k=k, s=s, ext=ext)(xs)

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


def _main():
    import sys
    plotter = ArpesPlotter.model_from_ebands(sys.argv[1]) #, aw, aw_meshes, tmesh)
    #ebands.plot()
    print(plotter.to_string(verbose=2))
    #plotter.plot_itemp(itemp=0, estep=0.05)
    plotter.plot_3dlines(itemp=0, estep=0.05, band_list=[1, 2, 3])
    #plotter.plot_surface(istep=0, estep=0.05)

if __name__ == "__main__":
    _main()

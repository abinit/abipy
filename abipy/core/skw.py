# coding: utf-8
"""
Shankland-Koelling-Wood Fourier interpolation scheme.
For the theoretical background see :cite:`Euwema1969,Koelling1986,Pickett1988,Madsen2006`.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import itertools
import pickle
import numpy as np
import scipy
import time

from collections import deque, OrderedDict
from monty.collections import dict2namedtuple
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.tools import gaussian
from abipy.core.kpoints import Ktables, Kpath
from abipy.core.symmetries import mati3inv


def n_fermi_dirac(enes, mu, temp):
    """
    Fermi-Dirac distribution.

    Args:
        enes: Energies in eV
        mu: Chemical potential in eV.
        temp: Temperature in Kelvin.

    Return:
        numpy array with occupations in [0, 1].
    """
    if temp > 1e-2:
        return 1. / (np.exp((enes - mu) / (kboltz * temp)) + 1)
    else:
        enes = np.asarray(enes)
        n = np.where(enes <= mu, 1, 0)
        return np.where(enes == mu, 0.5, enes)


def find_degs_sk(enesb, atol):
    """
    Return list of lists with the indices of the degenerated bands.

    Args:
        enesb: Iterable with energies for the different bands.
            Energies are assumed to be ordered.
        atol: Absolute tolerance. Two states are degerated if they differ by less than `atol`.

    Return:
        List of lists. The i-th item contains the indices of the degenerates states
            for the i-th degenerated set.

    :Examples:

    >>> find_degs_sk([1, 1, 2, 3.4, 3.401], atol=0.01)
    [[0, 1], [2], [3, 4]]
    """
    ndeg = 0
    degs = [[0]]
    e0 = enesb[0]
    for ib, ee in enumerate(enesb[1:]):
        ib += 1
        if abs(ee - e0) > atol:
            e0 = ee
            ndeg += 1
            degs.append([ib])
        else:
            degs[ndeg].append(ib)

    return degs


def map_bz2ibz(structure, ibz, ngkpt, has_timrev):
    ngkpt = np.asarray(ngkpt, dtype=dp.int)

    #ibz_grid = [k * ngkpt for k in ibz]
    #ibz_grid = ibz[..., :] * ngkpt
    #for i, k in enumerate(ibz_grid):
    #    ibz_grid[ik] = ibz_grid[ik] % ngkpt

    bzgrid2ibz = np.array(ngkpt.shape, dtype=np.int)
    #for i, gp in enumerate(itertools.product(range(ngkpt[0] + 1), range(ngkpt[1] + 1), range(ngkpt[2] + 1))):
    #      bzgrid2ibz[gp] = i

    symrec_fm = structure.spacegroup.symrec
    time_signs = (1, -1)

    for ik_ibz, kibz in enumerate(ibz):
        gp_ibz = np.rint(kibz * mesh)
        rotated_gps = np.dot(symref_fm, gp_ibz).T
        for tsign in tsigns:
            if tsign == -1: rotated_gps = tsign * rotated_gps
            for rot_gp in rotated_gps:
                gp_bz = rot_gp % ngkpt
                bzgrid2ibz[gp_bz] = ik_ibz

    #bz2ibz = -np.ones(len(bz), dtype=np.int)
    #from abipy.core.kpoints import issamek
    #for ik_bz, kbz in enumerate(bz):
    #    found = False
    #    for ik_ibz, kibz in enumerate(ibz):
    #        if found: break
    #        for symmop in structure.spacegroup:
    #            krot = symmop.rotate_k(kibz)
    #            if issamek(krot, kbz):
    #                bz2ibz[ik_bz] = ik_ibz
    #                found = True
    #                break

    return bz2ibz


class EDOS(object):
    def __init__(self, mesh, values, integral, is_shift, method, step, width):
        self.mesh, self.values, self.integral = mesh, values, integral
        self.is_shift, self.method, self.step, self.width = is_shift, method, step, width


class ElectronInterpolator(object):
    """
    """
    # Tolerances passed to spglib.
    symprec = 1e-5
    angle_tolerance = -1.0

    occtype = "insulator"
    # insulator: fermi level set to homo = nelect//2, occupation factors are either 0 or 1
    #   not suitable for metals, semi-metals, doped semiconductors)
    # fermi-dirac: metallic occupation scheme with physical temperature.
    # gaussian (metallic occupation scheme with gaussian broadening)

    # Disable cache
    use_cache = True

    @classmethod
    def pickle_load(cls, filepath):
        """Loads the object from a pickle file."""
        with open(filepath , "rb") as fh:
            return pickle.load(fh)

    def pickle_dump(self, filepath):
        """Save the status of the object in pickle format."""
        with open(filepath , "wb") as fh:
            pickle.dump(self, fh)

    def get_sampling(self, mesh, is_shift):
        """
        Use spglib to compute the k-points in the IBZ with the corresponding weights
        as well as the k-points in the full BZ.

        Args:
            mesh:
            is_shift:

        Return: named tuple with the following attributes:

            ibz:
            nibz
            weights:
            bz:
            nbz
            grid:
        """
        import spglib as spg
        mesh = np.array(mesh)
        mapping, grid = spg.get_ir_reciprocal_mesh(mesh, self.cell,
            is_shift=is_shift, is_time_reversal=self.has_timrev, symprec=self.symprec)

        uniq, weights = np.unique(mapping, return_counts=True)
        weights = np.asarray(weights, dtype=np.float) / len(grid)
        nkibz = len(uniq)
        ibz = grid[uniq] / mesh
        if self.verbose:
            print("Number of ir-kpoints: %d" % nkibz)

        kshift = 0.0 if is_shift is None else 0.5 * np.asarray(is_shift)
        bz = (grid + kshift) / mesh

        # All k-points and mapping to ir-grid points
        bz2ibz = np.empty(len(bz), dtype=np.int)
        for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
            inds = np.where(uniq == ir_gp_id)
            #print("inds", inds, "inds[0]", inds[0])
            assert len(inds) == 1
            bz2ibz[i] = inds[0]
            #print("%3d ->%3d %s" % (i, ir_gp_id, (gp + [0.5, 0.5, 0.5]) / mesh))
            #print("%3d ->%3d %s" % (i, ir_gp_id, (gp + kshift) / mesh))

        return dict2namedtuple(mesh=mesh, shift=kshift,
                               ibz=ibz, nibz=len(ibz), weights=weights,
                               bz=bz, nbz=len(bz), grid=grid, bz2ibz=bz2ibz)

        #return Ktables(structure, mesh, is_shift, has_timrev)

    #def recalc_fermie(self, kmesh, is_shift=None)
    #    # Compute DOS
    #    edos = _get_cached_edos(kmesh, is_shift)
    #    #if edos is None:
    #    #    edos = self.get_edos(kmesh, is_shift=is_shift, method="gaussian", step=0.1, width=0.2, wmesh=None)
    #    #    self._cache_edos(kmesh, is_shift, edos)
    #    self.fermie = fermie
    #    # Find number of electrons from new chemical potential from nelect.
    #    self.nelect = edos.integral.spline(fermie)
    #    return self.fermie

    @property
    def val_ib(self):
        """The index of the valence band."""
        if self.occtype != "insulator":
            print("Trying to access valence band index with occtype:", self.occtype)
        return int(self.nelect // 2 - 1)

    def set_fermie(self, fermie, kmesh, is_shift=None):
        """
        Change the Fermi level. Use the IDOS computed on the k-grid specifined by
        `kmesh` and `is_shift` to recompute the total number of electrons.

        Args:
            fermie: New Fermi level in eV
            kmesh: Three integers with the number of divisions along the reciprocal primitive axes.
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.

        Return:
            New value of `self.nelect`.
        """
        # Compute DOS
        edos = self._get_cached_edos(kmesh, is_shift)
        #if edos is None:
        #    edos = self.get_edos(kmesh, is_shift=is_shift, method="gaussian", step=0.1, width=0.2, wmesh=None)
        #    self._cache_edos(kmesh, is_shift, edos)

        self.fermie = fermie
        # Find number of electrons from new chemical potential from nelect.
        self.nelect = idos.spline(fermie)
        return self.nelect

    def set_nelect(self, nelect, kmesh, is_shift=None):
        """
        Change the total number of electrons. Use the IDOS computed on the k-grid specifined by
        `kmesh` and `is_shift` to recompute the new Fermi level

        Args:
            nelect: New numbre of electrons.
            kmesh: Three integers with the number of divisions along the reciprocal primitive axes.
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.

        Return:
            New value of `self.fermie`.
        """
        # Compute DOS
        edos = self._get_cached_edos(kmesh, is_shift)
        #if edos is None:
        #    edos = self.get_edos(kmesh, is_shift=is_shift, method="gaussian", step=0.1, width=0.2, wmesh=None)
        #    self._cache_edos(kmesh, is_shift, edos)
        self.nelect = nelect
        # Find new chemical potential from nelect.
        idos
        self.fermie = new_fermie
        return new_fermie

    #def set_charge_per_ucell(self, charge, kmesh, is_shift=None):

    def calc_occfacts(self, eigens, temp):
        """
        Compute occupation factors from the eigenvalues `eigens` and
        the temperature `temp` in K. occfacts in [0, 1].
        """
        if self.occtype == "insulator":
            occfacts = np.ones(eigens.shape)
            occfacts[:, :, self.val_ib + 1:] = 0.0
        else:
            return {
                #"gaussian": n_gaussian,
                "fermi-dirac": n_fermi_dirac,
            }[self.occtype](eigens, self.fermie, temp)

    def get_edos(self, kmesh, is_shift=None, method="gaussian", step=0.1, width=0.2, wmesh=None):
        """
        Compute the electron DOS on a linear mesh.

        Args:
            kmesh: Three integers with the number of divisions along the reciprocal primitive axes.
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            mesh: Frequency mesh to use. If None, the mesh is computed automatically from the eigenvalues.

        Returns:
            (mesh, values, integral)
        """
        k = self.get_sampling(kmesh, is_shift)

        # Interpolate eigenvalues in the IBZ.
        eigens = self._get_cached_eigens(kmesh, is_shift, "ibz")
        if eigens is None:
            eigens = self.interp_kpts(k.ibz).eigens
            self._cache_eigens(kmesh, is_shift, eigens, "ibz")

        # Compute the linear mesh.
        wmesh, step = self._get_wmesh_step(eigens, wmesh, step)
        nw = len(wmesh)
        values = np.zeros((self.nsppol, nw))

        # TODO: Write cython version.
        if method == "gaussian":
            for spin in range(self.nsppol):
                for ik, wtk in enumerate(k.weights):
                    for band in range(self.nband):
                        values[spin] += wtk * gaussian(wmesh, width, center=eigens[spin, ik, band])

            # Compute IDOS
            integral = scipy.integrate.cumtrapz(values, x=wmesh, initial=0.0)

        else:
            raise ValueError("Method %s is not supported" % method)

        return dict2namedtuple(mesh=wmesh, values=values, integral=integral)
        #return ElectronDos(wmesh, values, integral, is_shift, method, step, width)

    def get_jdos_q0(self, kmesh, is_shift=None, method="gaussian", step=0.1, width=0.2, wmesh=None):
        r"""
        Compute the join density of states at q==0

            :math:`\sum_{kbv} f_{vk} (1 - f_{ck}) \delta(\omega - E_{ck} + E_{vk})`

        Args:
            kmesh: Three integers with the number of divisions along the reciprocal primitive axes.
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            method: String defining the method.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            wmesh: Frequency mesh to use. If None, the mesh is computed automatically from the eigenvalues.

        Returns:
        """
        k = self.get_sampling(kmesh, is_shift)

        # Interpolate eigenvalues in the IBZ.
        eigens = self._get_cached_eigens(kmesh, is_shift, "ibz")
        if eigens is None:
            eigens = self.interp_kpts(k.ibz).eigens
            self._cache_eigens(kmesh, is_shift, eigens, "ibz")

        wmesh, step = self._get_w2mesh_step(eigens, wmesh, step)
        nw = len(wmesh)
        values = np.zeros((self.nsppol, nw))

        if self.occtype == "insulator":
            if method == "gaussian":
                for spin in range(self.nsppol):
                    for ik, wtk in enumerate(k.weights):
                        for icb in range(self.val_ib + 1, self.nband):
                            ec = eigens[spin, ik, icb]
                            for icv in range(self.val_ib):
                                ev = eigens[spin, ik, icv]
                                values[spin] += wtk * gaussian(wmesh, width, center=ec-ev)

            else:
                raise ValueError("Method %s is not supported" % method)

        else:
            #occfacts = self.
            if method == "gaussian":
                for spin in range(self.nsppol):
                    for ik, wtk in enumerate(k.weights):
                        for icb in conduction:
                            ec = eigens[spin, ik, icb]
                            fc = 1.0 - occfacts[spin, ik, icb]
                            for icv in valence:
                                ev = eigens[spin, ik, icv]
                                fv = occfacts[spin, ik, icv]
                                fact = wtk * fv * fc
                                values[spin] += fact * gaussian(wmesh, width, center=ec-ev)

            else:
                raise ValueError("Method %s is not supported" % method)

        if self.nsppol == 1: values *= 2.0
        integral = scipy.integrate.cumtrapz(values, x=wmesh, initial=0.0)

        return dict2namedtuple(mesh=wmesh, values=values, integral=integral)
        #return ElectronJointDos(wmesh, values, integral, is_shift, method, step, width)

    #def get_jdos_qpts(self, qpoints, kmesh, is_shift=None, method="gaussian", step=0.1, width=0.2, wmesh=None):
    #    qpoints = np.reshape(qpoints, (-1, 3))
    #    k = self.get_sampling(kmesh, is_shift)

    #    # Interpolate eigenvalues in the full BZ.
    #    eigens_kbz = self._get_cached_eigens(kmesh, is_shift, "bz")
    #    if eigens is None:
    #        eigens_kbz = self.interp_kpts(k.bz).eigens
    #        self._cache_eigens(kmesh, is_shift, eigens_kbz, "bz")
    #    g_skb = gaussian(eigens_kbz, width)

    #    # TODO: One could reduce the sum to IBZ(q) with appropriate weight.
    #    jdos_sqw = np.zeros((self.nsppol, len(qpoints), nw))
    #    for iq, qpt in enumerate(qpoints):
    #        kpq_bz = k.bz + qpt
    #        eigens_kqbz = self.interp_kpts(kpq_bz).eigens
    #        g_skqb = gaussian(eigens_kqbz, width)
    #        #vals = g_skb * g_skqb
    #        #jdos_sqw[:, iq] = vals.sum(axis=(1, 2))

    #    jdos_sqw *= 1. / k.nbz
    #    return jdos_sqw

    def get_nesting_at_e0(self, qpoints, kmesh, e0, width=0.2, is_shift=None):
        """
        Compute the nesting factor with gaussian broadening for an arbitrary list of q-points.

        Args:
            qpoints: List of q-points in reduced coordinates.
            kmesh: Three integers with the number of divisions along the reciprocal primitive axes.
            e0: Energy level in eV.
            width: Standard deviation (eV) of the gaussian.
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.

        Returns:
            numpy array of shape [self.nsppol, len(qpoints)]
        """
        qpoints = np.reshape(qpoints, (-1, 3))
        k = self.get_sampling(kmesh, is_shift)

        # Interpolate eigenvalues in the full BZ.
        eigens_kbz = self._get_cached_eigens(kmesh, is_shift, "bz")
        if eigens is None:
            eigens_kbz = self.interp_kpts(k.bz).eigens
            self._cache_eigens(kmesh, is_shift, eigens_kbz, "bz")

        eigens_kbz -= e0
        g_skb = gaussian(eigens_kbz, width)

        # TODO: One could reduce the sum to IBZ(q) with appropriate weight.
        nest_sq = np.zeros((self.nsppol, len(qpoints)))
        for iq, qpt in enumerate(qpoints):
            kpq_bz = kbz + qpt
            eigens_kqbz = self.interp_kpts(kpq_bz).eigens - e0
            g_skqb = gaussian(eigens_kqbz, width)
            vals = g_skb * g_skqb
            nest_sq[:, iq] = vals.sum(axis=(1, 2))

        nest_sq *= 1. / k.nbz
        return nest_sq

    #def get_unitcell_vals(self, kmesh):
    #    is_shift = None
    #    k = self.get_sampling(kmesh, is_shift=is_shift)

    #    # Interpolate eigenvalues in the IBZ.
    #    eigens_ibz = self._get_cached_eigens(kmesh, is_shift, "ibz")
    #    if eigens_ibz is None:
    #        eigens_ibz = self.interp_kpts(k.ibz).eigens
    #        self._cache_eigens(kmesh, is_shift, eigens_ibz, "ibz")

    #    egrid_sbuc = np.empty((self.nsppol, self.nband) + tuple(kmesh))
    #    egrid_sbuc[...] = np.inf
    #    print("shape eigens_ibz", eigens_ibz.shape)
    #    print("shape egrid_sbuc", egrid_sbuc.shape)
    #    print("nkibz: ", k.nibz, "nbz: ", k.nbz)
    #    for gp, ik_ibz in zip(k.grid, k.bz2ibz):
    #        ucgp = np.where(gp >= 0, gp, gp + k.mesh)
    #        print("gp", gp, "ucgp", ucgp, "ik_ibz", ik_ibz)
    #        egrid_sbuc[:, :, ucgp[0], ucgp[1], ucgp[2]] = eigens_ibz[:, ik_ibz, :]
    #    print(np.any(egrid_sbuc == np.inf))
    #    #print(egrid_sbuc[0, 0, ...])

    #    return egrid_sbuc

    def _get_wmesh_step(self, eigens, wmesh, step):
        if wmesh is not None:
            return wesh, wmesh[1] - wmesh[0]

        # Compute the linear mesh.
        epad = 1.0
        e_min = eigens.min() - epad
        e_max = eigens.max() + epad
        nw = int(1 + (e_max - e_min) / step)
        wmesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)

        return wmesh, step

    def _get_w2mesh_step(self, eigens, wmesh, step):
        if wmesh is not None:
            return wesh, wmesh[1] - wmesh[0]

        # Compute the linear mesh.
        cmin, cmax = +np.inf, -np.inf
        vmin, vmax = +np.inf, -np.inf
        for spin in range(self.nsppol):
            for c in range(self.val_ib + 1, self.nband):
                cmin = min(cmin, eigens[spin,:,c].min())
                cmax = max(cmax, eigens[spin,:,c].max())
            for v in range(self.val_ib + 1):
                vmin = min(vmin, eigens[spin,:,v].min())
                vmax = max(vmax, eigens[spin,:,v].max())
        #cmin = eigens.max()
        #cmax = cmin
        #vmax = eigens.min()
        #vmin = vmax

        e_min = cmin - vmax
        e_min -= 0.1 * abs(e_min)
        e_max = cmax - vmin
        e_max += 0.1 * abs(e_max)

        nw = int(1 + (e_max - e_min) / step)
        wmesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)

        return wmesh, step

    def _get_cached_eigens(self, kmesh, is_shift, kzone):
        """
        Return a copy of the interpolated eigenvalues associated to (kmesh, is_shift, kzone).
        Return None if eigens are not available.
        """
        if not self.use_cache: return None
        if not hasattr(self, "_cached_eigens"): self._cached_eigens = OrderedDict()
        kmesh = tuple(kmesh)
        if is_shift is not None: is_shift = tuple(is_shift)
        arr = self._cached_eigens.get((kmesh, is_shift, kzone))
        if arr is not None: arr = arr.copy()
        return arr

    def _cache_eigens(self, kmesh, is_shift, eigens, kzone):
        """
        Save interpolated eigenvalues associated to (kmesh, is_shift, kzone).
        """
        if not self.use_cache: return None
        if not hasattr(self, "_cached_eigens"): self._cached_eigens = OrderedDict()
        kmesh = tuple(kmesh)
        if is_shift is not None: is_shift = tuple(is_shift)
        self._cached_eigens[(kmesh, is_shift, kzone)] = eigens.copy()

    def _get_cached_edos(self, kmesh, is_shift):
        """
        Return the electron DOS obtained from the interpolated eigenvalues associated
        to (kmesh, is_shift). None if DOS is not available.
        """
        if not self.use_cache: return None
        if not hasattr(self, "_cached_edos"): self._cached_edos = OrderedDict()
        kmesh = tuple(kmesh)
        if is_shift is not None: is_shift = tuple(is_shift)
        return self._cached_edos.get((kmesh, is_shift))

    def _cache_edos(self, kmesh, is_shift, edos):
        """
        Save the electron DOS obtained from the interpolated eigenvalues associated to (kmesh, is_shift).
        """
        if not self.use_cache: return None
        if not hasattr(self, "_cached_edos"): self._cached_edos = OrderedDict()
        kmesh = tuple(kmesh)
        if is_shift is not None: is_shift = tuple(is_shift)
        self._cached_edos[(kmesh, is_shift)] = edos # .copy()

    @add_fig_kwargs
    def plot_dos_vs_kmeshes(self, kmeshes, is_shift=None, method="gaussian", step=0.1, width=0.2,
                            fontsize=12, ax=None, **kwargs):
        """
        Plot (interpolated) DOSes computed with different meshes.

        Args:
            kmeshes: List of kmeshes. Each item is given by three integers with the number of
                divisions along the reciprocal primitive axes.
            is_shift: three integers (spglib_ API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        kmeshes = np.reshape(np.asarray(kmeshes, dtype=np.int), (-1, 3))
        for kmesh in kmeshes:
            edos = self.get_edos(kmesh, is_shift=is_shift, method=method, step=step, width=width)
            for spin in range(self.nsppol):
                spin_sign = +1 if spin == 0 else -1
                ax.plot(edos.mesh, edos.values[spin] * spin_sign, label=str(kmesh) if spin == 0 else None)

        ax.grid(True)
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel('DOS (states/eV)')
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_jdosq0_vs_kmeshes(self, kmeshes, is_shift=None, method="gaussian", step=0.1, width=0.2,
                               ax=None, fontsize=12, **kwargs):
        """
        Plot (interpolated) Joint DOSes at q=0 computed with different meshes.

        Args:
            kmeshes: List of kmeshes. Each item is given by three integers with the number of
                divisions along the reciprocal primitive axes.
            is_shift: three integers (spglib API). If None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            ax: matplotlib `Axes` or None if a new figure should be created.
            fontsize: Legend and title fontsize.

        Returns:
            matplotlib figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        kmeshes = np.reshape(np.asarray(kmeshes, dtype=np.int), (-1, 3))
        for kmesh in kmeshes:
            jdos = self.get_jdos_q0(kmesh, is_shift=is_shift, method=method, step=step, width=width)
            for spin in range(self.nsppol):
               spin_sign = +1 if spin == 0 else -1
               ax.plot(jdos.mesh, jdos.values[spin] * spin_sign, label=str(kmesh) if spin == 0 else None)

        ax.grid(True)
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel('JDOS (states/eV)')
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_nesting_vs_widths(self, widths, kmesh, e0=None, qvertices_names=None,
                               line_density=20, is_shift=None, ax=None, fontsize=12, **kwargs):
        """
        Plot (interpolated) nesting factor computed with different broadening.

        Args:
            widths: List of standard deviations (eV) of the gaussian.
            kmeshes: List of kmeshes. Each item is given by three integers with the number of
                divisions along the reciprocal primitive axes.
            e0: Energy level in eV.
            qvertices_names
            line_density:
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            fontsize: Legend and title fontsize.

        Returns:
            matplotlib figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        qpoints = self._get_kpts_kticks_klabels(ax, qvertices_names, line_density)

        e0 = self.interpolated_fermie if e0 is None else e0
        for width in np.asarray(widths):
            nest_sq = self.get_nesting_at_e0(qpoints, kmesh, e0, width=width, is_shift=is_shift)
            for spin in range(self.nsppol):
                spin_sign = +1 if spin == 0 else -1
                ax.plot(nest_sq[spin] * spin_sign, label=str(kmesh) if spin == 0 else None)

        ax.grid(True)
        ax.set_ylabel('Nesting factor')
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_nesting_vs_kmeshes(self, width, kmeshes, e0=None, qvertices_names=None, line_density=20,
                                is_shift=None, ax=None, fontsize=12, **kwargs):
        """
        Plot (interpolated) nesting factor computed with different k-meshes.

        Args:
            width: Gaussian broadening (eV).
            kmeshes: List of kmeshes. Each item is given by three integers with the number of
                divisions along the reciprocal primitive axes.
            e0: Energy level in eV.
            qvertices_names
            line_density:
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            fontsize: legend and title fontsize.

        Returns:
            matplotlib figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        qpoints = self._get_kpts_kticks_klabels(ax, qvertices_names, line_density)

        kmeshes = np.reshape(np.asarray(kmeshes, dtype=np.int), (-1, 3))
        e0 = self.interpolated_fermie if e0 is None else e0
        for kmesh in kmeshes:
            nest_sq = self.get_nesting_at_e0(qpoints, kmesh, e0, width=width, is_shift=is_shift)
            for spin in range(self.nsppol):
                spin_sign = +1 if spin == 0 else -1
                ax.plot(nest_sq[spin] * spin_sign, label=str(kmesh) if spin == 0 else None)

        ax.grid(True)
        ax.set_ylabel('Nesting factor')
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_group_velocites(self, vertices_names=None, line_density=20, ax=None, **kwargs):
        """
        Plot (interpolated) group velocities computed along an arbitrary k-path.

        Args:
            vertices_names
            line_density:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            matplotlib figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        kfrac_coords = self._get_kpts_kticks_klabels(ax, vertices_names, line_density)

        #v_skb = self.interp_kpts(kfrac_coords, dk1=v_skb)
        for spin in range(self.nsppol):
            for band in range(self.nband):
                # plot |v|
                v_skb[spin, :, band]
                ax.plot(vals, color="k" if spin == 0 else "r")

        ax.grid(True)
        ax.set_ylabel('Group Velocities')

        return fig

    def _get_kpts_kticks_klabels(self, ax, vertices_names, line_density):
        if vertices_names is None:
            vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]

        kpath = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)

        #self._axset_ticks_labels(ax, xticks, xlabels)
        return kpath.frac_coords, list(range(len(kpath))), kpath.names


class SkwInterpolator(ElectronInterpolator):
    """
    This object implements the Shankland-Koelling-Wood Fourier interpolation scheme.
    It can be used to interpolate functions in k-space with the periodicity of the
    reciprocal lattice and satisfying F(k) = F(Sk) for each rotation S
    belonging to the point group of the crystal. For readability reason,
    the names of the variables are chosen assuming we are interpolating electronic eigenvalues
    but the same object can be used to interpolate other quantities. Just set the first dimension to 1.
    """

    def __init__(self, lpratio, kpts, eigens, fermie, nelect, cell, symrel, has_timrev,
                 filter_params=None, verbose=1):
        """
        Args:
            lpratio: Ratio between the number of star-functions and the number of ab-initio k-points.
                5-10 should be OK in many systems, larger values may be required for accurate derivatives.
            kpts: numpy array with the [nkpt, 3] ab-initio k-points in reduced coordinates.
            eigens: numpy array with the ab-initio energies. shape [nsppol, nkpt, nband].
            fermie: Fermi energy in eV.
            nelect: Number of electrons in the unit cell
            cell: (lattice, positions, numbers)
                lattice: numpy array with direct lattice vectors along the rows.
                positions: atomic positions in reduced coordinates.
                numbers: Atomic number for each atom in positions.
            symrel: [nsym, 3, 3] numpy array with the (ferromagnetic) symmetry operations of the direct lattice
                in reduced coordinates. anti-ferromagnetic symmetries (if any) should be removed by the caller.
            has_timrev: True is time-reversal can be used.
            filter_params: List with parameters used to filter high-frequency components (Eq 9 of PhysRevB.61.1639)
                First item gives rcut, second item sigma. Ignored if None.
            verbose: Verbosity level.
        """
        self.verbose = verbose
        self.cell = cell
        lattice = self.cell[0]
        self.original_fermie = fermie
        self.interpolated_fermie = self.original_fermie
        self.nelect = nelect
        self.has_timrev = has_timrev

        # iscomplexobj is used to handle lifetimes.
        eigens = np.atleast_3d(eigens)
        self.iscomplexobj = np.iscomplexobj(eigens)
        self.nsppol, self.nkpt, self.nband = eigens.shape

        if len(kpts) != self.nkpt:
            raise ValueError("Second dimension of eigens should be %d but got array of shape: %s" %
                (len(kpts), eigens.shape))
        if self.nkpt == 1:
            raise ValueError("Interpolation algorithm requires nkpt > 1")

        rprimd = np.asarray(lattice).T
        #print("rprimd", rprimd)
        self.rmet = np.matmul(rprimd.T, rprimd)

        # Find point group operations.
        symrel = np.reshape(symrel, (-1, 3, 3))
        self.ptg_symrel, self.ptg_symrec, has_inversion = extract_point_group(symrel, has_timrev)
        self.ptg_nsym = len(self.ptg_symrel)
        if self.verbose:
            print("Found", self.ptg_nsym, "symmetries in point group")

        # Find nrwant star points.
        self.lpratio = lpratio = int(lpratio)
        if lpratio <= 1:
            raise ValueError("lpratio must be > 1 but got %s" % lpratio)

        nrwant = lpratio * self.nkpt
        fact = 1/2 if has_inversion else 1
        rmax = int((1.0 + (lpratio * self.nkpt * self.ptg_nsym * fact) / 2.0) ** (1/3.)) * np.ones(3, dtype=np.int)
        #rmax = int((1.0 + (lpratio * self.nkpt) / 2.0) ** (1/3.)) * np.ones(3, dtype=np.int)

        while True:
            self.rpts, r2vals, ok = self._find_rstar_gen(nrwant, rmax)
            self.nr = len(self.rpts)
            if ok:
                break
            else:
                print("rmax: ", rmax," was not large enough to find", nrwant, "R-star points.")
                rmax *= 2
                print("Will try again with enlarged rmax:", rmax)

        print("Will use:", self.nr, "star-functions. nstars/nk:", self.nr / self.nkpt)

        # Compute (inverse) roughness function.
        c1, c2 = 0.25, 0.25
        r2min = r2vals[1]
        inv_rhor = 1.0 / ((1.0 - c1 * r2vals / r2min) ** 2 + c2 * (r2vals / r2min) ** 3)

        # Construct star functions for the ab-initio k-points.
        nsppol, nband, nkpt, nr = self.nsppol, self.nband, self.nkpt, self.nr
        self.skr = np.empty((nkpt, nr), dtype=np.complex)
        for ik, kpt in enumerate(kpts):
            self.skr[ik] = self.get_stark(kpt)

        # Build H(k,k') matrix (Hermitian)
        hmat = np.empty((nkpt-1, nkpt-1), dtype=np.complex)
        for jk in range(nkpt-1):
            v_jkr = self.skr[jk, 1:] - self.skr[nkpt-1, 1:]
            #for ik in range(jk + 1):
            for ik in range(nkpt-1):
                v_ikr = inv_rhor[1:] * (self.skr[ik, 1:] - self.skr[nkpt-1, 1:])
                hmat[ik, jk] = np.vdot(v_jkr, v_ikr)
                if ik == jk: hmat[ik, jk] = hmat[ik, jk].real

        # Solving system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)..."
        de_kbs = np.empty((nkpt-1, nband, nsppol), dtype=np.complex)
        for spin in range(nsppol):
            for ib in range(nband):
                de_kbs[:, ib, spin] = eigens[spin, 0:nkpt-1, ib] - eigens[spin, nkpt-1, ib]

        # Solve all bands and spins at once
        # FIXME: Portability problem with scipy 0.19 in which linalg.solve wraps the expert drivers
        # http://scipy.github.io/devdocs/release.0.19.0.html#foreign-function-interface-improvements
        if scipy.__version__ == "0.19.0":
            import warnings
            warnings.warn("linalg.solve in scipy 0.19.0 gives weird results. Use at your own risk!!!")

        try:
            lmb_kbs = scipy.linalg.solve(hmat, np.reshape(de_kbs, (-1, nband * nsppol)))
            #lmb_kbs = scipy.linalg.solve(hmat, np.reshape(de_kbs, (-1, nband * nsppol)),
                    #sym_pos=True, lower=False, overwrite_a=True, overwrite_b=True, check_finite=False)
                    #sym_pos=False, lower=False, overwrite_a=False, overwrite_b=False, check_finite=True)

        except scipy.linalg.LinAlgError as exc:
            print("Cannot solve system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)")
            print("This usually happens when there are symmetrical k-points passed to the interpolator.")
            raise exc

        lmb_kbs = np.reshape(lmb_kbs, (-1, nband, nsppol))

        # Compute coefficients.
        self.coefs = np.empty((nsppol, nband, nr), dtype=np.complex)
        for spin in range(nsppol):
            for ib in range(nband):
                for ir in range(1, nr):
                    self.coefs[spin, ib, ir] = inv_rhor[ir] \
                        * np.vdot(self.skr[:nkpt-1, ir] - self.skr[nkpt-1, ir], lmb_kbs[:nkpt-1, ib, spin])

                self.coefs[spin, ib, 0] = eigens[spin, nkpt-1, ib] \
                    - np.dot(self.coefs[spin, ib, 1:nr], self.skr[nkpt-1, 1:nr])

        # Filter high-frequency.
        self.rcut, self.rsigma = None, None
        if filter_params is not None:
            self.rcut = filter_params[0] * np.sqrt(r2vals[-1])
            self.rsigma = rsigma = filter_params[1]
            if self.verbose:
                print("Applying filter (Eq 9 of PhysRevB.61.1639) with rcut:", rcut, ", rsigma", self.rsigma)
            from scipy.special import erfc
            for ir in range(1, nr):
                self.coefs[:, :, ir] *= 0.5 * erfc((np.sqrt(r2vals[ir]) - self.rcut) / self.rsigma)

        # Prepare workspace arrays for star functions.
        self.cached_kpt = np.ones(3) * np.inf
        self.cached_kpt_dk1 = np.ones(3) * np.inf
        self.cached_kpt_dk2 = np.ones(3) * np.inf

        # Compare ab-initio data with interpolated results.
        mae = 0.0
        for spin in range(nsppol):
            for ik, kpt in enumerate(kpts):
                skw_eb = self.eval_sk(spin, kpt)
                mae += np.abs(eigens[spin, ik] - skw_eb).sum()
                if self.verbose >= 10:
                    # print interpolated eigenvales
                    for band in range(self.nband):
                        e0 = eigens[spin, ik, band]
                        eskw = skw_eb[band]
                        print("spin", spin, "band", band, "ikpt", ik, "e0", e0, "eskw", eskw, "diff", e0 - eskw)

        mae *= 1e3 / (nsppol * nkpt * nband)
        print("FIT vs input data: Mean Absolute Error=", mae, " [meV]")
        if mae > 10.0:
            # Issue warning if error too large.
            print("Large error in SKW interpolation!")
            print("MAE:", mae, "[meV]")

        if np.isnan(mae) or np.isinf(mae):
            raise RuntimeError("Interpolation went bananas! mae = %s" % mae)

        self.mae = mae

    def __str__(self):
        return self.to_string()

    def to_string(self, **kwargs):
        """String representation."""
        lines = []
        app = lines.append
        app("nsppol: %s, nband: %s" % (self.nsppol, self.nband))
        app("Number of operations in point-group: %s with time-reversal: %s" % (self.ptg_nsym, self.has_timrev))
        app("Number of ab-initio k-points: %s" % self.nkpt)
        app("Number of star functions: %s [nstars/nk = %s]" % (self.nr, (self.nr / self.nkpt)))
        if self.rcut is not None:
            app("Fourier filter (Eq 9 of PhysRevB.61.1639) with rcut: %s, rsigma: %s" % (self.rcut, self.rsigma))
        app("Comparison between ab-initio data and fit gave Mean Absolute Error: %s [meV]" % self.mae)

        return "\n".join(lines)

    def interp_kpts(self, kfrac_coords, dk1=False, dk2=False):
        """
        Interpolate energies on an arbitrary set of k-points. Optionally, compute
        gradients and Hessian matrices.

        Args:
            kfrac_coords: K-points in reduced coordinates.
            dk1 (bool): True if gradient is wanted.
            dk2 (bool): True to compute 2nd order derivatives.

        Return:
            namedtuple with:
            interpolated energies in eigens[nsppol, len(kfrac_coords), nband]
            gradient in dedk[self.nsppol, len(kfrac_coords), self.nband, 3))
            hessian in dedk2[self.nsppol, len(kfrac_coords), self.nband, 3, 3))

            gradient and hessian are set to None if not computed.
        """
        start = time.time()

        kfrac_coords = np.reshape(kfrac_coords, (-1, 3))
        new_nkpt = len(kfrac_coords)
        new_eigens = np.empty((self.nsppol, new_nkpt, self.nband))

        dedk = None if not dk1 else np.empty((self.nsppol, new_nkpt, self.nband, 3))
        dedk2 = None if not dk2 else np.empty((self.nsppol, new_nkpt, self.nband, 3, 3))

        der1, der2 = None, None
        for spin in range(self.nsppol):
            for ik, newk in enumerate(kfrac_coords):
                if dk1: der1 = dedk[spin, ik]
                if dk2: der2 = dedk2[spin, ik]
                new_eigens[spin, ik] = self.eval_sk(spin, newk, der1=der1, der2=der2)

        if self.verbose:
            print("Interpolation completed", time.time() - start)
        return dict2namedtuple(eigens=new_eigens, dedk=dedk, dedk2=dedk2)

    def interp_kpts_and_enforce_degs(self, kfrac_coords, ref_eigens, atol=1e-4):
        """
        Interpolate energies on an arbitrary set of k-points. Use `ref_eigens`
        to detect degeneracies and average the interpolated values in the degenerate subspace.
        """
        kfrac_coords = np.reshape(kfrac_coords, (-1, 3))
        new_nkpt = len(kfrac_coords)
        ref_eigens = np.reshape(ref_eigens, (self.nsppol, new_nkpt, self.nband))

        # Interpolate eigenvales.
        new_eigens = self.interp_kpts(kfrac_coords).eigens

        # Average interpolated values over degenerates bands.
        for spin in range(self.nsppol):
            for ik in range(new_nkpt):
                for dgbs in find_degs_sk(ref_eigens[spin, ik], atol):
                    if len(dgbs) == 1: continue
                    new_eigens[spin, ik, dgbs] = new_eigens[spin, ik, dgbs].sum() / len(dgbs)

        return dict2namedtuple(eigens=new_eigens, dedk=None, dedk2=None)

    def eval_sk(self, spin, kpt, der1=None, der2=None):
        """
        Interpolate eigenvalues for all bands at a given (spin, k-point).
        Optionally compute gradients and Hessian matrices.

        Args:
            spin: Spin index.
            kpt: K-point in reduced coordinates.
            der1: If not None, ouput gradient is stored in der1[nband, 3].
            der2: If not None, output Hessian is der2[nband, 3, 3].

        Return:
            oeigs[nband]
        """
        # [NB, NR] x [NR]
        oeigs = np.matmul(self.coefs[spin], self.get_stark(kpt))
        if not self.iscomplexobj: oeigs = oeigs.real

        if der1 is not None:
            skr_dk1 = self.get_stark_dk1(kpt)
            for ii in range(3):
                value = np.matmul(self.coefs[spin, :, :], skr_dk1[ii])
                if not self.iscomplexobj: value = value.real
                der1[:, ii] = value

        if der2 is not None:
            skr_dk2 = self.get_stark_dk2(kpt)
            for jj in range(3):
                for ii in range(jj + 1):
                    value = np.matmul(self.coefs[spin, :, :], skr_dk2[ii,jj])
                    if not self.iscomplexobj: value = value.real
                    der2[:, ii, jj] = value
                    if ii != jj: der2[jj, ii] = der2[ii, jj]

        return oeigs

    #def eval_skb(self, spin, kpt, band, der1=None, der2=None):
    #    """
    #    Interpolate eigenvalues for a given (spin, k-point, band).

    #    Args:

    #    Return:
    #    """
    #    # Compute star function for this k-point (if not already in memory)
    #    if np.any(kpt != self.cached_kpt):
    #        self.cached_skr, self.cached_kpt = self.get_stark(kpt), kpt

    #    # Do not take complex conjugate.
    #    oeig = np.dot(self.coefs[spin, band], self.cached_skr)
    #    if not self.iscomplexobj: oeigs = oeigs.real
    #    if der1 is None and der2 is None: return oeig

    #    # TODO: Test Derivatives
    #    if der1 is not None:
    #        # Compute first-order derivatives.
    #        if np.any(kpt != self.cached_kpt_dk1):
    #            self.cached_skr_dk1, self.cached_kpt_dk1 = self.get_stark_dk1(kpt), kpt

    #        for ii in range(3):
    #            value = np.dot(self.coefs[spin, band], self.cached_skr_dk1[ii])
    #            if not self.iscomplexobj: value = value.real
    #            der1[ii] = value

    #        if der2 is None: return oeig, der1

    #    if der2 is not None:
    #        # Compute second-order derivatives.
    #        if np.any(kpt != self.cached_kpt_dk2):
    #            self.cached_skr_dk2, self.cached_kpt_dk2 = self.get_stark_dk2(kpt), kpt

    #        der2 = zero
    #        for jj in range(3):
    #            for ii in range(jj + 1):
    #                value = np.dot(self.coefs[spin, band], self.cached_skr_dk2[ii,jj])
    #                if not self.iscomplexobj: value = value.real
    #                der2[ii, jj] = value
    #                if ii != jj: der2[jj, ii] = der2[ii, jj]

    #    return oeig, der1, der2

    def get_stark(self, kpt):
        """
        Return the star function for k-point `kpt`.

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            complex array of shape [self.nr]
        """
        two_pi = 2.0 * np.pi
        skr = np.zeros(self.nr, dtype=np.complex)
        _np_exp = np.exp
        for omat in self.ptg_symrel:
            sk = two_pi * np.matmul(omat.T, kpt)
            skr += _np_exp(1.j * np.matmul(self.rpts, sk))
        skr /= self.ptg_nsym

        return skr

    def get_stark_dk1(self, kpt):
        """
        Compute the 1st-order derivative of the star function wrt k

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            complex array [3, self.nr]  with the derivative of the
            star function wrt k in reduced coordinates.
        """
        srk_dk1 = np.zeros((3, self.nr), dtype=np.complex)
        two_pi = 2.0 * np.pi
        rpts_t = self.rpts.T

        for omat in self.ptg_symrel:
            sk = two_pi * np.matmul(omat.T, kpt)
            exp_skr = np.exp(1.j * np.matmul(self.rpts, sk))
            for ir, rr in enumerate(self.rpts):
                srk_dk1[:, ir] += exp_skr[ir] * np.matmul(omat, rr)
            #for ir, or in enumerate(np.matmul(omat, rpts_t).T):
            #    srk_dk1[:, ir] += exp_skr[ir] * or

        srk_dk1 *= 1.j / self.ptg_nsym
        return srk_dk1

    def get_stark_dk2(self, kpt):
        """
        Compute the 2nd-order derivatives of the star function wrt k.

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            Complex numpy array of shape [3, 3, self.nr] with the 2nd-order derivatives
            of the star function wrt k in reduced coordinates.
        """
        srk_dk2 = np.zeros((3, 3, self.nr), dtype=np.complex)
        raise NotImplementedError()
        #work = zero
        #do isym=1,self.ptg_nsym
        #   sk = two_pi * matmul(transpose(self%ptg_symrel(:,:,isym)), kpt)
        #   do ir=1,self%nr
        #     sr = matmul(self%ptg_symrel(:,:,isym), self%rpts(:,ir))
        #     eiskr = exp(j_dpc * dot_product(sk, self%rpts(:,ir)))
        #     do jj=1,3
        #        do ii=1,jj
        #            work(ii,jj,ir) = work(ii,jj,ir) + eiskr * sr(ii) * sr(jj)

        #work = - work / self.ptg_nsym

        #do jj=1,3
        #   do ii=1,jj
        #       srk_dk2(:, ii, jj) = work(ii, jj, :)
        #       if (ii /= jj) srk_dk2(:,jj,ii) = work(:,ii,jj)

        return srk_dk2

    #def find_stationary_points(self, kmesh, bstart=None, bstop=None, is_shift=None)
    #    k = self.get_sampling(kmesh, is_shift)
    #    if bstart is None: bstart = self.nelect // 2 - 1
    #    if bstop is None: bstop = self.nelect // 2
    #    nb = bstop - bstart + 1
    #    results = []
    #    for ik_ibz, kpt in enumerate(k.ibz):
    #        vk_b = self.eval_dk1(kpt, bstart, bstop)
    #        bands = []
    #        for ib, v in enumerate(vk_b):
    #            if v < atol: bands.append(ib + bstart)
    #        if bands:
    #            #results.append()
    #            for band in bands:
    #                d2k_b = self.eval_dk2(kpt, band)

    #    return results

    def _find_rstar_gen(self, nrwant, rmax):
        """
        Find all lattice points generating the stars inside the supercell defined by `rmax`

        Args:
            nrwant: Number of star-functions required.
            rmax: numpy array with the maximum number of cells along the 3 reduced directions.

        Returns:
            tuple: (rpts, r2vals, ok)
        """
        msize = (2 * rmax + 1).prod()
        rtmp = np.empty((msize, 3), dtype=np.int)
        r2tmp = np.empty(msize)
        if self.verbose: print("rmax", rmax, "msize:", msize)

        start = time.time()
        for cnt, l in enumerate(itertools.product(range(-rmax[0], rmax[0] + 1),
            range(-rmax[1], rmax[1] + 1), range(-rmax[2], rmax[2] + 1))):
              rtmp[cnt] = l
              r2tmp[cnt] = np.dot(l, np.matmul(self.rmet, l))
        if self.verbose: print("gen points", time.time() - start)

        start = time.time()
        # Sort r2tmp and rtmp
        iperm = np.argsort(r2tmp)
        r2tmp = r2tmp[iperm]
        rtmp = rtmp[iperm]
        #return rtmp, r2tmp, True

        # Find shells
        r2sh = np.empty(msize, dtype=np.int)    # Correspondence between R and shell index.
        shlim = np.empty(msize, dtype=np.int)   # For each shell, the index of the initial G-vector.
        nsh = 1
        r2sh[0] = 0
        shlim[0] = 0
        r2_prev = 0.0
        for ir in range(1, msize):
            if (abs(r2tmp[ir] - r2_prev) > r2tmp[ir] * 1e-8):
                r2_prev = r2tmp[ir]
                shlim[nsh] = ir
                nsh += 1
            r2sh[ir] = nsh
        shlim[nsh] = msize
        if self.verbose:
            print("nshells", nsh)
            print("shells", time.time() - start)

        # Find R-points generating the stars.
        # This is the most CPU critical part. I think it's difficult to do better than this without cython
        start = time.time()
        rgen = deque()

        for ish in range(nsh):
            ss, ee = shlim[ish], shlim[ish + 1]
            if ss + 1 == ee:
                rgen.append(rtmp[ss])
                continue

            if True:
                # This is faster
                rs = set([tuple(rtmp[ss])])
                for ir in range(ss + 1, ee):
                    rvec = rtmp[ir]
                    if all(tuple(np.matmul(rot, rvec)) not in rs for rot in self.ptg_symrel):
                        rs.add(tuple(rvec))
            else:
                rs = deque([rtmp[ss]])
                for ir in range(ss + 1, ee):
                    for rot in self.ptg_symrel:
                        rot_r = np.matmul(rot, rtmp[ir])
                        if any(np.all(rot_r == x) for x in rs): break
                    else:
                        # We have new R-point.
                        #print("Adding new point")
                        rs.append(rtmp[ir])

            #print(len(rs), rs)
            rgen.extend(rs)
        #print(rgen)
        if self.verbose: print("stars", time.time() - start)

        start = time.time()
        rgen = np.array(rgen, dtype=np.int)
        nstars = len(rgen)

        # Store rpts and compute ||R||**2.
        ok = nstars >= nrwant
        nr = min(nstars, nrwant)
        rpts = rgen[:nr].copy()
        r2vals = np.empty(nr)
        for ir in range(nr):
            r2vals[ir] = np.dot(rpts[ir], np.matmul(self.rmet, rpts[ir]))

        if self.verbose:
            print("r2max ", rpts[nr-1])
            print("end ", time.time() - start)
            if self.verbose > 10:
                print("nstars:", nstars)
                for r, r2 in zip(rpts, r2vals):
                    print(r, r2)

        return rpts, r2vals, ok


def extract_point_group(symrel, has_timrev):
    """
    Extract the point group rotations from the spacegroup. Add time-reversal
    if spatial inversion is not present and `has_timrev`.

    Return:
        (ptg_symrel, ptg_symrec) with rotations in real- and reciprocal-space.
    """
    nsym = len(symrel)
    tmp_nsym = 1
    work_symrel = np.empty((2*nsym, 3, 3), dtype=np.int)
    work_symrel[0] = symrel[0]

    for isym in range(1, nsym):
        found = any(np.all(symrel[isym] == work_symrel[search]) for search in range(tmp_nsym))
        if not found:
            work_symrel[tmp_nsym] = symrel[isym]
            tmp_nsym += 1

    inversion_3d = -np.eye(3, dtype=np.int)
    has_inversion = any(np.all(w == inversion_3d) for w in work_symrel[:tmp_nsym])

    # Now we know the symmetries of the point group.
    ptg_nsym = 2 * tmp_nsym if not has_inversion and has_timrev else tmp_nsym
    ptg_symrel = np.empty((ptg_nsym, 3, 3), dtype=np.int)
    ptg_symrec = np.empty((ptg_nsym, 3, 3), dtype=np.int)

    ptg_symrel[:tmp_nsym] = work_symrel[:tmp_nsym]
    for isym in range(tmp_nsym):
        ptg_symrec[isym] = mati3inv(ptg_symrel[isym])

    if not has_inversion and has_timrev:
        # Add inversion.
        ptg_symrel[tmp_nsym:] = -work_symrel[:tmp_nsym]
        for isym in range(tmp_nsym, ptg_nsym):
            ptg_symrec[isym] = mati3inv(ptg_symrel[isym])

    return ptg_symrel, ptg_symrec, has_inversion

# coding: utf-8
"""
Shankland-Koelling-Wood Fourier interpolation scheme.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import itertools
import pickle
import numpy as np
import scipy

from collections import deque, OrderedDict
from abipy.tools import gaussian

_np_matmul = np.matmul
_np_dot = np.dot
_np_vdot = np.vdot


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


class ElectronInterpolator(object):
    """
    """
    # Tolerances passed to spglib.
    symprec = 1e-5
    angle_tolerance = -1.0

    # insulator: fermi level set to homo = nelect//2, occupation factors are either 0 or 1
    #   not suitable for metals, semi-metals, doped semiconductors)
    # fermi-dirac: metallic occupation scheme with physical temperature.
    # gaussian (metallic occupation scheme with gaussian broadening)
    occtype = "insulator"

    # Disable cache
    no_cache = False

    @classmethod
    def pickle_load(cls, filepath):
        """Loads the object from a pickle file."""
        with open(filepath , "rb") as fh:
            return pickle.load(fh)

    def pickle_dump(self, filepath):
        """Save the status of the object in pickle format."""
        with open(filepath , "wb") as fh:
            pickle.dump(self, fh)

    def get_ibz_weighs_bz(self, mesh, is_shift):
        import spglib as spg
        mesh = np.asarray(mesh)
        mapping, grid = spg.get_ir_reciprocal_mesh(mesh, self.cell,
            is_shift=is_shift, is_time_reversal=self.is_time_reversal, symprec=self.symprec)

        uniq, weights = np.unique(mapping, return_counts=True)
        print("uniq", uniq, type(uniq))
        print("weights", weights, type(weights))
        print("grid", grid.shape, type(grid))
        weights = np.array(weights, dtype=np.float) / len(grid)
        print(weights.sum())
        nkibz = len(uniq)
        kibz = grid[uniq] / mesh
        print("Number of ir-kpoints: %d" % nkibz)
        print(kibz)

        kshift = 0.0 if is_shift is None else 0.5 * np.array(kshift)
        kbz = (grid + kshift) / mesh

        # All k-points and mapping to ir-grid points
        #for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
        #    print("%3d ->%3d %s" % (i, ir_gp_id, (gp + [0.5, 0.5, 0.5]) / mesh))

        return kibz, weights, kbz

    def set_chempot(self, chempot, kmesh, is_shift):
        # Compute DOS/IDOS
        wmesh, dos, idos = _get_cached_dos_idos(kmesh, is_shift)
        #if dos is None:
        #    wmesh, dos, idos = self.get_dos_idos(kmesh, is_shift=is_shift,
        #       method="gaussian", step=0.1, width=0.2, wmesh=None)
        #    self._cache_dos_idos(kmesh, is_shift, wmesh, dos, idos):

        self.chempot = chempot
        # Find number of electrons from new chemical potential from nelect.
        self.nelect = idos.spline(chempot)
        return self.nelect

    def set_nelect(self, nelect, kmesh, is_shift):
        # Compute DOS/IDOS
        wmesh, dos, idos = _get_cached_dos_idos(kmesh, is_shift)
        #if dos is None:
        #    wmesh, dos, idos = self.get_dos_idos(kmesh, is_shift=is_shift,
        #       method="gaussian", step=0.1, width=0.2, wmesh=None)
        #    self._cache_dos_idos(kmesh, is_shift, wmesh, dos, idos):
        self.nelect = nelect
        # Find new chemical potential from nelect.
        idos
        self.chempot = new_chempot
        return new_chempot

    def get_occfacts(self, eigens, kt):
        """occfacts in [0, 1]."""
        if self.occtype == "insulator":
            occfacts = np.ones(eigens.shape)
            occfacts[:, :, self.nelect // 2:] = 0.0
        else:
            return {
                "gaussian": theta_gauss,
                "fermi-dirac": theta_fd,
            }[self.occtype](eigens - self.chempot, kt) * full

    def get_dos_idos(self, kmesh, is_shift=None, method="gaussian", step=0.1, width=0.2, wmesh=None):
        """
        Compute the DOS/IDOS on a linear mesh.

        Args:
            kmesh:
            is_shift:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            mesh: Frequency mesh to use. If None, the mesh is computed automatically from the eigenvalues.

        Returns:
        """
        kibz, weights, kbz = self.get_ibz_weighs_bz(kmesh, is_shift)
        nkibz = len(kibz)

        # Interpolate eigenvalues in the IBZ.
        eigens = self._get_cached_eigens(kmesh, is_shift, "ibz")
        if eigens is None:
            eigens = self.eval_all(kibz)
            self._cache_eigens(kmesh, is_shift, eigens, "ibz")

        # Compute the linear mesh.
        wmesh, step = self._get_wmesh_step(eigens, wmesh, step)
        nw = len(wmesh)
        dos = np.zeros((self.nsppol, nw))

        # TODO: Write cython version.
        if method == "gaussian":
            for spin in range(self.nsppol):
                for ik, wtk in enumerate(weights):
                    for band in range(self.nband):
                        e = eigens[spin, ik, band]
                        dos[spin] += wtk * gaussian(wmesh, width, center=e)

        else:
            raise ValueError("Method %s is not supported" % method)

        idos = None
        return wmesh, dos, idos

    def get_jdos_q0(self, kmesh, is_shift=None, method="gaussian", step=0.1, width=0.2, wmesh=None):
        """
        Compute the join density of states at q==0
            :math:`\sum_{kbv} f_{vk} (1 - f_{ck}) \delta(\omega - E_{ck} + E_{vk})`

        Args:
            kmesh:
            is_shift:
            method: String defining the method.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            wmesh: Frequency mesh to use. If None, the mesh is computed automatically from the eigenvalues.

        Returns:
        """
        kibz, weights, kbz = self.get_ibz_weighs_bz(kmesh, is_shift)
        nkibz = len(kibz)

        # Interpolate eigenvalues in the IBZ.
        eigens = self._get_cached_eigens(kmesh, is_shift, "ibz")
        if eigens is None:
            eigens = self.eval_all(kibz)
            self._cache_eigens(kmesh, is_shift, eigens, "ibz")
        #occfacts = self.

        wmesh, step = self._get_w2mesh_step(eigens, wmesh, step)
        nw = len(wmesh)
        jdos = np.zeros((self.nsppol, nw))
        int_jdos = np.zeros((self.nsppol, nw))

        # Normalize the occupation factors.
        full = 2.0 if self.nsppol == 1 else 1.0

        if method == "gaussian":
            for ik, wtk in enumerate(weights):
                for c in conduction:
                    ec = eigens[spin, ik, c]
                    fc = 1.0 - occfacts[spin, ik, c] / full
                    for v in valence:
                        ev = eigens[spin, ik, v]
                        fv = occfacts[spin, ik, v] / full
                        fact = wtk * fv * fc
                        jdos += fact * gaussian(wmesh, width, center=ec-ev)

        else:
            raise ValueError("Method %s is not supported" % method)

        return wmesh, jdos, int_jdos

    def get_nesting_at_e0(self, qpoints, kmesh, e0, width=0.2, is_shift=None):
        """
        """
        qpoints = np.reshape(qpoints, (-1, 3))
        kibz, weights, kbz = self.get_ibz_weighs_bz(kmesh, is_shift)
        nkbz = len(kbz)

        # Interpolate eigenvalues in the full BZ.
        eigens_kbz = self._get_cached_eigens(kmesh, is_shift, "bz")
        if eigens is None:
            eigens_kbz = self.eval_all(kbz)
            self._cache_eigens(kmesh, is_shift, eigens_kbz, "bz")
        eigens_kbz -= e0
        g_skb = gaussian(eigens_kbz, width)

        # TODO: One could reduce the sum to IBZ(q) with appropriate weight.
        nest_sq = np.zeros((self.nsppol, len(qpoints)))
        for iq, qpt in enumerate(qpoints):
            kpq_bz = kbz + qpt
            eigens_kqbz = self.eval_all(kpq_bz) - e0
            g_skqb = gaussian(eigens_kqbz, width)
            vals = g_skb * g_skqb
            nest_sq[:, iq] = vals.sum(axis=(1, 2))

        nest_sq *= 1. / nkbz
        return nest_sq

    def _get_wmesh_step(self, eigens, wmesh, step):
        if wmesh is not None:
            return wesh, wmesh[1] - wmesh[0]

        # Compute the linear mesh.
        epad = 3.0
        e_min = eigens.min() - epad
        e_max = eigens.max() + epad
        nw = int(1 + (e_max - e_min) / step)
        wmesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)
        return wmesh, step

    def _get_w2mesh_step(self, eigens, wmesh):
        if wmesh is not None:
            return wesh, wmesh[1] - wmesh[0]

        # Compute the linear mesh.
        cmin, cmax = +np.inf, -np.inf
        vmin, vmax = +np.inf, -np.inf
        for c in conduction:
            cmin = min(cmin, eigens[spin,:,c].min())
            cmax = max(cmax, eigens[spin,:,c].max())
        for v in valence:
            vmin = min(vmin, eigens[spin,:,v].min())
            vmax = max(vmax, eigens[spin,:,v].max())

        e_min = cmin - vmax
        e_min -= 0.1 * abs(e_min)
        e_max = cmax - vmin
        e_max += 0.1 * abs(e_max)

        nw = int(1 + (e_max - e_min) / step)
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)

        return wmesh, step

    def _get_cached_eigens(self, kmesh, is_shift, kzone):
        """
        Return a copy of the interpolated eigenvalues associated to (kmesh, is_shift, kzone).
        None if eigens are not available.
        """
        if self.no_cache: return None
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
        if self.no_cache: return
        if not hasattr(self, "_cached_eigens"): self._cached_eigens = OrderedDict()
        kmesh = tuple(kmesh)
        if is_shift is not None: is_shift = tuple(is_shift)
        self._cached_eigens[(kmesh, is_shift, kzone)] = eigens.copy()

    def _get_cached_dos_idos(self, kmesh, is_shift):
        """
        Return the DOS/IDOS obtained from the interpolated eigenvalues associated to (kmesh, is_shift).
        (None, None) if DOS/IDOS is not available.
        """
        if self.no_cache: return None
        if not hasattr(self, "_cached_dos_idos"): self._cached_dos_idos = OrderedDict()
        kmesh = tuple(kmesh)
        if is_shift is not None: is_shift = tuple(is_shift)
        return self._cached_dos_idos.get((kmesh, is_shift))

    def _cache_dos_idos(self, kmesh, is_shift, wmesh, dos, idos):
        """
        Save DOS/IDOS obtained from the interpolated eigenvalues associated to (kmesh, is_shift).
        """
        if self.no_cache: return
        if not hasattr(self, "_cached_dos_idos"): self._cached_dos_idos = OrderedDict()
        kmesh = tuple(kmesh)
        if is_shift is not None: is_shift = tuple(is_shift)
        self._cached_dos_idos[(kmesh, is_shift)] = (wmesh.copy(), dos.dopy(), idos.copy())


class SkwInterpolator(ElectronInterpolator):
    """
    This object implements the Shankland-Koelling-Wood Fourier interpolation scheme.
    It can be used to interpolate functions in k-space with the periodicity of the
    reciprocal lattice and satisfying F(k) = F(Sk) for each rotation S
    belonging to the point group of the crystal. For readability reason,
    the names of the variables are chosen assuming we are interpolating electronic eigenvalues
    but the same object can be used to interpolate other quantities. Just set the first dimension to 1.
    """

    def __init__(self, lpratio, kpts, eigens, cell, symrel, is_time_reversal,
                 filter_params=None, verbose=1):
        """

        Args:
            self.cell = (lattice, positions, numbers)
            lpratio: Ratio between the number of star-functions and the number of ab-initio k-points.
                5 should be OK in many systems, larger values may be required for accurate derivatives.
            kpts: numpy array with the [nkpt, 3] ab-initio k-points in reduced coordinates.
            eigens: numpy array with the ab-initio energies. shape [nsppol, nkpt, nband].
            lattice: numpy array with direct lattice vectors along the rows.
            symrel: [nsym, 3, 3] numpy array with the (ferromagnetic) symmetry operations of the direct lattice
                in reduced coordinates. anti-ferromagnetic symmetries (if any) should be removed by the caller.
            is_time_reversal:
            filter_params:
        """
        self.verbose = verbose
        self.cell = cell
        lattice = self.cell[0]
        #self.original_chempot
        #self.chempot = self.chempot
        #self.nelect = nelect
        self.is_time_reversal = is_time_reversal

        eigens = np.atleast_3d(eigens)
        self.nsppol, self.nkpt, self.nband = eigens.shape
        if len(kpts) != self.nkpt:
            raise ValueError("Second dimension of eigens should be %d but got array of shape: %s" %
                (len(kpts), eigens.shape))
        if self.nkpt == 1:
            raise ValueError("Interpolation algorithm requires nkpt > 1")

        rprimd = np.array(lattice).T
        print("rprimd", rprimd)
        self.rmet = np.matmul(rprimd.T, rprimd)

        #self.get_ibz_weighs_bz([8, 8, 8], is_shift=None)
        #return

        # Find point group operations.
        #symrel = [s.T for s in symrel]
        symrel = np.reshape(symrel, (-1, 3, 3))
        self.ptg_symrel, self.ptg_symrec = self._get_point_group(symrel)
        self.ptg_nsym = len(self.ptg_symrel)
        if self.verbose:
            print("Found", self.ptg_nsym, "symmetries in point group")

        # Find nrwant star points.
        if lpratio <= 1:
            raise ValueError("lpratio must be > 1 but got %s" % lpratio)
        print("nkpt:", self.nkpt)
        nrwant = lpratio * self.nkpt
        rmax = int((1.0 + (lpratio * self.nkpt * self.ptg_nsym) / 2.0) ** (1/3.)) * np.ones(3, dtype=np.int)
        #rmax = 20 * np.ones(3, dtype=np.int)

        while True:
            self.rpts, r2vals, ok = self._find_rstar_gen(nrwant, rmax)
            self.nr = len(self.rpts)
            if ok:
                break
            else:
                print("rmax: ",rmax," was not large enough to find", nrwant, "R-star points.")
                rmax *= 2
                print("Will try again with enlarged rmax:",rmax)

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
            for ik in range(jk + 1):
            #for ik in range(nkpt-1):
                v_ikr = inv_rhor[1:] * (self.skr[ik, 1:] - self.skr[nkpt-1, 1:])
                hmat[ik, jk] = np.vdot(v_jkr, v_ikr)
                if (ik == jk): hmat[ik, jk] = hmat[ik, jk].real

        # Solving system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)..."
        de_kbs = np.empty((nkpt-1, nband, nsppol), dtype=np.complex)
        for spin in range(nsppol):
            for ib in range(nband):
                de_kbs[:, ib, spin] = eigens[spin, 0:nkpt-1, ib] - eigens[spin, nkpt-1, ib]

        # Solve all bands and spins at once
        #call np.linalg.zhesv("U", nkpt-1, nband*nsppol, hmat, nkpt-1, ipiv, lmbs, nkpt-1, work, lwork, ierr)
        try:
            lmb_kbs = scipy.linalg.solve(hmat, np.reshape(de_kbs, (-1, nband * nsppol)),
                    sym_pos=True, lower=False, overwrite_a=True, overwrite_b=True, debug=False, check_finite=False)
                    #sym_pos=False, lower=False, overwrite_a=False, overwrite_b=False, debug=False, check_finite=True)
        except scipy.linalg.LinAlgError as exc:
            print("Cannot solve system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)")
            print("This usually happens there are symmetrical k-points passed to the interpolator.")
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
            for ir in range(1, nr):
                self.coefs[:, :, ir] *= 0.5 * scipy.erfc((np.sqrt(r2vals(ir)) - self.rcut) / self.rsigma)

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
        print("MAE:", mae, " [meV]")
        if mae > 10.0:
            # Issue warning if error too large.
            print("Large error in SKW interpolation!")
            print("MAE:", mae, "[meV]")
        self.mae = mae

    def __str__(self):
        return self.to_string()

    def to_string(self, **kwargs):
        """String representation."""
        lines = []
        app = lines.append
        app("nsppol: %s, nband: %s" % (self.nsppol, self.nband))
        app("Number of operations in point-group: %s with time-reversal: %s" % (self.ptg_nsym, self.is_time_reversal))
        app("Number of ab-initio k-points: %s" % self.nkpt)
        app("Number of star functions: %s [nstars/nk = %s]" % (self.nr, (self.nr / self.nkpt)))
        if self.rcut is not None:
            app("Fourier filter (Eq 9 of PhysRevB.61.1639) with rcut: %s, rsigma: %s" % (self.rcut, self.rsigma))
        app("Comparison between ab-initio data and fit gave Mean Absolute Error: %s [meV]" % self.mae)

        return "\n".join(lines)

    def eval_all(self, kfrac_coords, der1=None, der2=None):
        # Interpolate energies.
        import time
        start = time.time()
        print("Begin interpolation...")

        kfrac_coords = np.reshape(kfrac_coords, (-1, 3))
        new_nkpt = len(kfrac_coords)
        new_eigens = np.empty((self.nsppol, new_nkpt, self.nband))

        for spin in range(self.nsppol):
            for ik, newk in enumerate(kfrac_coords):
                new_eigens[spin, ik] = self.eval_sk(spin, newk)

        if der1 is not None:
            new_der1 = np.empty((self.nsppol, new_nkpt, self.nband, 3))
        if der2 is not None:
            new_der2 = np.empty((self.nsppol, new_nkpt, self.nband, 3, 3))

        print("Interpolation completed", time.time() - start)
        return new_eigens

    def eval_all_and_enforce_degs(self, kfrac_coords, ref_eigens, atol=1e-4):
        kfrac_coords = np.reshape(kfrac_coords, (-1, 3))
        new_nkpt = len(kfrac_coords)
        ref_eigens = np.reshape(ref_eigens, (self.nsppol, new_nkpt, self.nband))

        # Interpolate eigenvales.
        new_eigens = self.eval_all(kfrac_coords)

        # Average interpolated values over degenerates bands.
        for spin in range(self.nsppol):
            for ik in range(new_nkpt):
                for dgbs in find_degs_sk(ref_eigens[spin, ik], atol):
                    if len(dgbs) == 1: continue
                    new_eigens[spin, ik, dgbs] = new_eigens[spin, ik, dgbs].sum() / len(dgbs)

        return new_eigens

    def eval_sk(self, spin, kpt, der1=None, der2=None):
        """
        Interpolate eigenvalues for all bands at a given (spin, k-point).

        Args:
            spin: Spin index
            kpt: K-point in reduced coordinates.

        Return:
            oeigs[nband]
            der1[nband, 3]
            der2[nband, 3, 3]
        """
        # [NB, NR] x [NR]
        oeigs = _np_matmul(self.coefs[spin], self.get_stark(kpt)).real

        if der1 is not None:
            skr_dk1 = self.get_stark_dk1(kpt)
            for ii in range(3):
                der1[:, ii] = np.matmul(self.coefs[spin, :, :], skr_dk1[ii]).real

        if der2 is not None:
            skr_dk2 = self.get_stark_dk2(kpt)
            for jj in range(3):
                for ii in range(jj + 1):
                    der2[:, ii, jj] = np.matmul(self.coefs[spin, :, :], skr_dk2[ii,jj]).real
                    if ii != jj: der2[jj, ii] = der2[ii, jj]

        return oeigs

    def eval_skb(self, spin, kpt, band, der1=None, der2=None):
        """
        Interpolate eigenvalues for a given (spin, k-point, band).

        Args:

        Return:
        """
        # Compute star function for this k-point (if not already in memory)
        if np.any(kpt != self.cached_kpt):
            self.cached_skr, self.cached_kpt = self.get_stark(kpt), kpt

        # Do not take complex conjugate.
        oeig = _np_dot(self.coefs[spin, band], self.cached_skr)
        if der1 is None and der2 is None: return oeig.real

        # TODO: Test Derivatives
        if der1 is not None:
            # Compute first-order derivatives.
            if np.any(kpt != self.cached_kpt_dk1):
                self.cached_skr_dk1, self.cached_kpt_dk1 = self.get_stark_dk1(kpt), kpt

            for ii in range(3):
                der1[ii] = np.dot(self.coefs[spin, band], self.cached_skr_dk1[ii])

            if der2 is None: return oeig.real, der1.real

        if der2 is not None:
            # Compute second-order derivatives.
            if np.any(kpt != self.cached_kpt_dk2):
                self.cached_skr_dk2, self.cached_kpt_dk2 = self.get_stark_dk2(kpt), kpt

            der2 = zero
            for jj in range(3):
                for ii in range(jj + 1):
                    der2[ii, jj] = np.dot(self.coefs[spin, band], self.cached_skr_dk2[ii,jj])
                    if ii != jj: der2[jj, ii] = der2[ii, jj]

        return oeig.real, der1.real, der2.real

    def get_stark(self, kpt):
        """
        Return the star function for a given k-point.

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            complex numpy arrays of shape [self.nr]
        """
        two_pi = 2.0 * np.pi
        skr = np.zeros(self.nr, dtype=np.complex)
        _np_exp = np.exp
        for omat in self.ptg_symrel:
            sk = two_pi * _np_matmul(omat.T, kpt)
            skr += _np_exp(1.j * _np_matmul(self.rpts, sk))
        skr /= self.ptg_nsym

        return skr

    def get_stark_dk1(self, kpt):
        """
        Compute the 1st derivative of the star function wrt k

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            complex numpy array [3, self.nr]  with the derivative of the
            star function wrt k in reduced coordinates.
        """
        srk_dk1 = np.zeros((3, self.nr), dtype=np.complex)
        two_pi = 2.0 * np.pi
        rpts_t = self.rpts.T

        for omat in self.ptg_symrel:
            sk = two_pi * _np_matmul(omat.T, kpt)
            exp_skr = _np_exp(1.j * _np_matmul(self.rpts, sk))
            for ir, rr in enumerate(self.rpts):
                srk_dk1[:, ir] += exp_skr[ir] * _np_matmul(omat, rr)
            #for ir, or in enumerate(_np_matmul(omat, rpts_t).T):
            #    srk_dk1[:, ir] += exp_skr[ir] * or

        srk_dk1 *= 1.j / self%ptg_nsym
        return srk_dk1

    def get_stark_dk2(self, kpt):
        """
        Compute the 2st derivatives of the star function wrt k.

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            Complex numpy array of shape [3, 3, self.nr] with the 2nd derivatives
            of the star function wrt k in reduced coordinates.
        """
        srk_dk2 = np.zeros((3, 3, self.nr), dtype=np.complex)

        #work = zero
        #do isym=1,self%ptg_nsym
        #  sk = two_pi * matmul(transpose(self%ptg_symrel(:,:,isym)), kpt)
        #  do ir=1,self%nr
        #    sr = matmul(self%ptg_symrel(:,:,isym), self%rpts(:,ir))
        #    eiskr = exp(j_dpc * dot_product(sk, self%rpts(:,ir)))
        #    do jj=1,3
        #      do ii=1,jj
        #        work(ii,jj,ir) = work(ii,jj,ir) + eiskr * sr(ii) * sr(jj)
        #      end do
        #    end do
        #  end do
        #end do
        #work = - work / self%ptg_nsym

        #do jj=1,3
        #  do ii=1,jj
        #    srk_dk2(:, ii, jj) = work(ii, jj, :)
        #    if (ii /= jj) srk_dk2(:,jj,ii) = work(:,ii,jj)
        #  end do
        #end do
        return srk_dk2

    def _find_rstar_gen(self, nrwant, rmax):
        """
        Find all lattice points generating the stars inside the supercell defined `rmax`

        Args:
            nrwant: Number of star-functions required.
            rmax: numpy array with the maximum number of cells along the 3 reduced directions.

        Returns:
            (rpts, r2vals, ok)
        """
        import time
        msize = (2 * rmax + 1).prod()
        rtmp = np.empty((msize, 3), dtype=np.int)
        r2tmp = np.empty(msize)
        print("rmax", rmax, "msize:", msize)

        start = time.time()
        for cnt, l in enumerate(itertools.product(range(-rmax[0], rmax[0] + 1),
            range(-rmax[1], rmax[1] + 1), range(-rmax[2], rmax[2] + 1))):
              rtmp[cnt] = l
              r2tmp[cnt] = np.dot(l, np.matmul(self.rmet, l))
        print("gen points", time.time() - start)

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

            rs = set([tuple(rtmp[ss])])
            for ir in range(ss + 1, ee):
                rvec = rtmp[ir]
                if all(tuple(_np_matmul(rot, rvec)) not in rs for rot in self.ptg_symrel):
                    rs.add(tuple(rvec))

            """
            rs = deque([rtmp[ss]])
            for ir in range(ss + 1, ee):
                for rot in self.ptg_symrel:
                    rot_r = np.matmul(rot, rtmp[ir])
                    if any(np.all(rot_r == x) for x in rs): break
                else:
                    # We have new R-point.
                    #print("Adding new point")
                    rs.append(rtmp[ir])
            """

            #print(len(rs), rs)
            rgen.extend(rs)
        #print(rgen)
        print("stars", time.time() - start)

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

        #if True:
        if self.verbose > 10:
            print("nstars:", nstars)
            for r, r2 in zip(rpts, r2vals):
                print(r, r2)
        print("end ", time.time() - start)

        return rpts, r2vals, ok

    def _get_point_group(self, symrel):
        """
        Find the symmetry of the point group. Return (ptg_symrel, ptg_symrec)
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
        ptg_nsym = 2 * tmp_nsym if not has_inversion and self.is_time_reversal else tmp_nsym
        ptg_symrel = np.empty((ptg_nsym, 3, 3), dtype=np.int)
        ptg_symrec = np.empty((ptg_nsym, 3, 3), dtype=np.int)

        from abipy.core.symmetries import mati3inv
        ptg_symrel[:tmp_nsym] = work_symrel[:tmp_nsym]
        for isym in range(tmp_nsym):
            ptg_symrec[isym] = mati3inv(ptg_symrel[isym])

        if not has_inversion and self.is_time_reversal:
            ptg_symrel[tmp_nsym:] = -work_symrel[tmp_nsym:]
            for isym in range(tmp_nsym, ptg_nsym):
                ptg_symrec[isym] = mati3inv(ptg_symrel[isym])

        return ptg_symrel, ptg_symrec

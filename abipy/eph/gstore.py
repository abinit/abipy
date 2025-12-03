"""
This module contains objects for postprocessing e-ph calculations
using the results stored in the GSTORE.nc file.
"""
from __future__ import annotations

import dataclasses
import itertools
import numpy as np
import pandas as pd
#import abipy.core.abinit_units as abu

from functools import cached_property
from monty.string import marquee #, list_strings
from monty.termcolor import cprint
from abipy.core.structure import Structure
from abipy.core.kpoints import kpoints_indices
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, Has_Header #, NotebookWriter
from abipy.tools.typing import PathLike
from abipy.tools.numtools import BzRegularGridInterpolator, nparr_to_df
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible, set_grid_legend,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles)
#from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.tools.typing import Figure
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader


def _allclose(arr_name, array1, array2, verbose: int, rtol=1e-5, atol=1e-8) -> bool:
    """
    Wraps numpy allclose.
    """
    if np.allclose(array1, array2, rtol=rtol, atol=atol):
        if verbose:
            cprint(f"The arrays for {arr_name} are almost equal within the tolerances {rtol=}, {atol=}", color="green")
        return True

    if verbose:
        cprint(f"The arrays for {arr_name} are not almost equal within the tolerances {rtol=}, {atol=}", color="red")

    #differing_indices = np.where(~np.isclose(array1, array2, atol=atol))
    #for index in zip(*differing_indices):
    #    print(f"Difference at index {index}: array1 = {array1[index]}, array2 = {array2[index]}, difference = {abs(array1[index] - array2[index])}")
    return False


class GstoreFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands):
    """
    This file stores the e-ph matrix elements produced by the EPH code of Abinit
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        from abipy.eph.gstore import GstoreFile
        with GstoreFile("out_GSTORE.nc") as gstore:
            print(gstore)

            for spin in range(gstore.nsppol):
                # Extract the object storing the g for this spin.
                gqk = gstore.gqk_spin[spin]
                print(gqk)

                # Get a Dataframe with g(k, q) for all modes and bands.
                df = gqk.get_gdf_at_qpt_kpt([1/2, 0, 0], [0, 0, 0])
                print(df)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreFile
    """
    @classmethod
    def from_file(cls, filepath: PathLike) -> GstoreFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = GstoreReader(filepath)

    @cached_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    def close(self) -> None:
        """Close the file."""
        self.r.close()

    @cached_property
    def gqk_spin(self) -> list:
        return [Gqk.from_gstore(self, spin) for spin in range(self.nsppol)]

    @cached_property
    def has_gwpt(self) -> bool:
        """True if GSTORE contains GWPT matrix elements."""
        return self.r.gtype == "gwpt"

    @cached_property
    def params(self) -> dict:
        """dict with the convergence parameters, e.g. ``nbsum``."""
        #od = OrderedDict([
        #    ("nbsum", self.nbsum),
        #    ("zcut", self.zcut),
        #    ("symsigma", self.symsigma),
        #    ("nqbz", self.r.nqbz),
        #    ("nqibz", self.r.nqibz),
        #])
        ## Add EPH parameters.
        #od.update(self.r.common_eph_params)

        od = {}
        return od

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        app(marquee("Gstore parameters", mark="="))
        app(f"nsppol: {self.r.nsppol}")
        app(f"gstore_completed: {bool(self.r.completed)}")
        app(f"kzone: {self.r.kzone}")
        app(f"kfilter: {self.r.kfilter}")
        app(f"gtype: {self.r.gtype}")
        app(f"qzone: {self.r.qzone}")
        app(f"with_vk: {self.r.with_vk}")
        app(f"kptopt: {self.r.kptopt}")
        app(f"qptopt: {self.r.qptopt}")
        app(f"use_lgk: {self.r.use_lgk}")
        app(f"use_lgq: {self.r.use_lgq}")

        for spin in range(self.r.nsppol):
            app(f"brange_k_spin[{spin}]: {self.r.brange_k_spin[spin]}")
            app(f"brange_kq_spin[{spin}]: {self.r.brange_kq_spin[spin]}")
            app(f"erange_spin[{spin}]: {self.r.erange_spin[spin]}")
            app(f"glob_spin_nq[{spin}]: {self.r.glob_spin_nq[spin]}")

        return "\n".join(lines)

    def check_unfilled_entries_in_gvals(self):
        """
        """
        r = self.r
        for spin in range(self.nsppol):
            # nctkarr_t("gvals", "dp", "gstore_cplex, nb_kq, nb_k, natom3, glob_nk, glob_nq)
            variable = r.read_variable("gvals", path=f"gqk_spin{spin+1}")
            fill_value = variable._FillValue
            # Read the data
            data = variable[:]
            missing_entries = np.where(data == fill_value)
            # Print the indices of missing entries
            print("Missing entries found at indices:", missing_entries)

            if self.r.kfilter == "none":
                raise ValueError("when kfilter == 'none' all the entries in gvals should have been written!")

    @add_fig_kwargs
    def plot_gwpt_vs_qpts(self,
                          kpoint,
                          band_k: int,
                          frac_bounds,
                          band_kq_range: range | list,
                          spin: int = 0,
                          dist_tol: float = 1e-12,
                          fontsize: int = 8,
                          **kwargs) -> Figure:
        """
        """
        if not self.has_gwpt:
            raise ValueError("GSTORE does not contain GWPT matrix elements.")

        gqk = self.gqk_spin[spin]
        in_k = band_k - gqk.bstart_k
        ik_glob, kpoint = self.r.find_ik_glob_kpoint(kpoint, spin)

        reciprocal_lattice = self.structure.reciprocal_lattice
        cart_bounds = reciprocal_lattice.get_cartesian_coords(frac_bounds)

        qpt_cart_coords = []
        for iq in range(gqk.glob_nq):
            iq_bz = self.r.qglob2bz[spin, iq]
            qpt_cart_coords.append(reciprocal_lattice.get_cartesian_coords(self.r.qbz[iq_bz]))

        from abipy.core.kpoints import find_points_along_path
        p = find_points_along_path(cart_bounds, qpt_cart_coords, dist_tol)
        if not len(p.ikfound):
            raise ValueError(f"Cannot find q-points with {dist_tol=}. Check input boundaries or try to increase dist_tol.")
        #p.ikfound, p.dist_list, p.path_ticks)

        natom = len(self.structure)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=natom, ncols=3,
                                               sharex=False, sharey=False, squeeze=False)

        # shape: (glob_nq, natom3, nb_kq)
        g_qpm = gqk.gvals[:,ik_glob,:,:,in_k]
        gks_qpm = gqk.gvals_ks[:,ik_glob,:,:,in_k]

        xs = list(range(len(p.ikfound)))
        for iat, idir in itertools.product(range(natom), range(3)):
            ax = ax_mat[iat, idir]
            ipert = iat + idir
            # (glob_nq, nb_kq) --> (nb_kq, glob_nq).
            gwpt_bq, ks_bq = g_qpm[:,ipert,:].T.copy(), gks_qpm[:,ipert,:].T.copy()

            gwpt_style = dict(marker="o", color="red")
            ks_style = dict(marker=".", color="blue")
            ratio_style = dict(ls="--", color="k")
            #band_kq_range = [0, 4]

            for ib_kq in range(gqk.nb_kq):
                if band_kq_range is not None and ib_kq not in band_kq_range: continue
                gwpt_ys = np.abs(gwpt_bq[ib_kq, p.ikfound])
                ks_ys = np.abs(ks_bq[ib_kq, p.ikfound])
                ax.plot(xs, gwpt_ys, **gwpt_style)
                ax.plot(xs, ks_ys, **ks_style)

                # Plot enhancement ratio.
                ratio_ax = ax.twinx()
                ratio = gwpt_ys / ks_ys
                ratio_ax.plot(xs, ratio, **ratio_style)

            #set_grid_legend(ax, fontsize, xlabel=r"band index (kq)")

        return fig

    @add_fig_kwargs
    def plot_gwpt_vs_ks(self,
                        kpoint,
                        band_k: int,
                        spin: int = 0,
                        fontsize: int = 8,
                        colormap: str = "jet",
                        **kwargs) -> Figure:
        """
        """
        if not self.has_gwpt:
            raise ValueError("GSTORE does not contain GWPT matrix elements.")

        natom = len(self.structure)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=natom, ncols=3,
                                               sharex=False, sharey=False, squeeze=False)

        gqk = self.gqk_spin[spin]
        in_k = band_k - gqk.bstart_k
        ik_glob, kpoint = self.r.find_ik_glob_kpoint(kpoint, spin)

        # shape: (glob_nq, natom3, nb_kq)
        g_qpm = gqk.gvals[:,ik_glob,:,:,in_k]
        gks_qpm = gqk.gvals_ks[:,ik_glob,:,:,in_k]

        cmap = plt.get_cmap(colormap)
        colors = [cmap(iq / gqk.glob_nq) for iq in range(gqk.glob_nq)]
        xs = list(range(gqk.nb_kq))

        for iat, idir in itertools.product(range(natom), range(3)):
            ax = ax_mat[iat, idir]
            ipert = iat + idir
            # (glob_nq, nb_kq)
            gwpt_ys, ks_ys = g_qpm[:,ipert,:], gks_qpm[:,ipert,:]

            for iq in range(gqk.glob_nq):
                numerator = np.abs(gwpt_ys[iq])
                denominator = np.abs(ks_ys[iq])
                #ratio = numerator / denominator
                #ratio = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=np.abs(denominator) > 1e-2)
                #ax.plot(xs, ratio, color=colors[iq])
                ax.plot(xs, np.abs(gwpt_ys[iq]), color=colors[iq], ls="-")
                ax.plot(xs, np.abs(ks_ys[iq]), color=colors[iq], ls="--")

            set_grid_legend(ax, fontsize, xlabel=r"band index (kq)")

        return fig


@dataclasses.dataclass(kw_only=True)
class Gqk:
    """
    This object stores the e-ph matrix elements (g or g^2) and the matrix elements
    of the velocity operator for a given spin.
    """
    spin: int          # Spin index.
    nb_k: int          # Number of bands at k.
    nb_kq: int         # Number of bands at k+q.
    bstart_k: int      # Initial band at k.
    bstart_kq: int     # Initial band at k+q.
    glob_nk: int       # Total number of k/q points in global matrix.
    glob_nq: int       # Note that k-points/q-points can be filtered.

    gstore: GstoreFile

    gvals: np.ndarray             # Array of shape (glob_nk, glob_nk, natom3, nb_kq, nb_k)
                                  # storing complex g(k,q) in the atom representation.

    gvals_ks: np.ndarray | None   # Same as gvals but for KS if we are in GWPT mode.

    vk_cart_ibz: np.ndarray | None
    vkmat_cart_ibz: np.ndarray | None

    def __post_init__(self):
        """Implement consistency check"""
        natom3 = len(self.structure) * 3
        expected_shape = (self.glob_nq, self.glob_nk, natom3, self.nb_kq, self.nb_k)
        if self.gvals.shape != expected_shape:
            raise ValueError(f"{self.gvals.shape=} != {expected_shape=}")
        if self.gvals_ks is not None and self.gvals_ks.shape != expected_shape:
            raise ValueError(f"{self.gvals_ks.shape=} != {expected_shape=}")

    @classmethod
    def from_gstore(cls, gstore: GstoreFile, spin: int) -> Gqk:
        """
        Build an instance from a GstoreFile and the spin index.
        """
        ncr = gstore.r
        path = f"gqk_spin{spin+1}"
        nb_k = ncr.read_dimvalue("nb_k", path=path)
        nb_kq = ncr.read_dimvalue("nb_kq", path=path)
        glob_nk = ncr.read_dimvalue("glob_nk", path=path)
        glob_nq = ncr.read_dimvalue("glob_nq", path=path)

        # Read e-ph matrix elements
        # nctkarr_t("gvals", "dp", "gstore_cplex, nb_kq, nb_k, natom3, glob_nk, glob_nq)
        # Have to transpose the (nb_kq, nb_k) submatrix written by Fortran.
        # Remember that gvals on disk are always complex, in Hartree units and in the atomic representation.

        gvals = ncr.read_value("gvals", path=path).transpose(0, 1, 2, 4, 3, 5).copy()
        gvals = gvals[...,0] + 1j*gvals[...,1]

        # Try to read KS gvals (produced by GWPT code)
        gvals_ks = None
        if gstore.has_gwpt:
            gvals_ks = ncr.read_value("gvals_ks", path=path).transpose(0, 1, 2, 4, 3, 5).copy()
            gvals_ks = gvals_ks[...,0] + 1j*gvals_ks[...,1]

        vk_cart_ibz, vkmat_cart_ibz = None, None
        if ncr.with_vk == 1:
            # nctk_def_arrays(spin_ncid, nctkarr_t("vk_cart_ibz", "dp", "three, nb_k, gstore_nkibz"))
            vk_cart_ibz = ncr.read_value("vk_cart_ibz", path=path)

        if ncr.with_vk == 2:
            # Full (nb_k x nb_k) matrix.
            # Have to transpose (nb_kq, nb_k) submatrix written by Fortran.
            # nctk_def_arrays(spin_ncid, nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb_k, nb_k, gstore_nkibz"))
            vkmat_cart_ibz = ncr.read_value("vkmat_cart_ibz", path=path).transpose(0, 1, 3, 2, 4).copy()
            vkmat_cart_ibz = vkmat_cart_ibz[...,0] + 1j*vkmat_cart_ibz[...,1]

        # Note conversion between Fortran and python indexing.
        bstart_k = ncr.read_value("bstart_k", path=path) - 1
        bstart_kq = ncr.read_value("bstart_kq", path=path) - 1

        data = locals()
        return cls(**{k: data[k] for k in [field.name for field in dataclasses.fields(Gqk)]})

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee(f"Gqk for spin: {self.spin}", mark="="))
        app(f"bstart_k: {self.bstart_k}")
        app(f"nb_k: {self.nb_k}")
        app(f"bstart_kq: {self.bstart_kq}")
        app(f"nb_kq: {self.nb_kq}")
        app(f"glob_nk: {self.glob_nk}")
        app(f"glob_nq: {self.glob_nq}")

        return "\n".join(lines)

    @property
    def structure(self) -> Structure:
        return self.gstore.structure

    @cached_property
    def g2(self) -> np.ndarray:
        """g2 in the atomic representation in Ha^2."""
        return np.abs(self.gvals) ** 2

    #@cached_property
    #def g_nu(self) -> np.ndarray:
    #   """g in the phonon representation in Ha."""
    #   self.gvals

    @cached_property
    def g2_ks(self) -> np.ndarray | None:
        """KS g2 in the atomic representation in Ha^2."""
        if self.gvals_ks is None:
            return None
        return np.abs(self.gvals_ks) ** 2

    #@cached_property
    #def g_ks_nu(self) -> np.ndarray:
    #   """g in the phonon representation in Ha."""
    #   if self.gvals_ks is None:
    #       return None

    def get_dataframe(self, what: str = "g2") -> pd.DataFrame:
        """
        Build and return a dataframe with all the |g(k,q)|^2 if what == "g2" or
        all |v_nk|^2 if what == "v2".
        """
        if what == "g2":
            df = nparr_to_df("g2", self.g2, ["iq", "ik", "imode", "m_kq", "n_k"])

        elif what == "v2":
            if self.vk_cart_ibz is None:
                raise ValueError("vk_cart_ibz is not available in GSTORE!")
            # Compute the squared norm of each vector
            v2 = np.sum(self.vk_cart_ibz ** 2, axis=2)
            df = nparr_to_df("v2", v2, ["ik", "n_k"])

        else:
            raise ValueError(f"Invalid {what=}")

        # Shift band indices.
        df["m_kq"] += self.bstart_kq
        df["n_k"] += self.bstart_k

        return df

    def get_g2q_interpolator_kpoint(self, kpoint, method="linear", check_mesh=1) -> BzRegularGridInterpolator:
        r"""
        Build and return an interpolator that can be used to interpolate g^2(q)

        NB: Invoking the interpolation with an arbitrary q-point returns a numpy array
        of shape (nb_kq, nb_k, natom3) with g_{m_kq n_k, \nu}(q)
        """
        r = self.gstore.r

        # Find the index of the kpoint.
        ik_g, kpoint = r.find_ik_glob_kpoint(kpoint, self.spin)

        # Compute indices of qpoints in the ngqpt mesh.
        ngqpt, shifts = r.ngqpt, [0, 0, 0]
        q_indices = kpoints_indices(r.qbz, ngqpt, shifts, check_mesh=check_mesh)

        natom3 = 3 * len(self.structure)
        nb_k = self.nb_k
        nb_kq = self.nb_kq
        nx, ny, nz = ngqpt
        assert nb_k == nb_kq

        # (glob_nq, glob_nk, natom3, m_kq, n_k)
        g2_qph_mn = self.g2[:,ik_g]

        # Insert g2 in g2_grid
        g2_grid = np.empty((nb_k, nb_kq, natom3, nx, ny, nz))
        for nu in range(natom3):
            for g2_mn, q_inds in zip(g2_qph_mn[:,nu], q_indices):
                ix, iy, iz = q_inds
                g2_grid[:, :, nu, ix, iy, iz] = g2_mn

        return BzRegularGridInterpolator(self.structure, shifts, g2_grid, method=method)

    def get_g_qpt_kpt(self, qpoint, kpoint, what) -> np.ndarray:
        """
        Return numpy array with e-ph matrix elements for the given (qpoint, kpoint) pair.

        Args:
            what="g2" for |g(k,q)|^2, "g" for g(k,q)
        """
        # Find the internal indices of (qpoint, kpoint)
        iq_g, qpoint = self.gstore.r.find_iq_glob_qpoint(qpoint, self.spin)
        ik_g, kpoint = self.gstore.r.find_ik_glob_kpoint(kpoint, self.spin)

        if what == "g2":
            return self.g2[iq_g, ik_g]
        if what == "g":
            return self.gvals[iq_g, ik_g]

        raise ValueError(f"Invalid {what=}")

    def get_gdf_at_qpt_kpt(self, qpoint, kpoint, what="g2") -> pd.DataFrame:
        """
        Build and return a dataframe with the |g(k,q)|^2 for the given (qpoint, kpoint) pair.

        Args:
            what="g2" for |g(k,q)|^2, "g" for g(k,q)
        """
        g2_slice = self.get_g_qpt_kpt(qpoint, kpoint, what)
        df = nparr_to_df(what, g2_slice, ["imode", "m_kq", "n_k"])

        # Shift band indices.
        df["m_kq"] += self.bstart_kq
        df["n_k"] += self.bstart_k

        return df

    def neq(self, other: Gqk, verbose: int) -> int:
        """
        Helper function to compare two GQK objects.
        """
        # This dimensions must agree in order to have a meaningful comparison.
        # so raise immediately if not equal.
        aname_list = ["spin", "nb_k", "nb_kq", "glob_nk", "glob_nq"]

        for aname in aname_list:
            val1, val2 = getattr(self, aname), getattr(other, aname)

            if isinstance(val1, (str, int, float)):
                eq = val1 == val2
            elif isinstance(val1, np.ndarray):
                eq = np.allclose(val1, val2)
            else:
                raise TypeError(f"Don't know how to handle comparison for type: {type(val1)}")

            if not eq:
                raise RuntimeError(f"Different values of {aname=}, {val1=}, {val2=}")

        ierr = 0
        kws = dict(verbose=verbose) # , atol= rtol)

        # Compare v_nk or v_mn_k.
        if self.vk_cart_ibz is not None:
            if not _allclose("vk_cart_ibz", self.vk_cart_ibz, other.vk_cart_ibz, **kws): ierr += 1

        if self.vkmat_cart_ibz is not None:
            if not _allclose("vkmat_cart_ibz", self.vkmat_cart_ibz, other.vkmat_cart_ibz, **kws): ierr += 1

        # Compare g or g^2.
        if not _allclose("g2", self.g2, other.g2, **kws): ierr += 1

        if self.gvals is not None:
            if not _allclose("gvals", self.gvals, other.gvals, **kws): ierr += 1

        return ierr

    #@add_fig_kwargs
    #def plot_g2_hist(self, ax_list=None, **kwargs) -> Figure:

    #    natom = len(self.structure)
    #    nrows, ncols, gridspec_kw = natom, 3, None
    #    ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
    #                                           sharex=True, sharey=True, squeeze=False, gridspec_kw=gridspec_kw)
    #    ax_list = ax_list.ravel()

    #    # (glob_nq, glob_nk, natom3, m_kq, n_k)
    #    for imode, ax in zip(range(natom * 3), ax_list):
    #        data = self.g2[:,:,imode,:,:].flatten()
    #        ax.hist(data)

    #    return fig


class GstoreReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreReader
    """
    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("number_of_spins")
        self.nkbz = self.read_dimvalue("gstore_nkbz")
        self.nkibz = self.read_dimvalue("gstore_nkibz")
        self.nqbz = self.read_dimvalue("gstore_nqbz")
        self.nqibz = self.read_dimvalue("gstore_nqibz")

        # Read important variables.
        self.completed = self.read_value("gstore_completed")
        self.done_spin_qbz = self.read_value("gstore_done_qbz_spin")
        self.with_vk = self.read_value("gstore_with_vk")
        self.qptopt = self.read_value("gstore_qptopt")
        self.kptopt = self.read_value("kptopt")
        self.kzone = self.read_string("gstore_kzone")
        self.qzone = self.read_string("gstore_qzone")
        self.kfilter = self.read_string("gstore_kfilter")
        self.gtype = self.read_string("gstore_gtype")

        # Note conversion Fortran --> C for the isym index.
        self.brange_k_spin = self.read_value("gstore_brange_k_spin")
        self.brange_k_spin[:,0] -= 1
        self.brange_kq_spin = self.read_value("gstore_brange_kq_spin")
        self.brange_kq_spin[:,0] -= 1

        self.erange_spin = self.read_value("gstore_erange_spin")
        # Total number of k/q points for each spin after filtering (if any)
        self.glob_spin_nq = self.read_value("gstore_glob_nq_spin")
        self.glob_nk_spin = self.read_value("gstore_glob_nk_spin")

        # K-points and q-points in the IBZ
        self.kibz = self.read_value("reduced_coordinates_of_kpoints")
        self.qibz = self.read_value("gstore_qibz")

        # K-points and q-points in the BZ
        self.kbz = self.read_value("gstore_kbz")
        self.qbz = self.read_value("gstore_qbz")
        self.ngqpt = self.read_value("gstore_ngqpt")

        # Mapping BZ --> IBZ. Note conversion Fortran --> C for the isym index.
        # nctkarr_t("gstore_kbz2ibz", "i", "six, gstore_nkbz"), &
        # nctkarr_t("gstore_qbz2ibz", "i", "six, gstore_nqbz"), &
        self.kbz2ibz = self.read_value("gstore_kbz2ibz")
        self.kbz2ibz[:,0] -= 1

        self.qbz2ibz = self.read_value("gstore_qbz2ibz")
        self.qbz2ibz[:,0] -= 1

        # Mapping q/k points in gqk --> BZ. Note conversion Fortran --> C for indexing.
        # nctkarr_t("gstore_qglob2bz", "i", "gstore_max_nq, number_of_spins"), &
        # nctkarr_t("gstore_kglob2bz", "i", "gstore_max_nk, number_of_spins") &
        self.qglob2bz = self.read_value("gstore_qglob2bz")
        self.qglob2bz -= 1
        self.kglob2bz = self.read_value("gstore_kglob2bz")
        self.kglob2bz -= 1

        self.use_lgk = self.read_value("gstore_use_lgk")
        self.use_lgq = self.read_value("gstore_use_lgq")

    def find_iq_glob_qpoint(self, qpoint, spin: int) -> tuple:
        """
        Find the internal index of the qpoint needed to access the gvals array.
        Return index of the q-point and qpoint as array.
        """
        qpoint = np.asarray(qpoint)
        for iq_g, iq_bz in enumerate(self.qglob2bz[spin]):
            if np.allclose(qpoint, self.qbz[iq_bz]):
                #print(f"Found {qpoint = } with index {iq_g = }")
                return iq_g, qpoint

        raise ValueError(f"Cannot find {qpoint=} in {self.path}")

    def find_ik_glob_kpoint(self, kpoint, spin: int) -> tuple:
        """
        Find the internal indices of the kpoint needed to access the gvals array.
        Return index of the k-point and kpoint as array.
        """
        kpoint = np.asarray(kpoint)
        for ik_g, ik_bz in enumerate(self.kglob2bz[spin]):
            if np.allclose(kpoint, self.kbz[ik_bz]):
                #print(f"Found {kpoint = } with index {ik_g = }")
                return ik_g, kpoint

        raise ValueError(f"Cannot find {kpoint=} in {self.path}")

    # TODO: This fix to read groups should be imported in pymatgen.
    @cached_property
    def path2group(self) -> dict:
        return self.rootgrp.groups


class GstoreRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple GSTORE.nc files.

    Usage example:

    .. code-block:: python

        robot = GstoreRobot.from_files([
            "t04o_GSTORE.nc",
            "t05o_GSTORE.nc",
            ])


    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreRobot
    """
    EXT = "GSTORE"

    def neq(self, ref_basename: str | None = None, verbose: int = 0) -> int:
        """
        Compare all GSTORE.nc files stored in the GstoreRobot
        """
        # Find reference gstore. By default the first file in the robot is used.
        ref_gstore = self._get_ref_abifile_from_basename(ref_basename)

        exc_list = []
        ierr = 0
        for other_gstore in self.abifiles:
            if ref_gstore.filepath == other_gstore.filepath:
                continue
            print("Comparing: ", ref_gstore.basename, " with: ", other_gstore.basename)
            try:
                ierr += self._neq_two_gstores(ref_gstore, other_gstore, verbose)
                cprint("EQUAL", color="green")
            except Exception as exc:
                exc_list.append(str(exc))

        for exc in exc_list:
            cprint(exc, color="red")

        return ierr

    @staticmethod
    def _neq_two_gstores(self: GstoreFile, gstore2: GstoreFile, verbose: int) -> int:
        """
        Helper function to compare two GSTORE files.
        """
        # These quantities must be the same to have a meaningful comparison.
        aname_list = ["structure", "nsppol", "nkbz", "nkibz",
                      "nqbz", "nqibz", "completed", "kzone", "qzone", "kfilter",
                      "brange_k_spin", "brange_kq_spin", "erange_spin", "glob_spin_nq", "glob_nk_spin",
                     ]

        for aname in aname_list:
            self._compare_attr_name(aname, self, gstore2)

        # Now compare the gkq objects for each spin.
        ierr = 0
        for spin in range(self.nsppol):
            gqk1, gqk2 = self.gqk_spin[spin], gstore2.gqk_spin[spin]
            ierr += gqk1.neq(gqk2, verbose)

        return ierr

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        #for fig in self.get_ebands_plotter().yield_figs(): yield fig

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.GstoreRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        # Mixins
        #nb.cells.extend(self.get_baserobot_code_cells())
        #nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

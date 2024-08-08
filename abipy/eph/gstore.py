"""
This module contains objects for postprocessing e-ph calculations
using the results stored in the GSTORE.nc file.

For a theoretical introduction see :cite:`Giustino2017`
"""
from __future__ import annotations

import dataclasses
import numpy as np
import pandas as pd
#import abipy.core.abinit_units as abu

from monty.string import marquee #, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.structure import Structure
from abipy.core.kpoints import kpoints_indices
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, Has_Header #, NotebookWriter
from abipy.tools.typing import PathLike
from abipy.tools.numtools import BzRegularGridInterpolator, nparr_to_df
#from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
#    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles)
#from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
#from abipy.tools.typing import Figure
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


class GstoreFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands): # , NotebookWriter):
    """
    This file stores the e-ph matrix elements produced by the EPH code of Abinit
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

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

    @lazy_property
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

    @lazy_property
    def gqk_spin(self) -> list:
        return [Gqk.from_gstore(self, spin) for spin in range(self.nsppol)]

    @lazy_property
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

    def to_string(self, verbose=0) -> str:
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

        app(f"nsppol: {self.r.nsppol}")
        app(f"gstore_completed: {bool(self.r.completed)}")
        app(f"gstore_cplex: {self.r.cplex}")
        app(f"gstore_kzone: {self.r.kzone}")
        app(f"gstore_kfilter: {self.r.kfilter}")
        app(f"gstore_gmode: {self.r.gmode}")
        app(f"gstore_qzone: {self.r.qzone}")
        app(f"gstore_with_vk: {self.r.with_vk}")
        app(f"gstore_kptopt: {self.r.kptopt}")
        app(f"gstore_qptopt: {self.r.qptopt}")
        for spin in range(self.r.nsppol):
            app(f"gstore_brange_spin[{spin}]: {self.r.brange_spin[spin]}")
            app(f"gstore_erange_spin[{spin}]: {self.r.erange_spin[spin]}")
            app(f"gstore_glob_spin_nq[{spin}]: {self.r.glob_spin_nq[spin]}")

        return "\n".join(lines)

    def check_unfilled_entries_in_gvals(self):
        """
        """
        r = self.r
        cplex = r.cplex
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


@dataclasses.dataclass(kw_only=True)
class Gqk:
    """
    This object stores the e-ph matrix elements (g or g^2) and the matrix elements
    of the velocity operator for a given spin.
    """
    cplex: int         # 1 if |g|^2 is stored
                       # 2 if complex valued g (mind the gauge)
    spin: int          # Spin index.
    nb: int            # Number of bands
    bstart: int
    #bstop: int

    glob_nk: int       # Total number of k/q points in global matrix.
    glob_nq: int       # Note that k-points/q-points can be filtered.
                       # Use kzone, qzone and kfilter to interpret these dimensions.

    gstore: GstoreFile

    gvals: np.ndarray | None
    g2: np.ndarray | None
    vk_cart_ibz: np.ndarray | None
    vkmat_cart_ibz: np.ndarray | None

    @classmethod
    def from_gstore(cls, gstore: GstoreFile, spin: int):
        """
        Build an istance from a GstoreFile and the spin index.
        """
        ncr = gstore.r
        path = f"gqk_spin{spin+1}"
        cplex = ncr.read_dimvalue("gstore_cplex")
        nb = ncr.read_dimvalue("nb", path=path)
        glob_nk = ncr.read_dimvalue("glob_nk", path=path)
        glob_nq = ncr.read_dimvalue("glob_nq", path=path)

        # Read e-ph matrix elements
        # nctkarr_t("gvals", "dp", "gstore_cplex, nb_kq, nb_k, natom3, glob_nk, glob_nq)
        # Have to transpose the (nb_kq, nb_k) submatrix written by Fortran.
        g2, gvals = None, None
        if cplex == 1:
            g2 = ncr.read_value("gvals", path=path).transpose(0, 1, 2, 4, 3, 5).copy()

        elif cplex == 2:
            gvals = ncr.read_value("gvals", path=path).transpose(0, 1, 2, 4, 3, 5).copy()
            gvals = gvals[...,0] + 1j*gvals[...,1]

        vk_cart_ibz, vkmat_cart_ibz = None, None
        if ncr.with_vk == 1:
            # nctk_def_arrays(spin_ncid, nctkarr_t("vk_cart_ibz", "dp", "three, nb, gstore_nkibz"))
            vk_cart_ibz = ncr.read_value("vk_cart_ibz", path=path)

        if ncr.with_vk == 2:
            # Full (nb x nb) matrix.
            # Have to transpose (nb_kq, nb_k) submatrix written by Fortran.
            # nctk_def_arrays(spin_ncid, nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb, nb, gstore_nkibz"))
            vkmat_cart_ibz = ncr.read_value("vkmat_cart_ibz", path=path).transpose(0, 1, 3, 2, 4).copy()
            vkmat_cart_ibz = vkmat_cart_ibz[...,0] + 1j*vkmat_cart_ibz[...,1]

        # Note conversion between Fortran and python indexing.
        bstart = ncr.read_value("bstart", path=path) - 1
        #bstop = ncr.read_value("stop", path=path)

        data = locals()
        return cls(**{k: data[k] for k in [field.name for field in dataclasses.fields(Gqk)]})

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee(f"Gqk for spin: {self.spin}", mark="="))
        app(f"cplex: {self.cplex}")
        app(f"nb: {self.nb}")
        app(f"bstart: {self.bstart}")
        app(f"glob_nk: {self.glob_nk}")
        app(f"glob_nq: {self.glob_nq}")

        return "\n".join(lines)

    @property
    def structure(self):
        return self.gstore.structure

    def get_dataframe(self, what: str = "g2") -> pd.DataFrame:
        """
        Build and return a dataframe with all the |g(k,q)|^2 if what == "g2" or
        all |v_nk|^2 if what == "v2".
        """
        if what == "g2":
            g2 = self.g2 if self.g2 is not None else np.abs(self.gvals) ** 2
            df = nparr_to_df("g2", g2, ["iq", "ik", "imode", "m_kq", "n_k"])

        elif what == "v2":
            if self.vk_cart_ibz is None:
                raise ValueError("vk_cart_ibz is not available in GSTORE!")
            # Compute the squared norm of each vector
            v2 = np.sum(self.vk_cart_ibz ** 2, axis=2)
            df = nparr_to_df("v2", v2, ["ik", "n_k"])

        else:
            raise ValueError(f"Invalid {what=}")

        #df["m_kq"] += bstart_mkq
        #df["n_k"] += bstart_nk

        return df

    def get_g2q_interpolator_kpoint(self, kpoint, method="linear", check_mesh=1):
        """
        """
        r = self.gstore.r

        # Find the index of the kpoint.
        ik_g, kpoint = r.find_ik_glob_kpoint(kpoint, self.spin)

        # Compute indices of qpoints in the ngqpt mesh.
        ngqpt, shifts = r.ngqpt, [0, 0, 0]
        q_indices = kpoints_indices(r.qbz, ngqpt, check_mesh=check_mesh)

        natom3 = 3 * len(self.structure)
        nb = self.nb
        nx, ny, nz = ngqpt

        # (glob_nq, glob_nk, natom3, m_kq, n_k)
        g2 = self.g2 if self.g2 is not None else np.abs(self.gvals) ** 2
        g2_qph_mn = g2[:,ik_g]

        # Insert g2 in g2_grid
        g2_grid = np.empty((nb, nb, natom3, nx, ny, nz))
        for nu in range(natom3):
            for g2_mn, q_inds in zip(g2_qph_mn[:,nu], q_indices):
                ix, iy, iz = q_inds
                g2_grid[:, :, nu, ix, iy, iz] = g2_mn

        return BzRegularGridInterpolator(self.structure, shifts, g2_grid, method=method)

    def get_g_qpt_kpt(self, qpoint, kpoint, what) -> np.ndarray:
        """
        Return numpy array with e-ph matrix elements the for the given (qpoint, kpoint) pair.

        Args:
            what="g2" for |g(k,q)|^2, "g" for g(k,q)
        """
        # Find the internal indices of (qpoint, kpoint)
        iq_g, qpoint = self.gstore.r.find_iq_glob_qpoint(qpoint, self.spin)
        ik_g, kpoint = self.gstore.r.find_ik_glob_kpoint(kpoint, self.spin)
        if what == "g2":
            g2 = self.g2 if self.g2 is not None else np.abs(self.gvals) ** 2
            return g2[iq_g, ik_g]
        if what == "g":
            if self.cplex != 2:
                raise ValueError("Gstore file stores g2 instead of complex g")
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
        #df["m_kq"] += bstart_mkq
        #df["n_k"] += bstart_nk

        return df

    def neq(self, other: Gqk, verbose: int) -> int:
        """
        Helper function to compare two GQK objects.
        """
        # This dimensions must agree in order to have a meaningfull comparison.
        # so raise immediately if not equal.
        aname_list = ["cplex", "spin", "nb", "glob_nk", "glob_nq"]

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
        if self.g2 is not None:
            if not _allclose("g2", self.g2, other.g2, **kws): ierr += 1

        if self.gvals is not None:
            if not _allclose("gvals", self.gvals, other.gvals, **kws): ierr += 1

        return ierr


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
        self.cplex = self.read_dimvalue("gstore_cplex")
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
        self.gmode = self.read_string("gstore_gmode")

        # Note conversion Fortran --> C for the isym index.
        self.brange_spin = self.read_value("gstore_brange_spin")
        self.brange_spin[:,0] -= 1
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

    def find_iq_glob_qpoint(self, qpoint, spin: int):
        """
        Find the internal index of the qpoint needed to access the gvals array.
        """
        qpoint = np.asarray(qpoint)
        for iq_g, iq_bz in enumerate(self.qglob2bz[spin]):
            if np.allclose(qpoint, self.qbz[iq_bz]):
                #print(f"Found {qpoint = } with index {iq_g = }")
                return iq_g, qpoint

        raise ValueError(f"Cannot find {qpoint = } in GSTORE.nc")

    def find_ik_glob_kpoint(self, kpoint, spin: int):
        """Find the internal indices of the kpoint needed to access the gvals array."""
        kpoint = np.asarray(kpoint)
        for ik_g, ik_bz in enumerate(self.kglob2bz[spin]):
            if np.allclose(kpoint, self.kbz[ik_bz]):
                #print(f"Found {kpoint = } with index {ik_g = }")
                return ik_g, kpoint

        raise ValueError(f"Cannot find {kpoint = } in GSTORE.nc")

    # TODO: This fix to read groups should be imported in pymatgen.
    @lazy_property
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

        robot.neq(verbose=1)

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
    def _neq_two_gstores(gstore1: GstoreFile, gstore2: GstoreFile, verbose: int) -> int:
        """
        Helper function to compare two GSTORE files.
        """
        # These quantities must be the same to have a meaningfull comparison.
        aname_list = ["structure", "nsppol", "cplex", "nkbz", "nkibz",
                      "nqbz", "nqibz", "completed", "kzone", "qzone", "kfilter", "gmode",
                      "brange_spin", "erange_spin", "glob_spin_nq", "glob_nk_spin",
                     ]

        for aname in aname_list:
            self._compare_attr_name(aname, gstore1, gstore2)

        # Now compare the gkq objects for each spin.
        ierr = 0
        for spin in range(gstore1.nsppol):
            gqk1, gqk2 = gstore1.gqk_spin[spin], gstore2.gqk_spin[spin]
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

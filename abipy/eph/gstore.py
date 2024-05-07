"""
This module contains objects for postprocessing e-ph calculations
using the results stored in the GSTORE.nc file.

For a theoretical introduction see :cite:`Giustino2017`
"""
from __future__ import annotations

#import tempfile
#import pickle
#import os
import dataclasses
import numpy as np
#import pandas as pd
#import abipy.core.abinit_units as abu

#from collections import OrderedDict, namedtuple
#from collections.abc import Iterable
#from tabulate import tabulate
from monty.string import marquee, list_strings
from monty.functools import lazy_property
#from monty.termcolor import cprint
from abipy.core.structure import Structure
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, Has_Header, NotebookWriter
#from abipy.core.kpoints import Kpoint, KpointList, Kpath, IrredZone, has_timrev_from_kptopt, find_points_along_path
#from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
#    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles)
#from abipy.tools import duck
#from abipy.tools.numtools import gaussian
from abipy.electrons.ebands import ElectronBands, ElectronDos, RobotWithEbands, ElectronBandsPlotter, ElectronDosPlotter
#from abipy.tools.typing import Figure
#from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader


#@dataclasses.dataclass
#class Gqk
#    cplex: int
#    spin: int
#    natom3: int
#    nb: int
#    bstart: int
#    bstop: int
#    glob_nk: int
#    glob_nq: int
#
#    #gvals: np.ndarray
#    #vk_cart_ib: np.ndarray
#    #vkmat_cart_ib: np.ndarray



class GstoreFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands): # , NotebookWriter):
    """
    This file contains the e-ph matrix elements
    Provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with GstoreFile("out_GSTORE.nc") as ncfile:
            print(ncfile)
            ncfile.ebands.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreFile
    """

    @classmethod
    def from_file(cls, filepath: str) -> SigEPhFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.r = r = GstoreReader(filepath)
        # Get important dimensions.

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

    def __str__(self):
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
        app(f"gstore_compled: {bool(self.r.completed)}")
        app(f"gstore_cplex: {self.r.cplex}")
        app(f"gstore_kzone: {self.r.kzone}")
        app(f"gstore_kfilter: {self.r.kfilter}")
        app(f"gstore_qzone: {self.r.qzone}")
        app(f"gstore_with_vk: {self.r.with_vk}")
        app(f"gstore_kptopt: {self.r.kptopt}")
        app(f"gstore_qptopt: {self.r.qptopt}")
        for spin in range(self.r.nsppol):
            app(f"gstore_brange_spin[{spin}]: {self.r.brange_spin[spin]}")
            app(f"gstore_erange_spin[{spin}]: {self.r.erange_spin[spin]}")
            app(f"gstore_glob_spin_nq[{spin}]: {self.r.glob_spin_nq[spin]}")

        return "\n".join(lines)


class GstoreReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreReader
    """
    def __init__(self, path: str):
        super().__init__(path)

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

        self.brange_spin = self.read_value("gstore_brange_spin")
        self.erange_spin = self.read_value("gstore_erange_spin")
        # Total number of k/q points for each spin after filtering (if any)
        self.glob_spin_nq = self.read_value("gstore_glob_nq_spin")
        self.glob_nk_spin = self.read_value("gstore_glob_nk_spin")

        #print("groups", self.rootgrp.groups)

        #print(self.path2group)
        #print(self.path2group.keys())

        # FIXME
        for spin in range(self.nsppol):
            path = f"gqk_spin{spin+1}"
            nb = self.read_dimvalue("nb", path=path)
            glob_nk = self.read_dimvalue("glob_nk", path=path)
            glob_nq = self.read_dimvalue("glob_nq", path=path)
            # nctkarr_t("gvals", "dp", "gstore_cplex, nb_kq, nb_k, natom3, glob_nk, glob_nq)
            gvals = self.read_value("gvals", path=path)
            # Have to transpose (nb_kq, nb_k) submatrix written by Fortran.

        ## [nsppol, nkcalc] arrays with index of KS bands computed.
        ## Note conversion between Fortran and python convention.
        #self.bstart_sk = self.read_value("bstart_ks") - 1
        #self.nbcalc_sk = self.read_value("nbcalc_ks")
        #self.bstop_sk = self.bstart_sk + self.nbcalc_sk
        #self.max_bstart = self.bstart_sk.max()
        #self.min_bstop = self.bstop_sk.min()

        ## Number of frequency points in Eliashberg functions
        ## This quantity is optional, 0 means *not available*
        #self.gfw_nomega = self.read_dimvalue("gfw_nomega", default=0)

    @lazy_property
    def path2group(self) -> dict:
        return self.rootgrp.groups


def compare_two_gstores(path1: str, path2: str):
    """Compare two GSTORE.nc files."""
    diffs = []
    app = diffs.append

    with GstoreFile(path1) as gstore1, GstoreFile(path2) as gstore2:

        def check_eq(aname: str) -> bool:
            """Helper function to compare two gstore attributes."""
            # Get attributes in gstore first, then in gstore.r, else raise.
            if hasattr(gstore1, aname):
                val1, val2 = getattr(gstore1, aname), getattr(gstore2, aname)
            elif hasattr(gstore1.r, aname):
                val1, val2 = getattr(gstore1.r, aname), getattr(gstore2.r, aname)
            else:
                raise AttributeError(f"Cannot find attribute `{aname}` neither in gstore not in gstore.r")
            #print("aname:", aname); print("val1", val1); print("val2", val2)

            # Now compare val1 and val2 taking into account the type.
            if isinstance(val1, (str, int, float, Structure)):
                eq = val1 == val2

            elif isinstance(val1, np.ndarray):
                eq = np.allclose(val1, val2)

            else:
                raise TypeError(f"Don't know how to handle comparison for type: {type(val1)}")

            if not eq:
                raise AssertionError(f"Different values of {aname=}, {val1=}, {val2=}")

        for aname in ["structure", "nsppol", "cplex", "nkbz", "nkibz",
                      "nqbz", "nqibz", "completed", "kzone", "qzone", "kfilter",
                      "brange_spin", "erange_spin", "glob_spin_nq", "glob_nk_spin",
                     ]:
            check_eq(aname)

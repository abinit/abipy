"""
This module contains objects for postprocessing polaron calculations
using the results stored in the VARPEQ.nc file.

For a theoretical introduction see ...
"""
from __future__ import annotations

import dataclasses
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from monty.string import marquee #, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.structure import Structure
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, Has_Header #, NotebookWriter
from abipy.tools.typing import PathLike
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles, Marker)
#from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.tools.typing import Figure
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader



class VarpeqFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands): # , NotebookWriter):
    """
    This file stores the e-ph matrix elements produced by the EPH code of Abinit
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        from abipy.eph.varpeq import VarpeqFile
        with VarpeqFile("out_VARPEQ.nc") as varpeq:
            print(varpeq)
            varpeq.plot_scf_cycle()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqFile
    """
    @classmethod
    def from_file(cls, filepath: PathLike) -> VarpeqFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = VarpeqReader(filepath)

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
    def polaron_spin(self) -> list:
        return [Polaron.from_varpeq(self, spin) for spin in range(self.r.nsppol)]

    @lazy_property
    def params(self) -> dict:
        """dict with the convergence parameters, e.g. ``nbsum``."""
        #od = OrderedDict([
        #    ("nbsum", self.nbsum),
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
        app("WARNING: Structure is missing")
        #app(self.structure.to_string(verbose=verbose, title="Structure"))

        app("WARNING: Ebands is missing")
        #app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        #if verbose > 1:
        #    app("")
        #    app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        app(f"nsppol: {self.r.nsppol}")
        #app(f"gstore_completed: {bool(self.r.completed)}")
        #app(f"gstore_cplex: {self.r.cplex}")
        #app(f"gstore_kzone: {self.r.kzone}")
        #app(f"gstore_kfilter: {self.r.kfilter}")
        #app(f"gstore_gmode: {self.r.gmode}")
        #app(f"gstore_qzone: {self.r.qzone}")
        #app(f"gstore_with_vk: {self.r.with_vk}")
        #app(f"gstore_kptopt: {self.r.kptopt}")
        #app(f"gstore_qptopt: {self.r.qptopt}")

        return "\n".join(lines)

    @add_fig_kwargs
    def plot_scf_cyle(self, ax_mat=None, fontsize=12, **kwargs) -> Figure:
        """
        Plot the SCF cycle.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        labels = [r'$E_{pol}$',
                  r'$E_{el}$',
                  r'$E_{ph}$',
                  r'$E_{elph}$',
                  r'$\varepsilon$'
        ]

        nsppol = self.r.nsppol
        nstep2cv_spin = self.r.read_value('nstep2cv')
        iter_rec_spin = self.r.read_value('iter_rec')

        # Build grid of plots.
        nrows, ncols = nsppol, 2
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        for spin in range(self.r.nsppol):
            nstep2cv = nstep2cv_spin[spin]
            iterations = iter_rec_spin[spin, :nstep2cv, :] * abu.Ha_eV
            xs = np.arange(1, nstep2cv + 1)

            for iax, ax in enumerate(ax_mat[spin]):
                for ilab, label in enumerate(labels):
                    ys = iterations[:,ilab]
                    if iax == 0:
                        # Plot energies in linear scale.
                        ax.plot(xs, ys, label=label)
                    else:
                        # Plot deltas in logscale.
                        ax.plot(xs, np.abs(ys - ys[-1]), label=label)
                        ax.set_yscale("log")

                ax.set_xlim(1, nstep2cv)
                ax.grid(True)
                ax.legend(loc="best", shadow=True, fontsize=fontsize)
                ax.set_xlabel("Iteration", fontsize=fontsize)
                ax.set_ylabel("Energy (eV)" if iax == 0 else r"$|\Delta|$ Energy (eV)", fontsize=fontsize)
                if nsppol == 2: ax.set_title(f"{spin=}", fontsize=fontsize)

        return fig


@dataclasses.dataclass(kw_only=True)
class Polaron:
    """
    This object stores the polaron coefficients A_kn for a given spin.
    """
    spin: int          # Spin index.
    nk: int            # Number of k-points
    nb: int            # Number of bands.

    #bstart: int        # First band starts at bstart
    #bstop: int
    #glob_nk: int       # Total number of k/q points in global matrix.
    #glob_nq: int       # Note that k-points/q-points can be filtered.
                        # Use kzone, qzone and kfilter to interpret these dimensions.

    #varpeq: VarpeqFile

    a_kn: np.ndarray
    kpoints: np.ndarray
    qpoints: np.ndarray

    @classmethod
    def from_varpeq(cls, varpeq: VarpeqFile, spin: int) -> Polaron:
        """
        Build an istance from a VarpeqFile and the spin index.
        """
        r = varpeq.r
        nk, nq, nb = r.nk_spin[spin], r.nq_spin[spin], r.nb_spin[spin]
        kpoints = r.read_value("kpts_spin")[spin, :nk]
        qpoints = r.read_value("qpts_spin")[spin, :nq]
        a_kn = r.read_value("a_spin", cmode="c")[spin, :nk, :nb]
        b_qnu = r.read_value("b_spin", cmode="c")[spin, :nq]

        data = locals()
        return cls(**{k: data[k] for k in [field.name for field in dataclasses.fields(Polaron)]})

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee(f"Ank for {self.spin=}", mark="="))
        #app(f"nb: {self.nb}")
        #app(f"bstart: {self.bstart}")
        #app(f"bstart: {self.bstop}")
        #app(f"glob_nk: {self.glob_nk}")
        #app(f"glob_nq: {self.glob_nq}")

        return "\n".join(lines)

    def interpolate_amods_n(self, kpoint):

        from abipy.tools.numtools import BlochRegularGridInterpolator

        BlochRegularGridInterpolator(structure, datar)

        from  scipy.interpolate import RegularGridInterpolator
        nx, ny, nz =
        ngkpt = np.array([nx, ny, nz])

        points = (
            np.linspace(0, 1, nx + 1),
            np.linspace(0, 1, ny + 1),
            np.linspace(0, 1, nz + 1),
        ]

        # Transforms x in its corresponding reduced number in the interval [0,1[
        k_indices = []
        for kpt in self.kpoints:
            kpt = kpt % 1
            k_indices.append(kpt * ngkpt)
        k_indices = np.array(k_indices, dtype=int)
        from abipy.tools.numtools import add_periodic_replicas

        interpol_band = []
        for ib in range(self.nb)
            values = np.zeros((nx+1, ny+1, nz+1), dtype=complex)
            for a, inds  in zip(self.a_kn[:,ib], k_indices):
                values[inds] = a
                # Have to fill other portions of the mesh

            rgi = RegularGridInterpolator(points, values, method='linear', bounds_error=True)
            interpol_band.append(rgi)

        return interpol_band


    #def interpolate_bmods_nu(self, qpoint):

    @add_fig_kwargs
    def plot_with_ebands(self, ebands, ax=None, **kwargs) -> Figure:
        """
        """
        ebands = ElectronBands.as_ebands(ebands)
        x, y, s = [], [], []
        for ik, kpoint in enumerate(ebands.kpoints):
            enes_n = ebands.eigens[self.spin, ik, self.bstart:self.bstop]
            amods_n = self.interpolate_amods_n(kpoint)
            for e, a in zip(enes_n, amods_n):
                x.append(ik); y.append(e); s.append(a)

        points = Marker(x, y, s)
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ebands.plot(ax=ax, points=points, show=False)

        return fig


class VarpeqReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqReader
    """

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        #char input_string(input_len) ;
        #int eph_task ;
        #int varpeq_nstep ;
        #int nkbz ;
        #int nqbz ;
        #double tolgrs ;
        #int nstep2cv(nsppol) ;
        #double iter_rec(nsppol, nstep, six) ;
        #int nk_spin(nsppol) ;
        #int nq_spin(nsppol) ;
        #int nb_spin(nsppol) ;
        #double kpts_spin(nsppol, max_nk, three) ;
        #double qpts_spin(nsppol, max_nq, three) ;
        #double a_spin(nsppol, max_nk, max_nb, two) ;
        #double b_spin(nsppol, max_nq, natom3, two) ;

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("nsppol")
        self.nk_spin = self.read_value("nk_spin")
        self.nq_spin = self.read_value("nq_spin")
        self.nb_spin = self.read_value("nb_spin")
        #self.nkbz = self.read_dimvalue("gstore_nkbz")
        #self.nkibz = self.read_dimvalue("gstore_nkibz")
        #self.nqbz = self.read_dimvalue("gstore_nqbz")
        #self.nqibz = self.read_dimvalue("gstore_nqibz")
        self.varpeq_pkind = self.read_string("varpeq_pkind")

        # Read important variables.
        #self.completed = self.read_value("gstore_completed")
        #self.done_spin_qbz = self.read_value("gstore_done_qbz_spin")
        #self.qptopt = self.read_value("gstore_qptopt")
        #self.kptopt = self.read_value("kptopt")
        #self.kzone = self.read_string("gstore_kzone")
        #self.qzone = self.read_string("gstore_qzone")
        #self.kfilter = self.read_string("gstore_kfilter")
        #self.gmode = self.read_string("gstore_gmode")

        #self.brange_spin = self.read_value("gstore_brange_spin")
        #self.erange_spin = self.read_value("gstore_erange_spin")
        # Total number of k/q points for each spin after filtering (if any)
        #self.glob_spin_nq = self.read_value("gstore_glob_nq_spin")
        #self.glob_nk_spin = self.read_value("gstore_glob_nk_spin")

        # K-points and q-points in the IBZ
        #self.kibz = self.read_value("reduced_coordinates_of_kpoints")
        #self.qibz = self.read_value("gstore_qibz")

        # K-points and q-points in the BZ
        #self.kbz = self.read_value("gstore_kbz")
        #self.qbz = self.read_value("gstore_qbz")

        # Mapping BZ --> IBZ. Note conversion Fortran --> C for the isym index.
        # nctkarr_t("gstore_kbz2ibz", "i", "six, gstore_nkbz"), &
        # nctkarr_t("gstore_qbz2ibz", "i", "six, gstore_nqbz"), &
        #self.kbz2ibz = self.read_value("gstore_kbz2ibz")
        #self.kbz2ibz[:,0] -= 1

        #self.qbz2ibz = self.read_value("gstore_qbz2ibz")
        #self.qbz2ibz[:,0] -= 1

        # Mapping q/k points in gqk --> BZ. Note conversion Fortran --> C for indexing.
        # nctkarr_t("gstore_qglob2bz", "i", "gstore_max_nq, number_of_spins"), &
        # nctkarr_t("gstore_kglob2bz", "i", "gstore_max_nk, number_of_spins") &
        #self.qglob2bz = self.read_value("gstore_qglob2bz")
        #self.qglob2bz -= 1

        #self.kglob2bz = self.read_value("gstore_kglob2bz")
        #self.kglob2bz -= 1

    #def find_iq_glob_qpoint(self, qpoint, spin: int):
    #    """Find the internal index of the qpoint needed to access the gvals array."""
    #    qpoint = np.asarray(qpoint)
    #    for iq_g, iq_bz in enumerate(self.qglob2bz[spin]):
    #        if np.allclose(qpoint, self.qbz[iq_bz]):
    #            return iq_g, qpoint

    #    raise ValueError(f"Cannot find {qpoint=} in GSTORE.nc")

    #def find_ik_glob_kpoint(self, kpoint, spin: int):
    #    """Find the internal indices of the kpoint needed to access the gvals array."""
    #    kpoint = np.asarray(kpoint)
    #    for ik_g, ik_bz in enumerate(self.kglob2bz[spin]):
    #        if np.allclose(kpoint, self.kbz[ik_bz]):
    #            return ik_g, kpoint

    #    raise ValueError(f"Cannot find {kpoint=} in GSTORE.nc")



#class VarpeqRobot(Robot, RobotWithEbands):
#    """
#    This robot analyzes the results contained in multiple VARPEQ.nc files.
#
#    Usage example:
#
#    .. code-block:: python
#
#        robot = VarpeqRobot.from_files([
#            "out1_VARPEQ.nc",
#            "out2_VARPEQ.nc",
#            ])
#
#        robot.neq(verbose=1)
#
#    .. rubric:: Inheritance Diagram
#    .. inheritance-diagram:: VarpeqRobot
#    """
#    EXT = "VARPEQ"
#
#    def neq(self, verbose: int, ref_basename: str | None) -> int:
#        """
#        Compare all GSTORE.nc files stored in the VarpeqRobot
#        """
#        exc_list = []
#
#        # Find reference gstore. By default the first file in the robot is used.
#        ref_gstore = self.abifiles[0]
#        if ref_basename is not None:
#            for i, gstore in enumerate(self.abifiles):
#                if gstore.basename == ref_basename:
#                    ref_gstore = gstore
#                    break
#            else:
#                raise ValueError(f"Cannot find {ref_basename=}")
#
#        ierr = 0
#        for other_gstore in self.abifiles:
#            if ref_gstore.filepath == other_gstore.filepath:
#                continue
#            print("Comparing ", ref_gstore.basename, " with: ", other_gstore.basename)
#            try:
#                ierr += self._neq_two_gstores(ref_gstore, other_gstore, verbose)
#                cprint("EQUAL", color="green")
#            except Exception as exc:
#                exc_list.append(str(exc))
#
#        for exc in exc_list:
#            cprint(exc, color="red")
#
#        return ierr
#
#    @staticmethod
#    def _neq_two_gstores(gstore1: GstoreFile, gstore2: GstoreFile, verbose: int) -> int:
#        """
#        Helper function to compare two GSTORE files.
#        """
#        # These quantities must be the same to have a meaningfull comparison.
#        aname_list = ["structure", "nsppol", "cplex", "nkbz", "nkibz",
#                      "nqbz", "nqibz", "completed", "kzone", "qzone", "kfilter", "gmode",
#                      "brange_spin", "erange_spin", "glob_spin_nq", "glob_nk_spin",
#                     ]
#
#        for aname in aname_list:
#            # Get attributes in gstore first, then in gstore.r, else raise.
#            if hasattr(gstore1, aname):
#                val1, val2 = getattr(gstore1, aname), getattr(gstore2, aname)
#            elif hasattr(gstore1.r, aname):
#                val1, val2 = getattr(gstore1.r, aname), getattr(gstore2.r, aname)
#            else:
#                raise AttributeError(f"Cannot find attribute `{aname=}` neither in gstore not in gstore.r")
#
#           # Now compare val1 and val2 taking into account the type.
#            if isinstance(val1, (str, int, float, Structure)):
#                eq = val1 == val2
#            elif isinstance(val1, np.ndarray):
#                eq = np.allclose(val1, val2)
#            else:
#                raise TypeError(f"Don't know how to handle comparison for type: {type(val1)}")
#
#            if not eq:
#                raise RuntimeError(f"Different values of {aname=}, {val1=}, {val2=}")
#
#        # Now compare the gkq objects for each spin.
#        ierr = 0
#        for spin in range(gstore1.nsppol):
#            gqk1, gqk2 = gstore1.gqk_spin[spin], gstore2.gqk_spin[spin]
#            ierr += gqk1.neq(gqk2, verbose)
#        return ierr
#
#    def yield_figs(self, **kwargs):  # pragma: no cover
#        """
#        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
#        Used in abiview.py to get a quick look at the results.
#        """
#        #yield self.plot_lattice_convergence(show=False)
#        #yield self.plot_gsr_convergence(show=False)
#        #for fig in self.get_ebands_plotter().yield_figs(): yield fig
#
#    def write_notebook(self, nbpath=None) -> str:
#        """
#        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
#        working directory is created. Return path to the notebook.
#        """
#        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)
#
#        args = [(l, f.filepath) for l, f in self.items()]
#        nb.cells.extend([
#            #nbv.new_markdown_cell("# This is a markdown cell"),
#            nbv.new_code_cell("robot = abilab.VarpeqRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
#            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
#        ])
#
#        # Mixins
#        #nb.cells.extend(self.get_baserobot_code_cells())
#        #nb.cells.extend(self.get_ebands_code_cells())
#
#        return self._write_nb_nbpath(nb, nbpath)

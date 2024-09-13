"""
This module contains objects for analyzing
the PATH.nc file with the e-ph matrix elements along a k/q path
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
from abipy.core.kpoints import Kpath
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.typing import PathLike
#from abipy.tools.numtools import BzRegularGridInterpolator, nparr_to_df
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles, Marker, set_grid_legend)
#from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.dfpt.phonons import PhononBands
from abipy.dfpt.phtk import NonAnalyticalPh
from abipy.tools.typing import Figure
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader


class GpathFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    This file stores the e-ph matrix elements along a k/q path
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with GpathFile("out_GPATH.nc") as gpath:
            print(gpath)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GpathFile
    """

    @classmethod
    def from_file(cls, filepath: PathLike) -> GpathFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = GpathReader(filepath)

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def ebands(self) -> ElectronBands:
        """Electron bands on the k+q path"""
        return self.r.read_ebands()

    @lazy_property
    def phbands(self) -> PhononBands:
        """Phonon bands on nq_path points."""
        return self.r.read_phbands()

    def close(self) -> None:
        """Close the file."""
        self.r.close()

    @lazy_property
    def params(self) -> dict:
        """dict with the convergence parameters, e.g. ``nbsum``."""
        #od = OrderedDict([
        #    ("nbsum", self.nbsum),
        #    ("nqibz", self.r.nqibz),
        #])
        ## Add EPH parameters.
        #od.update(self.r.common_eph_params)

        od = {}
        return od

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))

        #app(f"nsppol: {self.r.nsppol}")
        #app(f"gstore_cplex: {self.r.cplex}")
        #app(f"gstore_kptopt: {self.r.kptopt}")
        #app(f"gstore_qptopt: {self.r.qptopt}")

        return "\n".join(lines)

    @add_fig_kwargs
    def plot_g_qpath(self, band_range=None, which_g="sym", with_q=False, scale=1, ax_mat=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot ...

        Args:
            band_range:
            with_q
            average_mode:
            ax_mat: List of |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
        """
        nrows, ncols = 2, self.r.nsppol
        #nrows, ncols = 3, self.r.nsppol
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        qnorms = np.array([qpt.norm for qpt in self.phbands.qpoints]) if with_q else \
                 np.ones(len(self.phbands.qpoints))
        #band_range = (self.r.bstart, self.r.bstop) if band_range is None else band_range

        for spin in range(self.r.nsppol):
            g_nuq_sym, g_nuq_unsym = self.r.get_gnuq_average_spin(spin, band_range)
            g_nuq = g_nuq_sym if which_g == "sym" else g_nuq_unsym

            ax = ax_mat[0, spin]
            for mode in range(self.r.natom3):
                ax.plot(g_nuq[mode] * qnorms, label=f"{which_g} {mode=}")

            self.phbands.decorate_ax(ax, units="meV")
            set_grid_legend(ax, fontsize, xlabel=r"Wavevector $\mathbf{q}$", ylabel=r"$|g^{\text{avg}}| (meV)$")

            marker_color = "gold"
            x, y, s = [], [], []
            for iq, qpoint in enumerate(self.phbands.qpoints):
                omegas_nu = self.phbands.phfreqs[iq,:]
                for w, g2 in zip(omegas_nu, g_nuq[:,iq], strict=True):
                    x.append(iq); y.append(w); s.append(scale * g2)

            points = Marker(x, y, s, color=marker_color, edgecolors='gray', alpha=0.8, label=r'$|g^{\text{avg}}(\mathbf{q})|$ (meV)')

            ax = ax_mat[1, spin]
            self.phbands.plot(ax=ax, points=points, show=False)
            set_grid_legend(ax, fontsize, xlabel=r"Wavevector $\mathbf{q}$")

            #self.ebands.plot(ax=ax, points=points, show=False)
            #ax.plot(self.r.all_eigens_kq[spin]) # .transpose()) #, label=f"
            #set_grid_legend(ax, fontsize, xlabel=r"$\mathbf{k+q}$  xlabel=r"Wavevector $\mathbf{k}$")

        return fig

    #@add_fig_kwargs
    #def plot_g_kpath(self, band_range=None, which_g="sym", with_k=False, ax_mat=None, fontsize=8, **kwargs) -> Figure:

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.ebands.plot(show=False)
        if self.r.eph_fix_korq == "k":
            yield self.phbands.plot(show=False)
            yield self.plot_g_qpath()
        #if self.r.eph_fix_korq == "q":
        #    yield self.plot_g_kpath()

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("gpath = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(gpath)"),
            nbv.new_code_cell("gpath.ebands.plot();"),
            nbv.new_code_cell("gpath.phbands.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class GpathReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GpathReader
    """
    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("nsppol")
        self.natom = self.read_dimvalue("natom")
        self.natom3 = self.read_dimvalue("natom3")
        self.nb_in_g = self.read_dimvalue("nb_in_g")
        self.nk_path = self.read_dimvalue("nk_path")
        self.nq_path = self.read_dimvalue("nq_path")

        # eigens are in Ha, phfreq are in eV for historical reason
        self.phfreqs_ha = self.read_value("phfreqs")
        self.all_eigens_k = self.read_value("all_eigens_k")
        self.all_eigens_kq = self.read_value("all_eigens_kq")

        # Read important variables.
        self.eph_fix_korq = self.read_string("eph_fix_korq")
        self.eph_fix_wavec = self.read_value("eph_fix_wavevec")
        #self.completed = self.read_value("gstore_completed")

        # Note conversion Fortran --> C for the bstart index.
        nband = self.read_dimvalue("nband")
        self.bstart = self.read_value("bstart") - 1
        self.bstop = self.read_value("bstop")
        self.band_range = [self.bstart, self.bstop]

    def read_ebands(self):
        """
        Overrides method of superclass as we cannot rely of the etsf-io file format,
        and we have to build the ebands manually.
        """
        structure = self.read_structure()

        nspinor = self.read_dimvalue("nspinor")
        nspden = self.read_dimvalue("nspden")
        nelect = self.read_value("nelect")
        fermie = self.read_value("fermie") * abu.Ha_eV

        kpath_frac_coords = self.read_value("kpoints")
        qpath_frac_coords = self.read_value("qpoints")
        frac_coords = kpath_frac_coords + qpath_frac_coords
        path = Kpath(structure.lattice.reciprocal_lattice, frac_coords, weights=None, names=None, ksampling=None)

        # eigens are in Ha
        if self.nq_path > 1:
            all_eigens = self.read_value("all_eigens_kq") * abu.Ha_eV
        else:
            all_eigens = self.read_value("all_eigens_k") * abu.Ha_eV

        occfacts = np.zeros_like(all_eigens)

        return ElectronBands(structure, path, all_eigens, fermie, occfacts, nelect, nspinor, nspden)

    def read_phbands(self) -> PhononBands:
        """
        Read the phonon band structure along the q-path.
        """
        amu_list = self.read_value("atomic_mass_units")
        atomic_numbers = self.read_value("atomic_numbers")
        amu = {at: a for at, a in zip(atomic_numbers, amu_list)}

        # phfreqs are in eV for historical reason
        phfreqs = self.read_value("phfreqs")
        # Complex array with the Cartesian displacements in **Angstrom**
        phdispl_cart = self.read_value("phdispl_cart", cmode="c") * abu.Bohr_Ang

        structure = self.read_structure()
        qpath_frac_coords = self.read_value("qpoints")

        path_qq = Kpath(structure.lattice.reciprocal_lattice, qpath_frac_coords, weights=None, names=None, ksampling=None)

        non_anal_ph = NonAnalyticalPh.from_ncreader(self) if "non_analytical_directions" in self.rootgrp.variables else None

        return PhononBands(structure=structure, qpoints=path_qq, phfreqs=phfreqs, phdispl_cart=phdispl_cart, amu=amu,
                           non_anal_ph=non_anal_ph,
                           # TODO ?
                           #epsinf=epsinf,
                           #zcart=zcart,
                           )

    def get_gnuq_average_spin(self, spin: int, band_range: list|tuple, eps_mev: float=0.01):
        """
        Average ...

        Args:
            spin: Spin index
            band_range:
            eps_mev: Tolerance in meV used to detect degeneracies.
        """
        # Consistency check
        if self.nk_path != 1:
            raise ValueError(f"{self.nk_path=} != 1. In this case, one cannot ask for q-dependent g(k,q)!")

        eps_ha = eps_mev / abu.Ha_meV
        eps_ev = eps_ha * abu.Ha_eV

        nsppol, natom3 = self.nsppol, self.natom3
        # Number of m, n bands in g_mn, the first band starts at bstart.
        nb_in_g = self.nb_in_g
        bstart = self.bstart

        all_eigens_k = self.all_eigens_k    # eV units
        all_eigens_kq = self.all_eigens_kq  # eV units
        phfreqs_ha = self.phfreqs_ha        # Ha units

        # Now read the e-ph matrix elements. On disk we have
        #                                                  n-index, m-index
        # double gkq2_nu(nsppol, nk_path, nq_path, natom3, nb_in_g, nb_in_g) ;
		#  gkq2_nu:_FillValue = -1. ;
        #
        # in Ha^2 with nk_path == 1
        #                                      m-index, n-index
        # In memory we want: (nq_path, natom3, nb_in_g, nb_in_g)

        absg = np.sqrt(self.read_variable("gkq2_nu")[spin, 0][:].transpose(0, 1, 3, 2).copy()) * abu.Ha_meV
        absg_unsym = absg.copy()

        # Average over phonons.
        absg_sym = np.zeros_like(absg)

        for iq in range(self.nq_path):
            for nu in range(natom3):
                # Find all mu where |w_1 - w_nu| < eps_ha
                mask = np.abs(phfreqs_ha[iq, :] - phfreqs_ha[iq, nu]) < eps_ha
                # Sum the squared values of absg over the selected mu indices
                g2_mn = np.sum(absg[iq, mask, :, :]**2, axis=0)
                # Compute the symmetrized value (nn is the number of valid matches)
                nn = np.sum(mask)
                absg_sym[iq, nu, :, :] = np.sqrt(g2_mn / nn)

        # Average over k electrons.
        # MG FIXME: Note the difference with a similar function in gkq here I use absg and not absgk
        absg = absg_sym.copy()
        g2_nu = np.zeros((natom3), dtype=float)
        for iq in range(self.nq_path):
            for jbnd in range(nb_in_g):
                for ibnd in range(nb_in_g):
                    w_1 = all_eigens_k[spin, 0, ibnd]
                    g2_nu[:], nn = 0.0, 0
                    for pbnd in range(nb_in_g):
                        w_2 = all_eigens_k[spin, 0, pbnd]
                        if abs(w_2 - w_1) >= eps_ev: continue
                        nn += 1
                        g2_nu += absg[iq,:,jbnd,pbnd] ** 2
                    absg_sym[iq,:,jbnd,ibnd] = np.sqrt(g2_nu / nn)

        # Average over k electrons.
        # MG FIXME: Note the difference with a similar function in gkq here I use absg and not absgk
        #absg = absg_sym.copy()
        #for iq in range(self.nq_path):
        #    for jbnd in range(nb_in_g):
        #        for ibnd in range(nb_in_g):
        #            w_1 = all_eigens_k[spin, 0, ibnd]
        #            # Create mask to find all pbnd where |w_2 - w_1| < eps_ev
        #            mask = np.abs(all_eigens_k[spin, 0, :] - w_1) < eps_ev
        #            # Sum the squared values of absg over the selected pbnd indices (nn is the number of matching energy levels)
        #            nn = np.sum(mask)
        #            g2_nu = np.sum(absg[iq, :, jbnd, mask] ** 2, axis=-1)
        #            print(f"{mask=}")
        #            print(f"{absg[iq, :, jbnd, mask].shape=}")
        #            print(f"{absg.shape=}")
        #            print(f"{mask.shape=}")
        #            print(f"{g2_nu.shape=}")
        #            # Compute the symmetrized absg
        #            absg_sym[iq, :, jbnd, ibnd] = np.sqrt(g2_nu / nn)

        # Average over k+q electrons.
        absg = absg_sym.copy()
        for iq in range(self.nq_path):
          for ibnd in range(nb_in_g):
              for jbnd in range(nb_in_g):
                  w_1 = all_eigens_kq[spin, iq, jbnd]
                  g2_nu[:], nn = 0.0, 0
                  for pbnd in range(nb_in_g):
                      w_2 = all_eigens_kq[spin, iq, pbnd]
                      if abs(w_2 - w_1) >= eps_ev: continue
                      nn += 1
                      g2_nu += absg[iq,:,pbnd,ibnd] ** 2
                  absg_sym[iq,:,jbnd,ibnd] = np.sqrt(g2_nu / nn)

        # (nq_path, natom3, nb_in_g, nb_in_g) -> (natom3, nq_path, nb_in_g, nb_in_g)
        absg_sym = absg_sym.transpose(1, 0, 2, 3).copy()
        absg_unsym = absg_unsym.transpose(1, 0, 2, 3).copy()
        print(f"{absg_unsym.shape=}")

        # Average over bands 1/ n_b_in**2 sum_{mn}
        absg_sym = np.sum(absg_sym, axis=(-2, -1)) / nb_in_g**2
        absg_unsym = np.sum(absg_unsym, axis=(-2, -1)) / nb_in_g**2

        # [spin, nq_path]
        return absg_sym, absg_unsym


class GpathRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple GPATH.nc files.

    Usage example:

    .. code-block:: python

        robot = GpathRobot.from_files([
            "t04o_GPATH.nc",
            "t05o_GPATH.nc",
            ])


    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreRobot
    """
    EXT = "GPATH"

    #def neq(self, ref_basename: str | None = None, verbose: int = 0) -> int:
    #    """
    #    Compare all GPATHE.nc files stored in the robot
    #    """
    #    # Find reference gstore. By default the first file in the robot is used.
    #    ref_gstore = self._get_ref_abifile_from_basename(ref_basename)

    #    exc_list = []
    #    ierr = 0
    #    for other_gstore in self.abifiles:
    #        if ref_gstore.filepath == other_gstore.filepath:
    #            continue
    #        print("Comparing: ", ref_gstore.basename, " with: ", other_gstore.basename)
    #        try:
    #            ierr += self._neq_two_gstores(ref_gstore, other_gstore, verbose)
    #            cprint("EQUAL", color="green")
    #        except Exception as exc:
    #            exc_list.append(str(exc))

    #    for exc in exc_list:
    #        cprint(exc, color="red")

    #    return ierr

    #@staticmethod
    #def _neq_two_gstores(gstore1: GstoreFile, gstore2: GstoreFile, verbose: int) -> int:
    #    """
    #    Helper function to compare two GSTORE files.
    #    """
    #    # These quantities must be the same to have a meaningfull comparison.
    #    aname_list = ["structure", "nsppol", "cplex", "nkbz", "nkibz",
    #                  "nqbz", "nqibz", "completed", "kzone", "qzone", "kfilter", "gmode",
    #                  "brange_spin", "erange_spin", "glob_spin_nq", "glob_nk_spin",
    #                 ]

    #    for aname in aname_list:
    #        self._compare_attr_name(aname, gstore1, gstore2)

    #    # Now compare the gkq objects for each spin.
    #    ierr = 0
    #    for spin in range(gstore1.nsppol):
    #        gqk1, gqk2 = gstore1.gqk_spin[spin], gstore2.gqk_spin[spin]
    #        ierr += gqk1.neq(gqk2, verbose)

    #    return ierr

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

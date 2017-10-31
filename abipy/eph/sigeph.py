# coding: utf-8
"""
This module contains objects for the postprocessing of Sigma_eph calculations.

Warning:

    Work in progress, DO NOT USE THIS CODE
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pandas as pd
import pymatgen.core.units as units
import abipy.core.abinit_units as abu
from abipy.tools import duck

#from collections import OrderedDict
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import KpointList
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands
#from abipy.dfpt.phonons import PhononBands, RobotWithPhbands, factor_ev2units, unit_tag, dos_label_from_units
from abipy.abio.robots import Robot


# TODO QPState and QPList from electrons.gw (Define base abstract class?).


class SigEPhFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the Fan-Migdal self-energy, the ElectronBands on the k-mesh.
    Provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with SigEPhFile("out_SIGEPH.nc") as ncfile:
            print(ncfile)
            ncfile.ebands.plot()
    """
    # Markers used for up/down bands.
    marker_spin = {0: "^", 1: "v"}
    color_spin = {0: "k", 1: "r"}

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(SigEPhFile, self).__init__(filepath)
        self.reader = r = SigmaPhReader(filepath)

        # Get important dimensions.
        self.nkcalc = r.read_dimvalue("nkcalc")
        #self.max_nbcalc = r.read_dimvalue("max_nbcalc")
        self.ntemp = r.read_dimvalue("ntemp")
        #self.nqbz = r.read_dimvalue("nqbz")
        #self.nqibz = r.read_dimvalue("nqibz")
        self.ngqpt = r.read_value("ngqpt")

        self.symsigma = r.read_value("symsigma")
        # TODO zcut?
        self.zcut = r.read_value("eta")
        #self.nbsum = r.read_value("nbsum")

        self.ktmesh = r.read_value("kTmesh")
        self.tmesh = self.ktmesh / abu.kb_HaK

        # The K-points where QP corrections have been calculated.
        # gwkpoints in gw.py
        self.sigma_kpoints = KpointList(self.structure.reciprocal_lattice, r.read_value("kcalc"))

        # [nsppol, nkcalc] arrays with index of KS bands computed.
        # Note conversion between Fortran and python convention.
        self.bstart_sk = r.read_value("bstart_ks") - 1
        self.nbcalc_sk = r.read_value("nbcalc_ks")
        self.bstop_sk = self.bstart_sk + self.nbcalc_sk

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="KS Electronic Bands"))
        app("")
        # SigmaEPh section.
        app(marquee("SigmaEPh calculation", mark="="))
        app("Number of k-points computed: %d" % (self.nkcalc))
        #app("K-mesh: %s" % (str(self.ngqpt)))
        app("Q-mesh: %s" % (str(self.ngqpt)))
        app("zcut: %s [Ha]" % (self.zcut))
        app("Number of temperatures: %d, from %.1f to %.1f [K]" % (self.ntemp, self.tmesh[0], self.tmesh[-1]))
        #app("Number of bands in self-energy sum: %d" % (self.nbsum))

        app("K-points and bands for self-energy corrections:")
        for spin in range(self.nsppol):
            for ik, kpoint in enumerate(self.sigma_kpoints):
                post = "ik: %d" % ik if self.nsppol == 1 else "ik: %d, spin: %d" % (ik, spin)
                app("\t%s: bstart: %d, bstop: %d, %s" % (
                    repr(kpoint), self.bstart_sk[spin, ik], self.bstop_sk[spin, ik], post))

        return "\n".join(lines)

    @lazy_property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    def close(self):
         """Close the file."""
         self.reader.close()

    @lazy_property
    def ks_dirgaps(self):
        """
        Array of shape [nsppol, nkcalc] with the KS gaps in eV ordered as kcalc.
        """
        return self.reader.read_value("ks_gaps") * units.Ha_to_eV

    @lazy_property
    def qp_dirgaps_t(self):
        """
        Array of shape [nsppol, nkcalc, ntemp] with the QP gaps in eV ordered as kcalc.
        """
        return self.reader.read_value("qp_gaps") * units.Ha_to_eV

    #lazy_property
    #def mu_e(self):
    #    """mu_e[ntemp] chemical potential (eV) of electrons for the different temperatures."""
    #    return self.reader.read_value("mu_e") * units.Ha_to_eV

    #integer,allocatable :: kcalc2ibz(:,:)
    #!kcalc2ibz (nkcalc, 6))
    #! Mapping kcalc --> ibz as reported by listkk.

    def sigkpt2index(self, sig_kpoint):
        """
        Returns the index of the self-energy k-point in sigma_kpoints
        Used to access data in the arrays that are dimensioned [0:nkcalc]
        """
        return int(sig_kpoint) if duck.is_intlike(sig_kpoint) else self.sigma_kpoints.index(sig_kpoint)

    #@property
    #def params(self):
    #    """AttrDict dictionary with the convergence parameters, e.g. nbsum."""
    #    return self.reader.read_params()

    #@lazy_property
    #def qplist_spin(self):
    #    """Tuple of :class:`QPList` objects indexed by spin."""
    #    return self.reader.read_allqps()

    #def get_sigmaw(self, spin, kpoint, band):
    #    """"
    #    Read self-energy(t, w) for (spin, kpoint, band)
    #    Return :class:`Function1D` object
    #    """
    #    wmesh, sigxc_values = self.reader.read_sigmaw(spin, kpoint, band)
    #    wmesh, spf_values = self.reader.read_spfunc(spin, kpoint, band)
    #    return Sigmaw(spin, kpoint, band, wmesh, sigxc_values, spf_values)

    def get_dataframe(self, ignore_imag=False):
        """
        Returns pandas DataFrame with QP results for all k-points included in the calculation.

        Args:
            ignore_imag: Only real part is returned if `ignore_imag`.
        """
        df_list = []; app = df_list.append
        for spin in range(self.nsppol):
            for ik, kpoint in enumerate(self.sigma_kpoints):
                app(self.get_dataframe_sk(spin, ik, ignore_imag=ignore_imag))

        return pd.concat(df_list)

    def get_dataframe_sk(self, spin, sig_kpoint, index=None, ignore_imag=False):
        """
        Returns pandas DataFrame with QP results for the given (spin, k-point).

        Args:
            ignore_imag: Only real part is returned if `ignore_imag`.
        """
        ik = self.sigkpt2index(sig_kpoint)
        rows, bands = [], []
        for band in range(self.bstart_sk[spin, ik], self.bstop_sk[spin, ik]):
            bands.append(band)
            ## Build dictionary with the QP results.
            #qpstate = self.reader.read_qp(spin, kpoint, band, ignore_imag=ignore_imag)
            #d = qpstate.as_dict()
            ## Add other entries that may be useful when comparing different calculations.
            #d.update(self.params)
            #rows.append(d)

        #index = len(bands) * [index] if index is not None else bands
        #return pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))
        return pd.DataFrame()

    @add_fig_kwargs
    def plot_dirgaps_t(self, ax=None, **kwargs):
        """
        Plot the KS and the QP(T) direct gaps for all the k-points available on file.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for spin in range(self.nsppol):
            for ik, kpt in enumerate(self.sigma_kpoints):
                # Plot QP_(spin,k)(T)
                ax.plot(self.tmesh, self.qp_dirgaps_t[spin, ik], marker=self.marker_spin[spin],
                        label="QP gap @%s" % repr(kpt))
                # Add KS gaps (assumed at T=0)
                ax.scatter(0, self.ks_dirgaps[spin, ik], label="KS gap @%s" % repr(kpt))

        ax.grid(True)
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel("Direct gap [eV]")
        ax.legend(loc="best")

        return fig

    def write_notebook(self, nbpath=None, title=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=title)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.ebands.plot();"),
            nbv.new_code_cell("ncfile.plot_qpgaps();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class SigEPhRobot(Robot, RobotWithEbands, NotebookWriter):
    """
    This robot analyzes the results contained in multiple SIGEPH.nc files.
    """
    EXT = "SIGEPH"

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.SigEPhRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class SigmaPhReader(ElectronsReader):
    """
    Reads data from file and constructs objects.
    """
    #def read_params(self):
    #    """Read sigeph input parameters. Return OrderedDict"""
    #    od = OrderedDict()
    #    return od

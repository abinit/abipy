# coding: utf-8
"""Objects to analyze the results stored in the GRUNS.nc file produced by anaddb."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import abipy.core.abinit_units as abu
import scipy.constants as const

from collections import OrderedDict
from monty.string import marquee, list_strings
from monty.termcolor import cprint
from monty.collections import AttrDict
from monty.functools import lazy_property
from abipy.core.kpoints import Kpath, IrredZone, KSamplingInfo
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.dfpt.phonons import PhononBands, PhononBandsPlotter
from abipy.iotools import ETSF_Reader
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims
#from abipy.tools import duck
from abipy.core.func1d import Function1D


# DOS name --> meta-data
_ALL_DOS_NAMES = OrderedDict([
    ("gruns_wdos", dict(latex=r"$DOS$")),
    ("gruns_grdos", dict(latex=r"$DOS_{\gamma}$")),
    ("gruns_gr2dos", dict(latex=r"$DOS_{\gamma^2}$")),
    ("gruns_vdos", dict(latex=r"$DOS_v$")),
    ("gruns_v2dos", dict(latex=r"$DOS_{v^2}$")),
])


def get_cv_per_mode(t, w):
    """
    Returns the contribution to the constant-volume specific heat, in eV/K at temperature t for a mode with
    frequency w

    Args:
        t: temperature in K
        w: frequency in eV

    """




class GrunsNcFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    This object provides an interface to the ``GRUNS.nc`` file produced by abinit.
    This file contains Grunesein parameters computed via finite difference.

    Usage example:

    .. code-block:: python

        with GrunsNcFile("foo_GRUNS.nc") as ncfile:
            print(ncfile)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GrunsNcFile
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(GrunsNcFile, self).__init__(filepath)
        self.reader = GrunsReader(filepath)

    def close(self):
        """Close file."""
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        if self.phbands_qpath_vol:
            app(self.phbands_qpath_vol[self.iv0].to_string(with_structure=False, title="Phonon Bands at V0"))
        app("")
        app("Number of Volumes: %d" % self.reader.num_volumes)

        return "\n".join(lines)

    @property
    def structure(self):
        """|Structure| corresponding to the central point V0."""
        return self.reader.structure

    #@property
    #def volumes(self):
    #    """Volumes of the unit cell in Angstrom**3"""
    #    return self.reader.volumes

    @property
    def iv0(self):
        return self.reader.iv0

    @lazy_property
    def doses(self):
        """Dictionary with the phonon doses."""
        return self.reader.read_doses()

    @lazy_property
    def wvols_qibz(self):
        """Phonon frequencies on reagular grid for the different volumes in eV """
        return self.reader.read_value("gruns_wvols_qibz") * abu.Ha_eV

    @lazy_property
    def qibz(self):
        """q-points in the irreducible brillouin zone"""
        return self.reader.read_value("gruns_qibz")

    @lazy_property
    def gvals_qibz(self):
        """Gruneisen parameters in the irreducible brillouin zone"""
        if "gruns_gvals_qibz" not in self.reader.rootgrp.variables:
            raise RuntimeError("GRUNS.nc file does not contain `gruns_gvals_qibz`."
                               "Use prtdos in anaddb input file to compute these values in the IBZ")

        return self.reader.read_value("gruns_gvals_qibz")

    @lazy_property
    def phbands_qpath_vol(self):
        """List of |PhononBands| objects corresponding to the different volumes."""
        return self.reader.read_phbands_on_qpath()

    def to_dataframe(self):
        """
        Return a |pandas-DataFrame| with the following columns:

            ['qidx', 'mode', 'freq', 'qpoint']

        where:

        ==============  ==========================
        Column          Meaning
        ==============  ==========================
        qidx            q-point index.
        mode            phonon branch index.
        grun            Gruneisen parameter.
        groupv          Group velocity.
        freq            Phonon frequency in eV.
        ==============  ==========================
        """

        grun_vals = self.gvals_qibz
        nqibz, natom3 = grun_vals.shape
        phfreqs = self.reader.rootgrp.variables["gruns_wvols_qibz"][:, self.iv0, :] * abu.Ha_eV
        # TODO: units?
        dwdq = self.reader.read_value("gruns_dwdq_qibz")
        groupv = np.linalg.norm(dwdq, axis=-1)

        import pandas as pd
        rows = []
        for iq in range(nqibz):
            for nu in range(natom3):
                rows.append(OrderedDict([
                           ("qidx", iq),
                           ("mode", nu),
                           ("grun", grun_vals[iq, nu]),
                           ("groupv", groupv[iq, nu]),
                           ("freq", phfreqs[iq, nu]),
                           #("qpoint", self.qpoints[iq]),
                        ]))

        return pd.DataFrame(rows, columns=list(rows[0].keys()))

    @add_fig_kwargs
    def plot_doses(self, xlims=None, dos_names="all", with_idos=True, **kwargs):
        r"""
        Plot the different doses stored in the GRUNS.nc file.

        Args:
            xlims: Set the data limits for the x-axis in eV. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            dos_names: List of strings defining the DOSes to plot. Use "all" to plot all DOSes available.
            with_idos: True to display integrated doses.

        Return: |matplotlib-Figure|
        """
        if not self.doses: return None

        dos_names = _ALL_DOS_NAMES.keys() if dos_names == "all" else list_strings(dos_names)
        wmesh = self.doses["wmesh"]

        nrows, ncols = len(dos_names), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        for i, (name, ax) in enumerate(zip(dos_names, ax_list)):
            dos, idos = self.doses[name][0], self.doses[name][1]
            ax.plot(wmesh, dos, color="k")
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            ax.set_ylabel(_ALL_DOS_NAMES[name]["latex"])
            #ax.yaxis.set_ticks_position("right")

            if with_idos:
                other_ax = ax.twinx()
                other_ax.plot(wmesh, idos, color="k")
                other_ax.set_ylabel(_ALL_DOS_NAMES[name]["latex"].replace("DOS", "IDOS"))

            if i == len(dos_names) - 1:
                ax.set_xlabel(r"$\omega$ [eV]")
            #ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def get_plotter(self):
        """
        Return an instance of |PhononBandsPlotter| that can be use to plot
        multiple phonon bands or animate the bands
        """
        plotter = PhononBandsPlotter()
        for iv, phbands in enumerate(self.phbands_qpath_vol):
            label = "V=%.2f $A^3$ " % phbands.structure.volume
            plotter.add_phbands(label, phbands)

        return plotter

    @add_fig_kwargs
    def plot_phbands_with_gruns(self, fill_with="gruns", gamma_fact=1, alpha=0.6, with_doses="all", units="eV",
                                ylims=None, match_bands=False, **kwargs):
        """
        Plot the phonon bands corresponding to ``V0`` (the central point) with markers
        showing the value and the sign of the Grunesein parameters.

        Args:
            fill_with: Define the quantity used to plot stripes. "gruns" for Grunesein parameters,
                "groupv" for phonon group velocities.
            gamma_fact: Scaling factor for Grunesein parameters.
                Up triangle for positive values, down triangles for negative values.
            alpha: The alpha blending value for the markers between 0 (transparent) and 1 (opaque).
            with_doses: "all" to plot all DOSes available, ``None`` to disable DOS plotting,
                else list of strings with the name of the DOSes to plot.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            ylims: Set the data limits for the y-axis in eV. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            match_bands: if True tries to follow the band along the path based on the scalar product of the eigenvectors.

        Returns: |matplotlib-Figure|.
        """
        if not self.phbands_qpath_vol: return None
        phbands = self.phbands_qpath_vol[self.iv0]
        factor = phbands.phfactor_ev2units(units)
        gamma_fact *= factor

        # Build axes (ax_bands and ax_doses)
        if with_doses is None:
            ax_bands, fig, plt = get_ax_fig_plt(ax=None)
        else:
            import matplotlib.pyplot as plt
            from matplotlib.gridspec import GridSpec
            dos_names = list(_ALL_DOS_NAMES.keys()) if with_doses == "all" else list_strings(with_doses)
            ncols = 1 + len(dos_names)
            fig = plt.figure()
            width_ratios = [2] + len(dos_names) * [0.2]
            gspec = GridSpec(1, ncols, width_ratios=width_ratios, wspace=0.05)
            ax_bands = plt.subplot(gspec[0])
            ax_doses = []
            for i in range(len(dos_names)):
                ax_doses.append(plt.subplot(gspec[i + 1], sharey=ax_bands))

        # Plot phonon bands.
        phbands.plot(ax=ax_bands, units=units, match_bands=match_bands, show=False)

        if fill_with == "gruns":
            max_gamma = np.abs(phbands.grun_vals).max()
        elif fill_with == "groupv":
            # TODO: units?
            dwdq_qpath = self.reader.read_value("gruns_dwdq_qpath")
            groupv = np.linalg.norm(dwdq_qpath, axis=-1)
            max_gamma = np.abs(groupv).max()
        else:
            raise ValueError("Unsupported fill_with: `%s`" % fill_with)

        # Plot gruneisen markers on top of band structure.
        xvals = np.arange(len(phbands.phfreqs))
        max_omega = np.abs(phbands.phfreqs).max()

        for nu in phbands.branches:
            omegas = phbands.phfreqs[:, nu].copy() * factor

            if fill_with == "gruns":
                # Must handle positive-negative values
                sizes = phbands.grun_vals[:, nu].copy() * (gamma_fact * 0.02 * max_omega / max_gamma)
                yup = omegas + np.where(sizes >= 0, sizes, 0)
                ydown = omegas + np.where(sizes < 0, sizes, 0)
                ax_bands.fill_between(xvals, omegas, yup, alpha=alpha, facecolor="red")
                ax_bands.fill_between(xvals, ydown, omegas, alpha=alpha, facecolor="blue")

            elif fill_with == "groupv":
                sizes = groupv[:, nu].copy() * (gamma_fact * 0.04 * max_omega / max_gamma)
                ydown, yup = omegas - sizes / 2, omegas + sizes / 2
                ax_bands.fill_between(xvals, ydown, yup, alpha=alpha, facecolor="red")

        set_axlims(ax_bands, ylims, "x")

        if with_doses is None:
            return fig

        # Plot Doses.
        wmesh = self.doses["wmesh"] * factor
        for i, (name, ax) in enumerate(zip(dos_names, ax_doses)):
            dos, idos = self.doses[name][0], self.doses[name][1]
            ax.plot(dos, wmesh, label=name, color="k")
            set_axlims(ax, ylims, "x")
            ax.grid(True)
            ax.set_ylabel("")
            ax.tick_params(labelbottom='off')
            if i == len(dos_names) - 1:
                ax.yaxis.set_ticks_position("right")
            else:
                ax.tick_params(labelleft='off')
            ax.set_title(_ALL_DOS_NAMES[name]["latex"])

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.structure"),

            nbv.new_code_cell("ncfile.plot_doses();"),
            nbv.new_code_cell("ncfile.plot_phbands_with_gruns();"),

            #nbv.new_code_cell("phbands_qpath_v0.plot_fatbands(phdos_file=phdosfile);"),
            nbv.new_code_cell("plotter = ncfile.get_plotter()\nprint(plotter)"),
            nbv.new_code_cell("df_phbands = plotter.get_phbands_frame()\ndisplay(df_phbands)"),
            nbv.new_code_cell("plotter.ipw_select_plot()"),

            nbv.new_code_cell("gdata = ncfile.to_dataframe()\ngdata.describe()"),
            nbv.new_code_cell("""\
#import df_widgets.seabornw as snsw
#snsw.api_selector(gdata)"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class GrunsReader(ETSF_Reader):
    """
    This object reads the results stored in the GRUNS.nc file produced by ABINIT.
    It provides helper functions to access the most important quantities.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GrunsReader
    """
    # Fortran arrays (remember to transpose dimensions!)
    # Remember: Atomic units are used everywhere in this file.
    #nctkarr_t("gruns_qptrlatt", "int", "three, three"), &
    #nctkarr_t("gruns_shiftq", "dp", "three, gruns_nshiftq"), &
    #nctkarr_t("gruns_qibz", "dp", "three, gruns_nqibz"), &
    #nctkarr_t("gruns_wtq", "dp", "gruns_nqibz"), &
    #nctkarr_t("gruns_gvals_qibz", "dp", "number_of_phonon_modes, gruns_nqibz"), &
    #nctkarr_t("gruns_wvols_qibz", "dp", "number_of_phonon_modes, gruns_nvols, gruns_nqibz"), &
    #nctkarr_t("gruns_dwdq_qibz", "dp", "three, number_of_phonon_modes, gruns_nqibz"), &
    #nctkarr_t("gruns_omega_mesh", "dp", "gruns_nomega"), &
    #nctkarr_t("gruns_wdos", "dp", "gruns_nomega, two"), &
    #nctkarr_t("gruns_grdos", "dp", "gruns_nomega, two"), &
    #nctkarr_t("gruns_gr2dos", "dp", "gruns_nomega, two"), &
    #nctkarr_t("gruns_v2dos", "dp", "gruns_nomega, two"), &
    #nctkarr_t("gruns_vdos", "dp", "gruns_nomega, two") &
    #nctkarr_t("gruns_qpath", "dp", "three, gruns_nqpath")
    #nctkarr_t("gruns_gvals_qpath", "dp", "number_of_phonon_modes, gruns_nqpath")
    #nctkarr_t("gruns_wvols_qpath", "dp", "number_of_phonon_modes, gruns_nvols, gruns_nqpath")
    #nctkarr_t("gruns_dwdq_qpath", "dp", "three, number_of_phonon_modes, gruns_nqpath")
    #nctkarr_t("gruns_rprimd", "dp", "three, three, gruns_nvols"), &
    #nctkarr_t("gruns_xred", "dp", "three, number_of_atoms, gruns_nvols") &

    def __init__(self, filepath):
        super(GrunsReader, self).__init__(filepath)

        # Read and store important quantities.
        self.structure = self.read_structure()

        self.num_volumes = self.read_dimvalue("gruns_nvols")
        # The index of the volume used for the finite difference.
        self.iv0 = self.read_value("gruns_iv0") - 1  #  F --> C

    def read_doses(self):
        """
        Return a |AttrDict| with the DOSes available in the file. Empty dict if
        DOSes are not available.
        """
        if "gruns_nomega" not in self.rootgrp.dimensions:
            cprint("File `%s` does not contain ph-DOSes, returning empty dict" % self.path, "yellow")
            return {}

        # Read q-point sampling used to compute DOSes.
        qptrlatt = self.read_value("gruns_qptrlatt")
        shifts = self.read_value("gruns_shiftq")
        qsampling = KSamplingInfo.from_kptrlatt(qptrlatt, shifts, kptopt=1)

        frac_coords_ibz = self.read_value("gruns_qibz")
        weights = self.read_value("gruns_wtq")
        qpoints = IrredZone(self.structure.reciprocal_lattice, frac_coords_ibz,
                            weights=weights, names=None, ksampling=qsampling)

        # DOSes are in 1/Hartree.
        d = AttrDict(wmesh=self.read_value("gruns_omega_mesh") * abu.Ha_eV, qpoints=qpoints)

        for dos_name in _ALL_DOS_NAMES:
            dos_idos = self.read_value(dos_name)
            dos_idos[0] *= abu.eV_Ha  # Here we convert to eV. IDOS are not changed.
            d[dos_name] = dos_idos

        return d

    def read_phbands_on_qpath(self):
        """
        Return a list of |PhononBands| computed at the different volumes.
        The ``iv0`` entry corresponds to the central point used to compute Grunesein parameters
        with finite differences. This object stores the Grunesein parameters in ``grun_vals``.
        """
        if "gruns_qpath" not in self.rootgrp.variables:
            cprint("File `%s` does not contain phonon bands, returning empty list." % self.path, "yellow")
            return []

        qfrac_coords = self.read_value("gruns_qpath")
        grun_vals = self.read_value("gruns_gvals_qpath")
        freqs_vol = self.read_value("gruns_wvols_qpath") * abu.Ha_eV
        # TODO: Convert?
        dwdq_qpath = self.read_value("gruns_dwdq_qpath")

        amuz = self.read_amuz_dict()
        #print("amuz", amuz)

        # nctkarr_t("gruns_phdispl_cart_qpath", "dp", &
        # "two, number_of_phonon_modes, number_of_phonon_modes, gruns_nvols, gruns_nqpath") &
        # NOTE: in GRUNS the displacements are in Bohr. here we convert to Ang to be
        # consistent with the PhononBands API.
        phdispl_cart_qptsvol = self.read_value("gruns_phdispl_cart_qpath", cmode="c") * abu.Bohr_Ang

        lattices = self.read_value("gruns_rprimd") * abu.Bohr_Ang #, "dp", "three, three, gruns_nvols")
        gruns_xred = self.read_value("gruns_xred")                #, "dp", "three, number_of_atoms, gruns_nvols")

        phbands_qpath_vol = []
        for ivol in range(self.num_volumes):
            # TODO structure depends on:
            # volumes
            # non_anal_ph
            # amu: DONE
            # phdispl_cart DONE
            if ivol == self.iv0:
                structure = self.structure
            else:
                structure = self.structure.__class__(lattices[ivol], self.structure.species, gruns_xred[ivol])

            qpoints = Kpath(structure.reciprocal_lattice, qfrac_coords)
            phdispl_cart = phdispl_cart_qptsvol[:, ivol].copy()
            phb = PhononBands(structure, qpoints, freqs_vol[:, ivol], phdispl_cart, non_anal_ph=None, amu=amuz)
            # Add Grunesein parameters.
            if ivol == self.iv0: phb.grun_vals = grun_vals
            phbands_qpath_vol.append(phb)

        return phbands_qpath_vol

    def read_amuz_dict(self):
        """
        Dictionary that associates the atomic number to the values of the atomic
        mass units used for the calculation.
        """
        amu_typat = self.read_value("atomic_mass_units")
        znucl_typat = self.read_value("atomic_numbers")

        amuz = {}
        for symbol in self.chemical_symbols:
            type_idx = self.typeidx_from_symbol(symbol)
            amuz[znucl_typat[type_idx]] = amu_typat[type_idx]

        return amuz

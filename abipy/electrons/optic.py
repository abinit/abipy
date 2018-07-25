# coding: utf-8
"""
Objects to read and analyze optical properties stored in the optic.nc file produced by optic executable.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import abipy.core.abinit_units as abu

from collections import OrderedDict
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, data_from_cplx_mode
from abipy.abio.robots import Robot
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands

ALL_CHIS = OrderedDict([
    ("linopt", {
        "longname": "Dielectric function",
        "rank": 2,
        #"terms":
        #"latex": r"\chi(\omega)"
        }
    ),
    ("shg", {
        "longname": "Second Harmonic Generation",
        "rank": 3,
        "terms": ["shg_inter2w", "shg_inter1w", "shg_intra2w",
                  "shg_intra1w", "shg_intra1wS", "shg_chi2tot"],
        }
        #"latex": r"\chi(-2\omega, \omega, \omega)"
    ),
    ("leo", {
        "longname": "Linear Electro-optic",
        "rank": 3,
        "terms": ["leo_chi", "leo_eta", "leo_sigma", "leo_chi2tot"],
        }
        #"latex": r"\chi(-\omega, \omega, 0)"
    )
])

#LEO2_TERMS = OrderedDict([
#    ("leo2_chiw", None),
#    ("leo2_etaw", None),
#    ("leo2_chi2w", None),
#    ("leo2_eta2w", None),
#    ("leo2_sigmaw", None),
#    ("leo2_chi2tot", None),
#])


#####################################
# Helper functions for linear optic #
#####################################

def reflectivity(eps):
    """Reflectivity(w) from vacuum, at normal incidence"""
    return np.sqrt(0.5 * (np.abs(eps) + eps.real))


#def abs_coeff(eps):
#    """absorption coeff (in m-1) = omega Im(eps) / c n(eps)"""
#    if (abs(eps(iw)) + dble(eps(iw)) > zero) then
#       tmpabs = aimag(eps(iw))*ene / sqrt(half*( abs(eps(iw)) + dble(eps(iw)) )) / Sp_Lt / Bohr_meter
#    end if


def kappa(eps):
    """Im(refractive index(w)) aka kappa"""
    return np.sqrt(0.5 * (np.abs(eps) - eps.real))


def n(eps):
    """Re(refractive index(w)) aka n"""
    return np.sqrt(0.5 * (np.abs(eps) + eps.real))


#def eels(eps):
#    np.where(np.abs(eps)
#    return - (1 / eps).imag


LINEPS_WHAT2EFUNC = dict(
    n=n,
    reflectivity=reflectivity,
    kappa=kappa,
    re=lambda eps: eps.real,
    im=lambda eps: eps.imag,
    #abs: lambda: eps: np.abs(eps),
    #angle: lambda: eps: np.angle(eps, deg=False),
    #abs_coeff=abs_coeff
    #eels=lambda: eps /
)


class OpticNcFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the results produced by optic. Provides methods to plot optical
    properties and susceptibilty tensors. Used by :class:`OpticRobot` to analyze multiple files.

    Usage example:

    .. code-block:: python

        with OpticNcFile("out_optic.nc") as optic:
            optic.ebands.plot()
            optic.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: OpticNcFile
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(OpticNcFile, self).__init__(filepath)
        self.reader = OpticReader(filepath)

        # Read optic input variables and info on k-point sampling and store them in self.
        keys = [
            "kptopt",
            # optic input variables
            "broadening", "maxomega", "domega", "scissor", "tolerance", "do_antiresonant", "do_ep_renorm",
        ]
        for key in keys:
            setattr(self, key, self.reader.read_value(key))

    @lazy_property
    def wmesh(self):
        """
        Frequency mesh in eV. Note that the same frequency-mesh is used
        for the different optical properties.
        """
        return self.reader.read_value("wmesh")

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
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))

        app(marquee("Optic calculation", mark="="))
        # Show Optic variables.
        app("broadening: %s [Ha], %.3f (eV)" % (self.broadening, self.broadening * abu.Ha_eV))
        app("scissor: %s [Ha], %.3f (eV)" % (self.scissor, self.scissor * abu.Ha_eV))
        app("tolerance: %s [Ha], %.3f (eV)" % (self.tolerance, self.tolerance * abu.Ha_eV))
        app("maxomega: %s [Ha], %.3f (eV)" % (self.maxomega, self.maxomega * abu.Ha_eV))
        app("domega: %s [Ha], %.3f (eV)" % (self.domega, self.domega * abu.Ha_eV))
        app("do_antiresonant %s, do_ep_renorm %s" % (self.do_antiresonant, self.do_ep_renorm))
        app("Number of temperatures: %d" % self.reader.ntemp)

        # Show available tensors and computed components.
        for key, info in ALL_CHIS.items():
            if not self.reader.computed_components[key]: continue
            app("%s components computed: %s" % (
                info["longname"], ", ".join(self.reader.computed_components[key])))

        if verbose > 1:
            app(marquee("Abinit Header", mark="="))
            app(self.hdr.to_string(verbose=verbose))

        return "\n".join(lines)

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def has_linopt(self):
        """True if the ncfile contains Second Harmonic Generation tensor."""
        return "linopt" in self.reader.computed_components

    @lazy_property
    def has_shg(self):
        """True if the ncfile contains Second Harmonic Generation tensor."""
        return "shg" in self.reader.computed_components

    @lazy_property
    def has_leo(self):
        """True if the ncfile contains the Linear Electro-optic tensor"""
        return "leo" in self.reader.computed_components

    #@lazy_property
    #def xc(self):
    #    """:class:`XcFunc object with info on the exchange-correlation functional."""
    #    return self.reader.read_abinit_xcfunc()

    def close(self):
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @staticmethod
    def get_linopt_latex_label(what, comp):
        """
        Return latex label for linear-optic quantities. Used in plots.
        """
        return dict(
            n=r"$n_{%s}$" % comp,
            reflectivity=r"$R_{%s}$" % comp,
            kappa=r"$\kappa_{%s}$" % comp,
            re=r"$\Re(\epsilon_{%s})$" % comp,
            im=r"$\Im(\epsilon_{%s})$" % comp,
            #abs=r"$|\epsilon_{%s}|$" % comp,
            #abs_coeff=abs_coeff_{%s}} % comp,
            #eels:r"EELS_{%s}" % comp,
        )[what]

    def get_chi2_latex_label(self, key, what, comp):
        """
        Build latex label for chi2-related quantities. Used in plots.
        """
        symb = "{%s}" % key.capitalize()
        return dict(
            re=r"$\Re(%s_{%s})$" % (symb, comp),
            im=r"$\Im(%s_{%s})$" % (symb, comp),
            abs=r"$|%s_{%s}|$" % (symb, comp),
        )[what]

    @add_fig_kwargs
    def plot_linear_epsilon(self, components="all", what="im", itemp=0,
                            ax=None, xlims=None, with_xlabel=True, label=None, fontsize=12, **kwargs):
        """
        Plot components of the linear dielectric function.

        Args:
            components: List of cartesian tensor components to plot e.g. ["xx", "xy"].
                "all" if all components available on file should be plotted on the same ax.
            what: quantity to plot. "re" for real part, "im" for imaginary. Accepts also "abs", "angle".
            itemp: Temperature index.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            with_xlabel: True if x-label should be added.
            label: True to add legend label to each curve.
            fontsize: Legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        comp2eps = self.reader.read_lineps(components, itemp=itemp)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for comp, eps in comp2eps.items():
            values = LINEPS_WHAT2EFUNC[what](eps)
            # Note: I'm skipping the first point at w=0 because optic does not compute it!
            # The same trick is used in the other plots.
            ax.plot(self.wmesh[1:], values[1:],
                    label=self.get_linopt_latex_label(what, comp) if label is None else label)

        ax.grid(True)
        if with_xlabel: ax.set_xlabel('Photon Energy (eV)')
        set_axlims(ax, xlims, "x")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_linopt(self, select="all", itemp=0, xlims=None, **kwargs):
        """
        Subplots with all linear optic quantities selected by ``select`` at temperature ``itemp``.

        Args:
            select:
            itemp: Temperature index.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used.

        Returns: |matplotlib-Figure|
        """
        key = "linopt"
        if not self.reader.computed_components[key]: return None
        if select == "all": select = list(LINEPS_WHAT2EFUNC.keys())
        select = list_strings(select)

        nrows, ncols = len(select), 1
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=True)

        components = self.reader.computed_components[key]
        for i, (what, ax) in enumerate(zip(select, ax_mat)):
            self.plot_linear_epsilon(what=what, itemp=itemp, components=components,
                                     ax=ax, xlims=xlims, with_xlabel=(i == len(select) - 1),
                                     show=False)
        return fig

    @add_fig_kwargs
    def plot_chi2(self, key, components="all", what="abs", itemp=0, decompose=False,
                  ax=None, xlims=None, with_xlabel=True, label=None, fontsize=12, **kwargs):
        """
        Low-level function to plot chi2 tensor.

        Args:
            key:
            components: List of cartesian tensor components to plot e.g. ["xxx", "xyz"].
                "all" if all components available on file should be plotted on the same ax.
            what: quantity to plot. "re" for real part, "im" for imaginary, Accepts also "abs", "angle".
            itemp: Temperature index.
            decompose: True to plot individual contributions.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            with_xlabel: True to add x-label.
            label: True to add legend label to each curve.
            fontsize: Legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        if not self.reader.computed_components[key]: return None
        comp2terms = self.reader.read_tensor3_terms(key, components, itemp=itemp)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for comp, terms in comp2terms.items():
            for name, values in terms.items():
                if not decompose and not name.endswith("tot"): continue
                values = data_from_cplx_mode(what, values)
                ax.plot(self.wmesh[1:], values[1:],
                        label=self.get_chi2_latex_label(key, what, comp) if label is None else label,
                )

        ax.grid(True)
        if with_xlabel: ax.set_xlabel('Photon Energy (eV)')
        set_axlims(ax, xlims, "x")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_shg(self, **kwargs):
        """Plot Second Harmonic Generation. See plot_chi2 for args supported."""
        return self.plot_chi2(key="shg", show=False, **kwargs)

    @add_fig_kwargs
    def plot_leo(self, **kwargs):
        """Plot Linear Electro-optic. See plot_chi2 for args supported."""
        return self.plot_chi2(key="leo", show=False, **kwargs)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        if self.has_linopt:
            yield self.plot_linear_epsilon(what="re", show=False)
            yield self.plot_linear_epsilon(what="im", show=False)
            yield self.plot_linopt(show=False)
        if self.has_shg:
            yield self.plot_shg(show=False)
        if self.has_leo:
            yield self.plot_leo(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("optic = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(optic)"),
            nbv.new_code_cell("optic.ebands.plot();"),
        ])

        # Add plot calls if quantities have been computed.
        for key, info in ALL_CHIS.items():
            if not self.reader.computed_components[key]: continue
            pycall = "optic.plot_%s();" % key
            nb.cells.extend([
                nbv.new_code_cell(pycall),
            ])

        return self._write_nb_nbpath(nb, nbpath)


class OpticReader(ElectronsReader):
    """
    This object reads the results stored in the optic.nc file
    It provides helper function to access the most important quantities.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: OpticReader
    """
    def __init__(self, filepath):
        super(OpticReader, self).__init__(filepath)
        self.ntemp = self.read_dimvalue("ntemp")

        self.computed_components = OrderedDict()
        self.computed_ids = OrderedDict()
        for chiname, info in ALL_CHIS.items():
            comp_name = chiname + "_components"
            if comp_name not in self.rootgrp.variables:
                fi_comps, ids = [], []
            else:
                fi_comps = [str(i) for i in self.read_value(comp_name)]
                if info["rank"] == 2:
                    ids = [(int(i[0])-1, int(i[1])-1) for i in fi_comps]
                elif info["rank"] == 3:
                    ids = [(int(i[0])-1, int(i[1])-1, int(i[2])-1) for i in fi_comps]
                else:
                    raise NotImplementedError("rank %s" % info["rank"])

            self.computed_ids[chiname] = ids
            self.computed_components[chiname] = [abu.itup2s(it) for it in ids]

    def read_lineps(self, components, itemp=0):
        """
        Args:
            components: List of cartesian tensor components to plot e.g. ["xx", "xy"].
                "all" if all components available on file should be plotted on the same ax.
            itemp: Temperature index.
        """
        # linopt_epsilon has *Fortran* shape [two, nomega, num_comp, ntemp]
        key = "linopt"
        if components == "all": components = self.computed_components[key]
        if not (self.ntemp > itemp >= 0):
            raise ValueError("Invalid itemp: %s, ntemp: %s" % (itemp, self.ntemp))

        var = self.read_variable("linopt_epsilon")
        od = OrderedDict()
        for comp in list_strings(components):
            try:
                ijp = self.computed_components[key].index(comp)
            except ValueError:
                raise ValueError("epsilon_component %s was not computed" % comp)

            values = var[itemp, ijp]
            od[comp] = values[:, 0] + 1j * values[:, 1]
        return od

    def read_tensor3_terms(self, key, components, itemp=0):
        """
        Args:
            key: Name of the netcdf variable to read.
            components: List of cartesian tensor components to plot e.g. ["xxx", "xyz"].
                "all" if all components available on file should be plotted on the same ax.
            itemp: Temperature index.

        Return:
            :class:`OrderedDict` mapping cartesian components e.g. "xyz" to data dictionary.
            Individual entries are listed in ALL_CHIS[key]["terms"]
        """
        # arrays have Fortran shape [two, nomega, num_comp, ntemp]
        if components == "all": components = self.computed_components[key]
        components = list_strings(components)
        if not (self.ntemp > itemp >= 0):
            raise ValueError("Invalid itemp: %s, ntemp: %s" % (itemp, self.ntemp))

        od = OrderedDict([(comp, OrderedDict()) for comp in components])
        for chiname in ALL_CHIS[key]["terms"]:
            #print("About to read:", chiname)
            var = self.read_variable(chiname)
            for comp in components:
                try:
                    ijkp = self.computed_components[key].index(comp)
                except ValueError:
                    raise ValueError("%s component %s was not computed" % (key, comp))
                values = var[itemp, ijkp]
                od[comp][chiname] = values[:, 0] + 1j * values[:, 1]
        return od


class OpticRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple optic.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: OpticRobot
    """
    EXT = "OPTIC"

    @lazy_property
    def computed_components_intersection(self):
        """
        Dictionary with the list of cartesian tensor components
        available in each file. Use keys from ALL_CHIS.
        """
        od = OrderedDict()
        for ncfile in self.abifiles:
            for chiname in ALL_CHIS:
                comps = ncfile.reader.computed_components[chiname]
                if chiname not in od:
                    od[chiname] = comps
                else:
                    # Build intersection while preserving order.
                    od[chiname] = self.ordered_intersection(od[chiname], comps)
        return od

    @add_fig_kwargs
    def plot_linopt_convergence(self, components="all", what_list=("re", "im"),
                                sortby="nkpt", itemp=0, xlims=None, **kwargs):
        """
        Plot the convergence of the dielectric function tensor with respect to
        parameter defined by ``sortby``.

        Args:
            components: List of cartesian tensor components to plot e.g. ["xx", "xy"].
                "all" if all components available on file should be plotted on the same ax.
            what_list: List of quantities to plot. "re" for real part, "im" for imaginary.
                Accepts also "abs", "angle".
            sortby: Define the convergence parameter, sort files and produce plot labels. Can be None, string or function.
                If None, no sorting is performed.
                If string, it's assumed that the ncfile has an attribute with the same name and getattr is invoked.
                If callable, the output of callable(ncfile) is used.
            itemp: Temperature index.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.

        Returns: |matplotlib-Figure|
        """
        # Build grid plot: computed tensors along the rows, what_list along columns.
        key = "linopt"
        components = self.computed_components_intersection[key]

        nrows, ncols = len(components), len(what_list)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        label_ncfile_param = self.sortby(sortby)
        for i, comp in enumerate(components):
            for j, what in enumerate(what_list):
                ax = ax_mat[i, j]
                for ifile, (label, ncfile, param) in enumerate(label_ncfile_param):

                    ncfile.plot_linear_epsilon(components=comp, what=what, itemp=itemp, ax=ax,
                        xlims=xlims, with_xlabel=(i == len(components) - 1),
                        label="%s %s" % (sortby, param) if not callable(sortby) else str(param),
                        show=False)

                    if ifile == 0:
                        ax.set_title(ncfile.get_linopt_latex_label(what, comp))

                if (i, j) != (0, 0):
                    ax.legend().set_visible(False)

        return fig

    @add_fig_kwargs
    def plot_shg_convergence(self, **kwargs):
        """Plot Second Harmonic Generation. See plot_convergence_rank3 for args supported."""
        if "what_list" not in kwargs: kwargs["what_list"] = ["abs"]
        return self.plot_convergence_rank3(key="shg", **kwargs)

    @add_fig_kwargs
    def plot_leo_convergence(self, **kwargs):
        """Plot Linear electron-optic. See plot_convergence_rank3 for args supported."""
        if "what_list" not in kwargs: kwargs["what_list"] = ["abs"]
        return self.plot_convergence_rank3(key="leo", **kwargs)

    @add_fig_kwargs
    def plot_convergence_rank3(self, key, components="all", itemp=0, what_list=("abs",),
                               sortby="nkpt", decompose=False, xlims=None, **kwargs):
        """
        Plot convergence of arbitrary rank3 tensor. This is a low-level routine used in other plot methods.

        Args:
            key: Name of the quantity to analyze.
            components: List of cartesian tensor components to plot e.g. ["xxx", "xyz"].
                "all" if all components available on file should be plotted on the same ax.
            itemp: Temperature index.
            what_list: List of quantities to plot. "re" for real part, "im" for imaginary.
                Accepts also "abs", "angle".
            sortby: Define the convergence parameter, sort files and produce plot labels. Can be None, string or function.
                If None, no sorting is performed.
                If string, it's assumed that the ncfile has an attribute with the same name and ``getattr`` is invoked.
                If callable, the output of callable(ncfile) is used.
            decompose: True to plot individual contributions.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used.

        Returns: |matplotlib-Figure|
        """
        # Build grid plot: computed tensors along the rows, what_list along columns.
        components = self.computed_components_intersection[key]

        nrows, ncols = len(components), len(what_list)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        label_ncfile_param = self.sortby(sortby)
        for i, comp in enumerate(components):
            for j, what in enumerate(what_list):
                ax = ax_mat[i, j]
                for ifile, (label, ncfile, param) in enumerate(label_ncfile_param):

                    ncfile.plot_chi2(key=key, components=comp, what=what, itemp=itemp, decompose=decompose,
                        ax=ax, xlims=xlims, with_xlabel=(i == len(components) - 1),
                        label="%s %s" % (sortby, param) if not callable(sortby) else str(param),
                        show=False, **kwargs)

                    if ifile == 0:
                        ax.set_title(ncfile.get_chi2_latex_label(key, what, comp))

                if (i, j) != (0, 0):
                    ax.legend().set_visible(False)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        for key, comps in self.computed_components_intersection.items():
            if not comps: continue
            plot_fig = getattr(self, "plot_%s_convergence" % key)
            yield plot_fig(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.OpticRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
        ])

        for key, comps in self.computed_components_intersection.items():
            if not comps: continue
            pycall = "robot.plot_%s_convergence();" % key
            nb.cells.extend([
                nbv.new_code_cell(pycall),
            ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

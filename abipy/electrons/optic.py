# coding: utf-8
"""
Objects to read and analyze optical properties stored in the optic.nc file produced by optic executable.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons.ebands import ElectronsReader
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims, data_from_cplx_mode


def s2itup(comp):
    """
    Convert string in the form `xx`, `xyz` into tuple of two (three) indices
    that can be used to slice susceptibility tensors (numpy array).

    >>> assert s2itup("yy") == (0, 1)
    >>> assert s2itup("xyz") == (0, 1, 2)
    """
    d = {"x": 0, "y": 1, "z": 2}
    comp = str(comp).strip()
    if len(comp) == 2:
        return d[comp[0]], d[comp[1]]
    elif len(comp) == 3:
        return d[comp[0]], d[comp[1]], d[comp[2]]
    else:
        raise ValueError("Expecting component in the form `xy` or `xyz` but got `%s`" % comp)


def itup2s(t):
    """
    Convert tuple of 2 (3) integers into string in the form `xx` (`xyz`).
    Assume C-indexing e.g. 0 --> x

    >>> assert itup2s((0, 1)) == "yy"
    >>> assert itup2s((0, 1, 2)) == "xyz"
    """
    if not isinstance(t, tuple) and len(t) not in (2, 3):
        raise TypeError("Expecting tuple of len 2 or 3, got %s" % str(t))
    d = {0: "x", 1: "y", 2: "z"}
    return "".join(d[i] for i in t)

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


LO_WHAT2EFUNC = dict(
    n=n,
    reflectivity=reflectivity,
    kappa=kappa,
    re=lambda arr: arr.real,
    im=lambda arr: arr.imag,
    #abs: lambda: arr: np.abs(arr),
    #angle: lambda: arr: np.angle(arr, deg=False),
    #abs_coeff=abs_coeff
)


class OpticNcFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the results produced by optic. Provides methods to plot optical
    properties and susceptibilty tensors. Used by OpticRobot to analyze multiple file.

    Usage example:

    .. code-block:: python

        with OpticNcFile("out_optic.nc") as optic:
            optic.ebands.plot()
            optic.plot()
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(OpticNcFile, self).__init__(filepath)
        self.reader = OpticReader(filepath)

        # Read optic input variables, info on k-point sampling and store them in self.
        keys = [
            "kptopt",
            # optic input variables
            "broadening", "maxomega", "domega", "scissor", "tolerance", "do_antiresonant", "do_ep_renorm",
        ]
        for key in keys:
            setattr(self, key, self.reader.read_value(key))

        # Number of temperatures.
        self.ntemp = self.reader.ntemp
        # Frequency mesh in eV.
        self.wmesh = self.reader.wmesh

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(str(self.structure))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))

        app(marquee("Optic calculation", mark="="))
        # Show Optic variables.
        app("broadening: %s [Ha], %.3f [eV]" % (self.broadening, self.broadening * units.Ha_to_eV))
        app("scissor: %s [Ha], %.3f [eV]" % (self.scissor, self.scissor * units.Ha_to_eV))
        app("tolerance: %s [Ha], %.3f [eV]" % (self.tolerance, self.tolerance * units.Ha_to_eV))
        app("maxomega: %s [Ha], %.3f [eV]" % (self.maxomega, self.maxomega * units.Ha_to_eV))
        app("domega: %s [Ha], %.3f [eV]" % (self.domega, self.domega * units.Ha_to_eV))
        app("do_antiresonant %s, do_ep_renorm %s" % (self.do_antiresonant, self.do_ep_renorm))
        app("Number of temperatures: %d" % self.ntemp)

        # Show available quantities and computed components.
        for key, info in self.reader.ALL_CHIS.items():
            if not self.reader.computed_components[key]: continue
            app("%s components computed: %s" % (
                info["fullname"], ", ".join(self.reader.computed_components[key])))

        if verbose > 1:
            app(marquee("Abinit Header", mark="="))
            app(self.hdr.to_string(verbose=verbose))

        return "\n".join(lines)

    @lazy_property
    def hdr(self):
        """:class:`AttrDict` with the Abinit header e.g. hdr.ecut."""
        return self.reader.read_abinit_hdr()

    @lazy_property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    #@lazy_property
    #def xc(self):
    #    """:class:`XcFunc object with info on the exchange-correlation functional."""
    #    return self.reader.read_abinit_xcfunc()

    def close(self):
        """Close the file."""
        self.reader.close()

    @staticmethod
    def get_linopt_latex_label(what, comp):
        """
        Build latex label for linear-optic quantities. Used in plots.
        """
        return dict(
            n=r"$n_{%s}$" % comp,
            reflectivity=r"$R_{%s}$" % comp,
            kappa=r"$\kappa_{%s}$" % comp,
            re=r"$\Re(\epsilon_{%s})$" % comp,
            im=r"$\Im(\epsilon_{%s})$" % comp,
            #abs=r"$|\epsilon_{%s}|$" % comp,
            #abs_coeff=abs_coeff
        )[what]

    def get_chi2_latex_label(self, key, what, comp):
        symb = "{%s}" % key.capitalize()
        return dict(
            re=r"$\Re(%s_{%s})$" % (symb, comp),
            im=r"$\Im(%s_{%s})$" % (symb, comp),
            abs=r"$|%s_{%s}|$" % (symb, comp),
        )[what]

    @add_fig_kwargs
    def plot_linear_epsilon(self, components="all", what="im", itemp=0,
                            ax=None, xlims=None, with_xlabel=True, label=None, **kwargs):
        """
        Plot linear dielectric function tensor components.

        Args:
            components:
            what:
            itemp: Temperature index.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used.
            with_xlabel: True if x-label should be added.

        Returns:
            `matplotlib` figure
        """
        comp2eps = self.reader.read_epsilon(components, itemp=itemp)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for comp, eps in comp2eps.items():
             values = LO_WHAT2EFUNC[what](eps)
             # Note: I'm skipping the first point at w=0 because optic does not compute it!
             # The same trick is used in the other plots.
             ax.plot(self.wmesh[1:], values[1:],
                     label=self.get_linopt_latex_label(what, comp) if label is None else label)

        ax.grid(True)
        if with_xlabel: ax.set_xlabel('Photon Energy [eV]')
        set_axlims(ax, xlims, "x")
        ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_linopt(self, select="all", itemp=0, xlims=None, **kwargs):
        """
        Subplots with all linear optic quantities selected by `select` at temperature `itemp`

        Args:
            select:
            itemp: Temperature index.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used.

        Returns:
            `matplotlib` figure
        """
        key = "linopt"
        if not self.reader.computed_components[key]: return None
        if select == "all": select = list(LO_WHAT2EFUNC.keys())
        select = list_strings(select)

        import matplotlib.pyplot as plt
        fig, axmat = plt.subplots(nrows=len(select), ncols=1, sharex=True, sharey=False, squeeze=True)

        components = self.reader.computed_components[key]
        for i, (what, ax) in enumerate(zip(select, axmat)):
            self.plot_linear_epsilon(what=what, itemp=itemp, components=components,
                                     ax=ax, xlims=xlims, with_xlabel=(i == len(select) - 1),
                                     show=False)
        return fig

    @add_fig_kwargs
    def plot_chi2(self, key, components="all", what="abs", itemp=0, decompose=False,
                  ax=None, xlims=None, with_xlabel=True, label=None, **kwargs):
        """

        Args:
            key:
            components:
            select:
            itemp: Temperature index.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used.

        Returns:
            `matplotlib` figure
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
        if with_xlabel: ax.set_xlabel('Photon Energy [eV]')
        set_axlims(ax, xlims, "x")
        ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_shg(self, **kwargs):
        return self.plot_chi2(key="shg", **kwargs)

    @add_fig_kwargs
    def plot_leo(self, **kwargs):
        return self.plot_chi2(key="leo", **kwargs)

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("optic = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(optic)"),
            nbv.new_code_cell("optic.ebands.plot();"),
        ])

        # Add plot calls if quantities have been computed.
        for key, info in self.reader.ALL_CHIS.items():
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
    """

    ALL_CHIS = OrderedDict([
        ("linopt", {
            "fullname": "Dielectric function",
            "rank": 2,
            #"terms":
            #"latex": r"\chi(\omega)"
            }
        ),
        ("shg", {
            "fullname": "Second Harmonic Generation",
            "rank": 3,
            "terms": ["shg_inter2w", "shg_inter1w", "shg_intra2w",
                      "shg_intra1w", "shg_intra1wS", "shg_chi2tot"],
            }
            #"latex": r"\chi(-2\omega, \omega, \omega)"
        ),
        ("leo", {
            "fullname": "Linear Electro-optic",
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

    def __init__(self, filepath):
        super(OpticReader, self).__init__(filepath)
        self.ntemp = self.read_dimvalue("ntemp")

        self.computed_components = OrderedDict()
        self.computed_ids = OrderedDict()
        for vname, info in self.ALL_CHIS.items():
            comp_name = vname + "_components"
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

            self.computed_ids[vname] = ids
            self.computed_components[vname] = [itup2s(it) for it in ids]

    @lazy_property
    def wmesh(self):
        """
        Frequency mesh in eV. Note that the same frequency-mesh is used
        for the different optical properties.
        """
        return self.read_value("wmesh")

    def read_epsilon(self, components, itemp=0):
        """
        Args:
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
            itemp: Temperature index.
        """
        # arrays have Fortran shape [two, nomega, num_comp, ntemp]
        if components == "all": components = self.computed_components[key]
        components = list_strings(components)
        if not (self.ntemp > itemp >= 0):
            raise ValueError("Invalid itemp: %s, ntemp: %s" % (itemp, self.ntemp))

        od = OrderedDict([(comp, OrderedDict()) for comp in components])
        for vname in self.ALL_CHIS[key]["terms"]:
            #print("About to read:", vname)
            var = self.read_variable(vname)
            for comp in components:
                try:
                    ijkp = self.computed_components[key].index(comp)
                except ValueError:
                    raise ValueError("%s component %s was not computed" % (key, comp))
                values = var[itemp, ijkp]
                od[comp][vname] = values[:, 0] + 1j * values[:, 1]
        return od


from abipy.abio.robots import Robot, RobotWithEbands

def ordered_intersection(list_1, list_2):
    """Return ordered intersection of two lists. items must be hashable."""
    set_2 = frozenset(list_2)
    return [x for x in list_1 if x in set_2]


class OpticRobot(Robot, RobotWithEbands, NotebookWriter):
    """
    This robot analyzes the results contained in multiple optic.nc files.
    """
    EXT = "OPTIC"

    @lazy_property
    def computed_components_intersection(self):
        """
        Dictionary with the list tensor components available in each file.
        Use response function names as keys.
        """
        od = OrderedDict()
        for ncfile in self.ncfiles:
            for vname in ncfile.reader.ALL_CHIS:
                comps = ncfile.reader.computed_components[vname]
                if vname not in od:
                    od[vname] = comps
                else:
                    # Build intersection while maintain order.
                    od[vname] = ordered_intersection(od[vname], comps)
        return od

    @add_fig_kwargs
    def plot_linopt_convergence(self, components="all", what_list=("re", "im"), sortby="nkpt",
                                 itemp=0, xlims=None, **kwargs):
        """
        Plot convergence of the dielectric function tensor (Real and imaginary part)

        Args:
            components:
            what_list:
            sortby:
            itemp: Temperature index.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used.
        Returns:
            `matplotlib` figure
        """
        # Build grid plot: computed tensors along the rows, what_list along columns.
        key = "linopt"
        components = self.computed_components_intersection[key]
        import matplotlib.pyplot as plt
        fig, axmat = plt.subplots(nrows=len(components), ncols=len(what_list),
                                  sharex=True, sharey=False, squeeze=False)

        label_ncfile_param = self.sortby(sortby)
        for i, comp in enumerate(components):
            for j, what in enumerate(what_list):
                ax = axmat[i, j]
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
        if "what_list" not in kwargs: kwargs["what_list"] = ["abs"]
        return self.plot_convergence_rank3(key="shg", **kwargs)

    @add_fig_kwargs
    def plot_leo_convergence(self, **kwargs):
        if "what_list" not in kwargs: kwargs["what_list"] = ["abs"]
        return self.plot_convergence_rank3(key="leo", **kwargs)

    @add_fig_kwargs
    def plot_convergence_rank3(self, key, components="all", itemp=0, what_list=("abs",),
                               sortby="nkpt", decompose=False, xlims=None, **kwargs):
        """
        Plot convergence of chi2(-2w, w, w)

        Args:
            components:
            itemp: Temperature index.
            what_list:
            decompose:
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used.
        Returns:
            `matplotlib` figure
        """
        # Build grid plot: computed tensors along the rows, what_list along columns.
        components = self.computed_components_intersection[key]
        import matplotlib.pyplot as plt
        fig, axmat = plt.subplots(nrows=len(components), ncols=len(what_list),
                                  sharex=True, sharey=False, squeeze=False)

        label_ncfile_param = self.sortby(sortby)
        for i, comp in enumerate(components):
            for j, what in enumerate(what_list):
                ax = axmat[i, j]
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

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
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

        return self._write_nb_nbpath(nb, nbpath)

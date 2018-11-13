# coding: utf-8
"""
Interface to the GSR.nc_ file storing the Ground-state results and the electron band structure.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pandas as pd
import pymatgen.core.units as units
import abipy.core.abinit_units as abu

from collections import OrderedDict, Iterable, defaultdict
from tabulate import tabulate
from monty.string import is_string, list_strings, marquee
from monty.termcolor import cprint
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from pymatgen.core.units import EnergyArray, ArrayWithUnit
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.tools.plotting import add_fig_kwargs, get_axarray_fig_plt
from abipy.tools.tensors import Stress
from abipy.abio.robots import Robot
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "GsrFile",
]

_INVALID_STRESS_TENSOR = 9999999999e+99


class GsrFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    File containing the results of a ground-state calculation.

    Usage example:

    .. code-block:: python

        with GsrFile("foo_GSR.nc") as gsr:
            print("energy: ", gsr.energy)
            gsr.ebands.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsrFile
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(GsrFile, self).__init__(filepath)

        self.reader = r = GsrReader(filepath)

        # Initialize the electron bands from file
        self._ebands = r.read_ebands()

        # Add forces to structure
        if self.is_scf_run:
            self.structure.add_site_property("cartesian_forces", self.cart_forces)

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
        if self.is_scf_run:
            app("")
            app("Stress tensor (Cartesian coordinates in GPa):\n%s" % self.cart_stress_tensor)
            #if verbose:
            #    app("Stress tensor (Cartesian coordinates in Ha/Bohr**3):\n%s" % self.cart_stress_tensor / abu.HaBohr3_GPa)
            app("")
            app("Pressure: %.3f (GPa)" % self.pressure)
            app("Energy: %.8f (eV)" % self.energy)
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))

        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    @property
    def ebands(self):
        """|ElectronBands| object."""
        return self._ebands

    @lazy_property
    def is_scf_run(self):
        """True if the GSR has been produced by a SCF run."""
        # NOTE: We use kptopt to understand if we have a SCF/NSCF run
        # In principle one should use iscf but it's not available in the GSR.
        #return int(self.reader.read_value("kptopt")) >= 0
        return abs(self.cart_stress_tensor[0, 0] - _INVALID_STRESS_TENSOR) > 0.1

    @lazy_property
    def ecut(self):
        """Cutoff energy in Hartree (Abinit input variable)"""
        return units.Energy(self.reader.read_value("ecut"), "Ha")

    @lazy_property
    def pawecutdg(self):
        """Cutoff energy in Hartree for the PAW double grid (Abinit input variable)"""
        return units.Energy(self.reader.read_value("pawecutdg"), "Ha")

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def energy(self):
        """Total energy in eV."""
        return units.Energy(self.reader.read_value("etotal"), "Ha").to("eV")

    @lazy_property
    def energy_per_atom(self):
        """Total energy / number_of_atoms (eV units)"""
        return self.energy / len(self.structure)

    @lazy_property
    def cart_forces(self):
        """Cartesian forces in eV / Ang"""
        return self.reader.read_cart_forces()

    @lazy_property
    def max_force(self):
        fmods = np.sqrt([np.dot(force, force) for force in self.cart_forces])
        return fmods.max()

    def force_stats(self, **kwargs):
        """
        Return a string with information on the forces.
        """
        fmods = np.sqrt([np.dot(force, force) for force in self.cart_forces])
        imin, imax = fmods.argmin(), fmods.argmax()

        s = "\n".join([
            "fsum: %s" % self.cart_forces.sum(axis=0),
            "mean: %s, std %s" % (fmods.mean(), fmods.std()),
            "minimum at site %s, cart force: %s" % (self.structure.sites[imin], self.cart_forces[imin]),
            "maximum at site %s, cart force: %s" % (self.structure.sites[imax], self.cart_forces[imax]),
        ])

        table = [["Site", "Cartesian Force", "Length"]]
        for i, fmod in enumerate(fmods):
            table.append([self.structure.sites[i], self.cart_forces[i], fmod])
        s += "\n" + tabulate(table)

        return s

    @lazy_property
    def cart_stress_tensor(self):
        """Stress tensor in GPa."""
        return self.reader.read_cart_stress_tensor()

    @lazy_property
    def pressure(self):
        """Pressure in GPa."""
        pressure = - self.cart_stress_tensor.trace() / 3
        return units.FloatWithUnit(pressure, unit="GPa", unit_type="pressure")

    @lazy_property
    def residm(self):
        """Maximum of the residuals"""
        return self.reader.read_value("residm")

    @lazy_property
    def xc(self):
        """
        :class:`XcFunc` object with info on the exchange-correlation functional.
        Use libxc convention :cite:`Marques2012`.
        """
        return self.reader.read_abinit_xcfunc()

    @lazy_property
    def energy_terms(self):
        """:class:`EnergyTerms` with the different contributions to the total energy in eV."""
        return self.reader.read_energy_terms(unit="eV")

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        od["ecut"] = float(self.ecut)
        #if self.usepaw == 1
        #    od["pawecutdg"] = float(self.pawecutdg)
        return od

    def close(self):
        self.reader.close()

    # FIXME: This is deprecated. Must keep it to avoid breaking ScfTask.get_results
    def as_dict(self):
        return {}

    def get_computed_entry(self, inc_structure=True, parameters=None, data=None):
        """
        Returns a pymatgen :class:`ComputedStructureEntry` from the GSR file.
        Same API as the one used in vasp_output.get_computed_entry.

        Args:
            inc_structure (bool): Set to True if you want
                ComputedStructureEntries to be returned instead of ComputedEntries.
            parameters (list): Input parameters to include. It has to be one of
                the properties supported by the GSR object. If
                parameters is None, a default set of parameters that are
                necessary for typical post-processing will be set.
            data (list): Output data to include. Has to be one of the properties
                supported by the GSR object.

        Returns:
            ComputedStructureEntry/ComputedEntry
        """
        # TODO
        #param_names = {"is_hubbard", "hubbards", "potcar_symbols", "run_type"}
        if inc_structure:
            return ComputedStructureEntry(self.structure, self.energy,
                                          correction=0.0, parameters=parameters, data=data)
        else:
            return ComputedEntry(self.structure.composition, self.energy,
                                 parameters=parameters, data=data)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        for fig in self.yield_structure_figs(**kwargs): yield fig
        for fig in self.yield_ebands_figs(**kwargs): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("gsr = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(gsr)"),
            nbv.new_code_cell("gsr.ebands.plot();"),
            nbv.new_code_cell("gsr.ebands.kpoints.plot();"),
            nbv.new_code_cell("# gsr.ebands.plot_transitions(omega_ev=3.0, qpt=(0, 0, 0), atol_ev=0.1);"),
            nbv.new_code_cell("""\
if gsr.ebands.kpoints.is_ibz:
    gsr.ebands.get_edos().plot();"""),
            #nbv.new_code_cell("emass = gsr.ebands.effective_masses(spin=0, band=0, acc=4)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class EnergyTerms(AttrDict):
    """
    Contributions to the total GS energy. See energies_type in m_energies.F90.
    """

    _NAME2DOC = OrderedDict([
        # (Name, help)
        ("e_localpsp", "Local psp energy"),
        ("e_eigenvalues", "Sum of the eigenvalues - Band energy\n" +
                          "(valid for double-counting scheme dtset%optene == 1)"),
        ("e_ewald",  "Ewald energy, store also the ion/ion energy for free boundary conditions."),
        ("e_hartree", "Hartree part of the total energy"),
        ("e_corepsp", "psp core-core energy"),
        ("e_corepspdc", "psp core-core energy double-counting"),
        ("e_kinetic", "Kinetic energy part of total energy. (valid for direct scheme, dtset%optene == 0"),
        ("e_nonlocalpsp", "Nonlocal pseudopotential part of total energy."),
        ("e_entropy", "Entropy energy due to the occupation number smearing (if metal)\n" +
                      "Value is multiplied by dtset%tsmear, see %entropy for the entropy alone\n." +
                      "(valid for metals, dtset%occopt>=3 .and. dtset%occopt<=8)"),
        ("entropy", "Entropy term"),
        ("e_xc", "Exchange-correlation energy"),
        #("e_vdw_dftd2", "Dispersion energy from DFT-D2 Van der Waals correction"),
        ("e_xcdc", "enxcdc=exchange-correlation double-counting energy"),
        ("e_paw", "PAW spherical part energy"),
        ("e_pawdc", "PAW spherical part double-counting energy"),
        ("e_elecfield", "Electric enthalpy, by adding both ionic and electronic contributions"),
        ("e_magfield", "Orbital magnetic enthalpy, by adding orbital contribution"),
        ("e_fermie", "Fermie energy"),
        ("e_sicdc", "Self-interaction energy double-counting"),
        ("e_exactX", "Fock exact-exchange energy"),
        ("h0", "h0=e_kinetic+e_localpsp+e_nonlocalpsp"),
        ("e_electronpositron", "Electron-positron: electron-positron interaction energy"),
        ("edc_electronpositron", "Electron-positron: double-counting electron-positron interaction energy"),
        ("e0_electronpositron", "Electron-positron: energy only due to unchanged particles\n" +
                                "(if calctype=1, energy due to electrons only)\n" +
                                "(if calctype=2, energy due to positron only)\n"),
        ("e_monopole", "Monopole correction to the total energy for charged supercells"),
        # FIXME: Some names have been changed in Abinit8, I should recheck the code.
        #("e_xc_vdw", "vdW-DF correction to the XC energy"),
    ])

    ALL_KEYS = list(_NAME2DOC.keys())

    def __str__(self):
        return self.to_string(with_doc=False)

    __repr__ = __str__

    def to_string(self, verbose=0, with_doc=True):
        """String representation, with documentation if with_doc."""
        lines = [str(self.table)]
        if with_doc:
            for k, doc in self._NAME2DOC.items():
                lines.append("%s: %s" % (k, doc))

        return "\n".join(lines)

    @property
    def table(self):
        """string with results in tabular form."""
        table = [["Term", "Value"]]
        for k, doc in self._NAME2DOC.items():
            table.append([k, self[k]])
        return tabulate(table, tablefmt="plain")

    #def get_dataframe(self):
    #    """Return a |pandas-DataFrame|"""
    #    d = {k: float(self[k]) for k in self}
    #    return pd.DataFrame(d, index=[None], columns=list(d.keys()))
    #    #return pd.DataFrame(d, columns=list(d.keys()))


class GsrReader(ElectronsReader):
    """
    This object reads the results stored in the _GSR (Ground-State Results) file produced by ABINIT.
    It provides helper function to access the most important quantities.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsrReader
    """
    def read_cart_forces(self, unit="eV ang^-1"):
        """
        Read and return a |numpy-array| with the cartesian forces in unit ``unit``.
        Shape (natom, 3)
        """
        return ArrayWithUnit(self.read_value("cartesian_forces"), "Ha bohr^-1").to(unit)

    def read_cart_stress_tensor(self):
        """
        Return the stress tensor (3x3 matrix) in cartesian coordinates in GPa.
        If MaskedArray (i.e. tensor was not computed  e.g. Nscf run) set it to _INVALID_STRESS_TENSOR
        """
        # Abinit stores 6 unique components of this symmetric 3x3 tensor:
        # Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
        c = self.read_value("cartesian_stress_tensor")
        tensor = np.empty((3, 3), dtype=np.float)

        if np.ma.is_masked(c[()]):
            # NSCF
            tensor.fill(_INVALID_STRESS_TENSOR)
        else:
            for i in range(3): tensor[i, i] = c[i]
            for p, (i, j) in enumerate(((2, 1), (2, 0), (1, 0))):
                tensor[i, j] = c[3 + p]
                tensor[j, i] = c[3 + p]
            tensor *= abu.HaBohr3_GPa

        return Stress(tensor)

    def read_energy_terms(self, unit="eV"):
        """
        Return a dictionary with the different contributions to the total electronic energy.
        """
        convert = lambda e: units.Energy(e, unit="Ha").to(unit)
        d = OrderedDict()
        for k in EnergyTerms.ALL_KEYS:
            if k == "e_nonlocalpsp" and k not in self.rootgrp.variables:
                # Renamed in 8.9
                d[k] = convert(self.read_value("e_nlpsp_vfock"))
            else:
                d[k] = convert(self.read_value(k))

        return EnergyTerms(**d)


class GsrRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple GSR.nc_ files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsrRobot
    """
    EXT = "GSR"

    def get_dataframe(self, with_geo=True, abspath=False, funcs=None, **kwargs):
        """
        Return a |pandas-DataFrame| with the most important GS results.
        and the filenames as index.

        Args:
            with_geo: True if structure info should be added to the dataframe
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().

        kwargs:
            attrs:
                List of additional attributes of the |GsrFile| to add to the DataFrame.
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a |GsrFile| object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
        """
        # Add attributes specified by the users
        # TODO add more columns
        attrs = [
            "energy", "pressure", "max_force",
            "ecut", "pawecutdg",
            "tsmear", "nkpt",
            "nsppol", "nspinor", "nspden",
        ] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for label, gsr in self.items():
            row_names.append(label)
            d = OrderedDict()

            # Add info on structure.
            if with_geo:
                d.update(gsr.structure.get_dict4pandas(with_spglib=True))

            for aname in attrs:
                if aname == "nkpt":
                    value = len(gsr.ebands.kpoints)
                else:
                    value = getattr(gsr, aname, None)
                    if value is None: value = getattr(gsr.ebands, aname, None)
                d[aname] = value

            # Execute functions
            if funcs is not None: d.update(self._exec_funcs(funcs, gsr))
            rows.append(d)

        row_names = row_names if not abspath else self._to_relpaths(row_names)
        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    def get_eos_fits_dataframe(self, eos_names="murnaghan"):
        """
        Fit energy as function of volume to get the equation of state,
        equilibrium volume, bulk modulus and its derivative wrt to pressure.

        Args:
            eos_names: String or list of strings with EOS names.
                For the list of available models, see pymatgen.analysis.eos.

        Return:
            (fits, dataframe) namedtuple.
                fits is a list of ``EOSFit object``
                dataframe is a |pandas-DataFrame| with the final results.
        """
        # Read volumes and energies from the GSR files.
        energies, volumes = [], []
        for label, gsr in self.items():
            energies.append(float(gsr.energy))
            volumes.append(float(gsr.structure.volume))

        # Order data by volumes if needed.
        if np.any(np.diff(volumes) < 0):
            ves = sorted(zip(volumes, energies), key=lambda t: t[0])
            volumes = [t[0] for t in ves]
            energies = [t[1] for t in ves]

        # Note that eos.fit expects lengths in Angstrom, and energies in eV.
        # I'm also monkey-patching the plot method.
        from pymatgen.analysis.eos import EOS
        if eos_names == "all":
            # Use all the available models.
            eos_names = [n for n in EOS.MODELS if n not in ("deltafactor", "numerical_eos")]
        else:
            eos_names = list_strings(eos_names)

        fits, index, rows = [], [], []
        for eos_name in eos_names:
            try:
                fit = EOS(eos_name=eos_name).fit(volumes, energies)
            except Exception as exc:
                cprint("EOS %s raised exception:\n%s" % (eos_name, str(exc)))
                continue

            # Replace plot with plot_ax method
            fit.plot = fit.plot_ax
            fits.append(fit)
            index.append(eos_name)
            rows.append(OrderedDict([(aname, getattr(fit, aname)) for aname in
                ("v0", "e0", "b0_GPa", "b1")]))

        dataframe = pd.DataFrame(rows, index=index, columns=list(rows[0].keys()) if rows else None)
        return dict2namedtuple(fits=fits, dataframe=dataframe)

    def get_energyterms_dataframe(self, iref=None):
        """
        Build and return with the different contributions to the total energy in eV

        Args:
            iref: Index of the abifile used as reference: the energies of the
            ``iref`` gsrfile will be subtracted from the other rows. Ignored if None.

        Return: |pandas-Dataframe|
        """
        rows, row_names = [], []
        for label, gsr in self.items():
            row_names.append(label)
            # Add total energy.
            d = OrderedDict([("energy", gsr.energy)])
            d.update(gsr.energy_terms)
            d.update(gsr.params)
            rows.append(d)

        df = pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))
        if iref is not None:
            # Subtract iref row from the rest of the rows.
            iref_row = df.iloc[[iref]].values[0]
            df = df.apply(lambda row: row - iref_row, axis=1)

        return df

    @add_fig_kwargs
    def gridplot_eos(self, eos_names="all", fontsize=6, **kwargs):
        """
        Plot multiple EOS on a grid with captions showing the final results.

        Args:
            eos_names: String or list of strings with EOS names. See pymatgen.analysis.EOS
            fontsize: Fontsize used for caption text.

        Returns: |matplotlib-Figure|
        """
        r = self.get_eos_fits_dataframe(eos_names=eos_names)

        num_plots, ncols, nrows = len(r.fits), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        # Build grid of plots.
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        for i, (fit, ax) in enumerate(zip(r.fits, ax_list)):
            fit.plot_ax(ax=ax, fontsize=fontsize, label="", show=False)

        # Get around a bug in matplotlib
        if num_plots % ncols != 0:
            ax_list[-1].axis('off')

        return fig

    @add_fig_kwargs
    def plot_gsr_convergence(self, sortby=None, hue=None, fontsize=6,
                             items=("energy", "pressure", "max_force"), **kwargs):
        """
        Plot the convergence of the most important quantities available in the GSR file
        wrt to the ``sortby`` parameter. Values can optionally be grouped by ``hue``.

        Args:
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            items: List of GSR attributes (or callables) to be analyzed.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|

        Example:

             robot.plot_gsr_convergence(sortby="nkpt", hue="tsmear")
        """
        return self.plot_convergence_items(items, sortby=sortby, hue=hue, fontsize=fontsize, show=False, **kwargs)

    #def get_phasediagram_results(self):
    #    from abipy.core.restapi import PhaseDiagramResults
    #    entries = []
    #    for label, gsr in self.items():
    #        entries.append(gsr.get_computed_entry(inc_structure=True, parameters=None, data=None))
    #    return PhaseDiagramResults(entries)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.plot_lattice_convergence(show=False)
        yield self.plot_gsr_convergence(show=False)
        for fig in self.get_ebands_plotter().yield_figs(): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.GsrRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
            nbv.new_code_cell("df = ebands_plotter.get_ebands_frame()\ndisplay(df)"),
            nbv.new_code_cell("ebands_plotter.ipw_select_plot()"),
            nbv.new_code_cell("#anim = ebands_plotter.animate();"),
            nbv.new_code_cell("edos_plotter = robot.get_edos_plotter()"),
            nbv.new_code_cell("edos_plotter.ipw_select_plot()"),
            nbv.new_code_cell("#robot.gridplot_eos();"),
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

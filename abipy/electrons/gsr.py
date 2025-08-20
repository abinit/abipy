# coding: utf-8
"""
Interface to the GSR.nc_ file storing the Ground-state results and the electron band structure.
"""
from __future__ import annotations

import sys
import dataclasses
import numpy as np
import pandas as pd
import pymatgen.core.units as units
import abipy.core.abinit_units as abu

from collections import OrderedDict
from typing import Optional
from functools import cached_property
from tabulate import tabulate
from monty.string import list_strings, marquee
from monty.termcolor import cprint
from monty.collections import AttrDict, dict2namedtuple
from pymatgen.core.units import ArrayWithUnit
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.structure import Structure
from abipy.tools.plotting import add_fig_kwargs, get_axarray_fig_plt, get_ax_fig_plt, set_grid_legend, set_ax_xylabels
from abipy.tools.typing import Figure
from abipy.abio.robots import Robot
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands, ElectronBands


__all__ = [
    "GsrFile",
]

_INVALID_STRESS_TENSOR = 9999999999e+99


@dataclasses.dataclass(kw_only=True)
class MagneticData:
    spinat: np.ndarray
    use_gbt: int
    qgbt: np.ndarray | None
    #from pymatgen.electronic_structure.core import Magmom

    @classmethod
    def from_gsr(cls, gsr) -> MagneticData:
        """
        Build an instance from a GSR file.
        """
        spinat = gsr.r.read_value("spinat")
        intgden = gsr.r.read_value("intgden")
        nspden = intgden.shape[1]

        qgbt = None
        if (use_gbt := gsr.r.read_value("use_gbt", default=0)) != 0:
            qgbt = gsr.r.read_value("qgbt")

        energy_mev_pat = float(gsr.energy) * 1000 / len(gsr.structure)

        if nspden == 2:
            magmoms = Magmom(intgden[:, 1] - intgden[:, 0])
        elif nspden == 4:
            magmoms = [Magmom([intg_at[1], intg_at[2], intg_at[3]]) for intg_at in intgden]
        else:
            magmoms = None

        locs = locals()
        return cls**{locs[field] for field in dataclasses.fields(cls)}


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
    def from_file(cls, filepath: str) -> GsrFile:
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.r = self.reader = GsrReader(filepath)

        # Add forces to structure
        if self.is_scf_run:
            self.structure.add_site_property("cartesian_forces", self.cart_forces)

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
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

    @cached_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @cached_property
    def is_scf_run(self) -> bool:
        """True if the GSR has been produced by a SCF run."""
        # NOTE: We use kptopt to understand if we have a SCF/NSCF run
        # In principle one should use iscf but it's not available in the GSR.
        if "kptopt" in self.r.rootgrp.variables:
            return int(self.r.read_value("kptopt")) >= 0
        else:
            return abs(self.cart_stress_tensor[0, 0] - _INVALID_STRESS_TENSOR) > 0.1

    @cached_property
    def ecut(self):
        """Cutoff energy in Hartree (Abinit input variable)"""
        return units.Energy(self.r.read_value("ecut"), "Ha")

    @cached_property
    def pawecutdg(self):
        """Cutoff energy in Hartree for the PAW double grid (Abinit input variable)"""
        return units.Energy(self.r.read_value("pawecutdg"), "Ha")

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    @cached_property
    def energy(self):
        """Total energy in eV."""
        return units.Energy(self.r.read_value("etotal"), "Ha").to("eV")

    @cached_property
    def energy_per_atom(self):
        """Total energy / number_of_atoms (eV units)"""
        return self.energy / len(self.structure)

    @cached_property
    def cart_forces(self):
        """
        Cartesian forces in eV / Ang. None if forces are not available.
        """
        if self.is_scf_run:
            return self.r.read_cart_forces()
        return None

    @cached_property
    def max_force(self):
        """
        Max absolute cartesian force in eV/Ang. None if forces are not available.
        """
        cart_forces = self.cart_forces
        if cart_forces is None: return None

        fmods = np.sqrt([np.dot(force, force) for force in cart_forces])
        return fmods.max()

    def force_stats(self, **kwargs):
        """
        Return a string with information on the forces.
        Return None if forces are not available.
        """
        cart_forces = self.cart_forces
        if cart_forces is None: return None

        fmods = np.sqrt([np.dot(force, force) for force in cart_forces])
        imin, imax = fmods.argmin(), fmods.argmax()

        s = "\n".join([
            "fsum: %s" % cart_forces.sum(axis=0),
            "mean: %s, std %s" % (fmods.mean(), fmods.std()),
            "minimum at site %s, cart force: %s" % (self.structure.sites[imin], cart_forces[imin]),
            "maximum at site %s, cart force: %s" % (self.structure.sites[imax], cart_forces[imax]),
        ])

        table = [["Site", "Cartesian Force", "Length"]]
        for i, fmod in enumerate(fmods):
            #table.append([self.structure.sites[i], cart_forces[i], fmod])
            table.append([str(self.structure.sites[i]), str(cart_forces[i]), str(fmod)])
        s += "\n" + tabulate(table)

        return s

    @cached_property
    def cart_stress_tensor(self):
        """
        Stress tensor in GPa. Return None if not available e.g. if NSCF run.
        """
        if self.is_scf_run:
            return self.r.read_cart_stress_tensor()
        return None

    @cached_property
    def pressure(self):
        """
        Pressure in GPa. Return None if not available e.g. if NSCF run.
        """
        if self.is_scf_run:
            pressure = - self.cart_stress_tensor.trace() / 3
            return units.FloatWithUnit(pressure, unit="GPa", unit_type="pressure")
        return None

    @cached_property
    def residm(self):
        """
        Maximum of the residuals
        """
        return self.r.read_value("residm")

    @cached_property
    def xc(self):
        """
        :class:`XcFunc` object with info on the exchange-correlation functional.
        Use libxc convention :cite:`Marques2012`.
        """
        return self.r.read_abinit_xcfunc()

    @cached_property
    def energy_terms(self):
        """:class:`EnergyTerms` with the different contributions to the total energy in eV."""
        return self.r.read_energy_terms(unit="eV")

    def get_magnetization(self) -> np.ndarray:
        """
        Magnetization in Cartesian directions in atomic units.
        """
        rhomag = self.r.read_value("rhomag", default=None)
        if rhomag is None:
            raise RuntimeError("You GSR file does not contain the rhomag variable. Please use Abinit version >= 10.4")

        # rhomag(2, nspden) (Fortran array)
        #   in collinear case component 1 is total density and 2 is _magnetization_ up-down
        #   in non collinear case component 1 is total density, and 2:4 are the magnetization vector
        rhomag = rhomag[:,0]
        mag = np.zeros(3)
        if self.ebands.nspden == 2: mag[2] = rhomag[1]
        if self.ebands.nspden == 4: mag = rhomag[1:]
        return mag

    @cached_property
    def params(self) -> dict:
        """dict with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        od["ecut"] = float(self.ecut)
        #if self.hdr.usepaw == 1
        #    od["pawecutdg"] = float(self.pawecutdg)
        return od

    def close(self) -> None:
        self.r.close()

    # FIXME: This is deprecated. Must keep it to avoid breaking ScfTask.get_results
    def as_dict(self) -> dict:
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

    def print_efg_results(self, precision=4, file=sys.stdout) -> None:
        """
        @JOE: Can you please describe the goal of this method?

        Args:
            precision: print options precision.
            file: File handle for output.
        """
        # This code has been taken from efg_results.
        def _p(*args, **kwargs):
            return print(*args, file=file, **kwargs)

        scale_factor = 234.9599245

        if (efg := self.r.read_value("efg", default=None)) is None:
            raise ValueError(f"GSR file {self.filepath} does not contain EFG data!")

        if (quadmom := self.r.read_value("quadmom", default=None)) is None:
            _p("Found no quadrupole moment data, using 0.0 for all atoms")

        from numpy.linalg import eigvals
        atom_species = self.r.read_value('atom_species')
        atom_species_names = self.r.read_value('atom_species_names')

        #with np.set_printoptions(precision=precision):
        _p("Field gradient data")
        for iat in range(len(self.structure)):
            itypat = atom_species[iat]
            vpas = eigvals(efg[iat])
            vzz = vpas[np.argmax(np.abs(vpas))]
            vxx = vpas[np.argmin(np.abs(vpas))]
            vyy = -vzz -vxx

            eta = (vxx - vyy) / vzz if abs(vzz) > 1.0E-8 else 0.0

            cq = vzz * quadmom[itypat-1] * scale_factor
            _p('atom type '+ str(atom_species[iat]) + ' Cq(MHz): %7.3f   eta: %4.3f' % (cq, eta))

    def get_panel(self, **kwargs):
        """
        Build panel with widgets to interact with the |GsrFile| either in a notebook or in panel app.
        """
        from abipy.panels.gsr import GsrFilePanel
        return GsrFilePanel(self).get_panel(**kwargs)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        verbose = kwargs.get("verbose", 0)
        for fig in self.yield_ebands_figs(**kwargs): yield fig
        if verbose:
            for fig in self.yield_structure_figs(**kwargs): yield fig

    def yield_plotly_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of plotly figures with minimal input from the user.
        """
        verbose = kwargs.get("verbose", 0)
        for fig in self.yield_ebands_plotly_figs(**kwargs): yield fig
        if verbose:
            for fig in self.yield_structure_plotly_figs(**kwargs): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporaty file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)
        first_char = "" if self.has_panel() else "#"

        nb.cells.extend([
            nbv.new_code_cell("gsr = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(gsr)"),

            # Add panel GUI but comment the python code if panel is not available.
            nbv.new_markdown_cell("## Panel dashboard"),
            nbv.new_code_cell(f"""\
# Execute this cell to display the panel GUI (requires panel package).
# To display the dashboard inside the browser use `abiopen.py FILE --panel`.

{first_char}abilab.abipanel()
{first_char}gsr.get_panel()
"""),
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

    def to_string(self, verbose: int = 0, with_doc: bool = True) -> str:
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

    def get_dataframe(self) -> pd.DataFrame:
        """Return a |pandas-DataFrame|"""
        d = {k: float(self[k]) for k in self}
        return pd.DataFrame(d, index=[None], columns=list(d.keys()))


class EnergyTermsPlotter:

    @classmethod
    def from_label_file_dict(cls, label_file_dict: dict) -> EnergyTermsPlotter:
        """
        Build an object from a dictionary mapping labels to filepath.

        Usage example:

        .. code-block:: python

            plotter = EnergyTermsPlotter.from_label_file_dict({
               "label1": "out1_GSR.nc",
               "label2": "out2_GSR.nc",
            })

        """
        eterms_list = []
        for label, path in label_file_dict.items():
            if isinstance(path, GsrFile):
                eterms_list.append(path.eterms)
            else:
                with GsrFile(path) as gsr:
                    eterms_list.append(gsr.eterms)

        return cls(list(glabel_file_dict.keys()), eterms_list)

    def __init__(self, labels, eterms_list):
        self.labels = labels
        self.eterms_list = eterms_list

    #def add(self, label: str, gsr_path: str) -> None:
    #   self.labels.append(label)
    #    if isinstance(path, GsrFile):
    #        self.eterms_list.append(path.eterms)
    #    else:
    #    with GsrFile(gsr_path) as gsr:
    #        self.labels.append(

    def get_dataframe(self) -> pd.DataFrame:
        df_list = []
        for label, eterms in zip(self.labels, self.eterms_list):
            df = eterms.get_dataframe()
            df["label"] = label
            df_list.append(df)

        return pd.concat(df_list)

    #def plot(self, what_list=("foo", "bar"), fontsize=8, **kwargs) -> Figure:
    #    df = self.get_dataframe()
    #    #keys = [k for k in df.keys() if k != "label")
    #    #nkeys = len(keys)

    #    # Build grid of plots.
    #    ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=1,
    #                                            sharex=True, sharey=True, squeeze=False)
    #    ax_list = ax_list.ravel()

    #    for ax, what in zip(ax_list, what_list):

    #    return fig


class GsrReader(ElectronsReader):
    """
    This object reads the results stored in the _GSR (Ground-State Results) file produced by ABINIT.
    It provides helper function to access the most important quantities.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsrReader
    """
    def read_cart_forces(self, unit="eV ang^-1"):
        """
        Read and return a |numpy-array| with the cartesian forces in unit ``unit``. Shape (natom, 3)
        """
        return ArrayWithUnit(self.read_value("cartesian_forces"), "Ha bohr^-1").to(unit)

    def read_cart_stress_tensor(self, units: str = "GPa"):
        """
        Return the stress tensor (3x3 matrix) in cartesian coordinates in GPa.
        If MaskedArray (i.e. tensor was not computed  e.g. Nscf run) set it to _INVALID_STRESS_TENSOR

        Args:
            units: "GPa" for Gpa units or "au" for atomic units (Ha/Bohr^3)
        """
        # Abinit stores 6 unique components of this symmetric 3x3 tensor:
        # Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
        c = self.read_value("cartesian_stress_tensor")
        tensor = np.empty((3, 3), dtype=float)

        if np.ma.is_masked(c[()]):
            # NSCF run
            tensor.fill(_INVALID_STRESS_TENSOR)
        else:
            for i in range(3):
                tensor[i, i] = c[i]
            for p, (i, j) in enumerate(((2, 1), (2, 0), (1, 0))):
                tensor[i, j] = c[3 + p]
                tensor[j, i] = c[3 + p]

            if units == "GPa":
                tensor *= abu.HaBohr3_GPa
            elif units == "au":
                pass
            else:
                raise ValueError(f"Invalid {units=}")

        from abipy.tools.tensors import Stress
        return Stress(tensor)

    def read_energy_terms(self, unit: str ="eV") -> EnergyTerms:
        """
        Return a dictionary with the different contributions to the total electronic energy.
        """
        convert = lambda e: units.Energy(e, unit="Ha").to(unit)
        d = {}
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

    def get_dataframe(self, with_geo=True, abspath=False, with_paths=True, funcs=None, **kwargs) -> pd.DataFrame:
        """
        Return a |pandas-DataFrame| with the most important GS results and the filenames as index.

        Args:
            with_geo: True if structure info should be added to the dataframe
            abspath: True if paths in the index should be absolute. Default: Relative to getcwd().
            with_paths: False if filepaths should not be added

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
            "energy", "energy_per_atom", "pressure", "max_force",
            "ecut", "pawecutdg", "tsmear", "nkpt",
            "nsppol", "nspinor", "nspden",
        ] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for label, gsr in self.items():
            row_names.append(label)
            d = {}

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

        index = None
        if with_paths:
            index = row_names if not abspath else self._to_relpaths(row_names)

        return pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))

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

    def get_energyterms_dataframe(self, iref: Optional[int] = None) -> pd.DataFrame:
        """
        Build and return dataframe with the different contributions to the total energy in eV.

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
    def gridplot_eos(self, eos_names="all", fontsize=6, **kwargs) -> Figure:
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
    def plot_gsr_convergence(self, sortby=None, hue=None, fontsize=8,
                             items=("energy", "pressure", "max_force"), **kwargs) -> Figure:
        """
        Plot the convergence of the most important quantities available in the GSR file
        wrt to the ``sortby`` parameter. Values can be optionally grouped by ``hue``.

        Args:
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of `sortby(abifile)` is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of `hue(abifile)` is used.
            items: List of GSR attributes (or callables) to be analyzed.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|

        .. code-block::

             gsr.plot_gsr_convergence(sortby="ecut")
             gsr.plot_gsr_convergence(sortby="nkpt", hue="tsmear")
        """
        return self.plot_convergence_items(items, sortby=sortby, hue=hue,
                                           fontsize=fontsize, show=False, **kwargs)

    def get_spin_spiral_df(self,
                           with_params: bool = False,
                           with_geo: bool = False,
                           ) -> pd.DataFrame:
        """
        Build and return a dataframe with the atomic magnetization/charge for each
        site, the total energy, and the GBT q-point.

        Args:
            with_params: True if convergence params should be added to the dataframe
            with_geo: True if structure info should be added to the dataframe
        """
        # nctkarr_t("intgden", "dp", "number_of_components, number_of_atoms"), &
        # nctkarr_t("ratsph", "dp", "number_of_atom_species"), &
        # nctkarr_t("rhomag", "dp", "two, number_of_components") &

        # GBT calculations are done with nsym 1 but here we need a
        # structure with the spacegroup to find the equivalent q-points.
        structure0 = self.abifiles[0].structure.copy()
        structure0.spgset_abi_spacegroup(has_timerev=False, overwrite=True)

        rows = []
        for iq, (label, gsr) in enumerate(self.items()):
            #mag_data = MagneticData.from_gsr(gsr)
            spinat = gsr.r.read_value("spinat")
            intgden = gsr.r.read_value("intgden")
            nspden = intgden.shape[1]

            qgbt = None
            if (use_gbt := gsr.r.read_value("use_gbt", default=0)) != 0:
                qgbt = gsr.r.read_value("qgbt")
                qname = structure0.findname_in_hsym_stars(qgbt)

            # Convert energies to mev per atom.
            energy_mev_pat = float(gsr.energy) * 1000 / len(gsr.structure)

            for iat, site in enumerate(gsr.structure):
                magmom=intgden[iat, 1] - intgden[iat, 0] if nspden == 2 else intgden[iat, 1:]
                d = dict(
                    site_idx=iat,
                    symbol=site.specie.symbol,
                    frac_coords=site.frac_coords,
                    magmom=magmom,
                    magmom_norm=np.linalg.norm(magmom),
                    spinat=spinat[iat],
                    charge_isph=intgden[iat, 1] + intgden[iat, 0] if nspden == 2 else intgden[iat, 0],
                    energy_mev_pat=energy_mev_pat,
                )
                if qgbt is not None:
                    # Add qgbt and sequential index
                    d.update(qgbt=qgbt, qname=qname, iq=iq)

                if with_params:
                    d.update(gsr.params)
                if with_geo:
                    d.update(gsr.structure.get_dict4pandas(with_spglib=True))

                rows.append(d)

        return pd.DataFrame(rows)

    @add_fig_kwargs
    def plot_spin_spiral_magmom(self,
                                keys=("magmom_norm", "mx", "my", "mz", "charge_isph"),
                                symbols: str | list[str] | None = None,
                                site_inds: list[int] | None = None,
                                fontsize=8, **kwargs) -> Figure:
        """
        Plot the magnetic moments obtained with the generalized Bloch theorem
        as a function of the wave-vector q.

        Args:
            keys: List of quantities to plot.
            symbols: string or list of strings with chemical symbols to show.
                None means all chemical symbols in the structure.
            site_inds: List of site indices to show. None means all sites.
            fontsize: legend and label fontsize.
        """
        # Get dataframe with results
        df = self.get_spin_spiral_df()

        # Build grid of plots.
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=len(keys), ncols=1,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        structure0 = self.abifiles[0].structure
        if symbols is not None:
            symbols = list_strings(symbols)

        ticks, labels = None, None

        for site_idx, site in enumerate(structure0):
            # Filtering on symbol or site index.
            symbol = site.specie.symbol
            if symbols is not None and symbol not in symbols: continue
            if site_inds is not None and site_idx not in site_inds: continue

            # Select data for this site index.
            data = df[df["site_idx"] == site_idx]

            if ticks is None:
                # Get ticks and labels.
                ticks, labels = data["iq"].values, data["qname"].values
                # Filter and then unpack
                filtered_pairs = [(x, y) for x, y in zip(ticks, labels) if y is not None]
                ticks, labels = zip(*filtered_pairs)

            for ax, key in zip(ax_list, keys, strict=True):
                if key in ("mx", "my", "mz"):
                    # Convert magmom column to (nq, 3) array and select the Cartesian component.
                    idx = {"mx": 0, "my": 1, "mz": 2}[key]
                    ys = np.array([y for y in data["magmom"].values])[:,idx]
                else:
                    ys = data[key]

                ax.plot(ys, label="$%s_{%d}$" % (symbol, site_idx))

        for ix, (ax, key) in enumerate(zip(ax_list, keys, strict=True)):
            set_grid_legend(ax, fontsize=fontsize)
            ax.set_ylabel(key)
            if ix == len(ax_list) - 1:
                ax.set_xlabel("Wave Vector q")
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.get("qlabel_size", "large"))
                if len(ticks) > 1:
                    ax.set_xlim(ticks[0], ticks[-1])

        return fig

    @add_fig_kwargs
    def plot_spin_spiral(self, ax=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot spin spiral energy E(q) obtained with the generalized Bloch theorem.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: legend and label fontsize.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax, grid=True)

        # TODO: Should check that all structures are the same.
        # GBT calculations are done with nsym 1 but here we need a
        # structure with the spacegroup to find the equivalent q-points.
        structure0 = self.abifiles[0].structure
        structure0.spgset_abi_spacegroup(has_timerev=False, overwrite=True)
        natom = len(structure0)

        # Read GBT q-point and energies from the GSR files.
        qpoints, energies, spinat = [], [], None
        for label, gsr in self.items():
            if (use_gbt := gsr.r.read_value("use_gbt", default=0)) == 0:
                raise RuntimeError(f"{gsr.filepath=} has {use_gbt=}")

            if (spinat_ := gsr.r.read_value("spinat")) is None:
                spinat = spinat_
            if spinat is not None and not np.allclose(spinat, spinat_):
                cprint(f"spinat is not the same:\n {spinat_=}\n{spinat=}", color="yellow")

            qpoints.append(gsr.r.read_value("qgbt"))
            energies.append(float(gsr.energy))

        # Convert energies to mev per atom.
        energies_mev = np.array(energies) * 1000 / natom
        energies_mev -= energies_mev.min()

        # Find the q-point where we have the minimum.
        qmin = qpoints[np.argmin(energies_mev)]
        if (qmin_name := structure0.findname_in_hsym_stars(qmin)) is None:
            qmin_name = ""

        kw_color = kwargs.pop("color", "k")
        ax.plot(energies_mev, color=kw_color)

        ax.set_xlabel("Wave Vector q")
        ax.set_ylabel("Energy/atom (meV)")
        ax.set_title(f"Minimum at q: {qmin} {qmin_name}", fontsize=fontsize)

        ticks, labels = zip(*[(i, qname) for i, qpt in enumerate(qpoints)
            if (qname := structure0.findname_in_hsym_stars(qpt)) is not None])

        ax.set_xticks(ticks, minor=False)
        ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.get("qlabel_size", "large"))
        if len(ticks) > 1:
            ax.set_xlim(ticks[0], ticks[-1])

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.plot_lattice_convergence(show=False)
        yield self.plot_gsr_convergence(show=False)
        for fig in self.get_ebands_plotter().yield_figs(): yield fig

    def get_panel(self, **kwargs):
        """
        Build panel with widgets to interact with the |GsrRobot| either in a notebook or in panel app.
        """
        from abipy.panels.gsr import GsrRobotPanel
        return GsrRobotPanel(robot=self).get_panel(**kwargs)

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
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

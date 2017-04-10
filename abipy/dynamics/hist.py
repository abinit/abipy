# coding: utf-8
"""DDB File."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from monty.functools import lazy_property
from monty.collections import AttrDict
from monty.string import marquee # is_string, list_strings,
from pymatgen.core.units import EnergyArray, ArrayWithUnit
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.core.structure import Structure
from abipy.core.mixins import AbinitNcFile, NotebookWriter
from abipy.iotools import ETSF_Reader

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "HistFile",
]


class HistFile(AbinitNcFile, NotebookWriter):
    """
    File with the history of a structural relaxation or molecular dynamics calculation.

    Usage example:

    .. code-block:: python

        with HistFile("foo_HIST") as hist:
            hist.plot()
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(HistFile, self).__init__(filepath)
        self.reader = HistReader(filepath)

    def close(self):
        self.reader.close()

    def __str__(self):
        return self.to_string()

    def to_string(self):
        """Return string representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Final structure", mark="="))
        app("Number of relaxation steps performed: %d" % self.num_steps)
        app(str(self.final_structure))
        app("")

        cart_stress_tensors, pressures = self.reader.read_cart_stress_tensors()
        app("Stress tensor (Cartesian coordinates in Ha/Bohr**3):\n%s" % cart_stress_tensors[-1])
        app("Pressure: %.3f [GPa]" % pressures[-1])

        return "\n".join(lines)

    @property
    def num_steps(self):
        """Number of iterations performed."""
        return self.reader.num_steps

    @lazy_property
    def steps(self):
        """step indices."""
        return list(range(self.num_steps))

    @property
    def final_structure(self):
        """The structure of the last iteration."""
        return self.structures[-1]

    @lazy_property
    def structures(self):
        """List of :class:`Structure` objects at the different steps."""
        return self.reader.read_all_structures()

    @lazy_property
    def etotals(self):
        """numpy array with total energies in eV at the different steps."""
        return self.reader.read_eterms().etotals

    def export(self, filename, visu=None):
        """
        Export the crystalline structure on file filename.

        Args:
            filename: String specifying the file path and the file format.
                The format is defined by the file extension. filename="prefix.xsf", for example,
                will produce a file in XSF format. An *empty* prefix, e.g. ".xsf" makes the code use a temporary file.
            visu: `Visualizer` subclass. By default, this method returns the first available
                visualizer that supports the given file format. If visu is not None, an
                instance of visu is returned. See :class:`Visualizer` for the list of applications and formats supported.

        Returns: Instance of :class:`Visualizer`
        """
        print("Warning: work in progress")
        raise NotImplementedError("typat is missing in HIST --> wrong structures")

        if "." not in filename:
            raise ValueError("Cannot detect extension in filename %s: " % filename)

        from abipy.iotools.xsf import xsf_write_structure
        with open(filename, "w") as fh:
            xsf_write_structure(fh, self.structures)

    @add_fig_kwargs
    def plot(self, axlist=None, **kwargs):
        """
        Plot the evolution of structural parameters (lattice lengths, angles and volume)
        as well as pressure, info on forces and total energy.

        Args:
            axlist: List of matplotlib Axes. If None, a new figure is created.

        Returns:
            `matplotlib` figure
        """
        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=3, ncols=2, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()
        ax0, ax1, ax2, ax3, ax4, ax5 = ax_list
        for ax in ax_list: ax.grid(True)

        # Lattice parameters.
        for i, label in enumerate(["a", "b", "c"]):
            ax0.plot(self.steps, [s.lattice.abc[i] for s in self.structures], marker="o", label=label)
        ax0.set_ylabel('Lattice lengths [A]')
        ax0.legend(loc='best', shadow=True)

        # Lattice Angles
        for i, label in enumerate(["alpha", "beta", "gamma"]):
            ax1.plot(self.steps, [s.lattice.angles[i] for s in self.structures], marker="o", label=label)
        ax1.set_ylabel('Lattice Angles [degree]')
        ax1.legend(loc='best', shadow=True)

        ax2.plot(self.steps, [s.lattice.volume for s in self.structures], marker="o")
        ax2.set_ylabel('Lattice volume [A^3]')

        stress_cart_tensors, pressures = self.reader.read_cart_stress_tensors()
        ax3.plot(self.steps, pressures, marker="o", label="Pressure")
        ax3.set_ylabel('Pressure [GPa]')

        # Forces
        forces_hist = self.reader.read_cart_forces()
        fmin_steps, fmax_steps, fmean_steps, fstd_steps = [], [], [], []
        for step in range(self.num_steps):
            forces = forces_hist[step]
            fmods = np.sqrt([np.dot(force, force) for force in forces])
            fmean_steps.append(fmods.mean())
            fstd_steps.append(fmods.std())
            fmin_steps.append(fmods.min())
            fmax_steps.append(fmods.max())

        ax4.plot(self.steps, fmin_steps, marker="o", label="min |F|")
        ax4.plot(self.steps, fmax_steps, marker="o", label="max |F|")
        ax4.plot(self.steps, fmean_steps, marker="o", label="mean |F|")
        ax4.plot(self.steps, fstd_steps, marker="o", label="std |F|")
        ax4.set_ylabel('Force stats [eV/A]')
        ax4.legend(loc='best', shadow=True)
        ax4.set_xlabel('Step')

        # Total energy.
        ax5.plot(self.steps, self.etotals, marker="o", label="Energy")
        ax5.set_ylabel('Total energy [eV]')
        ax5.set_xlabel('Step')

        return fig

    @add_fig_kwargs
    def plot_energies(self, ax=None, **kwargs):
        """
        Plot the total energies as function of the iteration step.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        # TODO max force and pressure
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        terms = self.reader.read_eterms()
        for key, values in terms.items():
            if np.all(values == 0.0): continue
            ax.plot(self.steps, values, marker="o", label=key)

        ax.set_xlabel('Step')
        ax.set_ylabel('Energies [eV]')
        ax.grid(True)
        ax.legend(loc='best', shadow=True)

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("hist = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(hist)"),
            nbv.new_code_cell("fig = hist.plot_energies()"),
            nbv.new_code_cell("fig = hist.plot()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class HistReader(ETSF_Reader):
    """This object reads data from the HIST file."""

    @lazy_property
    def num_steps(self):
        """Number of iterations present in the HIST file."""
        return self.read_dimvalue("time")

    @lazy_property
    def natom(self):
        """Number of atoms un the unit cell"""
        return self.read_dimvalue("natom")

    def read_all_structures(self):
        """Return the list of structures at the different iteration steps."""
        rprimd_list = self.read_value("rprimd")
        xred_list = self.read_value("xred")

        # Alchemical mixing is not supported.
        num_pseudos = self.read_dimvalue("npsp")
        ntypat = self.read_dimvalue("ntypat")
        if num_pseudos != ntypat:
            raise NotImplementedError("Alchemical mixing is not supported, num_pseudos != ntypat")

        znucl, typat = self.read_value("znucl"), self.read_value("typat")
        #print(znucl.dtype, typat)
        cart_forces_step = self.read_cart_forces(unit="eV ang^-1")

        structures = []
        for step in range(self.num_steps):
            s = Structure.from_abivars(
                xred=xred_list[step],
                rprim=rprimd_list[step],
                acell=3 * [1.0],
                # FIXME ntypat, typat, znucl are missing!
                znucl=znucl,
                typat=typat,
            )
            s.add_site_property("cartesian_forces", cart_forces_step[step])
            structures.append(s)

        return structures

    def read_eterms(self, unit="eV"):
        return AttrDict(
            etotals=EnergyArray(self.read_value("etotal"), "Ha").to(unit),
            kinetic_terms=EnergyArray(self.read_value("ekin"), "Ha").to(unit),
            entropies=EnergyArray(self.read_value("entropy"), "Ha").to(unit),
        )

    def read_cart_forces(self, unit="eV ang^-1"):
        """
        Read and return a numpy array with the cartesian forces in unit `unit`.
        Shape (num_steps, natom, 3)
        """
        return ArrayWithUnit(self.read_value("fcart"), "Ha bohr^-1").to(unit)

    def read_reduced_forces(self):
        """
        Read and return a numpy array with the forces in reduced coordinates
        Shape (num_steps, natom, 3)
        """
        return self.read_value("fred")

    def read_cart_stress_tensors(self):
        """
        Return the stress tensors (nstep x 3 x 3) in cartesian coordinates (Hartree/Bohr^3)
        and the list of pressures in GPa unit.
        """
        # Abinit stores 6 unique components of this symmetric 3x3 tensor:
        # Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
        c = self.read_value("strten")
        tensors = np.empty((self.num_steps, 3, 3), dtype=np.float)

        for step in range(self.num_steps):
            for i in range(3): tensors[step, i,i] = c[step, i]
            for p, (i, j) in enumerate(((2,1), (2,0), (1,0))):
                tensors[step, i,j] = c[step, 3+p]
                tensors[step, j,i] = c[step, 3+p]

        HaBohr3_GPa = 29421.033 # 1 Ha/Bohr^3, in GPa
        pressures = np.empty(self.num_steps)
        for step, tensor in enumerate(tensors):
            pressures[step] = - (HaBohr3_GPa/3) * tensor.trace()

        return tensors, pressures

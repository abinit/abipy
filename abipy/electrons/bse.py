# coding: utf-8
"""Classes for the analysis of BSE calculations"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import itertools
import collections
import numpy as np

from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.string import marquee # is_string, list_strings,
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from abipy.core.func1d import Function1D
from abipy.core.kpoints import Kpoint, KpointList
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.core.tensor import SymmetricTensor
from abipy.iotools import ETSF_Reader


__all__ = [
    "DielectricTensor",
    "DielectricFunction",
    "MdfFile",
    "MdfReader",
    "MdfPlotter",
]


class DielectricTensor(object):
    """
    This object stores the frequency-dependent macroscopic dielectric tensor
    obtained from the dielectric functions for different q-directions.
    """
    def __init__(self, mdf, structure):
        nfreq = len(mdf.wmesh)

        self._wmesh = mdf.wmesh

        # Transform mdf emacros_q to numpy array
        all_emacros = []
        for emacro in mdf.emacros_q:
            all_emacros.append(emacro.values)

        all_emacros = np.array(all_emacros)

        # One tensor for each frequency
        all_tensors = []
        for (ifrq, freq) in enumerate(mdf.wmesh):
            tensor = SymmetricTensor.from_directions(mdf.qfrac_coords, all_emacros[:,ifrq],
                                                     structure.lattice.reciprocal_lattice, space="g")
            all_tensors.append(tensor)

        self._all_tensors = all_tensors

    def to_array(self, red_coords=True):

        table = []
        for tensor in self._all_tensors:
            if red_coords:
                table.append(tensor.reduced_tensor)
            else:
                table.append(tensor.cartesian_tensor)

        return np.array(table)

    def symmetrize(self, structure):

        for tensor in self._all_tensors:
            tensor.symmetrize(structure)

    def to_func1d(self, red_coords=True):

        table = self.to_array(red_coords)

        all_funcs = []

        for i in np.arange(3):
            for j in np.arange(3):
                all_funcs.append(Function1D(self._wmesh, table[:,i,j]))

        return all_funcs

    @add_fig_kwargs
    def plot(self, ax=None, *args, **kwargs):
        """
        Plot all the components of the tensor

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        red_coords      True to plot the reduced coordinate tensor (Default: True)
        ==============  ==============================================================

        Returns:
            matplotlib figure
        """
        red_coords = kwargs.pop("red_coords", True)
        ax, fig, plt = get_ax_fig_plt(ax)

        ax.grid(True)
        ax.set_xlabel('Frequency [eV]')
        ax.set_ylabel('Dielectric tensor')

        #if not kwargs:
        #    kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the 6 independent components
        for icomponent in [0,4,8,1,2,5]:
            self.plot_ax(ax, icomponent, red_coords, *args, **kwargs)

        return fig

    def plot_ax(self, ax, what, red_coords, *args, **kwargs):
        """
        Helper function to plot data on the axis ax.

        Args:
            ax: plot axis
            what: Sequential index of the tensor matrix element.
            args: Positional arguments passed to ax.plot
            kwargs: Keyword arguments passed to matplotlib. Accepts also:

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        cplx_mode:      string defining the data to print (case-insensitive).
                        Possible choices are:

                            - "re"  for real part
                            - "im" for imaginary part only.
                            - "abs' for the absolute value

                        Options can be concated with "-".
        ==============  ==============================================================
        """
        # Extract the function to plot according to qpoint.
        if isinstance(what, int):
            f = self.to_func1d(red_coords)[what]
        else:
            raise ValueError("Don't know how to handle %s" % str(what))

        return f.plot_ax(ax, *args, **kwargs)


class DielectricFunction(object):
    """
    This object stores the frequency-dependent macroscopic dielectric function
    computed for different q-directions in reciprocal space.

    .. note:
        Frequencies are in eV
    """

    def __init__(self, structure, qpoints, wmesh, emacros_q, info):
        """
        Args:
            structure: :class: Structure object.
            qpoints: :class:`KpointList` with the qpoints in reduced coordinates.
            wmesh: Array-like object with the frequency mesh (eV).
            emacros_q: Iterable with the macroscopic dielectric function for the different q-points.
            info: Dictionary containing info on the calculation that produced
                  the results (read from file). It must contain the following keywords:

                    - "lfe": True if local field effects are included.
                    - "calc_type": string defining the calculation type.

        """
        self.wmesh = np.array(wmesh)
        self.qpoints = qpoints
        assert len(self.qpoints) == len(emacros_q)
        self.info = info

        self.emacros_q, em_avg = [], np.zeros(len(wmesh), dtype=np.complex)
        for emq in emacros_q:
            em_avg += emq
            self.emacros_q.append(Function1D(wmesh, emq))
        self.emacros_q = tuple(self.emacros_q)

        # Compute the average value.
        # TODO: One should take into account the star of q, but I need the symops
        self.emacro_avg = Function1D(wmesh, em_avg / self.num_qpoints)

    def __str__(self):
        return self.__class__.__name__

    def __iter__(self):
        """Iterate over (q, em_q)."""
        return itertools.izip(self.qpoints, self.emacros_q)

    @property
    def num_qpoints(self):
        return len(self.qpoints)

    @property
    def qfrac_coords(self):
        """The fractional coordinates of the q-points as a ndarray."""
        return self.qpoints.frac_coords

    @property
    def has_lfe(self):
        """True if MDF includes local field effects."""
        return bool(self.info["lfe"])

    @property
    def calc_type(self):
        """String with the type of calculation."""
        return self.info["calc_type"]

    def show_info(self, stream=sys.stdout):
        """Pretty print of the info."""
        import pprint
        printer = pprint.PrettyPrinter(self, width=80, depth=None, stream=stream)
        printer.pprint(self.info)

    @add_fig_kwargs
    def plot(self, ax=None, **kwargs):
        """
        Plot the MDF.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        only_mean       True if only the averaged spectrum is wanted (default True)
        ==============  ==============================================================

        Returns:
            matplotlib figure
        """
        only_mean = kwargs.pop("only_mean", True)

        ax, fig, plt = get_ax_fig_plt(ax)

        ax.grid(True)
        ax.set_xlabel('Frequency [eV]')
        ax.set_ylabel('Macroscopic DF')

        #if not kwargs:
        #    kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the average value
        self.plot_ax(ax, qpoint=None, **kwargs)

        if not only_mean:
            # Plot the q-points
            for iq, qpoint in enumerate(self.qpoints):
                self.plot_ax(ax, iq, **kwargs)

        return fig

    def plot_ax(self, ax, qpoint=None, **kwargs):
        """
        Helper function to plot data on the axis ax.

        Args:
            ax: plot axis.
            qpoint: index of the q-point or Kpoint object or None) to plot emacro_avg.
            kwargs: Keyword arguments passed to matplotlib. Accepts also:

                cplx_mode:
                    string defining the data to print (case-insensitive).
                    Possible choices are

                        - "re"  for real part
                        - "im" for imaginary part only.
                        - "abs' for the absolute value

                    Options can be concated with "-".
        """
        # Extract the function to plot according to qpoint.
        if isinstance(qpoint, int):
            f = self.emacros_q[qpoint]

        elif isinstance(qpoint, Kpoint):
            iq = self.qpoints.index(qpoint)
            f = self.emacros_q[iq]

        elif qpoint is None:
            f = self.emacro_avg

        else:
            raise ValueError("Don't know how to handle %s" % str(qpoint))

        return f.plot_ax(ax, **kwargs)


class MdfFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    Usage example:

    .. code-block:: python

        with MdfFile("foo_MDF.nc") as mdf:
            mdf.plot_mdfs()
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(MdfFile, self).__init__(filepath)
        self.reader = MdfReader(filepath)

        # TODO Add electron Bands.
        #self._ebands = r.read_ebands()

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(str(self.structure))

        app(marquee("Q-points", mark="="))
        app(str(self.qpoints))

        return "\n".join(lines)

    def close(self):
        self.reader.close()

    @lazy_property
    def structure(self):
        """Returns the `Structure` object."""
        return self.reader.read_structure()

    @lazy_property
    def exc_mdf(self):
        "Excitonic macroscopic dieletric function."""
        return self.reader.read_exc_mdf()

    @lazy_property
    def rpanlf_mdf(self):
        """RPA dielectric function without local-field effects."""
        return self.reader.read_rpanlf_mdf()

    @lazy_property
    def gwnlf_mdf(self):
        """RPA-GW dielectric function without local-field effects."""
        return self.reader.read_gwnlf_mdf()

    @property
    def qpoints(self):
        return self.reader.qpoints

    @property
    def qfrac_coords(self):
        """The fractional coordinates of the q-points as a ndarray."""
        return self.qpoints.frac_coords

    @lazy_property
    def params(self):
        """
        Dictionary with the parameters that are usually tested for convergence.
        Used to build Pandas dataframes in Robots.
        """
        return self.reader.read_params()

    def get_mdf(self, mdf_type="exc"):
        """"Returns the macroscopic dielectric function."""
        d = {"exc": self.exc_mdf,
             "rpa": self.rpanlf_mdf,
             "gwrpa": self.gwnlf_mdf}

        try:
            return d[mdf_type.lower()]
        except KeyError:
            raise ValueError("Wrong value for mdf_type: %s" % mdf_type)

    def plot_mdfs(self, cplx_mode="Im", mdf_type="all", qpoint=None, **kwargs):
        """
        Plot the macroscopic dielectric function.

        Args:
            cplx_mode:
                string defining the data to print (case-insensitive).
                Possible choices are

                    - "re"  for real part
                    - "im" for imaginary part only.
                    - "abs' for the absolute value

                Options can be concated with "-".

            mdf_type:
                Select the type of macroscopic dielectric function.
                Possible choices are

                    - "exc" for the excitonic MDF.
                    - "rpa" for RPA MDF.
                    - "gwrpa" for GW-RPA MDF
                    - "all" if all types are wanted.

                Options can be concated with "-".

            qpoint:
                index of the q-point or Kpoint object or None to plot emacro_avg.
        """
        mdf_type, cplx_mode = mdf_type.lower(), cplx_mode.lower()

        plot_all = mdf_type == "all"
        mdf_type = mdf_type.split("-")

        # Build the plotter.
        plotter = MdfPlotter()

        # Excitonic MDF.
        if "exc" in mdf_type or plot_all:
            plotter.add_mdf("EXC", self.exc_mdf)

        # KS-RPA MDF
        if "rpa" in mdf_type or plot_all:
            plotter.add_mdf("KS-RPA", self.rpanlf_mdf)

        # GW-RPA MDF (obtained with the scissors operator).
        if "gwrpa" in mdf_type or plot_all:
            plotter.add_mdf("GW-RPA", self.gwnlf_mdf)

        # Plot spectra
        return plotter.plot(cplx_mode=cplx_mode, qpoint=qpoint, **kwargs)

    def get_tensor(self, mdf_type="exc"):
        """Get the macroscopic dielectric tensor from the MDF."""
        return DielectricTensor(self.get_mdf(mdf_type), self.structure)

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("mdf_file = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(mdf_file)"),
            nbv.new_code_cell("fig = mdf_file.plot_mdfs(cplx_mode='Re')"),
            nbv.new_code_cell("fig = mdf_file.plot_mdfs(cplx_mode='Im')"),
            # TODO:
            #nbv.new_code_cell("tensor_exc = mdf_file.get_tensor("exc")")
            #tensor_exc.symmetrize(mdf_file.structure)
            #tensor_exc.plot(title=title)
        ])

        return self._write_nb_nbpath(nb, nbpath)


# TODO Add band energies to MDF file.
#from abipy.electrons import ElectronsReader
class MdfReader(ETSF_Reader): #ElectronsReader
    """
    This object reads data from the MDF.nc file produced by ABINIT.
    """
    def __init__(self, path):
        """Initialize the object from a filename."""
        super(MdfReader, self).__init__(path)
        # Read the structure here to facilitate the creation of the other objects.
        self._structure = self.read_structure()

    @property
    def structure(self):
        return self._structure

    @lazy_property
    def qpoints(self):
        """List of q-points (ndarray)."""
        # Read the fractional coordinates and convert them to KpointList.
        return KpointList(self.structure.reciprocal_lattice, frac_coords=self.read_value("qpoints"))

    @lazy_property
    def wmesh(self):
        """The frequency mesh in eV."""
        return self.read_value("wmesh")

    def read_params(self):
        """Dictionary with the parameters of the run."""
        # TODO
        keys = [
            "nsppol", "ecutwfn", "ecuteps",
            "eps_inf", "soenergy", "broad", "nkibz", "nkbz", "nkibz_interp", "nkbz_interp",
            "wtype", "interp_mode", "nreh", "lomo_spin", "humo_spin"
        ]
        return self.read_keys(keys)

    def _read_mdf(self, mdf_type):
        """Read the MDF from file, returns numpy complex array."""
        return self.read_value(mdf_type, cmode="c")

    def read_exc_mdf(self):
        """Returns the excitonic MDF."""
        info = self.read_params()
        emacros_q = self._read_mdf("exc_mdf")
        return DielectricFunction(self.structure, self.qpoints, self.wmesh, emacros_q, info)

    def read_rpanlf_mdf(self):
        """Returns the KS-RPA MDF without LF effects."""
        info = self.read_params()
        emacros_q = self._read_mdf("rpanlf_mdf")
        return DielectricFunction(self.structure, self.qpoints, self.wmesh, emacros_q, info)

    def read_gwnlf_mdf(self):
        """Returns the GW-RPA MDF without LF effects."""
        info = self.read_params()
        emacros_q = self._read_mdf("gwnlf_mdf")
        return DielectricFunction(self.structure, self.qpoints, self.wmesh, emacros_q, info)


class MdfPlotter(object):
    """
    Class for plotting multiple MDFs.

    Usage example:

    .. code-block:: python

        plotter = MdfPlotter()
        plotter.add_mdf_from_file("foo_MDF.nc", label="foo mdf")
        plotter.add_mdf_from_file("bar_MDF.nc", label="bar mdf")
        plotter.plot()
    """
    def __init__(self):
        self._mdfs = collections.OrderedDict()

    def add_mdf(self, label, mdf):
        """
        Adds a :class:`DielectricFunction` for plotting.

        Args:
            name: name for the MDF. Must be unique.
            mdf: :class:`DielectricFunction` object.
        """
        if label in self._mdfs:
            raise ValueError("name %s is already in %s" % (label, self._mdfs.keys()))

        self._mdfs[label] = mdf

    def add_mdf_from_file(self, filepath, mdf_type="exc", label=None):
        """
        Adds a mdf for plotting. Reads data from file filepaths.

        Args:
            mdf_type: String defining the type of mdf.
            name: Optional string used to name the plot.
        """
        from abipy.abilab import abiopen
        with abiopen(filepath) as ncfile:
            mdf = ncfile.get_mdf(mdf_type=mdf_type)

        if label is None:
            label = mdf_type + ncfile.filepath
        self.add_mdf(label, mdf)

    @add_fig_kwargs
    def plot(self, ax=None, cplx_mode="Im", qpoint=None, **kwargs):
        """
        Get a matplotlib plot showing the MDFs.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            cplx_mode: string defining the data to print (case-insensitive).
                Possible choices are `re` for the real part, `im` for imaginary part only. `abs` for the absolute value.
                Options can be concated with "-".
            qpoint: index of the q-point or :class:`Kpoint` object or None to plot emacro_avg.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        xlim            x-axis limits. None (Default) for automatic determination.
        ylim            y-axis limits. None (Default) for automatic determination.
        ==============  ==============================================================
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        ax.grid(True)

        xlim = kwargs.pop("xlim", None)
        if xlim is not None: ax.set_xlim(xlim)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None: ax.set_ylim(ylim)

        ax.set_xlabel('Frequency [eV]')
        ax.set_ylabel('Macroscopic DF')

        cmodes = cplx_mode.split("-")
        qtag = "avg" if qpoint is None else repr(qpoint)

        lines, legends = [], []
        for label, mdf in self._mdfs.items():
            # Plot the q-points
            #for (iq, qpoint) in enumerate(self.qpoints):
            #    self.plot_ax(ax, iq, **kwargs)

            for cmode in cmodes:
                # Plot the average value
                l = mdf.plot_ax(ax, qpoint, cplx_mode=cmode, **kwargs)[0]
                lines.append(l)
                legends.append("%s: %s, %s $\,\\varepsilon$" % (cmode, qtag, label))

        # Set legends.
        ax.legend(lines, legends, loc='best', shadow=False)
        return fig


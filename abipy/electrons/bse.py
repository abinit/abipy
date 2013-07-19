from __future__ import print_function, division

import sys
import os
import itertools
import collections
import numpy as np

from abipy.core.func1d import Function1D
from abipy.iotools import ETSF_Reader, AbinitNcFile
from abipy.kpoints import Kpoint, kpoints_factory

__all__ = [
    "DielectricFunction",
    "MDF_File",
    "MDF_Reader",
    "MDF_Plotter",
]

#class DielectricTensor(object):
#    def __init__(self):
#    @classmethod
#    def from_file(cls, filobj):
#    def plot(self, *args, **kwargs):
#    def make_dielectric_function(self, qpoints):
#    def plot_ax(self, ax, what, *args, **kwargs):

#########################################################################################


class DielectricFunction(object):
    """
    This object stores the frequency-dependent macroscopic dielectric function
    computed for different q-directions in reciprocal space.

    .. note:
        Frequencies are in eV
    """

    def __init__(self, qpoints, wmesh, emacros_q, info):
        """
        Args:
            qpoints:
                List of qpoints in reduced coordinates.
            wmesh:
                Array-like object with the frequency mesh (eV).
            emacros_q:
                Iterable with the macroscopic dielectric function for the
                different q-points.
            info:
                Dictionary containing info on the calculation that produced
                the results (read from file). It must contain the following keywords:

                    - "lfe": True if local field effects are included.
                    - "calc_type": string defining the calculation type.

        """
        self.wmesh = np.array(wmesh)
        self.qpoints = np.reshape(qpoints, (-1, 3))
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

    def __iter__(self):
        """Iterate over (q, em_q)."""
        return itertools.izip(self.qpoints, self.emacros_q)

    @property
    def num_qpoints(self):
        return len(self.qpoints)

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

    def raw_print(self, stream=sys.stdout, fmt=None, delimiter=' '):
        """
        Write data on stream with format fmt. See also `numpy.savetxt`.

        Args:
            stream:
                filename or file handle. If the filename ends in .gz, the file is automatically
                saved in compressed gzip format.
            fmt: 
                str or sequence of strings, optional. A single format (%10.5f), a sequence of formats,
                or a multi-format string,
            delimiter:
                Character separating columns.
        """
        header = \
            """
            2 * (num_qpoints+1) columns representing num_qpoints+1 complex numbers (re, im).
            omega_re omega_im em_q[0]_re em_q[0]_im ... em_q[nq-1]_im
            """
        # Build table.
        table = []
        for (iw, omega) in enumerate(self.wmesh):
            line = [omega] + [em.values[iw] for em in self.emacros_q]
            table.append(line)

        if fmt is None: fmt = (1 + self.num_qpoints) * ['%.4f %.4f']

        np.savetxt(stream, table, fmt=fmt, delimiter=delimiter, header=header)

    def plot(self, *args, **kwargs):
        """
        Plot the MDF.

        args:
            Optional arguments passed to :mod:`matplotlib`.


        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            matplotlib figure
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)

        if title is not None:
            ax.set_title(title)

        ax.grid(True)
        ax.set_xlabel('Frequency [eV]')
        ax.set_ylabel('Macroscopic DF')

        #if not kwargs:
        #    kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the q-points
        for (iq, qpoint) in enumerate(self.qpoints):
            self.plot_ax(ax, iq, *args, **kwargs)

        # Plot the average value
        self.plot_ax(ax, "average", *args, **kwargs)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def plot_ax(self, ax, qpoint, *args, **kwargs):
        """
        Helper function to plot data on the axis ax.

        Args:
            ax:
                plot axis
            qpoint:
                index of the q-point or Kpoint object or "average" to plot emacro_avg.
            args:
                Positional arguments passed to ax.plot
            kwargs:
                Keyword arguments passed to matplotlib. Accepts also:

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

        elif isinstance(qpoint, str) and qpoint == "average":
            f = self.emacro_avg

        else:
            raise ValueError("Don't know how to handle %s" % str(qpoint))

        return f.plot_ax(ax, *args, **kwargs)

#########################################################################################

class MDF_File(AbinitNcFile):
    def __init__(self, filepath):
        super(MDF_File, self).__init__(filepath)

        with MDF_Reader(filepath) as r:
            self.structure = r.read_structure()

            self.exc_mdf = r.read_exc_mdf()
            self.rpanlf_mdf = r.read_rpanlf_mdf()
            self.gwnlf_mdf = r.read_gwnlf_mdf()

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def get_structure(self):
        """Returns the `Structure` object."""
        return self.structure

    def get_mdf(self, mdf_type="exc"):

        if mdf_type == "exc":
            return self.exc_mdf
        elif mdf_type == "rpa":
            return self.rpanlf_mdf
        elif mdf_type == "gwrpa":
            return self.gwnlf_mdf
        else:
            raise ValueError("Wrong value for mdf_type %s" % mdf_type)

    def plot_mdfs(self, cplx_mode="Im", mdf_select="all", **kwargs):
        """Plot the macroscopic dielectric function."""
        plot_all = mdf_select == "all"
        mdf_select = mdf_select.split("-")

        # Build the plotter.
        plotter = MDF_Plotter()

        # Excitonic MDF.
        if "exc" in mdf_select or plot_all:
            plotter.add_mdf("EXC", self.exc_mdf)

        # KS-RPA MDF
        if "rpa" in mdf_select or plot_all:
            plotter.add_mdf("KS-RPA", self.rpanlf_mdf)

        # GW-RPA MDF (obtained with the scissors operator).
        if "gwrpa" in mdf_select or plot_all:
            plotter.add_mdf("GW-RPA", self.gwnlf_mdf)

        # Plot spectra 
        plotter.plot(cplx_mode=cplx_mode, **kwargs)

#########################################################################################

class MDF_Reader(ETSF_Reader):
    """
    This object reads data from the MDF.nc file produced by ABINIT.
    """

    def __init__(self, path):
        """Initialize the object from a filename."""
        super(MDF_Reader, self).__init__(path)

    def _lazy_get(self, varname):
        """Helper function used to create lazy properties."""
        hiddename = "__" + varname
        try:
            return getattr(self, hiddename)
        except AttributeError:
            setattr(self, hiddename, self.read_value(varname))
            return self._lazy_get(varname)

    @property
    def qpoints(self):
        """List of q-points (ndarray)."""
        return self._lazy_get("qpoints")

    @property
    def wmesh(self):
        """The frequency mesh in eV."""
        return self._lazy_get("wmesh")

    def read_run_params(self):
        """Dictionary with the parameters of the run."""
        return {}
        # TODO
        #try:
        #    return self._run_params
        #except AttributeError:
        #    self._run_params = params = {}
        #    # Fill the dictionary with basic parameters of the run.
        #    return self._run_params.copy()

    def _read_mdf(self, mdf_type):
        """Read the MDF from file, returns numpy complex array."""
        return self.read_value(mdf_type, cmode="c")

    def read_exc_mdf(self):
        """Returns the excitonic MDF."""
        info = self.read_run_params()
        emacros_q = self._read_mdf("exc_mdf")
        return DielectricFunction(self.qpoints, self.wmesh, emacros_q, info)

    def read_rpanlf_mdf(self):
        """Returns the KS-RPA MDF without LF effects."""
        info = self.read_run_params()
        emacros_q = self._read_mdf("rpanlf_mdf")
        return DielectricFunction(self.qpoints, self.wmesh, emacros_q, info)

    def read_gwnlf_mdf(self):
        """Returns the GW-RPA MDF without LF effects."""
        info = self.read_run_params()
        emacros_q = self._read_mdf("gwnlf_mdf")
        return DielectricFunction(self.qpoints, self.wmesh, emacros_q, info)

#########################################################################################


class MDF_Plotter(object):
    """
    Class for plotting MDFs.
    """

    def __init__(self):
        self._mdfs = collections.OrderedDict()

    def add_mdf(self, label, mdf):
        """
        Adds a MDF for plotting.

        Args:
            label:
                label for the MDF. Must be unique.
            mdf:
                MacroscopicDielectricFunction object.
        """
        if label in self._mdfs:
            raise ValueError("label %s is already in %s" % (label, self._mdfs.keys()))

        self._mdfs[label] = mdf

    def add_mdf_from_file(self, filepath, mdf_type="exc", label=None):
        """
        Adds a mdf for plotting. Reads data from file filepaths.

        Args:
            mdf_type:
                String defining the type of mdf.
            label:
                Optional string used to label the plot.
        """
        from abipy import abiopen

        ncfile = abiopen(filepath)
        mdf = ncfile.get_mdf(mdf_type=mdf_type)
        if label is None:
            label = mdf_type + ncfile.filepath

        self.add_mdf(label, mdf)

    def plot(self, cplx_mode="Im", *args, **kwargs):
        """
        Get a matplotlib plot showing the MDFs.
        
        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig:        'abc.png' or 'abc.eps'* to save the figure to a file.
        xlim            x-axis limits. None (Default) for automatic determination.
        ylim            y-axis limits. None (Default) for automatic determination.
        ==============  ==============================================================
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        ax.grid(True)

        xlim = kwargs.pop("xlim", None)
        if xlim is not None: ax.set_xlim(xlim)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None: ax.set_ylim(ylim)

        ax.set_xlabel('Frequency [eV]')
        ax.set_ylabel('Macroscopic DF')

        if title is not None:
            ax.set_title(title)

        cmodes = cplx_mode.split("-")

        lines, legends = [], []
        for (label, mdf) in self._mdfs.items():

            # Plot the q-points
            #for (iq, qpoint) in enumerate(self.qpoints):
            #    self.plot_ax(ax, iq, *args, **kwargs)

            for cmode in cmodes:
                # Plot the average value
                l = mdf.plot_ax(ax, "average", *args, cplx_mode=cmode, **kwargs)[0]

                lines.append(l)
                legends.append("%s: %s $\,\\varepsilon$" % (cmode, label))

        # Set legends.
        ax.legend(lines, legends, 'best', shadow=True)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


class DIPME_File(object):
    """
    This object provides tools to analyze the dipole matrix elements produced by the BSE code.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

        # Save useful quantities
        with OME_Reader(self.path) as reader:
            self.structure = reader.structure
            self.nsppol = reader.nsppol
            self.kibz = reader.kibz
            self.minb_sk = reader.minb_sk
            self.maxb_sk = reader.maxb_sk

            # Dipole matrix elements.
            # opt_cvk(minb:maxb,minb:maxb,nkbz,Wfd%nsppol)
            self.dipme_scvk = reader.read_dipme_scvk()

    @classmethod
    def from_file(cls, path):
        return cls(path)

    def kpoint_index(self, kpoint):
        """The index of the kpoint"""
        return self.kibz.find(kpoint)

    def plot(self, qpoint=None, spin=None, kpoints=None, color_map=None, **kwargs):
        """
        Plot the dipole matrix elements.

        Args:
            qpoint:
                The qpoint for the optical limit.
                if qpoint is None, we plot  |<k|r|k>| else |<k|q.r|k>|
            spin:
                spin index. None if all spins are wanted
            kpoints:
                List of Kpoint objects, None if all k-points are wanted.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig:        'abc.png' or 'abc.eps'* to save the figure to a file.
        colormap        matplotlib colormap, see link below.
        ==============  ==============================================================

        Returns:
            matplotlib figure.

        .. see: 
            http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)
        color_map = kwargs.pop("color_map", None)

        if qpoint is not None:
            # Will compute scalar product with q
            qpoint = Kpoint.askpoint(qpoint, self.structure.reciprocal_lattice).versor()
        else:
            # Will plot |<psi|r|psi>|.
            qpoint = Kpoint((1, 1, 1), self.structure.reciprocal_lattice)

        if spin in None:
            spins = range(self.nsppol)
        else:
            spins = [spin]

        if kpoints is None:
            kpoints = self.ibz

        # Extract the matrix elements for the plot.
        from abipy.tools.plotting_utils import ArrayPlotter
        plotter = ArrayPlotter()
        for spin in spins:
            for kpoint in kpoints:
                ik = self.kpoint_index(kpoint)
                rme = self.dipme_scvk[spin, ik, :, :, :]

                #qrme = qpoint * rme
                label = "qpoint %s, spin %s, kpoint = %s" % (qpoint, spin, kpoint)
                plotter.add_array(label, rme)

        # Plot matrix elements and return matplotlib figure.
        return plotter.plot(title=title, color_map=color_map, show=show, savefig=savefig, **kwargs)


class DIPME_Reader(ETSF_Reader):
    """"
    This object reads the optical matrix elements from the OME.nc file.
    """

    def __init__(self, path):
        super(DIPME_Reader, self).__init__(path)

        # Read important dimensions and variables.
        frac_coords = self.read_value("kpoints_reduced_coordinates")
        self.kibz = kpoints_factory(self)

        # Minimum and maximum band index as function of [s,k]
        self.minb_ks = self.read_value("minb")
        self.maxb_ks = self.read_value("maxb")

        # Dipole matrix elements
        self.dipme_skvc = self.read_value("dipme", cmode="c")

    def kpoint_index(self, kpoint):
        return self.kibz.find(kpoint)

    def read_dipme(self, spin, kpoint):
        """Read the dipole matrix elements."""
        ik = self.kpoint_index(kpoint)
        return self.dipme_skvc[spin, ik, :, :, :]


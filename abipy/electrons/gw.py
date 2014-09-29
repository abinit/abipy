"""Classes for the analysis of GW calculations."""
from __future__ import print_function, division, unicode_literals

import sys
import copy
import collections
import warnings
import numpy as np

from monty.string import list_strings, is_string
from monty.collections import  AttrDict
from six.moves import cStringIO
from abipy.tools import find_le, find_ge, pprint_table
from abipy.core.func1d import Function1D
from abipy.core.kpoints import KpointList
from abipy.iotools import AbinitNcFile, ETSF_Reader, Has_Structure, Has_ElectronBands
from abipy.electrons.ebands import ElectronBands
from abipy.electrons.scissors import Scissors

__all__ = [
    "QPState",
    "SIGRES_File",
    "SIGRES_Plotter",
]


class QPState(collections.namedtuple("QPState", "spin kpoint band e0 qpe qpe_diago vxcme sigxme sigcmee0 vUme ze0")):
    """
    Quasi-particle result for given (spin, kpoint, band).

    .. Attributes:

        spin:
            spin index (C convention, i.e >= 0)
        kpoint:
            `Kpoint` object.
        band:
            band index. (C convention, i.e >= 0).
        e0:
            Initial KS energy.
        qpe:
            Quasiparticle energy (complex) computed with the perturbative approach.
        qpe_diago:
            Quasiparticle energy (real) computed by diagonalizing the self-energy.
        vxcme:
            Matrix element of vxc[n_val] with nval the valence charge density.
        sigxme:
            Matrix element of Sigma_x.
        sigcmee0:
            Matrix element of Sigma_c(e0) with e0 being the KS energy.
        vUme:
            Matrix element of the vU term of the LDA+U Hamiltonian.
        ze0:
            Renormalization factor computed at e=e0.

    .. note:: Energies are in eV.
    """
    @property
    def qpeme0(self):
        """E_QP - E_0"""
        return self.qpe - self.e0

    @property
    def skb(self):
        """Tuple with (spin, kpoint, band)"""
        return self.spin, self.kpoint, self.band

    def copy(self):
        d = {f: copy.copy(getattr(self, f)) for f in self._fields}
        return QPState(**d)

    @classmethod
    def get_fields(cls, exclude=()):
        fields = list(cls._fields) + ["qpeme0"]
        for e in exclude:
            fields.remove(e)
        return tuple(fields)

    def _asdict(self):
        od = super(QPState, self)._asdict()
        od["qpeme0"] = self.qpeme0
        return od

    def to_strdict(self, fmt=None):
        """Ordered dictionary mapping fields --> strings."""
        d = self._asdict()
        for (k, v) in d.items():
            if np.iscomplexobj(v):
                if abs(v.imag) < 1.e-3:
                    d[k] = "%.2f" % v.real
                else:
                    d[k] = "%.2f%+.2fj" % (v.real, v.imag)
            elif isinstance(v, int):
                d[k] = "%d" % v
            else:
                try:
                    d[k] = "%.2f" % v
                except TypeError as exc:
                    #print("k", k, str(exc))
                    d[k] = str(v)
        return d

    @property
    def tips(self):
        """Bound method of self that returns a dictionary with the description of the fields."""
        return self.__class__.TIPS()
    
    @classmethod
    def TIPS(cls):
        """
        Class method that returns a dictionary with the description of the fields.
        The string are extracted from the class doc string.
        """
        try:
            return cls._TIPS

        except AttributeError:
            # Parse the doc string.
            cls._TIPS = _TIPS = {}
            lines = cls.__doc__.splitlines()

            for i, line in enumerate(lines):
                if line.strip().startswith(".. Attributes"):
                    lines = lines[i+1:]
                    break

            def num_leadblanks(string):
                """Returns the number of the leading whitespaces."""
                return len(string) - len(string.lstrip())

            for field in cls._fields:
                for i, line in enumerate(lines):

                    if line.strip().startswith(field + ":"):
                        nblanks = num_leadblanks(line)
                        desc = []
                        for s in lines[i+1:]:
                            if nblanks == num_leadblanks(s) or not s.strip():
                                break
                            desc.append(s.lstrip())

                        _TIPS[field] = "\n".join(desc)

            diffset = set(cls._fields) - set(_TIPS.keys())
            if diffset:
                raise RuntimeError("The following fields are not documented: %s" % str(diffset))

            return _TIPS


class QPList(list):
    """A list of quasiparticle corrections."""

    def __init__(self, *args, **kwargs):
        super(QPList, self).__init__(*args)
        self.is_e0sorted = kwargs.get("is_e0sorted", False)

    def __repr__(self):
        return "<%s at %s, len=%d>" % (self.__class__.__name__, id(self), len(self))

    def __str__(self):
        """String representation."""
        table = self.to_table()

        strio = cStringIO()
        pprint_table(table, out=strio)
        strio.write("\n")
        strio.seek(0)
        return "".join(strio)

    def copy(self):
        """Copy of self."""
        return QPList([qp.copy() for qp in self], is_e0sorted=self.is_e0sorted)

    def sort_by_e0(self):
        """Return a new object with the E0 energies sorted in ascending order."""
        return QPList(sorted(self, key=lambda qp: qp.e0), is_e0sorted=True)

    def get_e0mesh(self):
        """Return the E0 energies."""
        if not self.is_e0sorted:
            raise ValueError("QPState corrections are not sorted. Use sort_by_e0")

        return np.array([qp.e0 for qp in self])

    def get_field(self, field):
        """ndarray containing the values of field."""
        return np.array([getattr(qp, field) for qp in self])

    def get_value(self, skb_tup, field):
        """ return the value of field for the given spin kp band tuple, None if not found"""
        for qp in self:
            if qp.skb == skb_tup:
                return getattr(qp, field)
        return None

    def get_qpenes(self):
        """Return an array with the QPState energies."""
        return self.get_field("qpe")

    def get_qpeme0(self):
        """Return an arrays with the QPState corrections."""
        return self.get_field("qpeme0")

    def to_table(self):
        """Return a table (list of list of strings)."""
        header = QPState.get_fields(exclude=["spin", "kpoint"])
        table = [header]

        for qp in self:
            d = qp.to_strdict(fmt=None)
            table.append([d[k] for k in header])

        return table

    def plot_qps_vs_e0(self, with_fields="all", exclude_fields=None, **kwargs):
        """
        Args:
            with_fields:
                The names of the qp attributes to plot as function of e0.
                Accepts:
                    List of strings or string with tokens separated by blanks.
                    See `QPState` for the list of available fields.
            args:
                Positional arguments passed to :mod:`matplotlib`.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure.
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)
        fermi = kwargs.pop("fermi", None)

        if is_string(with_fields):
            if with_fields == "all":
                fields = list(QPState.get_fields(exclude=["spin", "kpoint"]))
            else:
                fields = with_fields.split()

        if exclude_fields:
            if is_string(exclude_fields):
                exclude_fields = exclude_fields.split()
            for e in exclude_fields:
                fields.remove(e)

        # Build grid of plots.
        num_plots, ncols, nrows = len(fields), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        if title:
            fig.suptitle(title)

        if self.is_e0sorted:
            qps = self
        else:
            qps = self.sort_by_e0()

        e0mesh = qps.get_e0mesh()

        linestyle = kwargs.pop("linestyle", "o")
        for (field, ax) in zip(fields, ax_list):
            ax.grid(True)
            ax.set_xlabel('e0 [eV]')
            ax.set_ylabel(field)
            yy = qps.get_field(field)
            ax.plot(e0mesh, yy, linestyle, **kwargs)
            ax.plot(e0mesh, e0mesh)
            if fermi is not None:
                ax.plot(2*[fermi], [min(yy), max(yy)])

        # Get around a bug in matplotlib
        if (num_plots % ncols) != 0:
            ax_list[-1].plot([0,1], [0,1], lw=0)
            ax_list[-1].axis('off')

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def build_scissors(self, domains, bounds=None, plot=False, k=3, **kwargs):
        """
        Construct a scissors operator by interpolating the QPState corrections 
        as function of the initial energies E0.

        Args:
            domains:
                list in the form [ [start1, stop1], [start2, stop2]
                Domains should not overlap, cover e0mesh, and given in increasing order.
                Holes are permitted but the interpolation will raise an exception if the
                point is not in domains.
            bounds:
                Specify how to handle out-of-boundary conditions, i.e. how to treat
                energies that do not fall inside one of the domains (not used at present)
            plot:
                If true, use `matplolib` to compare input data  and fit.

        Return:
            instance of `Scissors`operator

        Example:

            # Build the scissors operator.
            scissors = qplist_spin[0].build_scissors(domains)

            # Compute list of interpolated QP energies.
            qp_enes = [scissors.apply(e0) for e0 in ks_energies]

        """
        # Sort QP corrections according to the initial KS energy.
        qps = self.sort_by_e0()
        e0mesh, qpcorrs = qps.get_e0mesh(), qps.get_qpeme0()

        # Check domains.
        domains = np.atleast_2d(domains)
        dsize, dflat = domains.size, domains.ravel()

        for idx, v in enumerate(dflat):
            if idx == 0 and v > e0mesh[0]:
                raise ValueError("min(e0mesh) %s is not included in domains" % e0mesh[0])

            if idx == dsize-1 and v < e0mesh[-1]:
                raise ValueError("max(e0mesh) %s is not included in domains" % e0mesh[-1])

            if idx != dsize-1 and dflat[idx] > dflat[idx+1]:
                raise ValueError("domain boundaries should be given in increasing order.")

            if idx == dsize-1 and dflat[idx] < dflat[idx-1]:
                raise ValueError("domain boundaries should be given in increasing order.")

        # Create the sub_domains and the spline functions in each subdomain.
        func_list = []
        for dom in domains[:]:
            low, high = dom[0], dom[1]
            start, stop = find_ge(e0mesh, low), find_le(e0mesh, high)

            dom_e0 = e0mesh[start:stop+1]
            dom_corr = qpcorrs[start:stop+1]

            # todo check if the number of non degenerate data points > k

            from scipy.interpolate import UnivariateSpline
            f = UnivariateSpline(dom_e0, dom_corr, w=None, bbox=[None, None], k=k, s=None)
            func_list.append(f)

        # Build the scissors operator.
        sciss = Scissors(func_list, domains, bounds)

        title = kwargs.pop("title", None)

        # Compare fit with input data.
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(e0mesh, qpcorrs, 'o', label="input data")
            if title:
                plt.suptitle(title)
            for dom in domains[:]:
                plt.plot(2*[dom[0]], [min(qpcorrs), max(qpcorrs)])
                plt.plot(2*[dom[1]], [min(qpcorrs), max(qpcorrs)])
            intp_qpc = [sciss.apply(e0) for e0 in e0mesh]
            plt.plot(e0mesh, intp_qpc, label="scissor")
            plt.legend()
            plt.show()

        # Return the object.
        return sciss

    def merge(self, other, copy=False):
        """
        Merge self with other. Return new QPList

        Raise:
            ValueError if merge cannot be done.
        """
        skb0_list = [qp.skb for qp in self]
        for qp in other:
            if qp.skb in skb0_list:
                raise ValueError("Found duplicated (s,b,k) indexes: %s" % str(qp.skb))

        if copy:
            qps = self.copy() + other.copy()
        else:
            qps = self + other

        return QPList(qps)


class Sigmaw(object):
    """This object stores the values of the self-energy as function of frequency"""

    def __init__(self, spin, kpoint, band, wmesh, sigmaxc_values, spfunc_values):
        self.spin, self.kpoint, self.band = spin, kpoint, band
        self.wmesh = np.array(wmesh)

        self.xc = Function1D(self.wmesh, sigmaxc_values)
        self.spfunc = Function1D(self.wmesh, spfunc_values)

    def plot_ax(self, ax, w="a", **kwargs):
        """Helper function to plot data on the axis ax."""
        #if not kwargs:
        #    kwargs = {"color": "black", "linewidth": 2.0}

        lines = []
        extend = lines.extend

        if w == "s":
            f = self.xc
            label = kwargs.get("label", "$\Sigma(\omega)$")
            extend(f.plot_ax(ax, cplx_mode="re", label="Re " + label))
            extend(f.plot_ax(ax, cplx_mode="im", label="Im " + label))
            ax.legend(loc="best")
            #ax.set_ylabel('Energy [eV]')

        elif w == "a":
            f = self.spfunc
            label = kwargs.get("label", "$A(\omega)$")
            extend(f.plot_ax(ax, label=label))
            # Plot I(w)
            #ax2 = ax.twinx()
            #extend(f.cumintegral().plot_ax(ax2, label="$I(\omega) = \int_{-\infty}^{\omega} A(\omega')d\omega'$"))
            #ax.set_ylabel('Energy [eV]')
            ax.legend(loc="best")

        else:
            raise ValueError("Don't know how to handle what option %s" % w)

        return lines

    def plot(self, what="sa", **kwargs):
        """
        Plot the self-energy and the spectral function

        Args:
            what:
                String specifying what to plot:
                    s for the self-energy
                    a for spectral function
                Characters can be concatenated.
            args:
                Positional arguments passed to :mod:`matplotlib`.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure.
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        nrows = len(what)
        fig, ax_list = plt.subplots(nrows=nrows, ncols=1, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        if title is None:
            title = 'spin %s, k-point %s, band %s' % (self.spin, self.kpoint, self.band)

        fig.suptitle(title)

        for i, w in enumerate(what):
            ax = ax_list[i]
            ax.grid(True)

            if i == len(what):
                ax.set_xlabel('Frequency [eV]')

            if not kwargs:
                kwargs = {"color": "black", "linewidth": 2.0}

            self.plot_ax(ax, w=w, **kwargs)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


def torange(obj):
    """
    Convert obj into a range. Accepts integer, slice object 
    or any object with an __iter__ method.
    Note that an integer is converted into range(int, int+1)

    >>> torange(1)
    [1]
    >>> torange(slice(0,4,2))
    [0, 2]
    >>> list(torange([1,4,2]))
    [1, 4, 2]
    """
    if isinstance(obj, int):
        return range(obj, obj+1)

    elif isinstance(obj, slice):
        start = obj.start if obj.start is not None else 0
        step = obj.step if obj.step is not None else 1
        return range(start, obj.stop, step)

    else:
        try:
            return obj.__iter__()
        except:
            raise TypeError("Don't know how to convert %s into a range object" % str(obj))


class SIGRES_Plotter(collections.Iterable):
    """
    This object receives a list of `SIGRES_File` objects and provides
    methods to inspect/analyze the GW results (useful for convergence studies)

    .. Attributes:
        
        nsppol:
            Number of spins (must be the same in each file)
        computed_gwkpoints:
            List of k-points where the QP energies have been evaluated.
            (must be the same in each file)

    Usage example:
                                                                  
    .. code-block:: python
        
        plotter = SIGRES_Plotter()
        plotter.add_file("foo_SIGRES.nc", label="foo bands")
        plotter.add_file("bar_SIGRES.nc", label="bar bands")
        plotter.plot_qpgaps()
    """
    def __init__(self):
        self._sigres_files = collections.OrderedDict()
        self._labels = []

    def __len__(self):
        return len(self._sigres_files)

    def __iter__(self):
        return iter(self._sigres_files.values())

    def __str__(self):
        s = ""
        for sigres in self:
            s += str(sigres) + "\n"
        return s

    def add_files(self, filepaths, labels=None):
        """Add a list of filenames to the plotter"""
        for i, filepath in enumerate(list_strings(filepaths)):
            label = None if labels is None else labels[i]
            self.add_file(filepath, label=label)

    def add_file(self, filepath, label=None):
        """Add a filename to the plotter"""
        from abipy.abilab import abiopen
        sigres = abiopen(filepath)
        self._sigres_files[sigres.filepath] = sigres
        # TODO: Not used 
        self._labels.append(label)

        # Initialize/check useful quantities.
        #
        # 1) Number of spins
        if not hasattr(self, "nsppol"): 
            self.nsppol = sigres.nsppol

        if self.nsppol != sigres.nsppol:
            raise ValueError("Found two SIGRES files with different nsppol")

        # The set of k-points where GW corrections have been computed.
        if not hasattr(self, "computed_gwkpoints"):
            self.computed_gwkpoints = sigres.gwkpoints

        if self.computed_gwkpoints != sigres.gwkpoints:
            raise ValueError("Found two SIGRES files with different list of GW k-points.")
            #self.computed_gwkpoints = (self.computed_gwkpoints + sigres.gwkpoints).remove_duplicated()

        if not hasattr(self, "max_gwbstart"):
            self.max_gwbstart = sigres.max_gwbstart
        else:
            self.max_gwbstart = max(self.max_gwbstart, sigres.max_gwbstart)

        if not hasattr(self, "min_gwbstop"):
            self.min_gwbstop = sigres.min_gwbstop
        else:
            self.min_gwbstop = min(self.min_gwbstop, sigres.min_gwbstop)

    @property
    def param_name(self):
        """
        The name of the parameter whose value is checked for convergence.
        This attribute is automatically individuated by inspecting the differences
        inf the sigres.params dictionaries of the files provided.
        """
        try: 
            return self._param_name
        except AttributeError:
            self.set_param_name(param_name=None)
            return self.param_name

    def _get_param_list(self):
        """Return a dictionary with the values of the parameters extracted from the SIGRES files."""
        param_list = collections.defaultdict(list)
                                                               
        for sigres in self:
            for pname in sigres.params.keys():
                param_list[pname].append(sigres.params[pname])

        return param_list

    def set_param_name(self, param_name):
        """
        Set the name of the parameter whose value is checked for convergence.
        if param_name is None, we try to find its name by inspecting 
        the values in the sigres.params dictionaries.
        """
        self._param_name = param_name

    def prepare_plot(self):
        """
        This method must be called before plotting data.
        It tries to figure the name of paramenter we are converging 
        by looking at the set of parameters used to compute the different SIGRES files.
        """
        param_list = self._get_param_list()

        param_name, problem = None, False
        for key, value_list in param_list.items():
            if any(v != value_list[0] for v in value_list):
                if param_name is None:
                    param_name = key
                else:
                    problem = True
                    warnings.warn("Cannot perform automatic detection of convergence parameter.\n" + 
                                  "Found multiple parameters with different values. Will use filepaths as plot labels.")

        self.set_param_name(param_name if not problem else None)

        if self.param_name is None:
            # Could not figure the name of the parameter.
            xvalues = range(len(self))
        else:
            xvalues = param_list[self.param_name]
                                                  
            # Sort xvalues and rearrange the files.
            items = sorted([iv for iv in enumerate(xvalues)], key=lambda item: item[1])
            indices = [item[0] for item in items]
                                                                                             
            files = self._sigres_files.values()
                                                                                             
            newd = collections.OrderedDict()
            for i in indices:
                sigres = files[i]
                newd[sigres.filepath] = sigres
                                                                                             
            self._sigres_files = newd

            # Use sorted xvalues for the plot.
            param_list = self._get_param_list()
            xvalues = param_list[self.param_name]

        self.set_xvalues(xvalues)

    @property
    def xvalues(self):
        """The values used for the X-axis."""
        return self._xvalues 

    def set_xvalues(self, xvalues):
        """xvalues setter."""
        assert len(xvalues) == len(self)
        self._xvalues = xvalues

    def decorate_ax(self, ax, **kwargs):
        ax.grid(True)
        if self.param_name is not None:
            ax.set_xlabel(self.param_name)
        ax.set_ylabel('Energy [eV]')
        ax.legend(loc="best")

        title = kwargs.pop("title", None)
        if title is not None:
            ax.set_title(title)
                                                                                 
        # Set ticks and labels. 
        if self.param_name is None:
            # Could not figure the name of the parameter ==> Use the basename of the files
            ticks, labels = range(len(self)), [f.basename for f in self]
        else:
            ticks, labels = self.xvalues, [f.params[self.param_name] for f in self]

        ax.set_xticks(ticks, minor=False)
        ax.set_xticklabels(labels, fontdict=None, minor=False)

    def extract_qpgaps(self, spin, kpoint):
        """
        Returns a `ndarray` with the QP gaps for the given spin, kpoint.
        Values are ordered with the list of SIGRES files in self.
        """
        qpgaps = []
        for sigres in self:
            k = sigres.ibz.index(kpoint)
            qpgaps.append(sigres.qpgaps[spin, k])
        
        return np.array(qpgaps)

    def extract_qpenes(self, spin, kpoint, band):
        """
        Returns a `ndarray` with the QP energies for the given spin, kpoint.
        Values are ordered with the list of SIGRES files in self.
        """
        qpenes = []
        for sigres in self:
            k = sigres.ibz.index(kpoint)
            qpenes.append(sigres.qpenes[spin,k,band])
        
        return np.array(qpenes)

    def plot_qpgaps(self, spin=None, kpoint=None, hspan=0.01, **kwargs):
        """
        Plot the QP gaps as function of the convergence parameter.

        Args:
            spin:
            kpoint: 
            hspan:
            kwargs:

        Returns:
            `matplotlib` figure
        """
        spin_range = range(self.nsppol) if spin is None else torange(spin)
        kpoints_for_plot = self.computed_gwkpoints  #if kpoint is None else KpointList.as_kpoints(kpoint)

        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        self.prepare_plot()

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        xx = self.xvalues
        for spin in spin_range:
            for kpoint in kpoints_for_plot:
                label = "spin %d, kpoint %s" % (spin, repr(kpoint))
                gaps = self.extract_qpgaps(spin, kpoint)
                ax.plot(xx, gaps, "o-", label=label, **kwargs)

                if hspan is not None:
                    last = gaps[-1]
                    ax.axhspan(last-hspan, last+hspan, facecolor='0.5', alpha=0.5)

        self.decorate_ax(ax)

        if title is not None:
            fig.suptitle(title)
                                 
        if show:
            plt.show()
                                 
        if savefig is not None:
            fig.savefig(savefig)
                                 
        return fig

    def plot_qpenes(self, spin=None, kpoint=None, band=None, hspan=0.01, **kwargs):
        """
        Plot the QP energies as function of the convergence parameter.

        Args:
            spin:
            kpoint: 
            band:
            hspan:
            kwargs:

        Returns:
            `matplotlib` figure
        """
        spin_range = range(self.nsppol) if spin is None else torange(spin)
        band_range = range(self.max_gwbstart, self.min_gwbstop) if band is None else torange(band)
        kpoints_for_plot = self.computed_gwkpoints #if kpoint is None else KpointList.askpoints(kpoint)

        self.prepare_plot()

        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)
                                              
        import matplotlib.pyplot as plt

        # Build grid of plots.
        num_plots, ncols, nrows = len(kpoints_for_plot), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, squeeze=False)
        ax_list = ax_list.ravel()

        if (num_plots % ncols) != 0:
            ax_list[-1].axis('off')

        xx = self.xvalues
        for kpoint, ax in zip(kpoints_for_plot, ax_list):
            
            for spin in spin_range:
                for band in band_range:
                    label = "spin %d, band %d" % (spin, band)
                    qpenes = self.extract_qpenes(spin, kpoint, band)
                    ax.plot(xx, qpenes, "o-", label=label, **kwargs)

                    if hspan is not None:
                        last = qpenes[-1]
                        ax.axhspan(last-hspan, last+hspan, facecolor='0.5', alpha=0.5)

            self.decorate_ax(ax, title="kpoint %s" % repr(kpoint))

        if title is not None:
            fig.suptitle(title)
                                 
        if show:
            plt.show()
                                 
        if savefig is not None:
            fig.savefig(savefig)
                                 
        return fig


class SIGRES_File(AbinitNcFile, Has_Structure, Has_ElectronBands):
    """
    Container storing the GW results reported in the SIGRES.nc file.

    Usage example:
                                                                  
    .. code-block:: python
        
        sigres = SIGRES_File("foo_SIGRES.nc")
        sigres.plot_qps_vs_e0()
        sigres.plot_ksbands_with_qpmarkers()
    """
    def __init__(self, filepath):
        """Read data from the netcdf file path."""
        super(SIGRES_File, self).__init__(filepath)

        # Keep a reference to the SIGRES_Reader.
        self.reader = reader = SIGRES_Reader(self.filepath)

        self._structure = reader.read_structure()
        self.gwcalctyp = reader.gwcalctyp
        self.ibz = reader.ibz
        self.gwkpoints = reader.gwkpoints

        self.gwbstart_sk = reader.gwbstart_sk 
        self.gwbstop_sk = reader.gwbstop_sk
        
        self.min_gwbstart = reader.min_gwbstart
        self.max_gwbstart = reader.max_gwbstart

        self.min_gwbstop = reader.min_gwbstop
        self.max_gwbstop = reader.max_gwbstop

        self._ebands = ebands = reader.ks_bands

        qplist_spin = self.qplist_spin

        # Add QPState markers to the KS band structure.
        # Each marker is a list of tuple(x,y,value)
        for qpattr in QPState.get_fields(exclude=("spin", "band", "kpoint",)):
            x, y, s = [], [], []

            for spin in range(self.nsppol):
                for qp in qplist_spin[spin]:
                    ik = self.ebands.kpoints.index(qp.kpoint)
                    x.append(ik)
                    y.append(qp.e0)
                    size = getattr(qp, qpattr)
                    # Handle complex quantities
                    if np.iscomplex(size): size = size.real
                    s.append(size)

            ebands.set_marker(qpattr, (x, y, s))

        # TODO handle the case in which nkptgw < nkibz
        self.qpgaps = reader.read_qpgaps()
        self.qpenes = reader.read_qpenes()

        self.params = reader.read_params()

    #def __del__(self):
    #    print("in %s __del__" % self.__class__.__name__)
    #    self.reader.close()
    #    super(SIGRES_File, self).__del__()

    @classmethod
    def from_file(cls, filepath):
        """Initialize an instance from file."""
        return cls(filepath)

    @property
    def nsppol(self):
        """Number of spins"""
        return self.ebands.nsppol

    def close(self):
        """Close the netcdf file."""
        self.reader.close()

    @property
    def structure(self):
        """Structure` instance."""
        return self._structure

    @property
    def ebands(self):
        """`ElectronBands` with the KS energies."""
        return self._ebands

    @property
    def qplist_spin(self):
        """Tuple of QPList objects indexed by spin."""
        try:
            return self._qplist_spin
        except AttributeError:
            self._qplist_spin = self.reader.read_allqps()
            return self._qplist_spin

    def get_qplist(self, spin, kpoint):
        qplist = self.reader.read_qplist_sk(spin, kpoint)
        return qplist

    def get_qpcorr(self, spin, kpoint, band):
        """Returns the `QPState` object for the given (s, k, b)"""
        return self.reader.read_qp(spin, kpoint, band)

    def get_qpgap(self, spin, kpoint):
        return self.qpgaps[spin, kpoint]

    def get_sigmaw(self, spin, kpoint, band):
        wmesh, sigxc_values = self.reader.read_sigmaw(spin, kpoint, band)
        wmesh, spf_values = self.reader.read_spfunc(spin, kpoint, band)

        return Sigmaw(spin, kpoint, band, wmesh, sigxc_values, spf_values)

    def get_spfunc(self, spin, kpoint, band):
        wmesh, spf_values = self.reader.read_spfunc(spin, kpoint, band)
        return Function1D(wmesh, spf_values)

    def plot_qps_vs_e0(self, with_fields="all", exclude_fields=None, **kwargs):
        """Plot QPState data as functio of the KS energy."""
        for spin in range(self.nsppol):
            qps = self.qplist_spin[spin].sort_by_e0()
            qps.plot_qps_vs_e0(with_fields=with_fields, exclude_fields=exclude_fields, **kwargs)

    def plot_spectral_functions(self, spin, kpoint, bands, **kwargs):
        """
        Args:
            spin:
                Spin index.
            kpoint:
                Required kpoint.
            bands:
                List of bands

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        if not isinstance(bands, collections.Iterable):
            bands = [bands]

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for band in bands:
            sigw = self.get_sigmaw(spin, kpoint, band)
            label = "skb = %s, %s, %s" % (spin, kpoint, band)
            sigw.plot_ax(ax, label="$A(\omega)$:" + label, **kwargs)

        if title is not None:
            fig.suptitle(title)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def plot_eigvec_qp(self, spin, kpoint, band=None, **kwargs):

        title = kwargs.pop("title", None)

        if kpoint is None:
            from abipy.tools.plotting_utils import ArrayPlotter
            plotter = ArrayPlotter()
            for kpoint in self.ibz:
                ksqp_arr = self.reader.read_eigvec_qp(spin, kpoint, band=band)
                plotter.add_array(str(kpoint), ksqp_arr)

            fig = plotter.plot(title=title)

        else:
            from abipy.tools.plotting_utils import plot_array
            ksqp_arr = self.reader.read_eigvec_qp(spin, kpoint, band=band)

            fig = plot_array(ksqp_arr)

        return fig

    def print_qps(self, spin=None, kpoint=None, bands=None, fmt=None, stream=sys.stdout):
        # TODO Is it still used?
        self.reader.print_qps(spin=spin, kpoint=kpoint, bands=bands, fmt=None, stream=stream)

    def plot_ksbands_with_qpmarkers(self, qpattr="qpeme0", fact=1, **kwargs):
        """
        Plot the KS energies as function an k and add markers 
        whose size is proportional to QPState attribute qpattr
        """
        with_marker = qpattr + ":" + str(fact)
        gwband_range = (self.min_gwbstart, self.max_gwbstop)

        fig = self.ebands.plot(marker=with_marker, band_range=gwband_range, **kwargs)

        return fig

    #def plot_matrix_elements(self, mel_name, spin, kpoint, *args, **kwargs):
    #   matrix = self.reader.read_mel(mel_name, spin, kpoint):
    #   return plot_matrix(matrix, *args, **kwargs)

    #def plot_mlda_to_qps(self, spin, kpoint, *args, **kwargs):
    #    matrix = self.reader.read_mlda_to_qps(spin, kpoint)
    #    return plot_matrix(matrix, *args, **kwargs)


# TODO  Write F90 routine to merge the SIGRES files.
#class SIGRES_Merger(object):
#    """This object merges multiple SIGRES files."""


class SIGRES_Reader(ETSF_Reader):
    """This object provides method to read data from the SIGRES file produced ABINIT.
    # See 70gw/m_sigma_results.F90
    # Name of the diagonal matrix elements stored in the file.
    # b1gw:b2gw,nkibz,nsppol*nsig_ab))
    #_DIAGO_MELS = [
    #    "sigxme",
    #    "vxcme",
    #    "vUme",
    #    "dsigmee0",
    #    "sigcmee0",
    #    "sigxcme",
    #    "ze0",
    #]

      integer :: b1gw,b2gw      ! min and Max gw band indeces over spin and k-points (used to dimension)
      integer :: gwcalctyp      ! Flag defining the calculation type.
      integer :: nkptgw         ! No. of points calculated
      integer :: nkibz          ! No. of irreducible k-points.
      integer :: nbnds          ! Total number of bands
      integer :: nomega_r       ! No. of real frequencies for the spectral function.
      integer :: nomega_i       ! No. of frequencies along the imaginary axis.
      integer :: nomega4sd      ! No. of real frequencies to evaluate the derivative of $\Sigma(E)$.
      integer :: nsig_ab        ! 1 if nspinor=1,4 for noncollinear case.
      integer :: nsppol         ! No. of spin polarizations.
      integer :: usepawu        ! 1 if we are using LDA+U as starting point (only for PAW)

      real(dp) :: deltae       ! Frequency step for the calculation of d\Sigma/dE
      real(dp) :: maxomega4sd  ! Max frequency around E_ks for d\Sigma/dE.
      real(dp) :: maxomega_r   ! Max frequency for spectral function.
      real(dp) :: scissor_ene  ! Scissor energy value. zero for None.

      integer,pointer :: maxbnd(:,:)   SET2NULL
      ! maxbnd(nkptgw,nsppol)
      ! Max band index considered in GW for this k-point.

      integer,pointer :: minbnd(:,:)   SET2NULL
      ! minbnd(nkptgw,nsppol)
      ! Min band index considered in GW for this k-point.

      real(dp),pointer :: degwgap(:,:)   SET2NULL
      ! degwgap(nkibz,nsppol)
      ! Difference btw the QPState and the KS optical gap.

      real(dp),pointer :: egwgap(:,:)   SET2NULL
      ! egwgap(nkibz,nsppol))
      ! QPState optical gap at each k-point and spin.

      real(dp),pointer :: en_qp_diago(:,:,:)   SET2NULL
      ! en_qp_diago(nbnds,nkibz,nsppol))
      ! QPState energies obtained from the diagonalization of the Hermitian approximation to Sigma (QPSCGW)

      real(dp),pointer :: e0(:,:,:)    SET2NULL
      ! e0(nbnds,nkibz,nsppol)
      ! KS eigenvalues for each band, k-point and spin. In case of self-consistent?

      real(dp),pointer :: e0gap(:,:)   SET2NULL
      ! e0gap(nkibz,nsppol),
      ! KS gap at each k-point, for each spin.

      real(dp),pointer :: omega_r(:)   SET2NULL
      ! omega_r(nomega_r)
      ! real frequencies used for the self energy.

      real(dp),pointer :: kptgw(:,:)  SET2NULL
      ! kptgw(3,nkptgw)
      ! ! TODO there is a similar array in sigma_parameters
      ! List of calculated k-points.

      real(dp),pointer :: sigxme(:,:,:)  SET2NULL
      ! sigxme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! Diagonal matrix elements of $\Sigma_x$ i.e $\<nks|\Sigma_x|nks\>$

      real(dp),pointer :: vxcme(:,:,:)  SET2NULL
      ! vxcme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! $\<nks|v_{xc}[n_val]|nks\>$ matrix elements of vxc (valence-only contribution).

      real(dp),pointer :: vUme(:,:,:)   SET2NULL
      ! vUme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! $\<nks|v_{U}|nks\>$ for LDA+U.

      complex(dpc),pointer :: degw(:,:,:)   SET2NULL
      ! degw(b1gw:b2gw,nkibz,nsppol))
      ! Difference between the QPState and the KS energies.

      complex(dpc),pointer :: dsigmee0(:,:,:)  SET2NULL
      ! dsigmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! Derivative of $\Sigma_c(E)$ calculated at the KS eigenvalue.

      complex(dpc),pointer :: egw(:,:,:)  SET2NULL
      ! egw(nbnds,nkibz,nsppol))
      ! QPState energies, $\epsilon_{nks}^{QPState}$.

      complex(dpc),pointer :: eigvec_qp(:,:,:,:)   SET2NULL
      ! eigvec_qp(nbnds,nbnds,nkibz,nsppol))
      ! Expansion of the QPState amplitude in the KS basis set.

      complex(dpc),pointer :: hhartree(:,:,:,:)   SET2NULL
      ! hhartree(b1gw:b2gw,b1gw:b2gw,nkibz,nsppol*nsig_ab)
      ! $\<nks|T+v_H+v_{loc}+v_{nl}|mks\>$

      complex(dpc),pointer :: sigcme(:,:,:,:)   SET2NULL
      ! sigcme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
      ! $\<nks|\Sigma_{c}(E)|nks\>$ at each nomega_r frequency

      complex(dpc),pointer :: sigmee(:,:,:)  SET2NULL
      ! sigmee(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! $\Sigma_{xc}E_{KS} + (E_{QPState}- E_{KS})*dSigma/dE_KS

      complex(dpc),pointer :: sigcmee0(:,:,:)   SET2NULL
      ! sigcmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! Diagonal mat. elements of $\Sigma_c(E)$ calculated at the KS energy $E_{KS}$

      complex(dpc),pointer :: sigcmesi(:,:,:,:)   SET2NULL
      ! sigcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
      ! Matrix elements of $\Sigma_c$ along the imaginary axis.
      ! Only used in case of analytical continuation.

      complex(dpc),pointer :: sigcme4sd(:,:,:,:)   SET2NULL
      ! sigcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
      ! Diagonal matrix elements of \Sigma_c around the zeroth order eigenvalue (usually KS).

      complex(dpc),pointer :: sigxcme(:,:,:,:)   SET2NULL
      ! sigxme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
      ! $\<nks|\Sigma_{xc}(E)|nks\>$ at each real frequency frequency.

      complex(dpc),pointer :: sigxcmesi(:,:,:,:)   SET2NULL
      ! sigxcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
      ! Matrix elements of $\Sigma_{xc}$ along the imaginary axis.
      ! Only used in case of analytical continuation.

      complex(dpc),pointer :: sigxcme4sd(:,:,:,:)   SET2NULL
      ! sigxcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
      ! Diagonal matrix elements of \Sigma_xc for frequencies around the zeroth order eigenvalues.

      complex(dpc),pointer :: ze0(:,:,:)   SET2NULL
      ! ze0(b1gw:b2gw,nkibz,nsppol))
      ! renormalization factor. $(1-\dfrac{\partial\Sigma_c} {\partial E_{KS}})^{-1}$

      complex(dpc),pointer :: omega_i(:)  SET2NULL
      ! omegasi(nomega_i)
      ! Frequencies along the imaginary axis used for the analytical continuation.

      complex(dpc),pointer :: omega4sd(:,:,:,:)  SET2NULL
      ! omega4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol).
      ! Frequencies used to evaluate the Derivative of Sigma.
    """
    def __init__(self, path):
        self.ks_bands = ElectronBands.from_file(path)
        self.nsppol = self.ks_bands.nsppol

        super(SIGRES_Reader, self).__init__(path)

        try:
            self.nomega_r = self.read_dimvalue("nomega_r")
        except self.Error:
            self.nomega_r = 0

        #self.nomega_i = self.read_dim("nomega_i")

        # Save important quantities needed to simplify the API.
        self.structure = self.read_structure()

        self.gwcalctyp = self.read_value("gwcalctyp")
        self.usepawu = self.read_value("usepawu")

        # 1) The K-points of the homogeneous mesh.
        self.ibz = self.ks_bands.kpoints

        # 2) The K-points where QPState corrections have been calculated.
        gwred_coords = self.read_redc_gwkpoints()
        self.gwkpoints = KpointList(self.structure.reciprocal_lattice, gwred_coords)

        # minbnd[nkptgw,nsppol] gives the minimum band index computed
        # Note conversion between Fortran and python convention.
        self.gwbstart_sk = self.read_value("minbnd") - 1
        self.gwbstop_sk = self.read_value("maxbnd")

        # min and Max band index for GW corrections.
        self.min_gwbstart = np.min(self.gwbstart_sk)
        self.max_gwbstart = np.max(self.gwbstart_sk)

        self.min_gwbstop = np.min(self.gwbstop_sk)
        self.max_gwbstop = np.max(self.gwbstop_sk)

        self._egw = self.read_value("egw", cmode="c")

        # Read and save important matrix elements.
        # All these arrays are dimensioned
        # vxcme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
        self._vxcme = self.read_value("vxcme")
        self._sigxme = self.read_value("sigxme")

        self._hhartree = self.read_value("hhartree", cmode="c")

        self._vUme = self.read_value("vUme")
        #if self.usepawu == 0: self._vUme.fill(0.0)

        # Complex arrays
        self._sigcmee0 = self.read_value("sigcmee0", cmode="c")
        self._ze0 = self.read_value("ze0", cmode="c")

        # Frequencies for the spectral function.
        if self.has_spfunc:
            self._omega_r = self.read_value("omega_r")

            self._sigcme = self.read_value("sigcme", cmode="c")
            self._sigxcme = self.read_value("sigxcme", cmode="c")

        # Self-consistent case
        self._en_qp_diago = self.read_value("en_qp_diago")

        # <KS|QPState>
        self._eigvec_qp = self.read_value("eigvec_qp", cmode="c")

        #self._mlda_to_qp

    #def is_selfconsistent(self, mode):
    #    return self.gwcalctyp

    @property
    def has_spfunc(self):
        """True if self contains the spectral function."""
        return self.nomega_r

    def kpt2fileindex(self, kpoint):
        """
        Helper function that returns the index of kpoint in the netcdf file.
        Accepts `Kpoint` instance of integer

        Raise:
            `KpointsError` if kpoint cannot be found.

        .. note:

            This tool is needed since arrays in the netcdf file are dimensioned
            with the total number of k-points in the IBZ.
        """
        if isinstance(kpoint, int):
            kpoint = self.gwkpoints[kpoint]

        try:
            return self.ibz.index(kpoint)
        except:
            raise

    def gwkpt2seqindex(self, gwkpoint):
        """
        This function returns the index of the GW k-point in (0:nkptgw)
        Used to access data in the arrays that are dimensioned [0:nkptgw] e.g. minbnd.
        """
        if isinstance(gwkpoint, int):
            return gwkpoint
        else:
            return self.gwkpoints.index(gwkpoint)

    def read_redc_gwkpoints(self):
        return self.read_value("kptgw")

    def read_allqps(self):
        qps_spin = self.nsppol * [None]

        for spin in range(self.nsppol):
            qps = []
            for gwkpoint in self.gwkpoints:
                ik = self.gwkpt2seqindex(gwkpoint)
                bands = range(self.gwbstart_sk[spin,ik], self.gwbstop_sk[spin,ik])
                for band in bands:
                    qps.append(self.read_qp(spin, gwkpoint, band))

            qps_spin[spin] = QPList(qps)

        return tuple(qps_spin)

    def read_qplist_sk(self, spin, kpoint):
        ik = self.gwkpt2seqindex(kpoint)
        bstart = self.gwbstart_sk[spin, ik]
        bstop = self.gwbstop_sk[spin, ik]

        qps = [self.read_qp(spin, kpoint, band) for band in range(bstart, bstop)]

        return QPList(qps)

    #def read_qpene(self, spin, kpoint, band)

    def read_qpenes(self):
        return self._egw[:, :, :]

    def read_qp(self, spin, kpoint, band):
        ik_file = self.kpt2fileindex(kpoint)
        ib_file = band - self.gwbstart_sk[spin, self.gwkpt2seqindex(kpoint)]

        return QPState(
            spin=spin,
            kpoint=kpoint,
            band=band,
            e0=self.read_e0(spin, ik_file, band),
            qpe=self._egw[spin, ik_file, band],
            qpe_diago=self._en_qp_diago[spin, ik_file, band],
            vxcme=self._vxcme[spin, ik_file, ib_file],
            sigxme=self._sigxme[spin, ik_file, ib_file],
            sigcmee0=self._sigcmee0[spin, ik_file, ib_file],
            vUme=self._vUme[spin, ik_file, ib_file],
            ze0=self._ze0[spin, ik_file, ib_file],
        )

    def read_qpgaps(self):
        """Read the QP gaps. Returns ndarray with shape [nsppol, nkibz] in eV"""
        return self.read_value("egwgap")

    def read_e0(self, spin, kfile, band):
        return self.ks_bands.eigens[spin, kfile, band]

    def read_sigmaw(self, spin, kpoint, band):
        """Returns the real and the imaginary part of the self energy."""
        if not self.has_spfunc:
            raise ValueError("%s does not contain spectral function data" % self.path)

        ik = self.kpt2fileindex(kpoint)

        return self._omega_r, self._sigxcme[spin,:,ik,band]

    def read_spfunc(self, spin, kpoint, band):
        """
        Returns the spectral function.

         one/pi * ABS(AIMAG(Sr%sigcme(ib,ikibz,io,is))) /
         ( (REAL(Sr%omega_r(io)-Sr%hhartree(ib,ib,ikibz,is)-Sr%sigxcme(ib,ikibz,io,is)))**2 &
        +(AIMAG(Sr%sigcme(ib,ikibz,io,is)))**2) / Ha_eV,&
        """
        if not self.has_spfunc:
            raise ValueError("%s does not contain spectral function data" % self.path)

        ik = self.kpt2fileindex(kpoint)
        ib = band - self.gwbstart_sk[spin, self.gwkpt2seqindex(kpoint)]

        aim_sigc = np.abs(self._sigcme[spin,:,ik,ib].imag)

        den = np.zeros(self.nomega_r)
        for (io, omega) in enumerate(self._omega_r):
            den[io] = (omega - self._hhartree[spin,ik,ib,ib] - self._sigxcme[spin,io,ik,ib].real) ** 2 + \
                np.imag(self._sigcme[spin,io,ik,ib]) ** 2

        return self._omega_r, 1./np.pi * (aim_sigc/den)

    def read_eigvec_qp(self, spin, kpoint, band=None):
        """
        Returns <KS|QPState> for the given spin, kpoint and band.

        If band is None, <KS_b|QP_{b'}> is returned.
        """
        ik = self.kpt2fileindex(kpoint)
        if band is not None:
            return self._eigvec_qp[spin,ik,:,band]
        else:
            return self._eigvec_qp[spin,ik,:,:]

    def read_params(self):
        """
        Read the parameters of the calculation. 
        Returns `AttrDict` instance with the value of the parameters.
        """
        param_names = [
            "ecutwfn",
            "ecuteps",
            "ecutsigx",
            "sigma_nband",
        ]

        params = {}
        for pname in param_names:
            try:
                params[pname] = self.read_value(pname)
            except self.Error:
                pass
        
        # Other quantities that might be subject to convergence studies.
        params["nkibz"] = len(self.ibz)

        return AttrDict(params)

    def print_qps(self, spin=None, kpoint=None, bands=None, fmt=None, stream=sys.stdout):
        spins = range(self.nsppol) if spin is None else [spin]
        kpoints = self.gwkpoints if kpoint is None else [kpoint]
        if bands is not None: bands = [bands]

        header = QPState.get_fields(exclude=["spin", "kpoint"])

        for spin in spins:
            for kpoint in kpoints:
                table_sk = [header]
                if bands is None:
                    ik = self.gwkpt2seqindex(kpoint)
                    bands = range(self.gwbstart_sk[spin,ik], self.gwbstop_sk[spin,ik])

                for band in bands:
                    qp = self.read_qp(spin, kpoint, band)
                    d = qp.to_strdict(fmt=fmt)
                    table_sk.append([d[k] for k in header])

                stream.write("\nkpoint: %s, spin: %s, energy units: eV (NB: bands start from zero)\n" % (kpoint, spin))
                pprint_table(table_sk, out=stream)
                stream.write("\n")

    #def read_mel(self, mel_name, spin, kpoint, band, band2=None):
    #    array = self.read_value(mel_name)
                                                                   
    #def read_mlda_to_qp(self, spin, kpoint, band=None):
    #    """Returns the unitary transformation KS-->QPS"""
    #    ik = self.kpt2fileindex(kpoint)
    #if band is not None:
    #    return self._mlda_to_qp[spin,ik,:,band]
    #else:
    #    return self._mlda_to_qp[spin,ik,:,:]
                                                                   
    #def read_qprhor(self):
    #    """Returns the QPState density in real space."""

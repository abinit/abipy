from __future__ import print_function, division

import sys
import collections
import copy
import cStringIO
import numpy as np

from pymatgen.io.abinitio.netcdf import ETSF_Reader
from abipy.tools import find_le, find_ge

from abipy.core.constants import Ha_eV
from abipy.core.func1d import Function1D

from abipy.tools.text import pprint_table
from abipy.kpoints.kpoints import askpoints, kpoints_factory
from abipy.electrons.ebands import ElectronBands
from abipy.iotools import AbinitNcFile

__all__ = [
    "QP",
    "SIGRES_File",
]

_QP_FIELDS = "spin kpoint band e0 qpe qpe_diago vxcme sigxme, sigcmee0 vUme ze0"



class QP(collections.namedtuple("QP", _QP_FIELDS)):
    """
    QP correction for given (spin, kpoint, band).

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
        return QP(**d)

    @classmethod
    def get_fields(cls, exclude=()):
        fields = list(cls._fields) + ["qpeme0"]
        for e in exclude:
            fields.remove(e)
        return tuple(fields)

    def _asdict(self):
        od = super(QP, self)._asdict()
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

#########################################################################################


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

        strio = cStringIO.StringIO()
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
            raise ValueError("QP corrections are not sorted. Use sort_by_e0")

        return np.array([qp.e0 for qp in self])

    def get_field(self, field):
        """ndarray containing the values of field."""
        return np.array([getattr(qp, field) for qp in self])

    def get_qpenes(self):
        """Return an array with the QP energies."""
        return self.get_field("qpe")

    def get_qpeme0(self):
        """Return an arrays with the QP corrections."""
        return self.get_field("qpeme0")

    def to_table(self):
        header = QP.get_fields(exclude=["spin", "kpoint"])
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
                    See `_QP_FIELDS` for the list of available fields.
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

        if isinstance(with_fields, str):
            if with_fields == "all":
                fields = list(QP.get_fields(exclude=["spin", "kpoint"]))
            else:
                fields = with_fields.split()

        if exclude_fields:
            if isinstance(exclude_fields, str):
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

        if (num_plots % ncols) != 0:
            ax_list[-1].axis('off')

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

        plt.tight_layout()

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def build_scissors(self, domains, bounds=None, plot=False, **kwargs):
        """
        Construct a scissors operator by interpolating the QP corrections as a function of E0.

        Args:
            domains:
                list in the form [ [start1, stop1], [start2, stop2]
                Domains should not overlap, cover e0mesh, and given in increasing order.
                Holes are permitted but the interpolation will raise an exception if the
                point is not in domains.
            plot:
                If true, use `matplolib` to compare input data  and fit.
        """
        qps = self.sort_by_e0()
        e0mesh, qpcorrs = qps.get_e0mesh(), qps.get_qpeme0()

        # Check domains
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

        # Create the sub_domains
        func_list = []
        for dom in domains[:]:
            low, high = dom[0], dom[1]
            start, stop = find_ge(e0mesh, low), find_le(e0mesh, high)

            dom_e0 = e0mesh[start:stop+1]
            dom_corr = qpcorrs[start:stop+1]

            from scipy.interpolate import UnivariateSpline
            f = UnivariateSpline(dom_e0, dom_corr, w=None, bbox=[None, None], k=3, s=None)
            func_list.append(f)

        # Build the scissors operator.
        from abipy.electrons import Scissors
        sciss = Scissors(func_list, domains, bounds)

        # Compare fit with input data.
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(e0mesh, qpcorrs, label="input data")
            intp_qpc = [sciss.apply(e0) for e0 in e0mesh]
            plt.plot(e0mesh, intp_qpc, label="scissor")
            plt.legend()
            plt.show()

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

#########################################################################################


class Sigmaw(object):
    """This object stores the values of the self-energy as function of frequency"""

    def __init__(self, spin, kpoint, band, wmesh, sigmaxc_values, spfunc_values):
        self.spin = spin
        self.kpoint = kpoint
        self.band = band
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
            label  = kwargs.get("label", "$A(\omega)$")
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
                ax.set_xlabel('Frequency [Ev]')

            if not kwargs:
                kwargs = {"color": "black", "linewidth": 2.0}

            self.plot_ax(ax, w=w, **kwargs)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

##########################################################################################


class SIGRES_File(AbinitNcFile):
    """Container storing the GW results reported in the SIGRES.nc file."""

    def __init__(self, filepath):
        """Read data from the netcdf file path."""
        super(SIGRES_File, self).__init__(filepath)

        # Keep a reference to the SIGRES_Reader.
        self.ncreader = ncreader = SIGRES_Reader(self.filepath)

        self.structure = ncreader.read_structure()
        self.gwcalctyp = ncreader.gwcalctyp
        self.kpoints = ncreader.kpoints
        self.gwkpoints = ncreader.gwkpoints

        self.gwbstart_sk = ncreader.gwbstart_sk 
        self.gwbstop_sk =  ncreader.gwbstop_sk
        
        self.min_gwbstart = ncreader.min_gwbstart
        self.max_gwbstop = ncreader.max_gwbstop

        # TODO
        # Add GW markers to ks_bands.
        self.ks_bands = ks_bands = ncreader.ks_bands

        # Each marker is a list of tuple(x,y,value)
        #qplist_spin = self.get_allqps()
        #x, y, s = [], [], []

        #min_band, max_band = np.inf, -np.inf
        #for spin in range(self.nsppol):
        #    for qp in qplist_spin[spin]:
        #        ik = ks_bands.kpoints.index(qp.kpoint)
        #        x.append(ik)
        #        y.append(qp.e0)
        #        size = getattr(qp, qpattr)
        #        # Handle complex quantities
        #        if np.iscomplex(size): size = size.real
        #        size = fact * size 
        #        s.append(size)

        #        # Plot only bands in this range.
        #        min_band = min(min_band, qp.band)
        #        max_band = max(max_band, qp.band)

        #ks_bands.set_markers(qpattr, (x, y, s))

    #def __del__(self):
    #    print("in %s __del__" % self.__class__.__name__)
    #    self.ncreader.close()
    #    super(SIGRES_File, self).__del__()

    @classmethod
    def from_file(cls, filepath):
        """Initialize an instance from file."""
        return cls(filepath)

    @property
    def nsppol(self):
        """Number of spins"""
        return self.ks_bands.nsppol

    def close(self):
        """Close the netcdf file."""
        self.ncreader.close()

    def get_structure(self):
        """Returns a `Structure` instance."""
        return self.structure

    def get_allqps(self):
        return self.ncreader.read_allqps()

    def get_qplist(self, spin, kpoint):
        qplist = self.ncreader.read_qplist_sk(spin, kpoint)
        return qplist

    def get_qpcorr(self, spin, kpoint, band):
        """Returns the `QP` object for the given (s, k, b)"""
        return self.ncreader.read_qp(spin, kpoint, band)

    def get_sigmaw(self, spin, kpoint, band):
        wmesh, sigxc_values = self.ncreader.read_sigmaw(spin, kpoint, band)
        wmesh, spf_values = self.ncreader.read_spfunc(spin, kpoint, band)
        return Sigmaw(spin, kpoint, band, wmesh, sigxc_values, spf_values)

    def get_spfunc(self, spin, kpoint, band):
        wmesh, spf_values = self.ncreader.read_spfunc(spin, kpoint, band)
        return Function1D(wmesh, spf_values)

    def plot_qps_vs_e0(self, with_fields="all", exclude_fields=None, **kwargs):
        """Plot QP data as functio of the KS energy."""
        qps_spin = self.get_allqps()

        for spin in range(self.nsppol):
            qps = qps_spin[spin].sort_by_e0()
            qps.plot_qps_vs_e0(with_fields=with_fields, exclude_fields=exclude_fields, **kwargs)

    def plot_spectral_functions(self, spin, kpoint, bands, *args, **kwargs):
        """
        Args:
            spin:
                Spin index.
            kpoint:
                Required kpoint.
            bands:
                List of bands
            args:
                Positional arguments passed to `matplotlib`.

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
            for kpoint in self.kpoints:
                ksqp_arr = self.ncreader.read_eigvec_qp(spin, kpoint, band=band)
                plotter.add_array(str(kpoint), ksqp_arr)
            plotter.plot(title=title)

        else:
            from abipy.tools.plotting_utils import plot_array
            ksqp_arr = self.ncreader.read_eigvec_qp(spin, kpoint, band=band)
            plot_array(ksqp_arr)

    def print_qps(self, spin=None, kpoint=None, bands=None, fmt=None, stream=sys.stdout):
        # TODO Is it still used?
        self.ncreader.print_qps(spin=spin, kpoint=kpoint, bands=bands, fmt=None, stream=stream)

    def plot_ksbands_with_qpmarkers(self, qpattr="qpeme0", fact=1000, **kwargs):
        """
        Plot the KS energies as function an k and add markers 
        whose size is proportional to QP attribute qpattr
        """
        # Each marker is a list of tuple(x,y,value)
        ks_bands = self.ks_bands
        qplist_spin = self.get_allqps()
        x, y, s = [], [], []

        min_band, max_band = np.inf, -np.inf
        for spin in range(self.nsppol):
            for qp in qplist_spin[spin]:
                ik = ks_bands.kpoints.index(qp.kpoint)
                x.append(ik)
                y.append(qp.e0)
                size = getattr(qp, qpattr)
                # Handle complex quantities
                if np.iscomplex(size): size = size.real
                size = fact * size 
                s.append(size)

                # Plot only bands in this range.
                min_band = min(min_band, qp.band)
                max_band = max(max_band, qp.band)

        ks_bands.set_markers(qpattr, (x, y, s))
        with_marker = qpattr + ":" + str(fact)
        ks_bands.plot(with_marker=with_marker, band_range=(min_band,max_band+1), **kwargs)
        ks_bands.del_markers(qpattr)

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
      ! Difference btw the QP and the KS optical gap.

      real(dp),pointer :: egwgap(:,:)   SET2NULL
      ! egwgap(nkibz,nsppol))
      ! QP optical gap at each k-point and spin.

      real(dp),pointer :: en_qp_diago(:,:,:)   SET2NULL
      ! en_qp_diago(nbnds,nkibz,nsppol))
      ! QP energies obtained from the diagonalization of the Hermitian approximation to Sigma (QPSCGW)

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
      ! Difference between the QP and the KS energies.

      complex(dpc),pointer :: dsigmee0(:,:,:)  SET2NULL
      ! dsigmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! Derivative of $\Sigma_c(E)$ calculated at the KS eigenvalue.

      complex(dpc),pointer :: egw(:,:,:)  SET2NULL
      ! degw(nbnds,nkibz,nsppol))
      ! QP energies, $\epsilon_{nks}^{QP}$.

      complex(dpc),pointer :: eigvec_qp(:,:,:,:)   SET2NULL
      ! eigvec_qp(nbnds,nbnds,nkibz,nsppol))
      ! Expansion of the QP amplitude in the KS basis set.

      complex(dpc),pointer :: hhartree(:,:,:,:)   SET2NULL
      ! hhartree(b1gw:b2gw,b1gw:b2gw,nkibz,nsppol*nsig_ab)
      ! $\<nks|T+v_H+v_{loc}+v_{nl}|mks\>$

      complex(dpc),pointer :: sigcme(:,:,:,:)   SET2NULL
      ! sigcme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
      ! $\<nks|\Sigma_{c}(E)|nks\>$ at each nomega_r frequency

      complex(dpc),pointer :: sigmee(:,:,:)  SET2NULL
      ! sigmee(b1gw:b2gw,nkibz,nsppol*nsig_ab))
      ! $\Sigma_{xc}E_{KS} + (E_{QP}- E_{KS})*dSigma/dE_KS

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
        self.kpoints = kpoints_factory(self)

        # 2) The K-points where QP corrections have been calculated.
        gwred_coords = self.read_redc_gwkpoints()
        self.gwkpoints = askpoints(gwred_coords, self.structure.reciprocal_lattice)

        # minbnd[nkptgw,nsppol] gives the minimum band index computed
        # Note conversion between Fortran and python convention.
        self.gwbstart_sk = self.read_value("minbnd") - 1
        self.gwbstop_sk = self.read_value("maxbnd")

        # min and Max band index for GW corrections.
        self.min_gwbstart = np.min(self.gwbstart_sk)
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

        # <KS|QP>
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
            return self.kpoints.index(kpoint)
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
            try:
                return self.gwkpoints.index(gwkpoint)
            except:
                raise

    def read_redc_gwkpoints(self):
        return self.read_value("kptgw")

    def read_allqps(self):
        qps_spin = self.nsppol * [None]
        for spin in range(self.nsppol):
            qps = []
            for gwkpoint in self.gwkpoints:
                i = self.gwkpt2seqindex(gwkpoint)
                bands = range(self.gwbstart_sk[spin,i], self.gwbstop_sk[spin,i])
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

    def read_qp(self, spin, kpoint, band):
        ik_file = self.kpt2fileindex(kpoint)
        ib_file = band - self.gwbstart_sk[spin, self.gwkpt2seqindex(kpoint)]

        d = dict(
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
        return QP(**d)

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
        Returns <KS|QP> for the given spin, kpoint and band.

        If band is None, <KS_b|QP_{b'}> is returned.
        """
        ik = self.kpt2fileindex(kpoint)
        if band is not None:
            return self._eigvec_qp[spin,ik,:,band]
        else:
            return self._eigvec_qp[spin,ik,:,:]

    def print_qps(self, spin=None, kpoint=None, bands=None, fmt=None, stream=sys.stdout):
        spins = range(self.nsppol) if spin is None else [spin]
        kpoints = self.gwkpoints if kpoint is None else [kpoint]
        if bands is not None: bands = [bands]

        header = QP.get_fields(exclude=["spin", "kpoint"])

        for spin in spins:
            for kpoint in kpoints:
                table_sk = [header]
                if bands is None:
                    i = self.gwkpt2seqindex(kpoint)
                    bands = range(self.gwbstart_sk[spin,i], self.gwbstop_sk[spin,i])

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
    #    """Returns the QP density in real space."""

from __future__ import division, print_function

import sys
import tempfile
import collections
import copy
import itertools
import warnings
import numpy as np

from abipy.core.constants import Ha_eV, eV_Ha, Bohr_Ang
from abipy.core.func1d import Function1D
from abipy.core.kpoints import Kpoint, Kpath, IrredZone, kpoints_factory
from abipy.tools import AttrDict
from abipy.iotools import ETSF_Reader, Visualizer, bxsf_write
from abipy.tools import gaussian
from .edos import ElectronDOS

__all__ = [
    "ElectronBands",
    "ElectronBandsPlotter",
    "ElectronDosPlotter",
]


class KSState(collections.namedtuple("KSState", "spin kpoint band eig occ")):
    """
    Kohn-Sham data for given (spin, kpoint, band).

    .. Attributes:

        spin:
            spin index (C convention, i.e >= 0)
        kpoint:
            `Kpoint` object.
        band:
            band index. (C convention, i.e >= 0).
        eig:
            KS eigenvalue.
        occ:
            Occupation factor

    .. note:: Energies are in eV.
    """
    @property
    def skb(self):
        """Tuple with (spin, kpoint, band)"""
        return self.spin, self.kpoint, self.band

    def copy(self):
        d = {f: copy.copy(getattr(self, f)) for f in self._fields}
        return KSState(**d)

    @classmethod
    def get_fields(cls, exclude=()):
        fields = list(cls._fields)
        for e in exclude:
            fields.remove(e)
        return tuple(fields)

    def _asdict(self):
        od = super(KSState, self)._asdict()
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


class Smearing(AttrDict):
    """Stores data and information about the smearing technique."""
    _MANDATORY_KEYS = [
        "scheme",
        "occopt",
        "tsmear_ev",
    ]

    def __init__(self, *args, **kwargs):
        super(Smearing, self).__init__(*args, **kwargs)
        for mkey in self._MANDATORY_KEYS:
            if mkey not in self:
                raise ValueError("Mandatory key %s must be provided" % str(mkey))


#class XcInfo(AttrDict):
#    """Stores data and information about the XC functional."""
#    _MANDATORY_KEYS = [
#    ]
#
#    def __init__(self, *args, **kwargs):
#        super(XcInfo, self).__init__(*args, **kwargs)
#        for mkey in self._MANDATORY_KEYS:
#            if mkey not in self:
#                raise ValueError("Mandatory key %s must be provided" % str(mkey))


class ElectronBands(object):
    """
    This object stores the electronic band structure.
    """
    def __init__(self, structure, kpoints, eigens, fermie, occfacts, nelect,
                 nband_sk=None, smearing=None, markers=None, widths=None):
        """
        Args:
            structure:
                pymatgen structure.
            kpoints:
                KpointList instance.
            eigens
                Array-like object with the eigenvalues (eV) stored as [s,k,b] where s: spin , k: kpoint, b: band index
            fermie:
                Fermi level in eV
            occfacts:
                Occupation factors (same shape as eigens)
            nelect:
                Number of valence electrons in the unit cell.
            smearing:
                `Smearing` object storing information on the smearing technique.
            nband_sk:
                Array-like object with the number of bands treated at each [spin,kpoint]
                If not given, nband_sk is initialized from eigens.
            markers:
                Optional dictionary containing markers labelled by a string.
                Each marker is a list of tuple(x,y,value)
                Used for plotting purpose e.g. QP data, energy derivatives...
            widths:
                Optional dictionary containing data used for the so-called fatbands
                Each entry is an array of shape [nsppol, nkpt, mband] giving the width
                of the band at that particular point.
                Used for plotting purpose e.g. L-projections.
        """
        self.structure = structure

        # Eigenvalues and occupancies are stored in ndarrays ordered by [spin,kpt,band]
        self._eigens = np.atleast_3d(eigens)
        self._occfacts = np.atleast_3d(occfacts)
        assert self.occfacts.shape == self.occfacts.shape
        self.nsppol, self.nkpt, self.mband = self.eigens.shape

        if nband_sk is not None:
            self.nband_sk = np.array(nband_sk)
        else:
            self.nband_sk = np.array(self.nsppol*self.nkpt*[self.mband])

        self.kpoints = kpoints
        assert self.nkpt == len(self.kpoints)

        self.fermie= fermie
        self.nelect = nelect

        self.smearing = {} if smearing is None else smearing

        # TODO: Make sure that we always instanciate the abipy structure
        # so that we have access to the symmetry operations.
        # 
        # Find the k-point names in the pymatgen database.
        # We'll use _auto_klabels to label the point in the matplotlib plot
        # if klabels are not specified by the user.
        self._auto_klabels = collections.OrderedDict()
        for idx, kpoint in enumerate(self.kpoints):
            name = self.structure.findname_in_hsym_stars(kpoint)
            if name is not None:
                self._auto_klabels[idx] = name

        if markers is not None:
            for key, xys in markers.items():
                self.set_marker(key, xys)

        if widths is not None:
            for key, width in widths.items():
                self.set_width(key, width)

    @classmethod
    def from_file(cls, filepath):
        """Initialize an instance of ElectronBands from a netCDF file."""
        with Ebands_Reader(filepath) as r:
            return cls(r.read_structure(),
                       r.read_kpoints(),
                       r.read_eigens(),
                       r.read_fermie(),
                       r.read_occfacts(),
                       r.read_nelect(),
                       smearing=r.read_smearing(),
                       nband_sk=r.read_nband_sk(),
                       )

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        """String representation."""
        lines = []
        app = lines.append
        for k in self._slots:
            try:
                value = self.__dict__[k]
                if prtvol == 0 and isinstance(value, np.ndarray):
                    continue
                app("%s = %s" % (k, value))
            except KeyError:
                pass
        return "\n".join(lines)

    # Handy variables used to loop
    @property
    def spins(self):
        return range(self.nsppol)

    @property
    def kidxs(self):
        return range(self.nkpt)

    @property
    def eigens(self):
        """Eigenvalues in eV. ndarray with shape (nspin, nkpt, mband)."""
        return self._eigens

    @property
    def occfacts(self):
        """Occupation factors. ndarray with shape (nspin, nkpt, mband)."""
        return self._occfacts

    @property
    def reciprocal_lattice(self):
        """Reciprocal lattice vectors in Angstrom."""
        return self.structure.reciprocal_lattice

    @property
    def shape(self):
        """Shape of the array with the eigenvalues."""
        return self.nsppol, self.nkpt, self.mband

    @property
    def markers(self):
        try:
            return self._markers
        except AttributeError:
            return {}

    def del_marker(self, key):
        """
        Delete the entry in self.markers with the specied key. 
        All markers are removed if key is None.
        """
        if key is not None:
            try:
                del self._markers[key]
            except AttributeError:
                pass
        else:
            try:
                del self._markers
            except AttributeError:
                pass

    def set_marker(self, key, xys, overwrite=False):
        """
        Set an entry in the markers dictionary.

        Args:
            key:
                string used to label the set of markers.
            xys:
                Three iterables x,y,s where x[i],y[i] gives the
                positions of the i-th markers in the plot and
                s[i] is the size of the marker.
            overwrite:
                True if key can overwrite a pre-existing entry.
        """
        if not hasattr(self, "_markers"):
            self._markers = collections.OrderedDict()

        if key in self.markers and not overwrite:
            raise ValueError("Cannot overwrite key %s in data" % key)

        for s in xys[-1]:
            if np.iscomplex(s):
                raise ValueError("Found ambiguous complex entry %s" % str(s))

        self._markers[key] = xys

    @property
    def widths(self):
        try:
            return self._widths
        except AttributeError:
            return {}

    def del_width(self, key):
        """
        Delete the entry in self.widths with the specified key.
        All keys are removed if key is None.
        """
        if key is not None:
            try:
                del self._widths[key]
            except AttributeError:
                pass
        else:
            try:
                del self._widths
            except AttributeError:
                pass

    def set_width(self, key, width, overwrite=False):
        """
        Set an entry in the widths dictionary.

        Args:
            key:
                string used to label the set of markers.
            width
                array-like of positive numbers, shape is [nsppol, nkpt, mband].
            overwrite:
                True if key can overwrite a pre-existing entry.
        """
        width = np.reshape(width, self.shape)

        if not hasattr(self, "_widths"):
            self._widths = collections.OrderedDict()

        if key in self.widths and not overwrite:
            raise ValueError("Cannot overwrite key %s in data" % key)

        if np.any(np.iscomplex(width)):
            raise ValueError("Found ambiguous complex entry %s" % str(width))

        if np.any(width < 0.0):
            raise ValueError("Found negative entry in width array %s" % str(width))

        self._widths[key] = width

    @property
    def has_bzmesh(self):
        """True if the k-point sampling is homogeneous."""
        return isinstance(self.kpoints, IrredZone)

    @property
    def has_bzpath(self):
        """True if the bands are computed on a k-path."""
        return isinstance(self.kpoints, Kpath)

    def iter_skb(self):
        """Iterator over (spin, band, kpt) indices."""
        for spin in self.spins:
            for k in self.kidxs:
                for band in self.nband_sk[spin,k]:
                    yield (spin, k, band)

    def show_bz(self):
        """Call `matplotlib` to show the Brillouin zone."""
        return self.structure.show_bz()

    def copy(self):
        """Deep copy of self."""
        return copy.deepcopy(self)

    def enemin(self, spin=None, band=None):
        """Compute the minimum of the eigenvalues."""
        my_spins = self.spins
        if spin is not None:
            assert isinstance(spin, int)
            my_spins = list(spin)

        my_kidxs = self.kidxs

        if band is not None:
            assert isinstance(band, int)
            my_bands = list(band)

        emin = np.inf
        for spin in my_spins:
            for k in my_kidxs:
                if band is None:
                    my_bands = range(self.nband_sk[spin,k])
                for band in my_bands:
                    e = self.eigens[spin,k,band]
                    emin = min(emin, e)
        return emin

    def enemax(self, spin=None, band=None):
        """Compute the maximum of the eigenvalues."""
        my_spins = self.spins
        if spin is not None:
            assert isinstance(spin, int)
            my_spins = list(spin)

        my_kidxs = self.kidxs

        if band is not None:
            assert isinstance(band, int)
            my_bands = list(band)

        emax = -np.inf
        for spin in my_spins:
            for k in my_kidxs:
                if band is None:
                    my_bands = range(self.nband_sk[spin,k])
                for band in my_bands:
                    e = self.eigens[spin,k,band]
                    emax = max(emax, e)
        return emax

    def raw_print(self, stream=sys.stdout):
        """Print k-points and energies on stream."""
        stream.write("# Band structure energies in Ev.\n")
        stream.write("# idx   kpt_red(1:3)  ene(b1) ene(b2) ...\n")

        fmt_k = lambda k: " %.6f" % k
        fmt_e = lambda e: " %.6f" % e
        for spin in self.spins:
            stream.write("# spin = " + str(spin) + "\n")
            for (k, kpoint) in enumerate(self.kpoints):
                nb = self.nband_sk[spin,k]
                ene_sk = self.eigens[spin,k,:nb]
                st = str(k+1)
                for c in kpoint: st += fmt_k(c)
                for e in ene_sk: st += fmt_e(e)
                stream.write(st+"\n")

        stream.flush()

    #@property
    #def homo_bands(self):
    #    try:
    #        return self._homo_bands

    #    except AttributeError:
    #        hband = self.nelect / 2
    #        if hband != int(hband):
    #             homo_bands = self.nsppol * [None]

    #        has_gap = self.nsppol * [True]
    #        for spin in self.spins:
    #            opt_gaps = self.eigens[spin,:,hband+1] - self.eigens[spin,:,hband]
    #            has_gap[spin] = (has_gap[spin] and np.all(opt_gaps > 0.0))

    #        hband if has_gap else None

    #def lomo_band(self):
    #    return None if self.homo_band is None else self.homo_band + 1

    #def _compute_hlstates(self):
    #    homos, lumos, lomos = [], [], []

    #    for spin in self.spins:
    #        homo_band = self.homo_bands(spin)
    #        lumo_band = homo_band + 1

    #        homo_kidx = self.eigens[spin,:,homo_band].argmax()
    #        lumo_kidx = self.eigens[spin,:,lumo_band].argmin()
    #        lomo_kidx = self.eigens[spin,:,1].argmin()

    #        homos.append(KSState(
    #            spin=spin,
    #            kpoint=self.kpoints[homo_kidx],
    #            band=homo_band,
    #            eig=self.eigens[spin,homo_kidx,homo_band],
    #            occ=self.occfacts[spin,homo_kidx,homo_band],
    #        ))

    #        lumos.append(KSState(
    #            spin=spin,
    #            kpoint=self.kpoints[lumo_kidx],
    #            band=lumo_band,
    #            eig=self.eigens[spin,lumo_kidx,lumo_band],
    #            occ=self.occfacts[spin,lumo_kidx,lumo_band],
    #        ))

    #        lomos.append(KSState(
    #            spin=spin,
    #            kpoint=self.kpoints[lomo_kidx],
    #            band=1,
    #            eig=self.eigens[spin,lomo_kidx,1],
    #            occ=self.occfacts[spin,lomo_kidx,1],
    #        ))

    #    return map(tuple, [homos, lumos, lomos])

    #@property
    #def is_metallic(self):
        #"""True if bands are metallic"""

    #@property
    #def bandwith(self):
    #    bandwiths = self*nsppol * [None]
    #    for spin in self.spins:
    #        bandwiths[spin] = self.homo_state[spin].eig - self.lomo_states[spin].eig
    #    return bandwiths

    #def direct_gap(self):
    #def fundamental_gap(self):

    #def info(self):
    #    for spin in self.spins:
    #        homo_band = self.homo_band(spin)
    #        opt_gaps = self.eigens[spin,:,homo_band+1] - self.eigens[spin,:,homo_band]
    #        kmax_idx = opt_gaps.argmax()
    #        kmix_idx = opt_gaps.argmin()
    #        kmax = self.kpoints[kmax_idx]
    #        kmin = self.kpoints[kmin_idx]
    #    np.minloc(opt_gaps)

    def get_dos(self, method="gaussian", step=0.1, width=0.2):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            method:
                String defining the method
            step:
                Energy step (eV) of the linear mesh.
            width:
                Standard deviation (eV) of the gaussian.

        Returns:
            `ElectronDOS` object.
        """
        if abs(self.kpoints.sum_weights - 1) > 1.e-6:
            raise ValueError("Kpoint weights should sum up to one")

        # Compute the linear mesh.
        e_min = self.enemin()
        e_min -= 0.1 * abs(e_min)

        e_max = self.enemax()
        e_max += 0.1 * abs(e_max)

        nw = 1 + (e_max - e_min) / step
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)

        dos = np.zeros((self.nsppol, nw))

        if method == "gaussian":
            for spin in self.spins:
                for (k, kpoint) in enumerate(self.kpoints):
                    weight = kpoint.weight
                    for band in range(self.nband_sk[spin,k]):
                        e = self.eigens[spin,k,band]
                        dos[spin] += weight * gaussian(mesh, width, center=e)

        else:
            raise ValueError("Method %s is not supported" % method)

        return ElectronDOS(mesh, dos)

    def get_jdos(self, spin, valence, conduction, method="gaussian", step=0.1, width=0.2,
                 mesh=None):
        """
        Compute the join density of states at q==0
            :math:`\sum_{kbv} f_{vk} (1 - f_{ck}) \delta(\omega - E_{ck} - E_{vk})`

        Args:
            spin:
                Spin index.
            valence:
                Int or iterable with the valence indices.
            conduction:
                Int or iterable with the conduction indices.
            method:
                String defining the method.
            step:
                Energy step (eV) of the linear mesh.
            width:
                Standard deviation (eV) of the gaussian.
            mesh: Frequency mesh to use. If None, the mesh is computed automatically from the
                  eigenvalues.

        Returns:
            `Function1D` object.
        """
        if abs(self.kpoints.sum_weights - 1) > 1.e-6:
            raise ValueError("Kpoint weights should sum up to one")

        if not isinstance(valence, collections.Iterable):
            valence = [valence]
        if not isinstance(conduction, collections.Iterable):
            conduction = [conduction]

        if mesh is None:
            # Compute the linear mesh.
            cmin, cmax = +np.inf, -np.inf
            vmin, vmax = +np.inf, -np.inf
            for c in conduction:
                cmin = min(cmin, self.eigens[spin,:,c].min())
                cmax = max(cmax, self.eigens[spin,:,c].max())
            for v in valence:
                vmin = min(vmin, self.eigens[spin,:,v].min())
                vmax = max(vmax, self.eigens[spin,:,v].max())

            e_min = cmin - vmax
            e_min -= 0.1 * abs(e_min)

            e_max = cmax - vmin
            e_max += 0.1 * abs(e_max)

            nw = 1 + (e_max - e_min) / step
            mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)
        else:
            nw = len(mesh)

        jdos = np.zeros(nw)

        # Normalize the occupation factors.
        full = 2.0 if self.nsppol == 1 else 1.0

        if method == "gaussian":
            for (k, kpoint) in enumerate(self.kpoints):
                weight = kpoint.weight
                for c in conduction:
                    ec = self.eigens[spin,k,c]
                    fc = 1 - self.occfacts[spin,k,c] / full
                    for v in valence:
                        ev = self.eigens[spin,k,v]
                        fv = self.occfacts[spin,k,v] / full
                        fact = weight * fv * fc
                        jdos += fact * gaussian(mesh, width, center=ec-ev)

        else:
            raise ValueError("Method %s is not supported" % method)

        return Function1D(mesh, jdos)

    def plot_jdosvc(self, vrange, crange, method="gaussian", step=0.1, width=0.2, cumulative=True,
                    **kwargs):
        """
        Plot the decomposition of the joint-density of States (JDOS).

        Args:
            vrange:
                Int or Iterable with the indices of the valence bands to consider.
            crange:
                Int or Iterable with the indices of the conduction bands to consider.
            method:
                String defining the method.
            step:
                Energy step (eV) of the linear mesh.
            width:
                Standard deviation (eV) of the gaussian.
            cumulative:
                True for cumulative plots (default).

            ================  ==============================================================
            kwargs            Meaning
            ================  ==============================================================
            title             Title of the plot (Default: None).
            show              True to show the figure (Default).
            savefig           'abc.png' or 'abc.eps'* to save the figure to a file.
            ================  ==============================================================

        Returns:
            `matplotlib` figure
        """
        if not isinstance(crange, collections.Iterable):
            crange = [crange]
        if not isinstance(vrange, collections.Iterable):
            vrange = [vrange]

        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        for s in self.spins:
            ax = fig.add_subplot(1, self.nsppol, s+1)
            # Get total JDOS
            tot_jdos = self.get_jdos(s, vrange, crange, method=method, step=step, width=width)

            jdos_vc = collections.OrderedDict()
            for v in vrange:
                for c in crange:
                    jd = self.get_jdos(s, v, c, method=method, step=step, width=width, mesh=tot_jdos.mesh)
                    jdos_vc[(v, c)] = jd

            # Plot data for this spin.
            if cumulative:
                cmap = plt.get_cmap("jet")
                cumulative = np.zeros(len(tot_jdos))
                num_plots, i = len(jdos_vc), 0

                for (v, c), jdos in jdos_vc.items():
                    label = "val=%s --> cond=%s, s=%s" % (v, c, s)
                    color = cmap(float(i)/(num_plots))
                    x, y = jdos.mesh, jdos.values
                    ax.plot(x, cumulative + y, lw=1.0, label=label, color=color)
                    ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=0.7)
                    cumulative += jdos.values
                    i += 1
                #tot_jdos.plot_ax(ax, color="k", lw=2, label=" Total JDOS: s=%s," % s, **kwargs)

            else:
                tot_jdos.plot_ax(ax, label="Total JDOS: s=%s," % s, **kwargs)
                for (v, c), jdos in jdos_vc.items():
                    jdos.plot_ax(ax, label="val=%s --> cond=%s, s=%s" % (v,c,s), **kwargs)

        plt.legend(loc="best")

        if title is not None:
            fig.suptitle(title)

        if show:
            plt.show()

        if savefig:
            fig.savefig(savefig)

        return fig

    def apply_scissors(self, scissors):
        """
        Modify the band structure with the scissors operator.

        Args:
            scissors:
                An instance of :class:`Scissors`.

        Returns:
            New instance of `ElectronBands` with modified energies.
        """
        if self.nsppol == 1 and not isinstance(scissors, collections.Iterable):
            scissors = [scissors]

        if self.nsppol == 2 and len(scissors) != 2:
            raise ValueError("Expecting two scissors operators for spin up and down")

        # Create new array with same shape as self.
        qp_energies = np.zeros(self.shape)

        # Calculate Quasi-particle energies with the scissors operator.
        for spin in self.spins:
            sciss = scissors[spin]
            for k in self.kidxs:
                for band in range(self.nband_sk[spin,k]):
                    e0 = self.eigens[spin,k,band]
                    try:
                        qp_ene = e0 + sciss.apply(e0)
                    except sciss.Error:
                        raise

                    # Update the energy.
                    qp_energies[spin,k,band] = qp_ene

        # Change the energies (NB: occupations and fermie are left unchanged).
        return ElectronBands(
            self.structure, self.kpoints, qp_energies, self.fermie, self.occfacts, self.nelect,
            nband_sk=self.nband_sk, smearing=self.smearing, markers=self.markers
        )

    def plot(self, klabels=None, marker=None, width=None, **kwargs):
        """
        Plot the band structure.

        Args:
            klabels:
                dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels.
                e.g. klabels = { (0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L" }.
            band_range:
                Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
            marker:
                String defining the marker to plot.
                accepts the syntax "markername:fact" where
                fact is a float used to scale the marker size.
            width:
                String defining the width to plot.
                accepts the syntax "widthname:fact" where
                fact is a float used to scale the stripe size.

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

        # Select the band range.
        band_range = kwargs.pop("band_range", None)
        if band_range is None:
            band_range = range(self.mband)
        else:
            band_range = range(band_range[0], band_range[1], 1)

        import matplotlib.pyplot as plt
        fig = plt.figure()

        ax = fig.add_subplot(1,1,1)

        # Set ticks and labels.
        self.decorate_ax(ax, klabels=klabels, title=title)

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the band energies.
        for spin in self.spins:
            for band in band_range:
                self.plot_ax(ax, spin=spin, band=band, **kwargs)

        # Add markers to the plot.
        if marker is not None:
            try:
                key, fact = marker.split(":")
            except ValueError:
                key = marker
                fact = 1
            fact = float(fact)

            self.plot_marker_ax(ax, key, fact=fact)

        # Plot fatbands.
        if width is not None:
            try:
                key, fact = width.split(":")
            except ValueError:
                key = width
                fact = 1

            for spin in self.spins:
                for band in band_range:
                    self.plot_width_ax(ax, key, fact=fact)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def plot_fatbands(self, **kwargs):
        """
        Plot the electronic fatbands.

        Args:
            klabels:
                dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels.
                e.g. klabels = { (0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L" }.
            band_range:
                Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
            width:
                String defining the width to plot.
                accepts the syntax "widthname:fact" where
                fact is a float used to scale the stripe size.

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

        # Build grid of plots.
        num_plots, ncols, nrows = len(self.widths), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        for ax, key in zip(ax_list, self.widths):
            # Decorate the axis
            self.decorate_ax(ax, title=key)
            # Plot the energies.
            self.plot_ax(ax)
            # Add width around each band.
            self.plot_width_ax(ax, key)

        if title is not None:
            fig.suptitle(title)

        if show:
            plt.show()
                                 
        if savefig is not None:
            fig.savefig(savefig)
                                 
        return fig

    def decorate_ax(self, ax, klabels=None, title=None):
        if title is not None:
            ax.set_title(title)

        ax.grid(True)
        ax.set_xlabel('k-point')
        ax.set_ylabel('Energy [eV]')
        ax.legend(loc="best")

        # Set ticks and labels.
        ticks, labels = self._make_ticks_and_labels(klabels)

        if ticks:
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False)

    def plot_ax(self, ax, spin=None, band=None, **kwargs):
        """Helper function to plot the energies for (spin,band) on the axis ax."""
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        lines, xx = [], range(self.nkpt)
        for spin in spin_range:
            for band in band_range:
                yy = self.eigens[spin,:,band]
                lines.extend(ax.plot(xx, yy, **kwargs))

        return lines

    def plot_width_ax(self, ax, key, fact=1.0, spin=None, band=None, **kwargs):
        """Helper function to plot fatbands for (spin,band) on the axis ax."""
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        facecolor = kwargs.pop("facecolor", "blue")
        alpha = kwargs.pop("alpha", 0.7)

        width = fact * self.widths[key]
        x = range(self.nkpt)
        for spin in spin_range:
            for band in band_range:
                y = self.eigens[spin,:,band]
                w = width[spin,:,band] * fact
                ax.fill_between(x, y-w/2, y+w/2, facecolor=facecolor, alpha=alpha)

    def plot_marker_ax(self, ax, key, fact=1.0):
        """Helper function to plot the markers for (spin,band) on the axis ax."""
        xvals, yvals, svals = self.markers[key]

        # Use different symbols depending on the value of s.
        # Cannot use negative s.
        pos_x, pos_y, pos_s = [], [], []
        neg_x, neg_y, neg_s = [], [], []
        for x, y, s in zip(xvals, yvals, svals):
            if s >= 0.0:
                pos_x.append(x)
                pos_y.append(y)
                pos_s.append(s)
            else:
                neg_x.append(x)
                neg_y.append(y)
                neg_s.append(s)

        if pos_s:
            ax.scatter(pos_x, pos_y, s=np.abs(pos_s)*fact, marker="^", label=key + " >0")

        if neg_s:
            ax.scatter(neg_x, neg_y, s=np.abs(neg_s)*fact, marker="v", label=key + " <0")

    def _make_ticks_and_labels(self, klabels):
        """Return ticks and labels from the mapping qlabels."""

        if klabels is not None:
            d = collections.OrderedDict()
            for (kcoord, kname) in klabels.items():
                # Build Kpoint instance.
                ktick = Kpoint(kcoord, self.reciprocal_lattice)
                for (idx, kpt) in enumerate(self.kpoints):
                    if ktick == kpt: 
                        d[idx] = kname

        else:
            d = self._auto_klabels

        # Return ticks, labels
        return d.keys(), d.values()

    def plot_with_dos(self, dos, klabels=None, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            dos:
                An instance of :class:`ElectronDOS`.
            klabels:
                dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels.
                e.g. klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.

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
        from matplotlib.gridspec import GridSpec

        gspec = GridSpec(1, 2, width_ratios=[2,1])
        ax1 = plt.subplot(gspec[0])
        # Align bands and DOS.
        ax2 = plt.subplot(gspec[1], sharey=ax1)

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the band structure
        for spin in self.spins:
            for band in range(self.mband):
                self.plot_ax(ax1, spin=spin, band=band, **kwargs)

        # Set ticks and labels.
        if klabels is not None:
            ticks, labels = self._make_ticks_and_labels(klabels)

            ax1.set_xticks(ticks, minor=False)
            ax1.set_xticklabels(labels, fontdict=None, minor=False)

        for ax in (ax1, ax2):
            ax.grid(True)

        if title:
            ax1.set_title(title)

        ax1.set_xlabel('k-point')
        ax1.set_ylabel('Energy [eV]')

        emin = np.min(self.eigens)
        emin -= 0.05 * abs(emin)

        emax = np.max(self.eigens)
        emax += 0.05 * abs(emax)

        ax1.yaxis.set_view_interval(emin, emax)

        # Plot the DOS
        dos.plot_ax(ax2, exchange_xy=True, **kwargs)

        ax2.yaxis.set_ticks_position("right")
        ax2.yaxis.set_label_position("right")

        if show:
            plt.show()

        fig = plt.gcf()
        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def export_bxsf(self, file, structure):
        """
        Export the full band structure on file
        Format is defined by the extension in file.
        """
        # Sanity check.
        errors = []
        eapp = errors.append
        if np.any(self.nband_sk != self.nband_sk[0,0]):
            eapp("nband must be constant")

        if not self.has_omesh:
            eapp("An omogeneous sampling is needed for the Fermi surface")

        if self.has_omesh and not np.allclose(self.shiftk, 0.0):
            eapp("shifted meshes are not supported by Xcryden")

        if errors:
            raise ValueError("\n".join(errors))

        if "." not in file:
            raise ValueError("Cannot detect file extension in path %s: " % file)

        tokens = file.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]:
            # fname == ".ext" ==> Create temporary file.
            path = tempfile.mkstemp(suffix="."+ext, text=True)[1]

        # Xcrysden uses C order and includes periodic images.
        ebands3d = EBands3D(self.kpoints, self.eigens, self.ngkpt, self.shiftk, structure,
                            pbc=3*[True], korder="c")

        with open(file, mode="w") as fh:
            bxsf_write(fh, ebands3d, structure, fermie=self.fermie)

        return Visualizer.from_file(file)

    def derivatives(self, spin, band, order=1, acc=4, asmarker=False):
        """Compute the derivative of the eigenvalues along the path."""
        if self.has_bzpath:
            # Extract the branch.
            branch = self.eigens[spin,:,band]

            # Simulate free-electron bands. This will produce all(effective masses == 1)
            #branch = [0.5 * Ha_eV * (k.norm * Bohr_Ang)**2 for k in self.kpoints]

            # Compute derivatives by finite differences.
            ders_onlines = self.kpoints.finite_diff(branch, order=order, acc=acc)

            if not asmarker:
                return ders_onlines
            else:
                x, y, s = [], [], []
                for i, line in enumerate(self.kpoints.lines):
                    #print(line)
                    x.extend(line)
                    y.extend(branch[line])
                    s.extend(ders_onlines[i])
                    assert len(x) == len(y) == len(s)

                return x, y, s

        else:
            raise ValueError("Derivatives on homogeneous k-meshes are not supported yet")

    def effective_masses(self, spin, band, acc=4):
        """Compute the effective masses."""
        ders2 = self.derivatives(spin, band, acc=acc) * eV_Ha / Bohr_Ang**2
        emasses = 1.0/ders2
        return emasses

#########################################################################################


class EBands3D(object):

    def __init__(self, structure, ibz, ene_ibz, ngkpt, shifts, pbc=3*[True], korder="c"):
        self.structure= structure
        self.ibz      = ibz
        self.ene_ibz  = np.atleast_3d(ene_ibz)
        shape = self.ene_ibz.shape
        self.nsppol, self.nkibz, self.nband = shape[0], shape[1], shape[2]
        self.ngkpt = ngkpt
        self.shifts = shifts
        self.pbc = pbc
        self.korder = korder
        raise NotImplementedError("This code must be tested!")

        # Generator for the K-mesh in the full Brillouin zone.
        from abipy.core.kpoints import KmeshGen, BZSymmetrizer
        self.kmesh_gen = KmeshGen(self.ngkpt, structure.lattice_vectors("g"),
                         shifts=self.shifts, pbc=pbc, korder=korder, wrap_tows=False)

        # Compute the mapping bz_mesh --> ibz
        try:
            self.ksymmetrizer = BZSymmetrizer(self.kmesh_gen, self.ibz, self.structure)
        except:
            raise

    def enebz(self, spin, band):
        """Return energies in the full BZone for given spin and band."""
        return self.ksymmetrizer(self.ene_ibz[spin,:,band])

#########################################################################################


class NestingFactor(object):

    def __init__(self, structure, bands):
        self.structure = structure
        self.bands = bands

        # check whether k-points form a homogeneous sampling.
        if not self.bands.has_bzmesh:
            msg = "The computation of the nesting factor requires a homogeneous k-point sampling"
            raise ValueError(msg)

    @classmethod
    def from_file(cls, filepath):
        """
        Initialize the object from a netcdf
        file containing an electronic band structure.
        """
        with ETSF_Reader(filepath) as r:
            return cls(r.get_structure(), r.get_bands())

    def compute_nesting(self, qpath):
        mesh, values = None, None
        return Function1D(mesh, values)

    def plot(self, qpath):
        nesting = self.compute_nesting(qpath)
        nesting.plot()

#########################################################################################


class ElectronBandsPlotter(object):
    """
    Class for plotting electronic bands structure and DOSes.
    Supports plots on the same graph or separated plots.
    """
    _LINE_COLORS = ["b", "r",]
    _LINE_STYLES = ["-",":","--","-.",]
    _LINE_WIDTHS = [2,]

    def __init__(self):
        self._bands = collections.OrderedDict()
        self._doses = collections.OrderedDict()

    def iter_lineopt(self):
        """Generates style options for lines."""
        for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def add_bands_from_file(self, filepath, label=None):
        """
        Adds a band structure for plotting. Reads data from a Netcdfile
        """
        from abipy import abiopen
        ncfile = abiopen(filepath)
        if label is None:
            label = ncfile.filepath

        self.add_bands(label, ncfile.get_bands())

    def add_bands(self, label, bands, dos=None):
        """
        Adds a band structure for plotting.

        Args:
            label:
                label for the bands. Must be unique.
            bands:
                `ElectronBands` object.
            dos:
                `ElectronDos` object.
        """
        if label in self._bands:
            raise ValueError("label %s is already in %s" % (label, self._bands.keys()))

        self._bands[label] = bands

        if dos is not None:
            self._doses[label] = dos

    def add_bands_list(self, labels, bands_list, dos_list=None):
        """
        Add a list of Bands and DOSes.

        Args:
            labels:
                List of labels.
            bands_list:
                List of `ElectronBands` objects.
            dos_list:
                List of `ElectronDos` objects.
        """
        assert len(labels) == len(bands_list)

        if dos_list is None:
            for label, bands in zip(labels, bands_list):
                self.add_bands(label, bands)
        else:
            assert len(dos_list) == len(bands_list)
            for label, bands, dos in zip(labels, bands_list, dos_list):
                self.add_bands(label, bands, dos=dos)

    def plot(self, klabels=None, *args, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            klabels:
                dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels.
                e.g. klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.
            args:
                Positional arguments passed to :mod:`matplotlib`.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        xlim            x-axis limits. None (default) for automatic determination.
        ylim            y-axis limits. None (default) for automatic determination.
        ==============  ==============================================================

        Returns:
            matplotlib figure.
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        # Build grid of plots.
        if self._doses:
            gspec = GridSpec(1, 2, width_ratios=[2, 1])
            ax1 = plt.subplot(gspec[0])
            # Align bands and DOS.
            ax2 = plt.subplot(gspec[1], sharey=ax1)
            ax_list = [ax1, ax2]
            fig = plt.gcf()
        else:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax_list = [ax1]

        if title is not None:
            fig.suptitle(title)

        for ax in ax_list:
            ax.grid(True)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None:
            [ax.set_ylim(ylim) for ax in ax_list]

        # Plot bands.
        ax = ax_list[0]
        ax.set_xlabel('k-point')
        ax.set_ylabel('Energy [eV]')

        lines, legends = [], []
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        for (label, bands), lineopt in zip(self._bands.items(), self.iter_lineopt()):
            i += 1
            my_kwargs.update(lineopt)
            opts_label[label] = my_kwargs.copy()
            l = bands.plot_ax(ax, spin=None, band=None, *args, **my_kwargs)
            lines.append(l[0])
            legends.append("%s" % label)

            # Set ticks and labels.
            if i == 0 and klabels is not None:
                ticks, labels = bands._make_ticks_and_labels(klabels)
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False)

        # Set legends.
        ax.legend(lines, legends, 'upper right', shadow=True)

        # Add DOSes
        if self._doses:
            ax = ax_list[1]
            for (label, dos) in self._doses.items():
                dos.plot_ax(ax, exchange_xy=True, **opts_label[label])

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


class ElectronDosPlotter(object):
    """
    Class for plotting electronic DOSes.
    """
    #_LINE_COLORS = ["b", "r",]
    #_LINE_STYLES = ["-",":","--","-.",]
    #_LINE_WIDTHS = [2,]

    def __init__(self):
        self._doses = collections.OrderedDict()

    #def iter_lineopt(self):
    #    """Generates style options for lines."""
    #    for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
    #        yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def add_dos_from_file(self, filepath, label=None, method="gaussian", step=0.1, width=0.2):
        """
        Adds a dos for plotting. Reads data from a Netcdfile
        """
        from abipy import abiopen
        bands = abiopen(filepath).get_bands()
        dos = bands.get_dos(method=method, step=step, width=width)
        if label is None: label = filepath

        self.add_dos(label, dos)

    def add_dos(self, label, dos):
        """
        Adds a DOS for plotting.

        Args:
            label:
                label for the DOS. Must be unique.
            dos:
                `ElectronDos` object.
        """
        if label in self._doses:
            raise ValueError("label %s is already in %s" % (label, self._doses.keys()))

        self._doses[label] = dos

    def plot(self, **kwargs):
        """
        Plot the band structure and the DOS.

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

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for (label, dos) in self._doses.items():
            dos.plot_ax(ax, label=label)

        ax.grid(True)
        ax.set_xlabel("Energy [eV]")
        ax.set_ylabel("DOS")
        ax.legend(loc="best")

        if title is not None:
            fig.suptitle(title)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


class Ebands_Reader(ETSF_Reader):
    """
    This object reads band structure data from a netcdf file written
    according to the ETSF-IO specifications.
    """
    #def read_structure(self):
    #    from abipy.core.structure import Structure
    #    return Structure.from_file(self.path)

    def read_kpoints(self):
        """Factory function. Returns KpointList instance."""
        return kpoints_factory(self)

    def read_kfrac_coords(self):
        """Returns a ndarray with the fractional coordinates of the k-points"""
        return self.read_value("reduced_coordinates_of_kpoints")

    def read_nband_sk(self):
        """Array with the number of bands indexed by [s,k]."""
        return self.read_value("number_of_states")

    def read_eigens(self):
        """Eigenvalues in eV."""
        return self.read_value("eigenvalues") * Ha_eV

    def read_occfacts(self):
        """Occupancies."""
        return self.read_value("occupations")

    def read_fermie(self):
        """Fermi level in eV."""
        return self.read_value("fermi_energy") * Ha_eV

    def read_nelect(self):
        """Number of valence electrons."""
        return self.read_value("number_of_electrons")

    def read_smearing(self):
        """Returns a `Smearing` instance with info on the smearing technique."""
        try:
            scheme = "".join(c for c in self.read_value("smearing_scheme"))
        except TypeError:
            scheme = None

        d = Smearing(
            scheme=scheme,
            occopt=self.read_value("occopt"),
            tsmear_ev=self.read_value("smearing_width") * Ha_eV
        )

        # FIXME there's a problem in smearing_scheme
        if scheme is None:
            warnings.warn("warning scheme is None, occopt %s" % d["occopt"])

        return d

    #def read_xcinfo(self):
    #   """Returns a dictionary with info on the XC functional."""
    #    return XcInfo.from_ixc(self.read_value("ixc"))


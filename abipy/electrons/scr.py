# coding: utf-8
"""Objects to analyze the screening file produced by ABINIT."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import six
import abc

from monty.string import marquee # is_string, list_strings,
from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.bisect import index as bs_index
from pymatgen.core.units import Ha_to_eV, eV_to_Ha
from abipy.core.func1d import Function1D
from abipy.core.kpoints import Kpoint, KpointList
from abipy.core.gsphere import GSphere
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.iotools import ETSF_Reader
from abipy.tools.plotting import ArrayPlotter, plot_array, data_from_cplx_mode, add_fig_kwargs, get_ax_fig_plt
from abipy.tools import duck

import logging
logger = logging.getLogger(__name__)


class ScrFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    ScrFile produced by the Abinit GW code.

    Usage example:

    .. code-block:: python

        with ScrFile("foo_SCR.nc") as scr:
            print(scr)
            em1 = scr.get_em1(kpoint=[0, 0, 0])
            em1.plot_w()
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(ScrFile, self).__init__(filepath)

        # Keep a reference to the reader.
        self.reader = ScrReader(filepath)

    def close(self):
        self.reader.close()

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(str(self.structure))
        app("")
        app(marquee("K-points for screening function", mark="="))
        app(str(self.kpoints))
        app("")
        #app("  Number of G-vectors: %d" % self.ng)
        #app("  Number of frequencies: %d (real:%d, imag%d)" % (self.nw, self.nrew, self.nimw))

        if verbose:
            app(str(self.params))

        return "\n".join(lines)

    @property
    def structure(self):
        """Crystalline structure."""
        return self.reader.structure

    @property
    def kpoints(self):
        """List of q-points for the dielectric function."""
        return self.reader.kpoints

    @property
    def wpoints(self):
        """
        Read the frequencies of the dielectric function in Ha.
        Returns numpy array of complex numbers.
        """
        return self.reader.wpoints

    @property
    def netcdf_name(self):
        """The etsf-io name associated to the data on disk."""
        return self.reader.netcdf_name

    @lazy_property
    def params(self):
        """
        :class:`AttrDict` with the most important parameters used to compute the screening
        keys can be accessed with the dot notation i.e. params.zcut
        """
        return self.reader.read_params()

    def get_em1(self, kpoint):
        return self.reader.read_wggmat(kpoint, InverseDielectricFunction)

    @add_fig_kwargs
    def plot_emacro(self, kpoint=(0, 0, 0), g0=0, g1=0, **kwargs):
        """
        Plot the macroscopic dielectric function with/without local-field effects.

        Returns:
            matplotlib figure.
        """
        #return self.reader.read_emacro_lf().plot(**kwargs)
        dataw = self.reader.read_wslice(kpoint, g0=g0, g1=g1)
        print(dataw)

    #@add_fig_kwargs
    #def plot_eelf(self, **kwargs):

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("scr = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(scr)"),
            #nbv.new_code_cell("fig = scr.plot_emacro_lf()"),
            #nbv.new_code_cell("fig = ncfile.phbands.get_phdos().plot()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class ScrReader(ETSF_Reader):
    """
    This object reads the results stored in the SCR (Screening) file produced by ABINIT.
    It provides helper functions to access the most important quantities.
    """
    def __init__(self, filepath):
        super(ScrReader, self).__init__(filepath)

        # Read and store important quantities.
        self.structure = self.read_structure()
        qred_coords = self.read_value("qpoints_dielectric_function")
        self.kpoints = KpointList(self.structure.reciprocal_lattice, qred_coords)
        self.wpoints = self.read_value("frequencies_dielectric_function", cmode="c")

        netcdf_names = ["polarizability", "dielectric_function", "inverse_dielectric_function"]
        nfound = 0
        for tryname in netcdf_names:
            if tryname in self.rootgrp.variables:
                self.netcdf_name = tryname
                nfound += 1

        if nfound == 0:
            raise RuntimeError("Cannot find %s in netcdf file" % str(netcdf_names))
        if nfound > 1:
            raise RuntimeError("Find multiple netcdf arrays (%s)to in netcdf file" % str(netcdf_names))

    def read_params(self):
        """
        Read the most important parameters used to compute the screening i.e.
        the parameters that may be subject to convergence studies.

        Returns:
            :class:`AttrDict` a dictionary whose keys can be accessed
            with the dot notation i.e. d.key.
        """
        # TODO: ecuteps is missing!
        keys = ["ikxc", "inclvkb", "gwcalctyp", "nbands_used", "npwwfn_used",
                "spmeth", "test_type", "tordering", "awtr", "icutcoul", "gwcomp",
                "gwgamma", "mbpt_sciss", "spsmear", "zcut", "gwencomp"]

        return AttrDict({k: self.read_value(k) for k in keys})

    def find_kpoint_fileindex(self, kpoint):
        """
        Returns the k-point and the index of in the netcdf file.
        Accepts `Kpoint` instance or integer

        Raise:
            `KpointsError` if kpoint cannot be found.
        """
        if duck.is_intlike(kpoint):
            ik = int(kpoint)
        else:
            ik = self.kpoints.index(kpoint)

        return self.kpoints[ik], ik

    def read_wslice(self, kpoint, g0=0, g1=0):
        kpoint, ik = self.find_kpoint_fileindex(kpoint)
	#double inverse_dielectric_function(number_of_qpoints_dielectric_function,
        # number_of_frequencies_dielectric_function, number_of_spins, number_of_spins,
        # number_of_coefficients_dielectric_function, number_of_coefficients_dielectric_function, complex) ;

        #var = self.read_value(cls.netcdf_name, cmode="c")
        ig0, ig1 = g0, g1
        var = self.rootgrp.variables[self.netcdf_name]
        values = var[ik, :, 0, 0, ig0, ig1, :]
        return values[:, 0] + 1j * values[:, 1]

    def read_emacro_lf(self, kpoint=(0, 0, 0)):
        """
        Read the macroscopic dielectric function *with* local field effects

            1/ em1_{0,0)(kpoint, omega).

        Return :class:`Function1D`
        """
        em1 = self.get_em1(kpoint=kpoint)
        emacro = 1 / em1.wggmat[:, 0, 0]
        return Function1D(em1.real_wpoints, emacro[:em1.nrew])

    def read_emacro_nlf(self, kpoint=(0, 0, 0)):
        """
        Read the macroscopic dielectric function *without* local field effects.

            e_{0,0)(kpoint, omega).

        Return: :class:`Function1D`

        .. warning::

            This function performs the inversion of e-1 to get e.
            that can be quite expensive and memory demanding for large matrices!
        """
        em1 = self.get_em1(kpoint=kpoint)
        e = np.linalg.inv(em1.wggmat[:em1.nrew, :, :])

        return Function1D(em1.real_wpoints, e[:, 0, 0])

    def read_wggmat(self, kpoint, cls):
        """
        Read data at the given q-point and return an instance
        of `cls` where `cls` is a subclass of `_AwggMatrix`
        """
        kpoint, ik = self.find_kpoint_fileindex(kpoint)
        # TODO: I don't remember how to slice in python-netcdf
        # ecuteps
        all_gvecs = self.read_value("reduced_coordinates_plane_waves_dielectric_function")
        ecuteps = 2
        # TODO: Gpshere.find is very slow if we don't take advantage of shells
        gsphere = GSphere(ecuteps, self.structure.reciprocal_lattice, kpoint, all_gvecs[ik])

        full_wggmat = self.read_value(cls.netcdf_name, cmode="c")
        wggmat = full_wggmat[ik]

        return cls(self.wpoints, gsphere, wggmat, inord="F")


class _AwggMatrix(object):
    r"""
    Base class for two-point functions expressed in reciprocal space
    i.e. a complex matrix $A_{G,G'}(\omega)$ where G, G' are reciprocal
    lattice vectors defines inside the sphere `gsphere`.

    This class is not supposed to be instanciated directly.

    .. attributes:

        wpoints:
        gpshere:
        nrew:
        nwim:
    """
    netcdf_name = "_AwggMatrix"
    latex_name = "Unknown"

    def __init__(self, wpoints, gsphere, wggmat, inord="C"):
        """"
        Args:
            gsphere: :class:`GSphere` with G-vectors and k-point object.
            wpoints: Complex frequency points in Hartree.
            wggmat: [nw, ng, ng] complex array.
            inord: storage order of wggmat. If inord == "F", wggmat in
                in Fortran column-major order. Default: "C" i.e. C row-major order
        """
        self.wpoints = np.array(wpoints, dtype=np.complex)
        self.gsphere = gsphere
        self.wggmat = np.reshape(wggmat, (self.nw, self.ng, self.ng))

        if inord.lower() == "f":
            # Fortran to C.
            for iw, _ in enumerate(wpoints):
                self.wggmat[iw] = self.wggmat[iw].T

        for i in (1, 2):
            assert len(gsphere) == wggmat.shape[-i]
        assert len(self.wpoints) == len(self.wggmat)

        # Find number of real/imaginary frequencies.
        self.nrew = self.nw; self.nimw = 0
        for i, w in enumerate(self.wpoints):
            if np.iscomplex(w):
                self.nrew = i
                break

        self.nimw = self.nw - self.nrew
        if self.nimw and not np.all(np.iscomplex(self.wpoints[self.nrew+1:])):
            raise ValueError("wpoints should contained real points packed in the first positions\n"
                "followed by imaginary points but got: %s" % str(self.wpoints))

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []
        app = lines.append

        #app(marquee("Structure", mark="="))
        #app(str(self.structure))
        #app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))

        app(self.netcdf_name)
        app("  k-point: %s" % self.kpoint)
        app("  Number of G-vectors: %d" % self.ng)
        app("  Total number of frequencies: %d (real: %s, imag: %s)" % (self.nw, self.nrew, self.nimw))
        if self.nrew:
            app("  real frequencies up to %.2f [eV]" % self.real_wpoints[-1].real)
        if self.nimw:
            app("  imaginary frequencies up to %.2f [eV]" % self.imag_wpoints[-1].imag)

        return "\n".join(lines)

    @property
    def kpoint(self):
        """Kpoint object."""
        return self.gsphere.kpoint

    @property
    def ng(self):
        """Number of G-vectors"""
        return len(self.gsphere)

    @property
    def nw(self):
        """Total number of frequencies."""
        return len(self.wpoints)

    @property
    def real_wpoints(self):
        """
        Real frequencies in Hartree. Empty list if not available.
        """
        if self.nrew > 0:
            return np.real(self.wpoints[:self.nrew])
        return []

    @property
    def imag_wpoints(self):
        """
        Imaginary frequencies in Hartree. Empty list if not available.
        """
        if self.nimw > 0:
            return self.wpoints[self.nrew:]
        return []

    @property
    def wggmat_realw(self):
        """The slice of ggmat along the real axis."""
        return self.wggmat[:self.nrew, :, :]

    @property
    def wggmat_imagw(self):
        """The slice of ggmat along the imaginary axis."""
        return self.wggmat[self.nrew:, :, :]

    def windex(self, w, atol=0.001):
        """
        Find the index of the **closest** frequency in `wpoints`.
        """
        if np.iscomplex(w):
            iw = bs_index(self.imag_wpoints.imag, w.imag, atol=atol)
            iw += self.nrew
        else:
            iw = bs_index(self.real_wpoints.real, w, atol=atol)

        return iw

    def gindex(self, gvec):
        """
        Find the index of gvec. If gvec is an int, gvec is returned.
        Raises `ValueError` if gvec is not found.
        """
        if duck.is_intlike(gvec): return int(gvec)
        return self.gsphere.index(gvec)

    def latex_label(self, cplx_mode):
        """Return a latex string to be passed to matplotlib."""
        if cplx_mode == "re": return r"$\Re(" + self.latex_name + ")$"
        if cplx_mode == "im": return r"$\Im(" + self.latex_name + ")$"
        if cplx_mode == "abs": return r"$\\abs(" + self.latex_name + ")$"
        if cplx_mode == "angle": return r"$Phase(" + self.latex_name + ")$"
        raise ValueError("Wrong value for cplx_mode: %s" % cplx_mode)

    @add_fig_kwargs
    def plot_w(self, gvec1, gvec2=None, waxis="real", cplx_mode="re-im", ax=None, **kwargs):
        """
        Plot the frequency dependence of W_{g1, g2}(omega)

        Args:
            gvec1, gvec2:
            waxis: "real" to plot along the real axis, "imag" for the imaginary axis.
            cplx_mode: string defining the data to print.
                       Possible choices are (case-insensitive): `re` for the real part
                       "im" for the imaginary part, "abs" for the absolute value.
                       "angle" will display the phase of the complex number in radians.
                       Options can be concatenated with "-" e.g. "re-im"
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            matplotlib figure.
        """
        # Select data to plot
        ig1 = self.gindex(gvec1)
        ig2 = ig1 if gvec2 is None else self.gindex(gvec2)

        ax, fig, plt = get_ax_fig_plt(ax)
        if waxis == "real":
            if self.nrew == 0: return fig
            xx = (self.real_wpoints * Ha_to_eV).real
            yy = self.wggmat_realw[:, ig1, ig2]

        elif waxis == "imag":
            if self.nimw == 0: return fig
            xx = (self.imag_wpoints * Ha_to_eV).imag
            yy = self.wggmat_imagw[:, ig1, ig2]

        else:
            raise ValueError("Wrong value for waxis: %s" % waxis)

        color_cmode = dict(re="red", im="blue", abs="black", angle="green")
        linewidth = kwargs.pop("linewidth", 2)
        linestyle = kwargs.pop("linestyle", "solid")

        lines = []
        for c in cplx_mode.lower().split("-"):
            l, = ax.plot(xx, data_from_cplx_mode(c, yy),
                         color=color_cmode[c], linewidth=linewidth,
                         linestyle=linestyle,
                         label=self.latex_label(c))
            lines.append(l)

        ax.grid(True)
        ax.set_xlabel(r"$\omega$ [eV]")
        ax.set_title("%s, kpoint: %s" % (self.netcdf_name, self.kpoint))
        ax.legend(loc="best")
        #ax.legend(loc="upper right")

        return fig

    @add_fig_kwargs
    def plot_ggmat(self, cplx_mode="abs", wpos=None, **kwargs):
        """
        Use imshow to plot W_{GG'} matrix

        Args:
            cplx_mode:
            wpos: List of frequency indices to plot. If None, the first frequency is used (usually w=0).
                If wpos == "all" all frequencies are shown (use it carefully)
                Other possible values: "real" if only real frequencies are wanted.
                "imag" for imaginary frequencies only.

        Returns:
            `matplotlib` figure.
        """
        # Get wpos indices.
        choice_wpos = {
            None: [0],
            "all": range(self.nw),
            "real": range(self.nrew),
            "imag": range(self.nrew, self.nw)}

        if any(wpos == k for k in choice_wpos):
            wpos = choice_wpos[wpos]
        else:
            if duck.is_intlike(wpos): wpos = [int(wpos)]
            wpos = np.array(wpos)

        # Build plotter.
        plotter = ArrayPlotter()
        for iw in wpos:
            label = cplx_mode + r" $\omega = %s" % self.wpoints[iw]
            data = data_from_cplx_mode(cplx_mode, self.wggmat[iw,:,:])
            plotter.add_array(label, data)

        return plotter.plot(show=False, **kwargs)


class Polarizability(_AwggMatrix):
    netcdf_name = "polarizability"
    latex_name = r"\tilde chi"


class DielectricFunction(_AwggMatrix):
    netcdf_name = "dielectric_function"
    latex_name = r"\epsilon"


class InverseDielectricFunction(_AwggMatrix):
    netcdf_name = "inverse_dielectric_function"
    latex_name = r"\epsilon^{-1}"

    def _add_ppmodel(self, ppm):
        """
        Add a :class:`PPModel` object to the internal list. Return ppm.
        """
        if not hasattr(self, "ppmodels"): self.ppmodels = []
        self.ppmodels.append(ppm)
        return ppm

    def build_godby_needs_ppmodel(self, wplasma):
        ppm = GodbyNeeds.from_em1(self, wplasma)
        return self._add_ppmodel(ppm)

    @add_fig_kwargs
    def plot_with_ppmodels(self, gvec1, gvec2=None, waxis="real", cplx_mode="re", zcut=0.1/Ha_to_eV, **kwargs):
        """
        Args:
            gvec1, gvec2:
            waxis: "real" to plot along the real axis, "imag" for the imaginary axis.
            cplx_mode: string defining the data to print.
                       Possible choices are (case-insensitive): `re` for the real part
                       "im" for the imaginary part, "abs" for the absolute value.
                       "angle" will display the phase of the complex number in radians.
                       Options can be concatenated with "-" e.g. "re-im"
            zcut: Small shift along the imaginary axis to avoid poles. 0.1 eV is the Abinit default.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            matplotlib figure.
        """
        # Select the (G,G') indices to plot.
        ig1 = self.gindex(gvec1)
        ig2 = ig1 if gvec2 is None else self.gindex(gvec2)

        ax, fig, plt = get_ax_fig_plt(None)

        self.plot_w(gvec1, gvec2=gvec2, waxis=waxis, cplx_mode=cplx_mode, ax=ax, show=False)

        # Compute em1 from the ppmodel on the same grid used for self.
        omegas = {"real": self.real_wpoints, "imag": self.imag_wpoints}[waxis]

        # Get y-limits of the ab-initio em1 to zoom-in the interesting region
        ymin_em1, ymax_em1 = ax.get_ylim()

        for ppm in self.ppmodels:
            em1_ppm = ppm.eval_em1(omegas, zcut)
            em1_ppm.plot_w(ig1, gvec2=ig2, waxis=waxis, cplx_mode=cplx_mode,
                           ax=ax, linestyle="--", show=False)

        ax.set_ylim(ymin_em1, ymax_em1)

        return fig


class PPModel(six.with_metaclass(abc.ABCMeta, object)):
    """
    Abstract base class for Plasmonpole models.
    """

    #@abc.abstractmethod
    #def from_em1(cls, em1):
    #    """Compute the plasmon-pole parameters from the inverse dielectric function."""

    @abc.abstractmethod
    def eval_em1(self, omegas, zcut):
        """Compute the plasmon-pole model at frequency omega (Ha units)."""


class GodbyNeeds(PPModel):

    def __init__(self, gsphere, omegatw, bigomegatwsq):
        r"""
        bigomegatwsq(:)
        Plasmon pole parameters $\tilde\Omega^2_{G Gp}(q)$.

        omegatw(:)
        omegatw(nqibz)%value(npwc,dm2_otq)
        Plasmon pole parameters $\tilde\omega_{G Gp}(q)$.
        """
        self.gsphere = gsphere
        self.kpoint = gsphere.kpoint
        self.omegatw = omegatw
        self.bigomegatwsq = bigomegatwsq

        self.ng = len(gsphere)
        assert len(omegatw) == len(bigomegatwsq)
        assert len(omegatw) == len(gsphere)

    @classmethod
    def from_em1(cls, em1, wplasma):
        # Find omega=0 and the second imaginary frequency to fit the ppm parameters.
        iw0 = -1; iw1 = -1
        for i, w in enumerate(em1.wpoints):
            if np.abs(w) <= 1e-6: iw0 = i
            if np.abs(w - 1j*wplasma) <= 1e-6: iw1 = i

        if iw0 == -1:
            raise ValueError("Cannot find omega=0 in em1")
        if iw1 == -1:
            raise ValueError("Cannot find second imaginary frequency at %s in em1!" % wplasma)

        w0gg = em1.wggmat[iw0, :, :]
        w1gg = em1.wggmat[iw1, :, :]

        aa = w0gg - np.eye(em1.ng)
        diff = w0gg - w1gg
        ratio = aa / diff
        omegatwsq = (ratio - 1.0) * (wplasma ** 2)

        # If omega-twiddle-squared is negative,set omega-twiddle-squared to 1.0 (a reasonable way of treating
        # such terms, in which epsilon**-1 was originally increasing along this part of the imaginary axis)
        # (note: originally these terms were ignored in Sigma; this was changed on 6 March 1990.)
        #if (REAL(omegatwsq) <= 0.0) omegatwsq=one
        omegatwsq[np.where(omegatwsq.real <= 0.0)] = 1.0
        #
        # Get omega-twiddle
        # * Neglect the imag part (if any) in omega-twiddle-squared
        #omegatw(ig,igp)=SQRT(REAL(omegatwsq))
        #omegatw = np.sqrt(omegatwsq)
        omegatw = np.sqrt(omegatwsq.real)
        #omegatw = omegatw + 7j * eV_to_Ha

        bigomegatwsq = -aa * omegatw**2

        return cls(em1.gsphere, omegatw, bigomegatwsq)

    def eval_em1(self, omegas, zcut):
        omegas = np.array(omegas)

        wggmat = np.empty((len(omegas), self.ng, self.ng), dtype=np.complex)
        for i, w in enumerate(omegas):
            # Add shift but only along the real axis.
            delta = 0.0 if w.imag != 0 else 1j * zcut
            den = w**2 - np.real((self.omegatw - delta)**2)
            em1gg = np.eye(self.ng) + self.bigomegatwsq / den
            #arg = (w - self.omegatw.real) / (0.01 * eV_to_Ha)
            #em1gg = em1gg * (1 - np.exp(-arg**2))
            wggmat[i] = em1gg

        return InverseDielectricFunction(self.kpoint, omegas, self.gsphere, wggmat)

    @add_fig_kwargs
    def plot_ggparams(self, **kwargs):
        plotter = ArrayPlotter(*[
            (r"$\\tilde\omega_{G G'}$", self.omegatw),
            (r"$\\tilde\Omega^2_{G, G'}$", self.bigomegatwsq)])

        return plotter.plot(show=False, **kwargs)

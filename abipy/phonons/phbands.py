from __future__ import print_function, division

import sys
import functools
import numpy as np

from abipy.iotools import ETSF_Reader, AbinitNcFile
from abipy.kpoints.kpoints import Kpoint
from abipy.tools import gaussian
from .phdos import PhononDOS

__all__ = [
    "PhononBands",
    "PHBST_File",
]

#########################################################################################

@functools.total_ordering
class PhononMode(object):
    """A phonon mode has a q-point, a frequency, a cartesian displacement and a structure."""
    __slots__ = [
        "qpoint",
        "freq",
        "displ_cart", # Cartesian displacement.
        "structure"
    ]

    def __init__(self, qpoint, freq, displ_cart, structure):
        """
        Args:
            qpoint:
                qpoint in reduced coordinates.
            freq:
                Phonon frequency in eV.
            displ:
                Displacement (Cartesian coordinates, Angstrom)
            structure:
                Pymatgen structure.
        """
        self.qpoint = Kpoint.askpoint(qpoint, structure.reciprocal_lattice)
        self.freq = freq
        self.displ_cart = displ_cart
        self.structure = structure

    #def __str__(self):

    # Rich comparison support (ordered is based on the frequency).
    # Note that missing operators are filled by total_ordering.
    def __eq__(self, other):
        return self.freq == other.freq

    def __lt__(self, other):
        return self.freq < other.freq

        #@property
        #def displ_red(self)
        #    return np.dot(self.xred, self.rprimd)

        #def make_supercell(self, delta):

        #def export(self, path):

        #def visualize(self, visualizer):

#########################################################################################


class PhononBands(object):
    """
    Container object storing the phonon band structure.

    .. Attributes:

      "phfreqs",       # (nqpt, 3*natom)
      "phdispl_cart",  # (nqpt, 3*natom, 3*natom)
                       # the last dimension stores the cartesian components.
      "qpoints",       # qpoints and wtq are replaced by self.ibz that is a list of Kpoints.
      "weights",

    .. note:
        Frequencies are in eV. Cartesian displacements are in Angstrom.
    """

    @classmethod
    def from_file(cls, path):
        """Create the object from a netCDF file."""
        with PHBST_Reader(path) as r:
            structure = r.read_structure()
            phfreqs = r.read_phfreqs()
            phdispl_cart = r.read_phdispl_cart()

            # Build list of q-points
            qcoords = r.read_qredcoords()
            qweights = r.read_qweights()

            #qpoints = kpoints_factory(self)
            qpoints = []
            for (qc, w) in zip(qcoords, qweights):
                qpoints.append(Kpoint(qc, structure.reciprocal_lattice, weight=w))

        return cls(structure, qpoints, phfreqs, phdispl_cart)

    def __init__(self, structure, qpoints, phfreqs, phdispl_cart):
        """
        Args:
            structure:
                Structure object
        """
        self.structure = structure
        self.qpoints = qpoints

        self.phfreqs = phfreqs
        self.phdispl_cart = phdispl_cart

        self.num_qpoints = len(self.qpoints)

        # Handy variables used to loop.
        self.num_atoms = structure.num_sites
        self.num_branches = 3 * self.num_atoms
        self.branches = range(self.num_branches)

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        """String representation."""
        lines = []
        app = lines.append
        for (key, value) in self.__dict__.items():
            if key.startswith("_"): continue
            if prtvol == 0 and isinstance(value, np.ndarray):
                continue
            app("%s = %s" % (key, value))
        return "\n".join(lines)

    def show_qpoints(self, stream=sys.stdout):
        lines = []
        for (i, qpoint) in enumerate(self.qpoints):
            lines.append("%d) %s\n" % (i, qpoint))
        stream.writelines("\n".join(lines))
        stream.flush()

    def displ_of_specie(self, specie):
        """Returns the displacement vectors for the given specie."""
        # TODO recheck the ordering
        # (nqpt, 3*natom, natom, 2) the last dimension stores the cartesian components.
        #raise NotImplementedError("")
        displ_specie = []
        for (i, site) in enumerate(self.structure):
            if site.specie == specie:
                displ_specie.append(self.phdispl_cart[:, :, i, :])
        return displ_specie

    @property
    def displ_shape(self):
        """The shape of phdispl_cart."""
        return self.phdispl_cart.shape

    @property
    def minfreq(self):
        """Minimum phonon frequency."""
        return self.get_minfreq_mode()

    @property
    def maxfreq(self):
        """Maximum phonon frequency."""
        return self.get_maxfreq_mode()

    def get_minfreq_mode(self, mode=None):
        """Compute the minimum of the frequencies."""
        if mode is None:
            return np.min(self.phfreqs)
        else:
            return np.min(self.phfreqs[:, mode])

    def get_maxfreq_mode(self, mode=None):
        """Compute the minimum of the frequencies."""
        if mode is None:
            return np.max(self.phfreqs)
        else:
            return np.max(self.phfreqs[:, mode])

    def raw_print(self, stream=sys.stdout, fmt=None, cvs=False):
        """Write data on stream with format fmt. Use CVS format if cvs."""
        raise NotImplementedError("")
        lines = []
        app = lines.append

        app("# Phonon band structure energies in Ev.")
        app("# idx   qpt_red(1:3)  freq(mode1) freq(mode2) ...")

        if fmt is None:
            significant_figures = 12
            format_str = "{{:.{0}f}}".format(significant_figures)
            fmt = format_str.format

        sep = ", " if cvs else " "
        for (q, qpoint) in enumerate(self.qpoints):
            freq_q = self.phfreqs[q, :]
            for c in qpoint: s += fmt(c)
            for w in freq_q: s += fmt(e)
            line = "%d " % q
            app(line)

        stream.writelines(sep.join(lines))
        stream.flush()

    def get_unstable_modes(self, below_mev=-5.0):
        """Return the list of unstable phonon modes."""
        raise NotImplemetedError("this is a stub")
        umodes = []
        for (q, qpoint) in enumerate(self.qpoints):
            for nu in self.branches:
                freq = self.phfreqs[q, nu]
                if freq < below_mev * 1000:
                    displ_cart = self.phdispl_cart[q, nu, :]
                    umodes.append(PhononMode(qpoint, freq, displ_cart, self.structure))
        return umodes

    def get_dos(self, method="gaussian", step=1.e-4, width=4.e-4):
        """
        Compute the phonon DOS on a linear mesh.

        Args:
            method:
                String defining the method
            step:
                Energy step (eV) of the linear mesh.
            width:
                Standard deviation (eV) of the gaussian.

        Returns:
            PhononDOS object.

        .. warning:
            The DOS computed here is reliable only if the set of q-points for
            a homogeneous sampling of the Brillouin zone.
        """
        # Compute the linear mesh for the DOS
        w_min = self.minfreq
        w_min -= 0.1 * abs(w_min)

        w_max = self.maxfreq
        w_max += 0.1 * abs(w_max)

        nw = 1 + (w_max - w_min) / step

        mesh, step = np.linspace(w_min, w_max, num=nw, endpoint=True, retstep=True)

        values = np.zeros(nw)
        if method == "gaussian":
            for (q, qpoint) in enumerate(self.qpoints):
                weight = qpoint.weight
                for nu in self.branches:
                    w = self.phfreqs[q, nu]
                    values += weight * gaussian(mesh, width, center=w)

        else:
            raise ValueError("Method %s is not supported" % method)

        return PhononDOS(mesh, values)

    def plot(self, qlabels=None, *args, **kwargs):
        """
        Plot the band structure.

        Args:
            qlabels:
                dictionary whose keys are tuple with the reduced
                coordinates of the q-points. The values are the labels.
                e.g. qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.
            args:
                Positional arguments passed to matplotlib.

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

        ax = fig.add_subplot(1, 1, 1)

        if title is not None:
            ax.set_title(title)

        ax.grid(True)
        ax.set_xlabel('q-point')
        ax.set_ylabel('Energy [eV]')

        args = [] if not args else args
        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Set ticks and labels.
        if qlabels is not None:
            ticks, labels = self._ticks_and_labels(qlabels)
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False)

        # Plot phonon branches.
        for nu in self.branches:
            self.plot_ax(ax, nu, *args, **kwargs)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def plot_ax(self, ax, mode, *args, **kwargs):
        """
        Plots the energies for a given mode as a function of the q index on axis ax.

        Return value is a list of lines that were added.
        """
        xx, yy = range(self.num_qpoints), self.phfreqs[:, mode]
        return ax.plot(xx, yy, *args, **kwargs)

    def _ticks_and_labels(self, qlabels):
        """Return ticks and labels from the mapping {qred: qstring} given in qlabels."""
        d = {}
        for (qcoord, qname) in qlabels.items():
            # Build Kpoint instancee
            qtick = Kpoint(qcoord, self.structure.reciprocal_lattice)
            for (q, qpoint) in enumerate(self.qpoints):
                if qtick == qpoint:
                    d[q] = qname
            # ticks, labels
        return d.keys(), d.values()

    def plot_fatbands(self, colormap="jet", max_stripe_width_mev=3.0, qlabels=None, **kwargs):
    #select_specie, select_red_dir
        """
        Plot phonon fatbands

        Args:
            title:
                Plot title.
            colormap
                Have a look at the colormaps here and decide which one you'd like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            max_stripe_width_mev:
                The maximum width of the stripe in meV.
            qlabels:
                dictionary whose keys are tuple with the reduced
                coordinates of the q-points. The values are the labels.
                e.g. qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.
            args:
                Positional arguments passed to matplotlib.

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
        # FIXME there's a bug in anaddb since we should orthogonalize
        # wrt the phonon displacement as done (correctly) here
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        structure = self.structure
        ntypat = structure.ntypesp

        # Grid with ntypat plots.
        nrows, ncols = (ntypat, 1)

        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True)
        xx = range(self.num_qpoints)

        # phonon_displacements are in cartesian coordinates and stored in an array with shape
        # (nqpt, 3*natom, 3*natom) where the last dimension stores the cartesian components.

        # Precompute normalization factor
        # d2(q,\nu) = \sum_{i=0}^{3*Nat-1) |d^{q\nu}_i|**2
        d2_qnu = np.zeros((self.num_qpoints, self.num_branches))
        for q in range(self.num_qpoints):
            for nu in self.branches:
                cvect = self.phdispl_cart[q, nu, :]
                d2_qnu[q, nu] = np.vdot(cvect, cvect).real

        # One plot per atom type.
        for (ax_idx, symbol) in enumerate(structure.symbol_set):
            ax = ax_list[ax_idx]
            if title and ax_idx == 0:
                ax.set_title(title)

            ax.grid(True)

            # dir_indices lists the coordinate indices for the atoms of the same type.
            atom_indices = structure.indices_from_symbol(symbol)
            dir_indices = []
            for aindx in atom_indices:
                start = 3 * aindx
                dir_indices.extend([start, start + 1, start + 2])

            for nu in self.branches:
                yy = self.phfreqs[:, nu]

                # Exctract the sub-vector associated to this atom type.
                displ_type = self.phdispl_cart[:, nu, dir_indices]
                d2_type = np.zeros(self.num_qpoints)
                for q in range(self.num_qpoints):
                    d2_type[q] = np.vdot(displ_type[q], displ_type[q]).real

                # Normalize and scale by max_stripe_width_mev.
                # The stripe is centered on the phonon branch hence the factor 2
                d2_type = max_stripe_width_mev * 1.e-3 * d2_type / (2. * d2_qnu[:, nu])

                # Plot the phonon branch and the stripe.
                color = plt.get_cmap(colormap)(float(ax_idx) / (ntypat - 1))
                if nu == 0:
                    ax.plot(xx, yy, lw=2, label=symbol, color=color)
                else:
                    ax.plot(xx, yy, lw=2, color=color)

                ax.fill_between(xx, yy + d2_type, yy - d2_type, facecolor=color, alpha=0.7, linewidth=0)

            ax.legend(loc="best")

            ylim = kwargs.pop("ylim", None)
            if ylim is not None:
                ax.set_ylim(ylim)

            ax.set_ylabel('Frequency [eV]')

            # Set ticks and labels.
            if qlabels is not None:
                ticks, labels = self._ticks_and_labels(qlabels)
                ax.set_xticks(ticks, minor=False)
                ax.set_xticklabels(labels, fontdict=None, minor=False)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def plot_with_dos(self, dos, qlabels=None, *args, **kwargs):
        """
        Plot the phonon band structure with the phonon DOS.

        Args:
            dos:
                An instance of :class:`PhononDOS`.
            qlabels:
                dictionary whose keys are tuple with the reduced
                coordinates of the q-points. The values are the labels.
                e.g. qlabels = {(0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L"}.
            args:
                Positional arguments passed to matplotlib.

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

        gspec = GridSpec(1, 2, width_ratios=[2, 1])

        ax1 = plt.subplot(gspec[0])
        # Align bands and DOS.
        ax2 = plt.subplot(gspec[1], sharey=ax1)

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the phonon band structure.
        for nu in self.branches:
            self.plot_ax(ax1, nu, *args, **kwargs)

        # Set ticks and labels.
        if qlabels is not None:
            ticks, labels = self._ticks_and_labels(qlabels)
            ax1.set_xticks(ticks, minor=False)
            ax1.set_xticklabels(labels, fontdict=None, minor=False)

        for ax in (ax1, ax2):
            ax.grid(True)

        if title:
            ax1.set_title(title)

        ax1.set_xlabel('q-point')
        ax1.set_ylabel('Energy [eV]')

        emin = np.min(self.minfreq)
        emin -= 0.05 * abs(emin)

        emax = np.max(self.maxfreq)
        emax += 0.05 * abs(emax)

        ax1.yaxis.set_view_interval(emin, emax)

        # Plot the DOS
        dos.plot_ax(ax2, what="d", exchange_xy=True, *args, **kwargs)

        ax2.yaxis.set_ticks_position("right")
        ax2.yaxis.set_label_position("right")

        if show:
            plt.show()

        fig = plt.gcf()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


class PHBST_Reader(ETSF_Reader):
    """This object reads data from PHBST.nc file produced by anaddb."""

    def read_qredcoords(self):
        """Array with the reduced coordinates of the q-points."""
        return self.read_value("qpoints")

    def read_qweights(self):
        """The weights of the q-points"""
        return self.read_value("qweights")

    def read_phfreqs(self):
        """Array with the phonon frequencies in eV."""
        return self.read_value("phfreqs")

    def read_phdispl_cart(self):
        """
        Complex array with the Cartesian displacements in Angstrom
        shape is (num_qpoints,  mu_mode,  cart_direction).
        """
        return self.read_value("phdispl_cart", cmode="c")


class PHBST_File(AbinitNcFile):
    def __init__(self, filepath):
        """
        Object used to access data stored in the PHBST file produced by ABINIT.

        Args:
            path:
                path to the file
        """
        super(PHBST_File, self).__init__(filepath)
        #
        # Initialize Phonon bands
        self.phbands = PhononBands.from_file(filepath)

    @property
    def structure(self):
        return self.phbands.structure

    def get_structure(self):
        return self.structure

    #def __str__(self):
    #    return self.tostring()

    #def tostring(self, prtvol=0):
    #    """
    #    String representation

    #    Args:
    #        prtvol:
    #            verbosity level.
    #    """
    #    return "\n".join(lines)

    def get_bands(self):
        """Return the electronic bands."""
        return self.bands

    def get_mode(self, qpoint, nu):
        """
        Returns he phonon with the given qpoint and branch nu.

        Args:
            qpoint:
                Either a vector with the reduced components of the q-point
                or an integer giving the sequential index (C-convention).
            nu:
                branch index (C-convention)

            returns:
                `PhononMode` instance.
        """
        if isinstance(qpoint, int):
            qidx = qpoint
        else:
            for (qidx, q) in enumerate(self.qpoints):
                if np.allclose(qpoint, q):
                    break
            else:
                raise ValueError("qpoint %s not found" % qpoint)

                #return PHMode

    def export_structure(self, path):
        """
        Export the structure structure on file filename.

        returns:
            Instance of :class:`Visualizer`
        """
        return self.structure.export(path)

    def visualize_structure_with(self, visualizer):
        """
        Visualize the crystalline structure with visualizer.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        extensions = Visualizer.exts_from_appname(visualizer)

        for ext in extensions:
            ext = "." + ext
            try:
                return self.export_structure(ext)
            except Visualizer.Error:
                pass
        else:
            raise Visualizer.Error("Don't know how to export data for visualizer " % visualizer)

            #def plot_bands(self, title=None, klabels=None, *args, **kwargs):
            #    self.bands.plot(title, klabels, *args, **kwargs)


#########################################################################################

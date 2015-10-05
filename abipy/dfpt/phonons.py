# coding: utf-8
from __future__ import print_function, division, unicode_literals

import sys
import functools
import numpy as np

from collections import OrderedDict
from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.core.units import Ha_to_eV, eV_to_Ha
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from abipy.core.func1d import Function1D
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_PhononBands
from abipy.core.kpoints import Kpoint, KpointList
from abipy.iotools import ETSF_Reader
from abipy.tools import gaussian
from abipy.tools.plotting_utils import Marker


__all__ = [
    "PhononBands",
    #"PhononBandsPlotter",
    "PhbstFile",
    "PhononDos",
    "PhdosReader",
    "PhdosFile",
]


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
            qpoint: qpoint in reduced coordinates.
            freq: Phonon frequency in eV.
            displ: Displacement (Cartesian coordinates, Angstrom)
            structure: :class:`Structure` object.
        """
        self.qpoint = Kpoint.as_kpoint(qpoint, structure.reciprocal_lattice)
        self.freq = freq
        self.displ_cart = displ_cart
        self.structure = structure

    # Rich comparison support (ordered is based on the frequency).
    # Missing operators are automatically filled by total_ordering.
    def __eq__(self, other):
        return self.freq == other.freq

    def __lt__(self, other):
        return self.freq < other.freq

    def __str__(self):
        return self.to_string()

    def to_string(self, with_displ=False):
        lines = ["%s: q-point %s, frequency %.5f [eV]" % (self.__class__.__name__, self.qpoint, self.freq)]
        app = lines.append

        if with_displ:
            app("Phonon displacement in cartesian coordinates [Angstrom]")
            app(str(self.displ_cart))

        return "\n".join(lines)

    #@property
    #def displ_red(self)
    #    return np.dot(self.xred, self.rprimd)

    #def export(self, path):
    #def visualize(self, visualizer):
    #def build_supercell(self):


class PhononBands(object):
    """
    Container object storing the phonon band structure.

    .. note::

        Frequencies are in eV. Cartesian displacements are in Angstrom.
    """
    @classmethod
    def from_file(cls, filepath):
        """Create the object from a netCDF file."""
        with PHBST_Reader(filepath) as r:
            structure = r.read_structure()

            # Build the list of q-points
            qpoints = KpointList(structure.reciprocal_lattice, 
                                 frac_coords=r.read_qredcoords(),
                                 weights=r.read_qweights(),
                                 names=None)
                                                                                   
            return cls(structure=structure,
                       qpoints=qpoints, 
                       phfreqs=r.read_phfreqs(),
                       phdispl_cart=r.read_phdispl_cart())

    def __init__(self, structure, qpoints, phfreqs, phdispl_cart, markers=None, widths=None,
                 non_anal_directions=None, non_anal_phfreqs=None):
        """
        Args:
            structure: :class:`Structure` object.
            qpoints: :class:`KpointList` instance.
            phfreqs: Phonon frequencies in eV.
            phdispl_cart: Displacement in Cartesian coordinates.
            markers: Optional dictionary containing markers labelled by a string.
                Each marker is a list of tuple(x, y, s) where x,and y are the position 
                in the graph and s is the size of the marker.
                Used for plotting purpose e.g. QP data, energy derivatives...
            widths: Optional dictionary containing data used for the so-called fatbands
                Each entry is an array of shape [nsppol, nkpt, mband] giving the width
                of the band at that particular point. Used for plotting purpose e.g. fatbands.
        """
        self.structure = structure

        #: :class:`KpointList` with the q-points
        self.qpoints = qpoints
        self.num_qpoints = len(self.qpoints)

        #: numpy array with phonon frequencies. Shape=(nqpt, 3*natom)
        self.phfreqs = phfreqs

        #: phonon displacements in Cartesian coordinates.
        #: `ndarray` of shape (nqpt, 3*natom, 3*natom).
        #: The last dimension stores the cartesian components.
        self.phdispl_cart = phdispl_cart

        # Handy variables used to loop.
        self.num_atoms = structure.num_sites
        self.num_branches = 3 * self.num_atoms
        self.branches = range(self.num_branches)

        if markers is not None:
            for key, xys in markers.items():
                self.set_marker(key, xys)
                                                                            
        if widths is not None:
            for key, width in widths.items():
                self.set_width(key, width)

        # numpy arrays with the non analytical directions and corresponding frequences calculated at Gamma
        self.non_anal_directions = non_anal_directions
        self.non_anal_phfreqs = non_anal_phfreqs

    def __str__(self):
        return self.to_string()

    def to_string(self, prtvol=0):
        """String representation."""
        lines = []
        app = lines.append

        for key, value in self.__dict__.items():
            if key.startswith("_"): continue
            if prtvol == 0 and isinstance(value, np.ndarray):
                continue
            app("%s = %s" % (key, value))

        return "\n".join(lines)

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

    @lazy_property
    def _auto_qlabels(self):
        # Find the q-point names in the pymatgen database.
        # We'll use _auto_klabels to label the point in the matplotlib plot
        # if qlabels are not specified by the user.
        _auto_qlabels = OrderedDict()
        for idx, qpoint in enumerate(self.qpoints):
            name = self.structure.findname_in_hsym_stars(qpoint)
            if name is not None:
                _auto_qlabels[idx] = name

        return _auto_qlabels

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

    @property
    def shape(self):
        """Shape of the array with the eigenvalues."""
        return self.num_qpoints, self.num_branches

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

    def set_marker(self, key, xys, extend=False):
        """
        Set an entry in the markers dictionary.

        Args:
            key: string used to label the set of markers.
            xys: Three iterables x,y,s where x[i],y[i] gives the
                positions of the i-th markers in the plot and
                s[i] is the size of the marker.
            extend:
                True if the values xys should be added to a pre-existing marker.
        """
        if not hasattr(self, "_markers"):
            self._markers = OrderedDict()

        if extend:
            if key not in self.markers:
                self._markers[key] = Marker(*xys)
            else:
                # Add xys to the previous marker set.
                self._markers[key].extend(*xys)
        else:
            if key in self.markers:
                raise ValueError("Cannot overwrite key %s in data" % key)

            self._markers[key] = Marker(*xys)

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

    def set_width(self, key, width):
        """
        Set an entry in the widths dictionary.

        Args:
            key: string used to label the set of markers.
            width: array-like of positive numbers, shape is [nqpt, num_modes].
        """
        width = np.reshape(width, self.shape)

        if not hasattr(self, "_widths"):
            self._widths = OrderedDict()

        if key in self.widths:
            raise ValueError("Cannot overwrite key %s in data" % key)

        if np.any(np.iscomplex(width)):
            raise ValueError("Found ambiguous complex entry %s" % str(width))

        if np.any(width < 0.0):
            raise ValueError("Found negative entry in width array %s" % str(width))

        self._widths[key] = width

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

    def qindex(self, qpoint):
        """Returns the index of the qpoint. Accepts integer or reduced coordinates."""
        if isinstance(qpoint, int):
            return qpoint
        else:
            return self.qpoints.index(qpoint)

    def qindex_qpoint(self, qpoint):
        """Returns (qindex, qpoint) from an integer or a qpoint."""
        qindex = self.qindex(qpoint)
        qpoint = self.qpoints(qindex)
        return qindex, qpoint

    def get_unstable_modes(self, below_mev=-5.0):
        """
        Return a list of :class:`PhononMode` objects with the unstable modes.
        A mode is unstable if its frequency is < below_mev. Output list is sorted
        and modes with lowest frequency come first.
        """
        umodes = []

        for iq, qpoint in enumerate(self.qpoints):
            for nu in self.branches:
                freq = self.phfreqs[iq, nu]
                if freq < below_mev * 1000:
                    displ_cart = self.phdispl_cart[iq, nu, :]
                    umodes.append(PhononMode(qpoint, freq, displ_cart, self.structure))

        return sorted(umodes)

    #def find_irreps(self, qpoint, tolerance):
    #    """
    #    Find the irreducible representation at this q-point
    #    Raise: QIrrepsError if algorithm fails
    #    """
    #    qindex, qpoint = self.qindex_qpoint(qpoint)

    def get_phdos(self, method="gaussian", step=1.e-4, width=4.e-4):
        """
        Compute the phonon DOS on a linear mesh.

        Args:
            method: String defining the method
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns:
            :class:`PhononDos` object.

        .. warning::

            Requires a homogeneous sampling of the Brillouin zone.
        """
        if abs(self.qpoints.sum_weights() - 1) > 1.e-6:
            raise ValueError("Qpoint weights should sum up to one")

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

        return PhononDos(mesh, values)

    def create_xyz_vib(self, iqpt, filename, pre_factor=200, do_real=True, scale_matrix=None, max_supercell=None):
        """
        Create vibration XYZ file for visualization of phonons.

        Args:
            iqpt: index of qpoint in self
            filename: name of the XYZ file that will be created
            pre_factor: Multiplication factor of the eigendisplacements
            do_real: True if we want only real part of the displacement, False means imaginary part
            scale_matrix: Scaling matrix of the supercell
            max_supercell: Maximum size of the supercell with respect to primitive cell
        """
        if scale_matrix is None:
            if max_supercell is None:
                raise ValueError("If scale_matrix is not provided, please provide max_supercell !")

            scale_matrix = self.structure.get_smallest_supercell(self.qpoints[iqpt].frac_coords, max_supercell=max_supercell)

        natoms = int(np.round(len(self.structure)*np.linalg.det(scale_matrix)))
        with open(filename, "w") as xyz_file:
            for imode in np.arange(self.num_branches):
                xyz_file.write(str(natoms) + "\n")
                xyz_file.write("Mode " + str(imode) + " : " + str(self.phfreqs[iqpt, imode]) + "\n")
                self.structure.write_vib_file(
                    xyz_file, self.qpoints[iqpt].frac_coords, pre_factor * np.reshape(self.phdispl_cart[iqpt, imode,:],(-1,3)), 
                    do_real=True, frac_coords=False, max_supercell=max_supercell, scale_matrix=scale_matrix)

    def decorate_ax(self, ax, **kwargs):
        title = kwargs.pop("title", None)
        if title is not None: ax.set_title(title)
                                                                   
        ax.grid(True)
        ax.set_xlabel('q-point')
        ax.set_ylabel('Energy [eV]')
        ax.legend(loc="best")
                                                                   
        # Set ticks and labels.
        ticks, labels = self._make_ticks_and_labels(kwargs.pop("qlabels", None))
                                                                   
        if ticks:
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False)

    @add_fig_kwargs
    def plot(self, ax=None, qlabels=None, branch_range=None, marker=None, width=None, **kwargs):
        """
        Plot the phonon band structure.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            qlabels: dictionary whose keys are tuple with the reduced coordinates of the q-points. 
                The values are the labels. e.g. qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}
            branch_range: Tuple specifying the minimum and maximum branch index to plot (default: all branches are plotted).
            marker: String defining the marker to plot. Syntax `markername:fact` where fact is a float used to scale the marker size.
            width: String defining the width to plot. Syntax `widthname:fact` where fact is a float used to scale the stripe size.

        Returns:
            `matplotlib` figure.
        """
        # Select the band range.
        if branch_range is None:
            branch_range = range(self.num_branches)
        else:
            branch_range = range(branch_range[0], branch_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, qlabels=qlabels)

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the phonon branches.
        for nu in branch_range:
            self.plot_ax(ax, nu, **kwargs)

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

            self.plot_width_ax(ax, key, fact=fact)

        return fig

    def plot_ax(self, ax, branch, **kwargs):
        """
        Plots the frequencies for the given branch index as a function of the q index on axis ax.
        If branch is None, all phonon branches are plotted.

        Return:
            The list of 'matplotlib' lines added.
        """
        branch_range = range(self.num_branches) if branch is None else [branch]

        first_xx = 0
        lines = []

        for pf in self.split_phfreqs:
            xx = range(first_xx, first_xx+len(pf))
            for branch in branch_range:
                lines.extend(ax.plot(xx, pf[:, branch], **kwargs))
            first_xx = xx[-1]

        return lines

    @property
    def split_qpoints(self):
        try:
            return self._split_qpoints
        except AttributeError:
            self._set_split_intervals()
            return self._split_qpoints

    @property
    def split_phfreqs(self):
        try:
            return self._split_phfreqs
        except AttributeError:
            self._set_split_intervals()
            return self._split_phfreqs

    def _set_split_intervals(self):
        # Calculations available for LO-TO splitting
        # Split the lines at each Gamma to handle possible discontinuities
        if self.non_anal_phfreqs is not None and self.non_anal_directions is not None:
            end_points_indeces = [0]

            end_points_indeces.extend(
                [i for i in range(1, self.num_qpoints - 1) if np.array_equal(self.qpoints.frac_coords[i], [0, 0, 0])])
            end_points_indeces.append(self.num_qpoints - 1)

            # split the list of qpoints and frequencies at each end point. The end points are in both the segments.
            # Lists since the array contained have different shapes
            split_qpoints = [np.array(self.qpoints.frac_coords[end_points_indeces[i]:end_points_indeces[i + 1] + 1])
                             for i in range(len(end_points_indeces) - 1)]
            split_phfreqs = [np.array(self.phfreqs[end_points_indeces[i]:end_points_indeces[i + 1] + 1])
                             for i in range(len(end_points_indeces) - 1)]

            for i, q in enumerate(split_qpoints):
                if np.array_equal(q[0], (0, 0, 0)):
                    split_phfreqs[i][0, :] = self._get_non_anal_freqs(q[1])
                if np.array_equal(q[-1], (0, 0, 0)):
                    split_phfreqs[i][-1, :] = self._get_non_anal_freqs(q[-2])
        else:
            split_qpoints = [self.qpoints]
            split_phfreqs = [self.phfreqs]

        self._split_qpoints = split_qpoints
        self._split_phfreqs = split_phfreqs
        return split_phfreqs, split_qpoints

    def _get_non_anal_freqs(self, direction):
        # directions for the qph2l in anaddb are given in cartesian coordinates
        direction = self.structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(direction)
        direction = direction / np.linalg.norm(direction)

        for i, d in enumerate(self.non_anal_directions):
            d = d / np.linalg.norm(d)
            if np.allclose(direction, d):
                return self.non_anal_phfreqs[i]

        raise ValueError("Non analytical contribution has not been calcolated for direction {} ".format(direction))

    def plot_width_ax(self, ax, key, branch=None, fact=1.0, **kwargs):
        """Helper function to plot fatbands for given branch on the axis ax."""
        branch_range = range(self.num_branches) if branch is None else [branch]

        facecolor = kwargs.pop("facecolor", "blue")
        alpha = kwargs.pop("alpha", 0.7)

        x, width = range(self.num_qpoints), fact * self.widths[key]

        for branch in branch_range:
           y, w = self.phfreq[:, branch], width[:,branch] * fact
           ax.fill_between(x, y-w/2, y+w/2, facecolor=facecolor, alpha=alpha)

    def plot_marker_ax(self, ax, key, fact=1.0):
        """Helper function to plot the markers for (spin,band) on the axis ax."""
        pos, neg = self.markers[key].posneg_marker()

        if pos:
            ax.scatter(pos.x, pos.y, s=np.abs(pos.s)*fact, marker="^", label=key + " >0")

        if neg:
            ax.scatter(neg.x, neg.y, s=np.abs(neg.s)*fact, marker="v", label=key + " <0")

    def _make_ticks_and_labels(self, qlabels):
        """Return ticks and labels from the mapping {qred: qstring} given in qlabels."""
        if qlabels is not None:
            d = OrderedDict()

            for qcoord, qname in qlabels.items():
                # Build Kpoint instancee
                qtick = Kpoint(qcoord, self.structure.reciprocal_lattice)
                for q, qpoint in enumerate(self.qpoints):
                    if qtick == qpoint:
                        d[q] = qname
        else:
            d = self._auto_qlabels

        # Return ticks, labels
        return list(d.keys()), list(d.values())

    @add_fig_kwargs
    def plot_fatbands(self, colormap="jet", max_stripe_width_mev=3.0, qlabels=None, **kwargs):
                      #select_specie, select_red_dir
        """
        Plot phonon fatbands

        Args:
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            max_stripe_width_mev: The maximum width of the stripe in meV.
            qlabels: dictionary whose keys are tuple with the reduced coordinates of the q-points. 
                The values are the labels. e.g. qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.

        Returns:
            `matplotlib` figure.
        """
        # FIXME there's a bug in anaddb since we should orthogonalize
        # wrt the phonon displacement as done (correctly) here
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

            self.decorate_ax(ax, qlabels=qlabels)

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

            ylim = kwargs.pop("ylim", None)
            if ylim is not None:
                ax.set_ylim(ylim)

        return fig

    @add_fig_kwargs
    def plot_with_phdos(self, dos, qlabels=None, **kwargs):
        """
        Plot the phonon band structure with the phonon DOS.

        Args:
            dos: An instance of :class:`PhononDos`.
            qlabels: dictionary whose keys are tuple with the reduced coordinates of the q-points. 
                The values are the labels e.g. qlabels = {(0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L"}.

        Returns:
            `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        gspec = GridSpec(1, 2, width_ratios=[2, 1])

        ax1 = plt.subplot(gspec[0])
        # Align bands and DOS.
        ax2 = plt.subplot(gspec[1], sharey=ax1)

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the phonon band structure.
        self.plot_ax(ax1, branch=None, **kwargs)

        self.decorate_ax(ax1, qlabels=qlabels)

        emin = np.min(self.minfreq)
        emin -= 0.05 * abs(emin)

        emax = np.max(self.maxfreq)
        emax += 0.05 * abs(emax)

        ax1.yaxis.set_view_interval(emin, emax)

        # Plot the DOS
        dos.plot_ax(ax2, what="d", exchange_xy=True, **kwargs)

        ax2.grid(True)
        ax2.yaxis.set_ticks_position("right")
        #ax2.yaxis.set_label_position("right")

        fig = plt.gcf()
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


class PhbstFile(AbinitNcFile, Has_Structure, Has_PhononBands):

    def __init__(self, filepath):
        """
        Object used to access data stored in the PHBST file produced by ABINIT.

        Args:
            path: path to the file
        """
        super(PhbstFile, self).__init__(filepath)
        self.reader = PHBST_Reader(filepath)

        # Initialize Phonon bands
        self._phbands = PhononBands.from_file(filepath)

    @property
    def structure(self):
        """:class:`Structure` object"""
        return self.phbands.structure

    @property
    def qpoints(self):
        return self.phbands.qpoints

    @property
    def phbands(self):
        """:class:`PhononBands` object"""
        return self._phbands

    def close(self):
        self.reader.close()

    def qindex(self, qpoint):
        """Returns the index of the qpoint. Accepts integer or reduced coordinates."""
        if isinstance(qpoint, int):
            return qpoint
        else:
            return self.qpoints.index(qpoint)

    def qindex_qpoint(self, qpoint):
        """Returns (qindex, qpoint) from an intege or a qpoint."""
        qindex = self.qindex(qpoint)
        qpoint = self.qpoints[qindex]
        return qindex, qpoint

    def get_phframe(self, qpoint):
        """
        Return a pandas :class:`DataFrame` with the phonon frequencies at the given q-point and 
        information on the crystal structure (used for convergence studies).

        Args:
            qpoint: integer, vector of reduced coordinates or :class:`Kpoint` object.
        """
        qindex, qpoint = self.qindex_qpoint(qpoint)
        phfreqs = self.phbands.phfreqs

        d = dict(
            omega=phfreqs[qindex, :],
            branch=list(range(3 * len(self.structure))),
        )

        # Add geo information
        structure = self.structure
        abc, angles = structure.lattice.abc, structure.lattice.angles
        d.update(dict(
            a=abc[0], b=abc[1], c=abc[2], volume=structure.volume,
            angle0=angles[0], angle1=angles[1], angle2=angles[2],
            formula=structure.formula,
        ))

        # Build the pandas Frame and add the q-point as attribute.
        import pandas as pd
        frame = pd.DataFrame(d, columns=d.keys())
        frame.qpoint = qpoint

        return frame

    def get_phmode(self, qpoint, branch):
        """
        Returns the :class:`PhononMode` with the given qpoint and branch nu.

        Args:
            qpoint: Either a vector with the reduced components of the q-point
                or an integer giving the sequential index (C-convention).
            branch: branch index (C-convention)

        Returns:
            :class:`PhononMode` instance.
        """
        qindex, qpoint = self.qindex_qpoint(qpoint)

        return PhononMode(qpoint=qpoint, 
                          freq=self.phbands.phfreqs[qindex, branch], 
                          displ_cart=self.phbands.phdispl_cart[qindex, branch, :],
                          structure=self.structure)


class PhononDos(Function1D):
    """
    This object stores the phonon density of states.
    An instance of `PhononDos` has a `mesh` (numpy array with the points of the mesh)
    and another numpy array, `values`, with the DOS on the mesh..

    .. note::

        mesh is given in eV, values are in states/eV.
    """
    #def __init__(self, mesh, values, qmesh):
    #    super(PhononDos, self).__init__(mesh, values)
    #    self.qmesh = qmesh

    @lazy_property
    def idos(self):
        """Integrated DOS."""
        return self.integral()

    #@lazy_property
    #def zpm(self)
    #    """Zero point motion"""
    #    return 0.5 * hbar * (self.dos.values * self.dos.mesh).integral[-1]

    def plot_ax(self, ax, what="d", exchange_xy=False, *args, **kwargs):
        """
        Helper function to plot the data on the axis ax.

        Args:
            ax: matplotlib axis
            what: string selecting the quantity to plot:
                "d" for DOS, "i" for IDOS. chars can be concatenated
                hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy: True to exchange exis
            args, kwargs:
                Options passes to matplotlib.

        Return:
            list of lines added to the plot
        """
        opts = [c.lower() for c in what]

        # Use super because we are overwriting the plot_ax provided by Func1D
        cases = {"d": super(PhononDos, self),  
                 "i": self.idos}

        lines = []
        for c in opts:
            f = cases[c]
            ls = f.plot_ax(ax, exchange_xy=exchange_xy, *args, **kwargs)
            lines.extend(ls)

        return lines

    @add_fig_kwargs
    def plot(self, *args, **kwargs):
        """
        Plot DOS and IDOS.

        Args:
            args: Positional arguments passed to :mod:`matplotlib`.

        Returns:
            `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        gspec = GridSpec(2, 1, height_ratios=[1, 2])
        ax1 = plt.subplot(gspec[0])
        ax2 = plt.subplot(gspec[1])

        for ax in (ax1, ax2):
            ax.grid(True)

        ax2.set_xlabel('Energy [eV]')
        ax1.set_ylabel("IDOS")
        ax2.set_ylabel("DOS")

        self.plot_ax(ax1, what="i", *args, **kwargs)
        self.plot_ax(ax2, what="d", *args, **kwargs)

        fig = plt.gcf()
        return fig

    def get_harmonic_thermo(self, tstart, tstop, num=50):
        """
        Compute thermodinamic properties from the phonon DOS within the harmonic approximation.

        tstart: The starting value (in Kelvin) of the temperature mesh. 
        tstop: The end value (in Kelvin) of the mesh.
        num: int, optional Number of samples to generate. Default is 50.
        """
        tmesh = np.linspace(tstart, tstop, num=num)

        coth = lambda x: 1.0 / np.tanh(x)
        csch2 = lambda x: 1.0 / (np.sinh(x) ** 2)

        # Boltzmann constant in Ha/K
        kb_HaK = 8.617343e-5 / Ha_to_eV 

        for i, gw in enumerate(self.values):
            if gw > 0: break
        #i = self.dos.find_mesh_index(0.0)

        # Use atomic units
        w =  self.mesh[i:] * eV_to_Ha
        gw = self.values[i:] * Ha_to_eV

        from scipy.integrate import cumtrapz
        def integrate(values):
            return cumtrapz(values, x=w)[-1]

        #w, gw = w[i:], gw[i:]
        # TODO
        # Check for possible numerical instabilities when w ~ 0 or negative
        # Prefactors are missing!
        df, de, cv, s = map(np.empty, 4 * (len(tmesh),))
        for i, temp in enumerate(tmesh):
            # Equations in Xavier's paper.
            kt = kb_HaK * temp
            wd2kt = w / (2 * kt)
            df[i] = kt * integrate(np.log(2 * np.sinh(wd2kt)) * gw)
            de[i] = integrate(w * coth(wd2kt) * gw)
            cv[i] = integrate(wd2kt * wd2kt * csch2(wd2kt) * gw)
            s[i] = integrate((wd2kt * coth(wd2kt) - np.log(2 * np.sinh(wd2kt))) * gw)

        locvars = locals()
        return HarmonicThermo(**{name: Function1D(tmesh, locvars[name]) for name in ("df", "de", "cv", "s")})


class HarmonicThermo(AttrDict):
    """
    This object gather the thermodinamic properties computed within the harmonic approximation.
    It also provides methods to plot the data and to compare two calculations.
    """
    LATEX_LABELS = dict(
        df="$\Delta F(T)$",
        de="$\Delta E(T)$",
        cv="$C_{v}(T)$",
        s="S(T)",
    )

    @add_fig_kwargs
    def plot(self, **kwargs):
        """Plot thermodinamic properties with matplotlib."""
        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=2, ncols=2, sharex=True, squeeze=True)
        ax_list = ax_list.ravel()

        for (name, func), ax in zip(self.items(), ax_list):
            func.plot_ax(ax)
            ax.grid(True)
            ax.set_xlabel('Temperature')
            ax.set_ylabel(self.LATEX_LABELS[name])

        return fig

    #def compare(self, other, tol=0.1, ret_diffs=False):
    #    """
    #    Compare two :class:`HarmonicThermo` objects

    #    Returns: 
    #        (iseq, diffs) 
    #        where iseq is True if all values are converged withing tol 
    #        and diffs is a :class:`AttrDict`that gives, for each property the  l2 norm of (f1 - f2)
    #    """
    #    # TODO: Relative diff or absolute diff?
    #    diffs, iseq = AttrDict(), True
    #    for name, self_func in self.items():
    #        other_func = other[name]
    #        diffs[name] = (self_func - other_func).l2_norm
    #        if diffs[name] > tol: iseq = False

    #    if ret_diffs: return iseq, diffs
    #    return iseq


class PhdosReader(ETSF_Reader):
    """
    This object reads data from the PHDOS.nc file produced by anaddb.

    .. note::
            Frequencies are in eV, DOSes are in states/eV.
    """
    @lazy_property
    def wmesh(self):
        """The frequency mesh in eV."""
        return self.read_value("wmesh")

    @lazy_property
    def pjdos_type(self):
        """DOS projected over atom types e.g. pjdos_type(ntypat,nomega)."""
        return self.read_value("pjdos_type")

    @lazy_property
    def pjdos_rc_type(self):
        """DOS projected over atom types and reduced directions e.g. pjdos_type(3,ntypat,nomega)."""
        return self.read_value("pjdos_rc_type")

    @lazy_property
    def pjdos(self):
        """DOS projected over atoms and reduced directions pjdos(natom,3,nomega)."""
        return self.read_value("pjdos")

    @lazy_property
    def structure(self):
        """The crystalline structure."""
        return self.read_structure()

    def read_phdos(self):
        """Return the :class:`PhononDOS`."""
        return PhononDos(self.wmesh, self.read_value("phdos"))

    def read_pjdos_type(self, symbol):
        """
        The contribution to the DOS due to the atoms of given chemical symbol.
        pjdos_type(ntypat,nomega)
        """
        type_idx = self.typeidx_from_symbol(symbol)
        return PhononDos(self.wmesh, self.pjdos_type[type_idx])

    # def read_pjdos(self, atom_idx=None):
    #     """
    #     projected DOS (over atoms)
    #     """
    #     return self.read_value("phonon_frequencies")

    # def read_pjdos_rc_type(self, symbol=None):
    #     """
    #     phdos(3,ntypat,nomega)
    #     phonon DOS contribution arising from a particular atom-type
    #     decomposed along the three reduced directions.
    #     """
    #     return self.read_value("phonon_frequencies")


class PhdosFile(AbinitNcFile, Has_Structure):
    """
    Container object storing the different DOSes stored in the
    PHDOS.nc file produced by anaddb. Provides helper function
    to visualize/extract data.
    """

    def __init__(self, filepath):
        # Open the file, read data and create objects.
        super(PhdosFile, self).__init__(filepath)

        self.reader = r = PhdosReader(filepath)
        self.wmesh = r.wmesh

    def close(self):
        self.reader.close()

    @lazy_property
    def structure(self):
        """Returns the :class:`Structure` object."""
        return self.reader.structure

    @lazy_property
    def phdos(self):
        return self.reader.read_phdos()

    @lazy_property
    def pjdos_type_dict(self):
        pjdos_type_dict = OrderedDict()
        for symbol in self.reader.chemical_symbols:
            #print(symbol, ncdata.typeidx_from_symbol(symbol))
            pjdos_type_dict[symbol] = self.reader.read_pjdos_type(symbol)

        return pjdos_type_dict

    @add_fig_kwargs
    def plot_pjdos_type(self, ax=None, colormap="jet", **kwargs):
        """
        Stacked Plot of the projected DOS (projection is for atom types)

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            colormap: Have a look at the colormaps 
                `here <http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html>`_
                and decide which one you'd like:

        Returns:
            matplotlib figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        ax.grid(True)

        xlim = kwargs.pop("xlim", None)
        if xlim is not None: ax.set_xlim(xlim)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None: ax.set_ylim(ylim)

        ax.set_xlabel('Frequency [eV]')
        ax.set_ylabel('PJDOS [states/eV]')

        # Type projected DOSes.
        num_plots = len(self.pjdos_type_dict)
        cumulative = np.zeros(len(self.wmesh))
        for i, (symbol, pjdos) in enumerate(self.pjdos_type_dict.items()):
            x, y = pjdos.mesh, pjdos.values
            color = plt.get_cmap(colormap)(float(i) / (num_plots - 1))
            ax.plot(x, cumulative + y, lw=2, label=symbol, color=color)
            ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=0.7)
            cumulative += y

        # Total PHDOS
        ax.plot(self.phdos.mesh, self.phdos.values, lw=2, label="Total PHDOS", color='black')

        ax.legend(loc="best")

        return fig


class PhononDosPlotter(object):
    """
    Class for plotting multiple phonon DOSes.
    """
    def __init__(self, *args):
        self._phdoses_dict = OrderedDict()
        for label, phdos in args:
            self.add_dos(label, phdos)

    def add_phdos(self, label, phdos):
        """
        Adds a DOS for plotting.

        Args:
            label: label for the phonon DOS. Must be unique.
            phdos: :class:`PhononDos` object.
        """
        if label in self._phdoses_dict:
            raise ValueError("label %s is already in %s" % (label, self._phdoses_dict.keys()))

        self._phdoses_dict[label] = phdos

    @add_fig_kwargs
    def plot(self, ax=None, *args, **kwargs):
        """
        Get a matplotlib plot showing the DOSes.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        xlim            x-axis limits. None (default) for automatic determination.
        ylim            y-axis limits.  None (default) for automatic determination.
        ==============  ==============================================================
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        ax.grid(True)

        xlim = kwargs.pop("xlim", None)
        if xlim is not None: ax.set_xlim(xlim)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None: ax.set_ylim(ylim)

        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('DOS [states/eV]')

        lines, legends = [], []
        for (label, dos) in self._phdoses_dict.items():
            l = dos.plot_ax(ax, *args, **kwargs)[0]

            lines.append(l)
            legends.append("DOS: %s" % label)

        # Set legends.
        ax.legend(lines, legends, loc='best', shadow=True)

        return fig

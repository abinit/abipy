# coding: utf-8
"""This module contains the class describing a planewave wavefunction."""
from __future__ import print_function, division, unicode_literals, absolute_import

import tempfile
import copy
import six
import itertools
import numpy as np

from monty.termcolor import cprint
from abipy.core import Mesh3D
from abipy.core.kpoints import Kpoint
from abipy.iotools import Visualizer
from abipy.iotools.xsf import xsf_write_structure, xsf_write_data
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt


__all__ = [
    "PWWaveFunction",
]

def latex_label_ispinor(ispinor, nspinor):
    if nspinor == 1:
        return ""
    elif nspinor == 2:
        return {k: v.replace("myuparrow", "uparrow") for k, v in
            {0: r"$\sigma=\myuparrow$", 1: r"$\sigma=\downarrow$"}.items()}[ispden]
    else:
        raise ValueError("Wrong value for nspinor: %s" % nspinor)


class WaveFunction(object):
    """
    Abstract class defining base and abstract methods for wavefunction objects.
    """
    def __eq__(self, other):
        if other is None: return False
        if self.gsphere != other.gsphere: return False
        return np.allclose(self.ug, other.ug)

    def __ne__(self, other):
        return not (self == other)

    def __iter__(self):
        """Yields G, ug[0:nspinor, G]"""
        if six.PY2:
            return itertools.izip(self.gvecs, self.ug.T)
        else:
            return zip(self.gvecs, self.ug.T)

    def __getitem__(self, slice):
        return self.gvecs[slice], self.ug[:, slice]

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.to_string()

    @property
    def shape(self):
        """Shape of ug i.e. (nspinor, npw)"""
        return self.nspinor, self.npw

    @property
    def gsphere(self):
        """:class:`GSphere` object"""
        return self._gsphere

    @property
    def kpoint(self):
        """|Kpoint| object"""
        return self.gsphere.kpoint

    @property
    def gvecs(self):
        """G-vectors in reduced coordinates."""
        return self.gsphere.gvecs

    @property
    def npw(self):
        """Number of G-vectors."""
        return len(self.gsphere)

    @property
    def ecut(self):
        """Cutoff energy in Hartree."""
        return self.gsphere.ecut

    @property
    def isnc(self):
        """True if norm-conserving wavefunction."""
        return isinstance(self, PWWaveFunction)

    @property
    def ispaw(self):
        """True if PAW wavefunction."""
        return isinstance(self, PAW_WaveFunction)

    @property
    def ug(self):
        """Periodic part of the wavefunctions in G-space."""
        return self._ug

    def set_ug(self, ug):
        """Set the value of the u(nspinor, G) array."""
        if ug.shape != self.shape:
            raise ValueError("Input ug shape %s differs from the one stored in self %s" % (ug.shape, self.shape))
        self._ug = ug
        self.delete_ur()

    @property
    def ur(self):
        """Periodic part of the wavefunction in real space."""
        try:
            return self._ur
        except AttributeError:
            self._ur = self.fft_ug()
            return self._ur

    def delete_ur(self):
        """Delete _u(r) (if it has been computed)."""
        try:
            del self._ur
        except AttributeError:
            pass

    #@property
    #def ur_xyz(self):
    #    """
    #    Returns a copy with ur[nspinor, nx, ny, nz]. Mainly used for post-processing.
    #    """
    #    return self.mesh.reshape(self.ur).copy()

    #@property
    #def ur2_xyz(self):
    #    """
    #    Returns ur2[nx, ny, nz]. Mainly used for post-processing.
    #    """
    #    return self.mesh.reshape(self.ur2)

    @property
    def mesh(self):
        """The mesh used for the FFT."""
        return self._mesh

    def set_mesh(self, mesh):
        """Change the FFT mesh. `u(r)` will be computed on this box."""
        assert isinstance(mesh, Mesh3D)
        self._mesh = mesh
        self.delete_ur()

    #def deepcopy(self):
    #    """Deep copy of self."""
    #    return copy.deepcopy(self)

    def get_ug_mesh(self, mesh=None):
        """
        Returns u(G) on the FFT mesh

        Args:
            mesh: |Mesh3d| object. If mesh is None, the internal mesh is used.
        """
        mesh = self.mesh if mesh is None else mesh
        return self.gsphere.tofftmesh(mesh, self.ug)

    def get_ur_mesh(self, mesh, copy=True):
        """
        Returns u(r) on the FFT mesh. This routine is mainly used if we need u(r) on a
        differen mesh. Data on the initial mesh is already available in `self.ur`

        Args:
            mesh: |Mesh3d| object.
            copy: By default, we return a copy of ur if mesh == self.mesh.
        """
        if mesh == self.mesh:
            return self.ur.copy() if copy else self.ur
        else:
            return self.fft_ug(mesh=mesh)

    def fft_ug(self, mesh=None):
        """
        Performs the FFT transform of :math:`u(g)` on mesh.

        Args:
            mesh: |Mesh3d| object. If mesh is None, self.mesh is used.

        Returns:
            :math:`u(r)` on the real space FFT box.
        """
        mesh = self.mesh if mesh is None else mesh
        ug_mesh = self.get_ug_mesh(mesh=mesh)
        return mesh.fft_g2r(ug_mesh, fg_ishifted=False)

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app("%s: nspinor: %d, spin: %d, band: %d " % (
            self.__class__.__name__, self.nspinor, self.spin, self.band))

        if hasattr(self, "gsphere"):
            app(self.gsphere.to_string(verbose=verbose))
        if hasattr(self, "mesh"):
            app(self.mesh.to_string(verbose=verbose))

        return "\n".join(lines)

    # TODO: get_ur2?
    @property
    def ur2(self):
        """
        [nx, ny, nz] array with :math:`||u(r)||^2` in real space.
        """
        ur2 = (self.ur.conj() * self.ur).real.copy()
        #if self.nspinor == 2: ur2 = ur2.sum(axis=3)
        return ur2


class PWWaveFunction(WaveFunction):
    """
    This object describes a wavefunction expressed in a plane-wave basis set.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PWWaveFunction
    """
    def __init__(self, structure, nspinor, spin, band, gsphere, ug):
        """
        Creation method.

        Args:
            structure: |Structure| object.
            nspinor: number of spinorial components.
            spin: spin index (only used if collinear-magnetism).
            band: band index (>=0)
            gsphere |GSphere| instance.
            ug: 2D array containing u[nspinor,G] for G in gsphere.
        """
        self.structure = structure
        self.nspinor, self.spin, self.band = nspinor, spin, band
        # Sanity check.
        assert ug.ndim == 2
        assert ug.shape[0] == nspinor
        assert ug.shape[1] == gsphere.npw

        self._gsphere = gsphere
        self._ug = np.array(ug)

    #def kinetic_energy(self):
    #    """Computes the matrix element of the kinetic operator in reciprocal space."""
    #    tug = -0.5 * self.gsphere.kpg2 * self.ug
    #    return np.vdot(self.ug, tug).sum()

    def norm2(self, space="g"):
        r"""
        Return :math:`||\psi||^2` computed in G- or r-space.

        space:  Integration space. Possible values ["g", "gsphere", "r"]
            if "g" or "r" the scalar product is computed in G- or R-space on the FFT box.
            if "gsphere" the integration is done on the G-sphere.
        """
        space = space.lower()

        if space == "g":
            return np.real(np.vdot(self.ug, self.ug))
        elif space == "gsphere":
            return np.real(np.vdot(self.ug, self.ug))
        elif space == "r":
            return np.vdot(self.ur, self.ur) / self.mesh.size
        else:
            raise ValueError("Wrong space: %s" % str(space))

    def braket(self, other, space="g"):
        """
        Returns the scalar product <u1|u2> of the periodic part of two wavefunctions
        computed in G-space or r-space, depending on the value of space.
        Note that selection rules introduced by k-points is not taken into accout.

        Args:
            other: Other wave (right-hand side)
            space:  Integration space. Possible values ["g", "gsphere", "r"]
                if "g" or "r" the scalar product is computed in G- or R-space on the FFT box.
                if "gsphere" the integration is done on the G-sphere. Note that
                this option assumes that self and other have the same list of G-vectors.
        """
        space = space.lower()

        if space == "g":
            ug1_mesh = self.gsphere.tofftmesh(self.mesh, self.ug)
            ug2_mesh = other.gsphere.tofftmesh(self.mesh, other.ug) if other is not self else ug1_mesh
            return np.vdot(ug1_mesh, ug2_mesh)
        elif space == "gsphere":
            return np.vdot(self.ug, other.ug)
        elif space == "r":
            return np.vdot(self.ur, other.ur) / self.mesh.size
        else:
            raise ValueError("Wrong space: %s" % str(space))

    def get_interpolator(self):
        """
        Return an interpolator object that interpolates periodic functions in real space.
        """
        from abipy.tools.numtools import BlochRegularGridInterpolator
        return BlochRegularGridInterpolator(self.structure, self.ur)

    #def pww_translation(self, gvector, rprimd):
    #    """Returns the pwwave of the kpoint translated by one gvector."""
    #    gsph = self.gsphere.copy()
    #    wpww = PWWaveFunction(self.structure, self.nspinor, self.spin, self.band, gsph, self.ug.copy())
    #    wpww.mesh = self.mesh
    #    wpww.pww_translation_inplace(gvector, rprimd)
    #    return wpww

    #def pww_translation_inplace(self, gvector, rprimd):
    #    """Translates the pwwave from 1 kpoint by one gvector."""
    #    self.gsphere.kpoint = self.gsphere.kpoint + gvector
    #    self.gsphere.gvecs = self.gsphere.gvecs + gvector
    #    fft_ndivs = (self.mesh.shape[0] + 2, self.mesh.shape[1] + 2, self.mesh.shape[2] + 2)
    #    newmesh = Mesh3D(fft_ndivs, rprimd, pbc=True)
    #    self.mesh = newmesh

    #def pwwtows_inplace(self):
    #    """Wrap the kpoint to the interval ]-1/2,1/2] and update pwwave accordingly."""
    #    kpoint = Kpoint(self.gsphere.kpoint, self.gsphere.gprimd)
    #    wkpt = kpoint.wrap_to_ws()

    #    if np.allclose(wkpt.rcoord, kpoint.rcoord):
    #        return

    #    #@David FIXME this is wrong
    #    gvector = np.array(kpoint.rcoord - wkpt.rcoord, np.int)
    #    self.gsphere.gvecs = self.gsphere.gvecs + gvector
    #    self.gsphere.kpoint = wkpt.rcoord

    #def pwwtows(self):
    #    """Return the pwwave of the kpoint wrapped to the interval ]-1/2,1/2]."""
    #    gsph = self.gsphere.copy()
    #    wpww = PWWaveFunction(self.structure, self.nspinor, self.spin, self.band, gsph, self.ug.copy())
    #    wpww.pwwtows_inplace()
    #    wpww.mesh = self.mesh
    #    return wpww

    #def rotate(self, symmop, mesh=None):
    #    """
    #    Apply the symmetry operation `symmop` to the periodic part.

    #    Args:
    #        symmop: :class:`Symmetry` operation
    #        mesh: mesh for the FFT, if None the mesh of self is used.

    #    Returns:
    #        New wavefunction object.
    #    """
    #    if self.nspinor != 1:
    #        raise ValueError("Spinor rotation not available yet.")

    #    rot_gsphere = self.gsphere.rotate(symmop)
    #    #rot_istwfk = istwfk(rot_kpt)

    #    if not np.allclose(symmop.tau, np.zeros(3)):
    #        rot_ug = np.empty_like(self.ug)
    #        rot_gvecs = rot_gsphere.gvecs
    #        rot_kpt = rot_gsphere.kpoint.frac_coords

    #        ug = self._ug
    #        #phase = np.exp(-2j * np.pi * (np.dot(rot_gvecs + rot_kpt, symmop.tau)))
    #        for ig in range(self.npw):
    #            rot_ug[:, ig] = ug[:, ig] * np.exp(-2j * np.pi * (np.dot(rot_gvecs[ig] + rot_kpt, symmop.tau)))
    #    else:
    #        rot_ug = self.ug.copy()

    #    # Invert the collinear spin if we have an AFM operation.
    #    rot_spin = self.spin
    #    if self.nspinor == 1:
    #        rot_spin = self.spin if symmop.is_fm else (self.spin + 1) % 2

    #    # Build new wave and set the mesh.
    #    new = self.__class__(self.nspinor, rot_spin, self.band, rot_gsphere, rot_ug)
    #    new.set_mesh(mesh if mesh is not None else self.mesh)
    #    return new

    @add_fig_kwargs
    def plot_line(self, point1, point2, num=200, with_krphase=False, cartesian=False,
                  ax=None, fontsize=12, **kwargs):
        """
        Plot (interpolated) wavefunction in real space along a line defined by ``point1`` and ``point2``.

        Args:
            point1: First point of the line. Accepts 3d vector or integer.
                The vector is in reduced coordinates unless ``cartesian`` is True.
                If integer, the first point of the line is given by the i-th site of the structure
                e.g. ``point1=0, point2=1`` gives the line passing through the first two atoms.
            point2: Second point of the line. Same API as ``point1``.
            num: Number of points sampled along the line.
            with_krphase: True to include the :math:`e^{ikr}` phase-factor.
            cartesian: By default, ``point1`` and ``point1`` are interpreted as points in fractional
                coordinates (if not integers). Use True to pass points in cartesian coordinates.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: legend and title fontsize.

        Return: |matplotlib-Figure|
        """
        # Interpolate along line.
        interpolator = self.get_interpolator()
        r = interpolator.eval_line(point1, point2, num=num, cartesian=cartesian,
                                   kpoint=None if not with_krphase else self.kpoint)
        # Plot data.
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        which = r"\psi(r)" if with_krphase else "u(r)"
        for ispinor in range(self.nspinor):
            spinor_label = latex_label_ispinor(ispinor, self.nspinor)
            ur = r.values[ispinor]
            ax.plot(r.dist, ur.real, label=r"$\Re %s$ %s" % (which, spinor_label))
            ax.plot(r.dist, ur.imag, label=r"$\Im %s$ %s" % (which, spinor_label))
            ax.plot(r.dist, ur.real**2 + ur.imag**2, label=r"$|\psi(r)|^2$ %s" % spinor_label)

        ax.grid(True)
        ax.set_xlabel("Distance from site1 [Angstrom]")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_line_neighbors(self, site_index, radius, num=200, with_krphase=False, max_nn=10, fontsize=12, **kwargs):
        """
        Plot (interpolated) density/potential in real space along the lines connecting
        an atom specified by ``site_index`` and all neighbors within a sphere of given ``radius``.

        .. warning:

            This routine can produce lots of plots! Be careful with the value of ``radius``.
            See also ``max_nn``.

        Args:
            site_index: Index of the atom in the structure.
            radius: Radius of the sphere in Angstrom.
            num: Number of points sampled along the line.
            with_krphase: True to include the :math:`e^{ikr}` phase-factor.
            max_nn: By default, only the first ``max_nn`` neighbors are showed.
            fontsize: legend and label fontsize.

        Return: |matplotlib-Figure|
        """
        site = self.structure[site_index]
        nn_list = self.structure.get_neighbors(site, radius, include_index=True)
        if not nn_list:
            cprint("Zero neighbors found for radius %s Ang. Returning None." % radius, "yellow")
            return None

        # Sort sites by distance.
        nn_list = list(sorted(nn_list, key=lambda t: t[1]))
        if max_nn is not None and len(nn_list) > max_nn:
            cprint("For radius %s, found %s neighbors but only max_nn %s sites are show." %
                    (radius, len(nn_list), max_nn), "yellow")
            nn_list = nn_list[:max_nn]

        # Get grid of axes (one row for neighbor)
        nrows, ncols = len(nn_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=True)
        ax_list = ax_list.ravel()

        interpolator = self.get_interpolator()
        kpoint = None if not with_krphase else self.kpoint
        which = r"\psi(r)" if with_krphase else "u(r)"

        # For each neighbor, plot psi along the line connecting site to nn.
        for i, (nn, ax) in enumerate(zip(nn_list, ax_list)):
            nn_site, nn_dist, nn_sc_index  = nn
            title = "%s, %s, dist=%.3f A" % (nn_site.species_string, str(nn_site.frac_coords), nn_dist)

            r = interpolator.eval_line(site.frac_coords, nn_site.frac_coords, num=num, kpoint=kpoint)

            for ispinor in range(self.nspinor):
                spinor_label = latex_label_ispinor(ispinor, self.nspinor)
                ur = r.values[ispinor]
                ax.plot(r.dist, ur.real, label=r"$\Re %s$ %s" % (which, spinor_label))
                ax.plot(r.dist, ur.imag, label=r"$\Im %s$ %s" % (which, spinor_label))
                ax.plot(r.dist, ur.real**2 + ur.imag**2, label=r"$|\psi(r)|^2$ %s" % spinor_label)

            ax.set_title(title, fontsize=fontsize)
            ax.grid(True)

            if i == nrows - 1:
                ax.set_xlabel("Distance from site_index %s [Angstrom]" % site_index)
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def export_ur2(self, filename, visu=None):
        """
        Export :math:`|u(r)|^2` to file ``filename``.

        Args:
            filename: String specifying the file path and the file format.
                The format is defined by the file extension. filename="prefix.xsf", for example,
                will produce a file in XSF format. An *empty* prefix, e.g. ".xsf" makes the code use a temporary file.
            visu: :class:`Visualizer` subclass. By default, this method returns the first available
                visualizer that supports the given file format. If visu is not None, an
                instance of visu is returned. See :class:`Visualizer` for the list of
                applications and formats supported.

        Returns:
            Instance of :class:`Visualizer`
        """
        if "." not in filename:
            raise ValueError("Cannot detect file extension in: %s" % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]:
            # fname == ".ext" ==> Create temporary file.
            # dir = os.getcwd() is needed when we invoke the method from a notebook.
            from abipy.core.globals import abinb_mkstemp
            _, filename = abinb_mkstemp(suffix="." + ext, text=True)
            print("Creating temporary file: %s" % filename)

        # Compute |u(r)|2 and write data according to ext.
        ur2 = np.reshape(self.ur2, (1,) + self.ur2.shape)

        with open(filename, mode="wt") as fh:
            if ext == "xsf":
                # xcrysden
                xsf_write_structure(fh, structures=self.structure)
                xsf_write_data(fh, self.structure, ur2, add_replicas=True)
            else:
                raise NotImplementedError("extension %s is not supported." % ext)

        if visu is None:
            return Visualizer.from_file(filename)
        else:
            return visu(filename)

    def visualize_ur2(self, appname="vesta"):
        """
        Visualize :math:`|u(r)|^2|`.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        # Get the Visualizer subclass from the string.
        visu = Visualizer.from_name(appname)

        # Try to export data to one of the formats supported by the visualizer
        # Use a temporary file (note "." + ext)
        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                return self.export_ur2(ext, visu=visu)
            except visu.Error:
                pass
        else:
            raise visu.Error("Don't know how to export data for %s" % str(appname))

    #def mvplot_cutplanes(self, show=True):
    #    data = self.ur2
    #    from abipy.display import mvtk
    #    figure, mlab = mvtk.get_fig_mlab(figure=None)
    #    source = mlab.pipeline.scalar_field(data)
    #    mlab.pipeline.image_plane_widget(source, plane_orientation='x_axes', slice_index=data.shape[0]//2)
    #    mlab.pipeline.image_plane_widget(source, plane_orientation='y_axes', slice_index=data.shape[1]//2)
    #    mlab.pipeline.image_plane_widget(source, plane_orientation='z_axes', slice_index=data.shape[2]//2)
    #    #mlab.pipeline.iso_surface(source, contours=contours) #, opacity=0.1)
    #    #mlab.pipeline.iso_surface(source, contours=[data.min()+ 0.1 * data.ptp()], opacity=0.1)
    #    mlab.outline()
    #    if show: mlab.show()
    #    return figure

    #def mvplot_volume(self, figure=None, vmin=0.65, vmax=0.9, show=True):
    #    from abipy.display import mvtk
    #    figure, mlab = mvtk.get_fig_mlab(figure=figure)
    #    data = self.ur2
    #    source = mlab.pipeline.scalar_field(data)
    #    data_min, data_max = data.min(), data.max()
    #    mlab.pipeline.volume(source,
    #                         vmin=data_min + vmin * (data_max - data_min),
    #                         vmax=data_min + vmax * (data_max - data_min))
    #    if show: mlab.show()
    #    return figure


class PAW_WaveFunction(WaveFunction):
    """
    All the methods that are related to the all-electron representation should start with ae.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PaW_WaveFunction
    """

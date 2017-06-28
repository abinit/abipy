# coding: utf-8
"""Fold2Bloch netcdf file."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from collections import OrderedDict
from monty.string import marquee # is_string, list_strings,
from monty.functools import lazy_property
import pymatgen.core.units as units
from pymatgen.core.lattice import Lattice
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import KpointList
from abipy.tools.plotting import set_axlims, add_fig_kwargs, get_ax_fig_plt
from abipy.electrons.ebands import ElectronsReader
from abipy.tools.numtools import gaussian


def dist_point_from_line(x0, x1, x2):
    """
    Return distance from point x0 to line x1 - x2. Cartesian coordinates are used.
    See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    """
    denom = x2 - x1
    denomabs = np.sqrt(np.dot(denom, denom))
    numer = np.cross(x0 - x1, x0 - x2)
    numerabs = np.sqrt(np.dot(numer, numer))
    return numerabs / denomabs


def find_points_along_path(cart_bounds, cart_coords, dist_tol, frac_coords):
    """
    Find points in `cart_coords` lying on the path defined by `cart_bounds`.

    Args:
        cart_bounds: [N, 3] array with the boundaries of the path in Cartesian coordinates.
        cart_coords: [M, 3] array with the points in Cartesian coordinate
        dist_tol: A point is considered to be on the path if its distance from the line
            is less than dist_tol.

    Return:
        (klines, dist_list, ticks)

        klines is a numpy array with the indices of the points lying on the path. Empty if no point is found.
        dist_list: numpy array with the distance of the points along the line.
        ticks:
    """
    klines, dist_list, ticks = [], [], [0]

    dl = 0  # cumulative length of the path
    for ibound, x0 in enumerate(cart_bounds[:-1]):
        x1 = cart_bounds[ibound + 1]
        B = x0 - x1
        #B = x1 - x0
        dk = np.sqrt(np.dot(B,B))
        #print("x0", x0, "x1", x1)
        ticks.append(ticks[ibound] + dk)
        for ik, k in enumerate(cart_coords):
            dist = dist_point_from_line(k, x0, x1)
            print(frac_coords[ik], dist)
            if dist > dist_tol: continue
            # k-point is on the cart_bounds
            A = x0 - k
            #A = k - x0
            x = np.dot(A,B)/dk
            #print("k-x0", A, "B", B)
            print(frac_coords[ik], x, x > 0 and x < dist_tol + dk)
            if dist_tol + dk >= x >= 0:
                # k-point is within the cart_bounds range
                # append k-point coordinate along the cart_bounds
                klines.append(ik)
                dist_list.append(x + dl)

        dl = dl + dk

    return np.array(klines), np.array(dist_list), np.array(ticks)


class Fold2BlochNcfile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    Netcdf file with output data produced by Fold2Bloch.

    Usage example:

    .. code-block:: python

        with Fold2BlochNcfile("foo_FOLD2BLOCH.nc") as fb:
            fb.plot_unfolded()
    """
    @classmethod
    def from_wfkpath(cls, wfkpath, folds, mpi_procs=1, workdir=None, manager=None, verbose=0):
        # Usage: $fold2Bloch file_WFK x:y:z (folds)
        # Build a simple manager to run the job in a shell subprocess
        manager = TaskManager.as_manager(manager).to_shell_manager(mpi_procs=mpi_procs)

        import tempfile
        workdir = tempfile.mkdtemp() if workdir is None else workdir

        # Run task in workdir
        fold2bloch = wrappers.Fold2Bloch(manager=manager, verbose=verbose)
        #fold2bloch.merge(self.outdir.path, ddb_files, out_ddb=out_ddb, description=desc)
        filepaths = [f for f in os.listdir(workdir) if f.endswith("_FOLD2BLOCH.nc")]
        if len(filepaths) != 1:
            raise RuntimeError("Cannot find *_FOLD2BLOCH.nc file in directory: %s" % os.listdir(workdir))

        return cls(os.path.join(workdir, filepaths[0]))

    def __init__(self, filepath):
        super(Fold2BlochNcfile, self).__init__(filepath)
        self.reader = ElectronsReader(filepath)

        # Initialize the electron bands from file.
        # Spectral weights are dimensioned with `nss`
	# Fortran arrays.
	# nctkarr_t("reduced_coordinates_of_unfolded_kpoints", "dp", "number_of_reduced_dimensions, nk_unfolded")
	# nctkarr_t("unfolded_eigenvalues", "dp", "max_number_of_states, nk_unfolded, number_of_spins")
	# nctkarr_t("spectral_weights", "dp", "max_number_of_states, nk_unfolded, nsppol_times_nspinor")
        self._ebands = self.reader.read_ebands()
        self.nss = max(self.nsppol, self.nspinor)
        self.folds = self.reader.read_value("folds")
        # Direct lattice of the primitive cell.
        self.pc_lattice = Lattice((self.structure.lattice.matrix.T * (1.0 / self.folds)).T.copy())

        # Read fold2bloch output data.
        self.uf_kfrac_coords = self.reader.read_value("reduced_coordinates_of_unfolded_kpoints")
        self.uf_kpoints = KpointList(self.pc_lattice.reciprocal_lattice, self.uf_kfrac_coords)
        self.uf_nkpt = len(self.uf_kpoints)

    def __str__(self):
        return self.to_string()

    def to_string(self, func=str, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(func(self.structure))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))

        app("Direct lattice of the primitive cell:")
        to_s = lambda x: "%0.6f" % x
        app("abc   : " + " ".join([to_s(i).rjust(10) for i in self.pc_lattice.abc]))
        app("angles: " + " ".join([to_s(i).rjust(10) for i in self.pc_lattice.angles]))
        app("Folds: %s" % (str(self.folds)))

        if verbose:
            app("Unfolded k-points:")
            app(self.uf_kpoints.to_string(func=func, verbose=verbose))

        return "\n".join(lines)

    @property
    def ebands(self):
        """:class:`ElectronBands` object with folded band energies."""
        return self._ebands

    @property
    def structure(self):
        """:class:`Structure` object defining the supercell."""
        return self.ebands.structure

    def close(self):
        self.reader.close()

    @lazy_property
    def uf_eigens(self):
        """[nsppol, nk_unfolded, nband] array with unfolded eigenvalues in eV."""
        # nctkarr_t("unfolded_eigenvalues", "dp", "max_number_of_states, nk_unfolded, number_of_spins")
        return self.reader.read_value("unfolded_eigenvalues") * units.Ha_to_eV

    @lazy_property
    def uf_weights(self):
        """[nss, nk_unfolded, nband] array with spectral weights. nss = max(nspinor, nsppol)."""
        # nctkarr_t("spectral_weights", "dp", "max_number_of_states, nk_unfolded, nsppol_times_nspinor")
        return self.reader.read_value("spectral_weights")

    def get_spectral_functions(self, step=0.01, width=0.02):
        """
        Args:
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Return:
            mesh, sfw, int_sfw
        """
        # Compute linear mesh.
        epad = 3.0 * width
        e_min = self.uf_eigens.min() - epad
        e_max = self.uf_eigens.max() + epad
        nw = int(1 + (e_max - e_min) / step)
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)

        sfw = np.zeros((self.nss, self.uf_nkpt, nw))
        for spin in range(self.nss):
            for ik in range(self.uf_nkpt):
                for band in range(self.nband):
                    e = self.uf_eigens[spin, ik, band]
                    sfw[spin, ik] += self.uf_weights[spin, ik, band] * gaussian(mesh, width, center=e)

        from scipy.integrate import cumtrapz
        int_sfw = cumtrapz(sfw, x=mesh, initial=0.0)

        return mesh, sfw, int_sfw

    @add_fig_kwargs
    def plot_unfolded(self, klabels=None, ylims=None, dist_tol=1e-12, verbose=0,
                      colormap="afmhot", facecolor="black", ax=None, **kwargs):
        """
        Plot unfolded band structure with spectral weights.

        Args:
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. `klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}`.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            dist_tol: A point is considered to be on the path if its distance from the line
                is less than dist_tol.
            verbose: Verbosity level.
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            facecolor:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
	"""
        uf_frac_coords = np.reshape([k.frac_coords for k in self.uf_kpoints], (-1, 3))
        cart_bounds = np.reshape([0.0, 1/2, 0, 0, 0, 0, 0, 0, 1/2], (-1, 3))
        cart_bounds = [self.pc_lattice.reciprocal_lattice.get_cartesian_coords(c) for c in cart_bounds]
        bound_labels = ["Y", "Gamma", "X"]
        uf_cart = np.reshape([k.cart_coords for k in self.uf_kpoints], (-1, 3))

        klines, xs, ticks = find_points_along_path(cart_bounds, uf_cart, dist_tol, uf_frac_coords)
        if len(klines) == 0: return None
        if verbose:
            fcoords = uf_frac_coords[klines]
            print("Found %d points along input k-path" % len(fcoords))
            print("k-points of path in reduced coordinates:")
            print(fcoords)

        fact = 8.0
        e0 = self.ebands.fermie
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.set_facecolor(facecolor)

        xs = np.tile(xs, self.nband)
        marker_spin = {0: "^", 1: "v"} if self.nss == 2 else {0: "o"}
        for spin in range(self.nss):
            ys = self.uf_eigens[spin, klines, :] - e0
            ws = self.uf_weights[spin, klines, :]
            s = ax.scatter(xs, ys.T, s=fact * ws.T, c=ws.T,
                           marker=marker_spin[spin], label=None if self.nss == 1 else "spin %s" % spin,
                           linewidth=1, edgecolors='none', cmap=plt.get_cmap(colormap))
            plt.colorbar(s, ax=ax, orientation='vertical')

        ax.set_xticks(ticks, minor=False)
        ax.set_xticklabels(bound_labels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))
        ax.grid(True)
        ax.set_ylabel('Energy [eV]')
        set_axlims(ax, ylims, "y")
        if self.nss == 2: ax.legend(loc="best")

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If `nbpath` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("f2b = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(f2b)"),
            nbv.new_code_cell("fig = f2b.ebands.plot()"),
            nbv.new_code_cell("# fig = f2b.unfolded_kpoints.plot()"),
            nbv.new_code_cell("fig = f2b.plot_unfolded()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

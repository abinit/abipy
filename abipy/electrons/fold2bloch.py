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


def dp2l(X0, X1, X2):
    # distance from point {X0} to line {X1}-{X2}
    # see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    denom = X2 - X1
    denomabs = np.sqrt(np.dot(denom, denom))
    numer = np.cross(X0-X1 , X0-X2)
    numerabs = np.sqrt(np.dot(numer, numer))
    return numerabs / denomabs


def find_points_onlines(cart_bounds, cart_coords, dist_tol, frac_coords):
    klines, dist_list, ticks = [], [], [0]

    dl = 0  # cumulative length of the path
    for ibound, x0 in enumerate(cart_bounds[:-1]):
        x1 = cart_bounds[ibound + 1]
        B = x0 - x1
        #B = x1 - x0
        dk = np.sqrt(np.dot(B,B))
        #print("x0", x0, "x1", x1)
        ticks.append(ticks[ibound] +dk)
        for ik, k in enumerate(cart_coords):
            dist = dp2l(k, x0, x1)
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
                #if x <= dist_tol and bound_labels is not None:
                #    key = round(x + dl, 12)
                #    ticks[key] = x + dl
                #    labels[key] =

        dl = dl + dk

    return np.array(klines), np.array(dist_list), np.array(ticks)


class Fold2BlochNcfile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    File containing the results of a ground-state calculation.

    Usage example:

    .. code-block:: python

        with Fold2BlochNcfile("foo_FOLD2BLOCH.nc") as fb:
            fb.ebands.plot()
    """
    #@classmethod
    #def from_wfkpath(cls, wfkpath, folds, mpi_procs=1, workdir=None, manager=None, verbose=0):
    #    # Run task in workdir
    #    # Usage: $fold2Bloch file_WFK x:y:z (folds)
    #    #manager = TaskManager.as_manager(manager).to_shell_manager(mpi_procs=mpi_procs)

    #    # Build a simple manager to run the job in a shell subprocess
    #    #import tempfile
    #    #workdir = tempfile.mkdtemp() if workdir is None else workdir

    #    #mrgddb = wrappers.Mrgddb(manager=manager, verbose=verbose)
    #    #mrgddb.merge(self.outdir.path, ddb_files, out_ddb=out_ddb, description=desc)

    #    return cls.from_directory(cls, workdir=workdir)

    def __init__(self, filepath):
        super(Fold2BlochNcfile, self).__init__(filepath)
        self.reader = Fold2BlochReader(filepath)

        # Initialize the electron bands from file
        self._ebands = self.reader.read_ebands()
        self.folds = self.reader.read_value("folds")
        # Direct lattice of the primitive cell
        self.pc_lattice = Lattice((self.structure.lattice.matrix.T * (1.0 / self.folds)).T.copy())

        self.uf_kfrac_coords = self.reader.read_value("reduced_coordinates_of_unfolded_kpoints")
        self.uf_kpoints = KpointList(self.pc_lattice.reciprocal_lattice, self.uf_kfrac_coords)
        self.uf_nkpt = len(self.uf_kpoints)
        return

        # Names of output files (depend on nsppol, nspinor)
        datafiles = []
        if self.nsppol == 1:
            if self.nspinor == 1:
                datafiles.append(seedname + ".f2b")
            else:
                # Weights for spin up/down spinor component.
                datafiles.append(seedname + "_SPOR_1.f2b")
                datafiles.append(seedname + "_SPOR_2.f2b")
        elif self.nsppol == 2:
            # Weights for spin up/down
            datafiles.append(seedname + "_UP.f2b")
            datafiles.append(seedname + "_DOWN.f2b")

        datafiles = [os.path.join(self.workdir, p) for p in datafiles]

        # Now Read output files.
	# Columns 1-3 correspond to kx, ky and kz of the unfolded bands;
	# the 4th column is the energy eigenvalue in [Ha] and the 5th column corresponds
	# to a spectral weight of the k-point after unfolding.
        data = []
        for isp in range(len(datafiles)):
            od = OrderedDict()
            with open(datafiles[isp], "rt") as fh:
                for line in fh:
                    tokens = line.split()
                    kstr, eig, w = " ".join(tokens[:3]), float(tokens[3]), float(tokens[4])
                    if kstr not in od: od[kstr] = []
                    od[kstr].append((eig, w))
            data.append(od)

	# Build numpy arrays.
        for spin, od in enumerate(data):
            kpoints = np.reshape([tuple(map(float, t.split())) for k in od], (-1, 3))
            if spin == 1:
                self.kpoints = kpoints
                self.nkpt = len(self.kpoints)
                self.uf_eigens = np.empty((self.nsppol, self.nkpt, self.nband))
                self.sf_weights = np.empty((self.nsppol * self.nspinor, self.nkpt, self.nband))
            else:
                if not np.all(self.kpoints == kpoints):
                    print("self.kpoints with shape:", self.kpoints.shape, "\nfrac_coords:\n", self.kpoints)
                    print("kpoints with shape:", kpoints.shape, "\nfrac_coords:\n", kpoints)
                    raise RuntimeError("Something wrong in the k-point parser")

            for ik, (kstr, tuple_list) in enumerate(od.items()):
                self.sf_weights[spin, ik] = [t[1] for t in tuple_list]
                if spin == 2 and self.nspinor == 2: continue
                self.uf_eigens[spin, ik] = [t[0] for t in tuple_list]

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(str(self.structure))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))

        app("Direct lattice of the primitive cell:")
        to_s = lambda x: "%0.6f" % x
        app("abc   : " + " ".join([to_s(i).rjust(10) for i in self.pc_lattice.abc]))
        app("angles: " + " ".join([to_s(i).rjust(10) for i in self.pc_lattice.angles]))
        app("Folds: %s" % (str(self.folds)))

        if verbose:
            app("Unfolded k-points:")
            app(self.uf_kpoints.to_string(verbose=verbose))

        return "\n".join(lines)

    @property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self._ebands

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    def close(self):
        self.reader.close()

    @property
    def uf_eigens(self):
        # nctkarr_t("unfolded_eigenvalues", "dp", "max_number_of_states, nk_unfolded, number_of_spins")
        return self.reader.read_value("unfolded_eigenvalues") * units.Ha_to_eV

    @property
    def uf_weights(self):
        # nctkarr_t("spectral_weights", "dp", "max_number_of_states, nk_unfolded, nsppol_times_nspinor")
        return self.reader.read_value("spectral_weights")

    def get_spectral_functions(self, step=0.1, width=0.2):
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

        nss = max(self.nsppol, self.nspinor)
        sfw = np.zeros((nss, self.uf_nkpt, nw))
        for spin in range(nss):
            for ik in range(self.uf_nkpt):
                for band in range(self.nband):
                    e = self.uf_eigens[spin, ik, band]
                    sfw[spin, ik] += self.uf_weights[spin, ik, band] * gaussian(mesh, width, center=e)

        from scipy.integrate import cumtrapz
        int_sfw = cumtrapz(sfw, x=mesh, initial=0.0)

        return mesh, sfw, int_sfw

    @add_fig_kwargs
    def plot_unfolded(self, klabels=None, ylims=None, ax=None, dist_tol=1e-12,
                      colormap="afmhot", facecolor="black", **kwargs):
        """
        Args:
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. `klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}`.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html

        Returns:
            `matplotlib` figure
	"""
        uf_frac_coords = np.reshape([k.frac_coords for k in self.uf_kpoints], (-1, 3))
        cart_bounds = np.reshape([0.0, 1/2, 0, 0, 0, 0, 0, 0, 1/2], (-1, 3))
        cart_bounds = [self.pc_lattice.reciprocal_lattice.get_cartesian_coords(c) for c in cart_bounds]
        bound_labels = ["Y", "Gamma", "X"]
        uf_cart = np.reshape([k.cart_coords for k in self.uf_kpoints], (-1, 3))

        klines, xs, ticks = find_points_onlines(cart_bounds, uf_cart, dist_tol, uf_frac_coords)
        if len(klines) == 0: return None

        print("k points of path")
        print(uf_frac_coords[klines])

        fact = 8.0
        e0 = self.ebands.fermie
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.set_facecolor(facecolor)

        xs = np.tile(xs, self.nband)
        nss = max(self.nsppol, self.nspinor)
        marker_spin = {0: "^", 1: "v"} if nss == 2 else {0: "o"}
        for spin in range(nss):
            ys = self.uf_eigens[spin, klines, :] - e0
            ws = self.uf_weights[spin, klines, :]
            s = ax.scatter(xs, ys.T, s=fact * ws.T, c=ws.T,
                           marker=marker_spin[spin], label=None if nss == 1 else "spin %s" % spin,
                           linewidth=1, edgecolors='none', cmap=plt.get_cmap(colormap))
                           #,vmin=0, vmax=1);
            plt.colorbar(s, ax=ax, orientation='vertical')

        ax.set_xticks(ticks, minor=False)
        ax.set_xticklabels(bound_labels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))
        ax.grid(True)
        ax.set_ylabel('Energy [eV]')
        set_axlims(ax, ylims, "y")
        if nss == 2: ax.legend(loc="best")

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("f2b = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(f2b)"),
            nbv.new_code_cell("fig = f2b.ebands.plot()"),
            #nbv.new_code_cell("fig = f2b.unfolded_kpoints.plot()"),
            nbv.new_code_cell("fig = f2b.plot_unfolded()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class Fold2BlochReader(ElectronsReader):
    """
    This object reads the results stored in the _GSR (Ground-State Results) file produced by ABINIT.
    It provides helper function to access the most important quantities.
    """

    # Fortran arrays.
    # nctkarr_t("reduced_coordinates_of_unfolded_kpoints", "dp", "number_of_reduced_dimensions, nk_unfolded")
    # nctkarr_t("unfolded_eigenvalues", "dp", "max_number_of_states, nk_unfolded, number_of_spins")
    # nctkarr_t("spectral_weights", "dp", "max_number_of_states, nk_unfolded, nsppol_times_nspinor")

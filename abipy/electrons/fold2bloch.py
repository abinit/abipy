# coding: utf-8
"""Fold2Bloch netcdf file."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from collections import OrderedDict
from monty.string import marquee # is_string, list_strings,
from monty.functools import lazy_property
from monty.collections import dict2namedtuple
from monty.termcolor import cprint
import pymatgen.core.units as units
from pymatgen.core.lattice import Lattice
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import KpointList, is_diagonal, find_points_along_path
from abipy.tools.plotting import set_axlims, add_fig_kwargs, get_ax_fig_plt
from abipy.electrons.ebands import ElectronsReader
from abipy.tools.numtools import gaussian


class Fold2BlochNcfile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    Netcdf file with output data produced by Fold2Bloch.

    Usage example:

    .. code-block:: python

        with Fold2BlochNcfile("foo_FOLD2BLOCH.nc") as fb:
            fb.plot_unfolded()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Fold2BlochNcfile
    """
    @classmethod
    def from_wfkpath(cls, wfkpath, folds, workdir=None, manager=None, mpi_procs=1, verbose=0):
        """
        Run fold2bloch in workdir.

        Args:
            wfkpath:
            folds:
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: :class:`TaskManager` of the task. If None, the manager is initialized from the config file.
            mpi_procs: Number of MPI processors to use.
            verbose: Verbosity level.
        """
        # Build a simple manager to run the job in a shell subprocess
        from abipy import flowtk
        manager = flowtk.TaskManager.as_manager(manager).to_shell_manager(mpi_procs=mpi_procs)
        fold2bloch = flowtk.Fold2Bloch(manager=manager, verbose=verbose)

        # Create temporary directory and link to the WFK file
        import tempfile
        workdir = tempfile.mkdtemp() if workdir is None else workdir
        wfkpath = os.path.abspath(wfkpath)
        link = os.path.join(workdir, os.path.basename(wfkpath))
        os.symlink(wfkpath, link)

        # Run fold2bloch
        ncpath = fold2bloch.unfold(link, folds, workdir=workdir)

        return cls(ncpath)

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
        self.fold_matrix = self.reader.read_value("fold_matrix")

        # Compute direct lattice of the primitive cell from fold_matrix.
        if is_diagonal(self.fold_matrix):
            folds = np.diagonal(self.fold_matrix)
            self.pc_lattice = Lattice((self.structure.lattice.matrix.T * (1.0 / folds)).T.copy())
        else:
            raise NotImplementedError("non diagonal fold_matrix: %s" % str(self.fold_matrix))

        # Read fold2bloch output data.
        self.uf_kfrac_coords = self.reader.read_value("reduced_coordinates_of_unfolded_kpoints")
        self.uf_kpoints = KpointList(self.pc_lattice.reciprocal_lattice, self.uf_kfrac_coords)
        self.uf_nkpt = len(self.uf_kpoints)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))

        app("Direct lattice of the primitive cell:")
        to_s = lambda x: "%0.6f" % x
        app("abc   : " + " ".join([to_s(i).rjust(10) for i in self.pc_lattice.abc]))
        app("angles: " + " ".join([to_s(i).rjust(10) for i in self.pc_lattice.angles]))
        if is_diagonal(self.fold_matrix):
            app("Diagonal folding: %s" % (str(np.diagonal(self.fold_matrix))))
        else:
            app("Folding matrix: %s" % str(self.fold_matrix))

        if verbose:
            app(self.uf_kpoints.to_string(verbose=verbose, title="Unfolded k-points"))

        if verbose > 2:
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    @property
    def ebands(self):
        """|ElectronBands| object with folded band energies."""
        return self._ebands

    @property
    def structure(self):
        """|Structure| object defining the supercell."""
        return self.ebands.structure

    def close(self):
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @lazy_property
    def uf_eigens(self):
        """[nsppol, nk_unfolded, nband] |numpy-array| with unfolded eigenvalues in eV."""
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

        return dict2namedtuple(mesh=mesh, sfw=sfw, int_sfw=int_sfw)

    @add_fig_kwargs
    def plot_unfolded(self, kbounds, klabels, ylims=None, dist_tol=1e-12, verbose=0,
                      colormap="afmhot", facecolor="black", ax=None, fontsize=12, **kwargs):
        r"""
        Plot unfolded band structure with spectral weights.

        Args:
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            dist_tol: A point is considered to be on the path if its distance from the line
                is less than dist_tol.
            verbose: Verbosity level.
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            facecolor:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
	"""
        cart_bounds = [self.pc_lattice.reciprocal_lattice.get_cartesian_coords(c)
                       for c in np.reshape(kbounds, (-1, 3))]
        uf_cart = self.uf_kpoints.get_cart_coords()

        p = find_points_along_path(cart_bounds, uf_cart, dist_tol)
        if len(p.ikfound) == 0:
            cprint("Warning: find_points_along_path returned zero points along the path. Try to increase dist_tol.", "yellow")
            return None
        if verbose:
            uf_frac_coords = np.reshape([k.frac_coords for k in self.uf_kpoints], (-1, 3))
            fcoords = uf_frac_coords[p.ikfound]
            print("Found %d points along input k-path" % len(fcoords))
            print("k-points of path in reduced coordinates:")
            print(fcoords)

        fact = 8.0
        e0 = self.ebands.fermie
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.set_facecolor(facecolor)

        xs = np.tile(p.dist_list, self.nband)
        marker_spin = {0: "^", 1: "v"} if self.nss == 2 else {0: "o"}
        for spin in range(self.nss):
            ys = self.uf_eigens[spin, p.ikfound, :] - e0
            ws = self.uf_weights[spin, p.ikfound, :]
            s = ax.scatter(xs, ys.T, s=fact * ws.T, c=ws.T,
                           marker=marker_spin[spin], label=None if self.nss == 1 else "spin %s" % spin,
                           linewidth=1, edgecolors='none', cmap=plt.get_cmap(colormap))
            plt.colorbar(s, ax=ax, orientation='vertical')

        ax.set_xticks(p.path_ticks, minor=False)
        ax.set_xticklabels(klabels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))
        ax.grid(True)
        ax.set_ylabel('Energy (eV)')
        set_axlims(ax, ylims, "y")
        if self.nss == 2: ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        print("TODO: Add call to plot_unfolded")
        yield self.ebands.plot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If `nbpath` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("f2b = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(f2b)"),
            nbv.new_code_cell("f2b.ebands.plot();"),
            nbv.new_code_cell("#f2b.unfolded_kpoints.plot();"),
            nbv.new_code_cell(r"""\
# kbounds = [0, 1/2, 0, 0, 0, 0, 0, 0, 1/2]
# klabels = ["Y", "$\Gamma$", "X"]
# f2b.plot_unfolded(kbounds, klabels);"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)

# coding: utf-8
"""
Interface to the ABIWAN netcdf file produced by abinit when calling wannier90 in library mode.
Inspired to the Fortran version of wannier90.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pandas as pd
import time

from collections import OrderedDict
from tabulate import tabulate
from monty.string import marquee
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import Kpath, IrredZone
from abipy.abio.robots import Robot
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt #, get_axarray_fig_plt
from abipy.electrons.ebands import ElectronBands, ElectronsReader, ElectronBandsPlotter, RobotWithEbands


class AbiwanFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    File produced by Abinit with the unitary matrices obtained by
    calling wannier90 in library mode.

    Usage example:

    .. code-block:: python

        with abilab.abiopen("foo_ABIWAN.nc") as abiwan:
            print(abiwan)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbiwanFile
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(AbiwanFile, self).__init__(filepath)

        self.reader = AbiwanReader(filepath)

        # Number of bands actually used to construct the Wannier functions
        self.num_bands_spin = self.reader.read_value("num_bands")

    @lazy_property
    def nwan_spin(self):
        """Number of Wannier functions for each spin."""
        return self.reader.read_value("nwan")

    @lazy_property
    def mwan(self):
        """
        Max number of Wannier functions over spins, i.e max(nwan_spin)
        Used to dimension arrays.
        """
        return self.reader.read_dimvalue("mwan")

    @lazy_property
    def nntot(self):
        """Number of k-point neighbours."""
        return int(self.reader.read_value("nntot"))

    @lazy_property
    def bands_in(self):
        """
        [nsppol, mband] logical array. Set to True if (spin, band) is included
        in the calculation. Set by exclude_bands
        """
        return self.reader.read_value("band_in_int").astype(np.bool)

    @lazy_property
    def lwindow(self):
        """
        [nsppol, nkpt, max_num_bands] array. Only if disentanglement.
        True if this band at this k-point lies within the outer window
        """
        return self.reader.read_value("lwindow_int").astype(np.bool)

    #@lazy_property
    #def ndimwin(self):
    #    """
    #    [nsppol, nkpt] array giving the number of bands inside the outer window for each k-point and spin.
    #    """
    #    return self.reader.read_value("ndimwin")

    @lazy_property
    def have_disentangled_spin(self):
        """[nsppol] bool array. Whether disentanglement has been performed."""
        #return self.reader.read_value("have_disentangled_spin").astype(np.bool)
        # TODO: Exclude bands
        return self.nwan_spin != self.num_bands_spin

    @lazy_property
    def wann_centers(self):
        """[nsppol, mwan, 3] array with Wannier centers in Ang."""
        return self.reader.read_value("wann_centres")

    @lazy_property
    def wann_spreads(self):
        """[nsppol, mwan] array with spreads in Ang^2"""
        return self.reader.read_value("wann_spreads")

    #@lazy_property
    #def spread_spin(self):
    #    """[nsppol, 3] array with spread in Ang^2"""
    #    return self.reader.read_value("spread")

    @lazy_property
    def irvec(self):
        """
        [nrpts, 3] array with the lattice vectors in the Wigner-Seitz cell
        in the basis of the lattice vectors defining the unit cell
        """
        return self.reader.read_value("irvec")

    @lazy_property
    def ndegen(self):
        """
        [nrpts] array with the degeneracy of each point.
        It will be weighted using 1 / ndegen[ir]
        """
        return self.reader.read_value("ndegen")

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        # TODO
        return od

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosity level verbose."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(title="Electronic Bands", with_structure=False, with_kpoints=True, verbose=verbose))
        app("")
        app(marquee("Wannier90 Results", mark="="))

        for spin in range(self.nsppol):
            if self.nsppol == 2: app("For spin: %d" % spin)
            app("No of Wannier functions: %d, No bands: %d, Number of k-point neighbours: %d" %
                (self.nwan_spin[spin], self.num_bands_spin[spin], self.nntot))
            app("Disentanglement: %s, exclude_bands: %s" %
                (self.have_disentangled_spin[spin], "no" if np.all(self.bands_in[spin]) else "yes"))
            app("")
            table = [["WF_index", "Center", "Spread"]]
            for iwan in range(self.nwan_spin[spin]):
                table.append([iwan,
                             "%s" % np.array2string(self.wann_centers[spin, iwan], precision=5),
                             "%.3f" % self.wann_spreads[spin, iwan]])
            app(tabulate(table, tablefmt="plain"))
            app("")

        if verbose and np.any(self.have_disentangled_spin):
            app(marquee("Lwindow", mark="="))
            app("[nsppol, nkpt, max_num_bands] array. True if state lies within the outer window.\n")
            for spin in range(self.nsppol):
                if self.nsppol == 2: app("For spin: %d" % spin)
                for ik in range(self.nkpt):
                    app("For ik: %d, %s" % (ik, self.lwindow[spin, ik]))
                app("")

        if verbose and np.any(self.bands_in):
            app(marquee("Bands_in", mark="="))
            app("[nsppol, mband] array. True if (spin, band) is included in the calculation. Set by exclude_bands.\n")
            for spin in range(self.nsppol):
                if self.nsppol == 2: app("For spin: %d" % spin)
                app("%s" % str(self.bands_in[spin]))
                app("")

        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))
            #app("irvec and ndegen")
            #for r, n in zip(self.irvec, self.ndegen):
            #    app("%s %s" % (r, n))

        return "\n".join(lines)

    def close(self):
        """Close file."""
        self.reader.close()

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def hwan(self):
        """
        Construct the matrix elements of the KS Hamiltonian in real space
        """
        start = time.time()

        nrpts, num_kpts = len(self.irvec), self.ebands.nkpt
        kfrac_coords = self.ebands.kpoints.frac_coords
        # Init datastructures needed by HWanR
        spin_rmn = [None] * self.nsppol
        spin_vmatrix = np.empty((self.nsppol, num_kpts), dtype=object)

        #kptopt = self.read.read_value("kptopt")
        #has_timrev =

        # Read unitary matrices from file.
        # Here be very careful with F --> C because we have to transpose.
        # complex U_matrix[nsppol, nkpt, mwan, mwan]
        u_matrix = self.reader.read_value("U_matrix", cmode="c")

        # complex U_matrix_opt[nsppol, mkpt, mwan, mband]
        if np.any(self.have_disentangled_spin):
            u_matrix_opt = self.reader.read_value("U_matrix_opt", cmode="c")

        for spin in range(self.nsppol):
            num_wan = self.nwan_spin[spin]

            # Real-space Hamiltonian H(R) is calculated by Fourier
            # transforming H(q) defined on the ab-initio reciprocal mesh
            HH_q = np.zeros((num_kpts, num_wan, num_wan), dtype=np.complex)

            for ik in range(num_kpts):
                eigs_k = self.ebands.eigens[spin, ik]
                # May have num_wan != mwan
                uk = u_matrix[spin, ik][:num_wan, :num_wan].transpose().copy()

                # Calculate the matrix that describes the combined effect of
                # disentanglement and maximal localization. This is the combination
                # that is most often needed for interpolation purposes
                # FIXME: problem with netcdf file
                if not self.have_disentangled_spin[spin]:
                    # [num_wann, num_wann] matrices, bands_in needed if exclude_bands
                    hks = np.diag(eigs_k[self.bands_in[spin]])
                    v_matrix = uk
                else:
                    # Select bands within the outer window
                    # TODO: Test if bands_in?
                    mask = self.lwindow[spin, ik]
                    hks = np.diag(eigs_k[mask])
                    #v_matrix = u_matrix_opt[spin, ik][:num_wan, mask].transpose() @ uk
                    v_matrix = np.matmul(u_matrix_opt[spin, ik][:num_wan, mask].transpose(), uk)

                #HH_q[ik] = v_matrix.transpose().conjugate() @ hks @ v_matrix
                HH_q[ik] = np.matmul(v_matrix.transpose().conjugate(), np.matmul(hks, v_matrix))
                spin_vmatrix[spin, ik] = v_matrix

            # Fourier transform Hamiltonian in the wannier-gauge representation.
            # O_ij(R) = (1/N_kpts) sum_q e^{-iqR} O_ij(q)
            rmn = np.zeros((nrpts, num_wan, num_wan), dtype=np.complex)
            j2pi = 2.0j * np.pi

            #for ir in range(nrpts):
            #   for ik, kfcs in enumerate(kfrac_coords):
            #      jqr = j2pi * np.dot(kfcs, self.irvec[ir])
            #      rmn[ir] += np.exp(-jqr) * HH_q[ik]
            #rmn *= (1.0 / num_kpts)

            for ik, kfcs in enumerate(kfrac_coords):
                jqr = j2pi * np.dot(self.irvec, kfcs)
                phases = np.exp(-jqr)
                rmn += phases[:, None, None] * HH_q[ik]
            rmn *= (1.0 / num_kpts)

            # Save results
            spin_rmn[spin] = rmn

        print("HWanR built in %.3f (s)" % (time.time() - start))
        return HWanR(self.structure, self.nwan_spin, spin_vmatrix, spin_rmn, self.irvec, self.ndegen)

    def interpolate_ebands(self, vertices_names=None, line_density=20, ngkpt=None, shiftk=(0, 0, 0), kpoints=None):
        """
        Build new |ElectronBands| object by interpolating the KS Hamiltonian with Wannier functions.
        Supports k-path via (vertices_names, line_density), IBZ mesh defined by ngkpt and shiftk
        or input list of kpoints.

        Args:
            vertices_names: Used to specify the k-path for the interpolated QP band structure
                List of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            line_density: Number of points in the smallest segment of the k-path. Used with ``vertices_names``.
            ngkpt: Mesh divisions. Used if bands should be interpolated in the IBZ.
            shiftk: Shifts for k-meshs. Used with ngkpt.
            kpoints: |KpointList| object taken e.g from a previous ElectronBands.
                Has precedence over vertices_names and line_density.

        Returns: |ElectronBands| object with Wannier-interpolated energies.
        """
        # Need KpointList object.
        if kpoints is None:
            if ngkpt is not None:
                # IBZ sampling
                kpoints = IrredZone.from_ngkpt(self.structure, ngkpt, shiftk, kptopt=1, verbose=0)
            else:
                # K-Path
                if vertices_names is None:
                    vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]
                kpoints = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)

        nk = len(kpoints)
        eigens = np.zeros((self.nsppol, nk, self.mwan))

        # Interpolate Hamiltonian for each kpoint and spin.
        start = time.time()
        write_warning = True
        for spin in range(self.nsppol):
            num_wan = self.nwan_spin[spin]
            for ik, kpt in enumerate(kpoints):
                oeigs = self.hwan.eval_sk(spin, kpt.frac_coords)
                eigens[spin, ik, :num_wan] = oeigs
                if num_wan < self.mwan:
                    # May have different number of wannier functions if nsppol == 2.
                    # Here I use the last value to fill eigens matrix (not very clean but oh well).
                    eigens[spin, ik, num_wan:self.mwan] = oeigs[-1]
                    if write_warning:
                        cprint("Different number of wannier functions for spin. Filling last bands with oeigs[-1]",
                               "yellow")
                        write_warning = False

        print("Interpolation completed in %.3f [s]" % (time.time() - start))
        occfacts = np.zeros_like(eigens)

        return ElectronBands(self.structure, kpoints, eigens, self.ebands.fermie,
                             occfacts, self.ebands.nelect, self.nspinor, self.nspden,
                             smearing=self.ebands.smearing)

    def get_plotter_from_ebands(self, ebands):
        """
        Interpolate energies using the k-points given in input |ElectronBands| ebands.

        Return: |ElectronBandsPlotter| object.
        """
        ebands = ElectronBands.as_ebands(ebands)
        wan_ebands = self.interpolate_ebands(kpoints=ebands.kpoints)

        return ElectronBandsPlotter(key_ebands=[("Ab-initio", ebands), ("Interpolated", wan_ebands)])

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.interpolate_ebands().plot(show=False)
        yield self.hwan.plot(show=False)
        if kwargs.get("verbose"):
            linestyle_dict = {"Interpolated": dict(linewidth=0, color="red", marker="o")}
            yield self.get_plotter_from_ebands(self.ebands).combiplot(linestyle_dict=linestyle_dict, show=False)
        else:
            cprint("Use verbose option to compare ab-initio points with interpolated values", "yellow")

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abiwan = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(abiwan.to_string(verbose=0))"),
            nbv.new_code_cell("abiwan.ebands.plot();"),
            nbv.new_code_cell("abiwan.ebands.kpoints.plot();"),
            nbv.new_code_cell("abiwan.hwan.plot();"),
            nbv.new_code_cell("#abiwan.get_plotter_from_ebands(;file_with_abinitio_ebands').combiplot();"),
            nbv.new_code_cell("ebands_kpath = abiwan.interpolate_ebands()"),
            nbv.new_code_cell("ebands_kpath.plot();"),
            nbv.new_code_cell("ebands_kmesh = abiwan.interpolate_ebands(ngkpt=[8, 8, 8])"),
            nbv.new_code_cell("ebands_kpath.plot_with_edos(ebands_kmesh.get_edos());"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


from abipy.core.skw import ElectronInterpolator
class HWanR(ElectronInterpolator):
    """
    This object represents the KS Hamiltonian in the wannier-gauge representation.
    It provides low-level methods to interpolate the KS eigenvalues, and a high-level API
    to interpolate bandstructures and plot the decay of the matrix elements in real space.

    # <0n|H|Rm>
    """

    def __init__(self, structure, nwan_spin, spin_vmatrix, spin_rmn, irvec, ndegen):
        self.structure = structure
        self.nwan_spin = nwan_spin
        self.spin_vmatrix = spin_vmatrix
        self.spin_rmn = spin_rmn
        self.irvec = irvec
        self.ndegen = ndegen
        self.nrpts = len(ndegen)
        self.nsppol = len(nwan_spin)
        assert self.nsppol == len(self.spin_rmn)
        assert len(self.irvec) == self.nrpts
        for spin in range(self.nsppol):
            assert len(self.spin_rmn[spin]) == self.nrpts

        # To call spglib
        self.cell = (self.structure.lattice.matrix, self.structure.frac_coords, self.structure.atomic_numbers)
        self.has_timrev = True
        self.verbose = 0
        self.nband = nwan_spin[0]
        #self.nelect

    def eval_sk(self, spin, kpt, der1=None, der2=None):
        """
        Interpolate eigenvalues for all bands at a given (spin, k-point).
        Optionally compute gradients and Hessian matrices.

        Args:
            spin: Spin index.
            kpt: K-point in reduced coordinates.
            der1: If not None, ouput gradient is stored in der1[nband, 3].
            der2: If not None, output Hessian is der2[nband, 3, 3].

        Return:
            oeigs[nband]
        """
        if der1 is not None or der2 is not None:
            raise NotImplementedError("Derivatives")

        # O_ij(k) = sum_R e^{+ik.R}*O_ij(R)
        j2pi = 2.0j * np.pi
        """
        num_wan = self.nwan_spin[spin]
        hk_ij = np.zeros((num_wan, num_wan), dtype=np.complex)
        for ir in range(self.nrpts):
            jrk = j2pi * np.dot(kpt, self.irvec[ir])
            hk_ij += self.spin_rmn[spin][ir] * (np.exp(jrk) / self.ndegen[ir])
        oeigs, _ = np.linalg.eigh(hk_ij)
        """

        # This is a bit faster.
        jrk = j2pi * np.dot(self.irvec, kpt)
        phases = np.exp(jrk) / self.ndegen
        hk_ij = (self.spin_rmn[spin] * phases[:, None, None]).sum(axis=0)
        oeigs, _ = np.linalg.eigh(hk_ij)

        return oeigs

    # TODO
    #def interpolate_omat(self, omat, kpoints):
    #def interpolate_sigres(self, sigres):
    #def interpolate_sigeph(self, sigeph):

    @add_fig_kwargs
    def plot(self, ax=None, fontsize=12, yscale="log", **kwargs):
        """
        Plot the matrix elements of the KS Hamiltonian in real space in the Wannier Gauge.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
            yscale: Define scale for y-axis. Passed to ax.set_yscale
            kwargs: options passed to ``ax.plot``.

        Return: |matplotlib-Figure|
        """
        # Sort R-points by length and build sortmap.
        irs = [ir for ir in enumerate(self.structure.lattice.norm(self.irvec))]
        items = sorted(irs, key=lambda t: t[1])
        sortmap = np.array([item[0] for item in items])
        rvals = np.array([item[1] for item in items])

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)
        ax.set_xlabel(r"$|R|$ (Ang)")
        ax.set_ylabel(r"$Max |H^W_{ij}(R)|$")

        marker_spin = {0: "^", 1: "v"}
        needs_legend = False
        for spin in range(self.nsppol):
            amax_r = [np.abs(self.spin_rmn[spin][ir]).max() for ir in range(self.nrpts)]
            amax_r = [amax_r[i] for i in sortmap]
            label = kwargs.get("label", None)
            if label is not None:
                label = "spin: %d" % spin if self.nsppol == 2 else None
            if label: needs_legend = True
            ax.plot(rvals, amax_r, marker=marker_spin[spin],
                    lw=kwargs.get("lw", 2),
                    color=kwargs.get("color", "k"),
                    markeredgecolor="r",
                    markerfacecolor="r",
                    label=label)

            ax.set_yscale(yscale)

        if needs_legend:
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig


class AbiwanReader(ElectronsReader):
    """
    This object reads the results stored in the ABIWAN file produced by ABINIT after having
    called wannier90 in library mode.
    It provides helper function to access the most important quantities.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbiwanReader
    """


class AbiwanRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple ABIWAN.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbiwanRobot
    """
    EXT = "ABIWAN"

    def get_dataframe(self, with_geo=True, abspath=False, funcs=None, **kwargs):
        """
        Return a |pandas-DataFrame| with the most important Wannier90 results.
        and the filenames as index.

        Args:
            with_geo: True if structure info should be added to the dataframe
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().

        kwargs:
            attrs:
                List of additional attributes of the |GsrFile| to add to the DataFrame.
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a |GsrFile| object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
        """
        # TODO
        # Add attributes specified by the users
        #attrs = [
        #    "energy", "pressure", "max_force",
        #    "ecut", "pawecutdg",
        #    "tsmear", "nkpt",
        #    "nsppol", "nspinor", "nspden",
        #] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for label, abiwan in self.items():
            row_names.append(label)
            d = OrderedDict()

            # Add info on structure.
            if with_geo:
                d.update(abiwan.structure.get_dict4pandas(with_spglib=True))

            #for aname in attrs:
            #    if aname == "nkpt":
            #        value = len(abiwan.ebands.kpoints)
            #    else:
            #        value = getattr(abiwan, aname, None)
            #        if value is None: value = getattr(abiwan.ebands, aname, None)
            #    d[aname] = value

            # Execute functions
            if funcs is not None: d.update(self._exec_funcs(funcs, abiwan))
            rows.append(d)

        row_names = row_names if not abspath else self._to_relpaths(row_names)
        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    @add_fig_kwargs
    def plot_hwanr(self, ax=None, colormap="jet", fontsize=8, **kwargs):
        """
        Plot the matrix elements of the KS Hamiltonian in real space in the Wannier Gauge.
        on the same Axes.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: matplotlib color map.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        for i, abiwan in enumerate(self.abifiles):
            abiwan.hwan.plot(ax=ax, fontsize=fontsize, color=cmap(i / len(self)), show=False)

        return fig

    def get_interpolated_ebands_plotter(self, vertices_names=None, knames=None, line_density=20,
            ngkpt=None, shiftk=(0, 0, 0), kpoints=None, **kwargs):
        """
        Args:
            vertices_names: Used to specify the k-path for the interpolated QP band structure
                It's a list of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            knames: List of strings with the k-point labels defining the k-path. It has precedence over `vertices_names`.
            line_density: Number of points in the smallest segment of the k-path. Used with ``vertices_names``.
            ngkpt: Mesh divisions. Used if bands should be interpolated in the IBZ.
            shiftk: Shifts for k-meshs. Used with ngkpt.
            kpoints: |KpointList| object taken e.g from a previous ElectronBands.
                Has precedence over vertices_names and line_density.

        Return: |ElectronBandsPlotter| object.
        """
        diff_str = self.has_different_structures()
        if diff_str: cprint(diff_str, "yellow")

        # Need KpointList object (assume same structures in Robot)
        nc0 = self.abifiles[0]
        if kpoints is None:
            if ngkpt is not None:
                # IBZ sampling
                kpoints = IrredZone.from_ngkpt(nc0.structure, ngkpt, shiftk, kptopt=1, verbose=0)
            else:
                # K-Path
                if knames is not None:
                    kpoints = Kpath.from_names(nc0.structure, knames, line_density=line_density)
                else:
                    if vertices_names is None:
                        vertices_names = [(k.frac_coords, k.name) for k in nc0.structure.hsym_kpoints]
                    kpoints = Kpath.from_vertices_and_names(nc0.structure, vertices_names, line_density=line_density)

        plotter = ElectronBandsPlotter()
        for label, abiwan in self.items():
            plotter.add_ebands(label, abiwan.interpolate_ebands(kpoints=kpoints))

        return plotter

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        p = self.get_eband_plotter()
        yield p.combiplot(show=False)
        yield p.gridplot(show=False)
        yield self.plot_hwanr(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.AbiwanRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.get_dataframe()"),
            nbv.new_code_cell("robot.plot_hwanr();"),
            nbv.new_code_cell("ebands_plotter = robot.get_interpolated_ebands_plotter()"),
            nbv.new_code_cell("ebands_plotter.ipw_select_plot()"),
        ])

        # Mixins
        #nb.cells.extend(self.get_baserobot_code_cells())
        #nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

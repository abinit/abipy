# coding: utf-8
"""Post-processing tools to analyze orbital magnetism."""
from __future__ import annotations

import numpy as np

from numpy.linalg import inv, det, eigvals
#from monty.termcolor import cprint
from monty.functools import lazy_property
from monty.string import list_strings, marquee
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.structure import Structure
from abipy.tools.numtools import BzRegularGridInterpolator
from abipy.electrons.ebands import ElectronBands, ElectronsReader
from abipy.core.kpoints import kpoints_indices, kmesh_from_mpdivs, map_grid2ibz
#from abipy.tools.numtools import gaussian
from abipy.tools.typing import Figure
from abipy.tools.plotting import set_axlims, get_ax_fig_plt, get_axarray_fig_plt, add_fig_kwargs, Marker


def print_options_decorator(**kwargs):
    """
    A decorator to temporarily set np.printoptions inside a function.
    Use the decorator to apply print settings inside the function.

    @print_options_decorator(precision=2, suppress=True)
    def print_array():
        print(np.array([1.234567, 7.891234, 3.456789]))
    """
    def decorator(func):
        import functools
        @functools.wraps(func)
        def wrapper(*args, **kwargs_inner):
            with np.printoptions(**kwargs):
                return func(*args, **kwargs_inner)
        return wrapper
    return decorator


class OrbmagAnalyzer:
    """
    This object gathers three ORBMAG.nc files, post-processes the data and
    provides tools to analyze/plot the results.

    Usage example:

    .. code-block:: python

        from abipy.electrons.orbmag import OrbmagAnalyzer
        orban = OrbmagAnalyzer(["gso_DS1_ORBMAG.nc", "gso_DS2_ORBMAG.nc", "gso_DS3_ORBMAG.nc"])
        print(orban)
        orban.report_eigvals(report_type="S")
        orban.plot_fatbands("bands_GSR.nc")
    """

    def __init__(self, filepaths: list):
        """
        Args:
            filepaths: List of filepaths to ORBMAG.nc files
        """
        if not isinstance(filepaths, (list, tuple)):
            raise TypeError(f"Expecting list or tuple with paths but got {type(filepaths)=}")

        if len(filepaths) != 3:
            raise ValueError(f"{len(filepaths)=} != 3")

        # @JOE TODO: One should store the direction in the netcdf file
        # so that we can check that the files are given in the right order.
        self.orb_files = [OrbmagFile(path) for path in filepaths]
        self.verbose = 0

        # This piece of code is taken from merge_orbmag_mesh. The main difference
        # is that here ncroots[0] is replaced by the reader instance of the first OrbmagFile.
        r0 = self.orb_files[0].r

        self.mband = mband = r0.read_dimvalue('mband')
        self.nkpt = nkpt = r0.read_dimvalue('nkpt')
        self.nsppol = nsppol = r0.read_dimvalue('nsppol')
        self.ndir = ndir = r0.read_dimvalue('ndir')
        orbmag_nterms = r0.read_dimvalue('orbmag_nterms')

        rprimd = r0.read_value('primitive_vectors')
        gprimd = inv(rprimd)
        ucvol = det(rprimd)

        # merging here means combine the 3 files: each delivers a 3-vector (ndir = 3),
        # output is a 3x3 matrix (ndir x ndir) for each term, sppol, kpt, band
        self.orbmag_merge_mesh = np.zeros((orbmag_nterms, ndir, ndir, nsppol, nkpt, mband))

        # here ncroots have been replaced by a list of reader instances.
        readers = [orb.r for orb in self.orb_files]

        for iband in range(mband):
            for isppol in range(nsppol):
                for ikpt in range(nkpt):
                    for iterm in range(orbmag_nterms):
                        for idir, r in enumerate(readers):
                            # convert terms to Cart coords; formulae differ depending on term. First
                            # four were in k space, remaining two already in real space
                            if iterm < 4:
                                omtmp = ucvol * np.matmul(gprimd, r.read_variable('orbmag_mesh')[iterm,0:ndir,isppol,ikpt,iband])
                            else:
                                omtmp = np.matmul(rprimd, r.read_variable('orbmag_mesh')[iterm,0:ndir,isppol,ikpt,iband])

                            self.orbmag_merge_mesh[iterm,idir,0:ndir,isppol,ikpt,iband] = omtmp

        # This piece of code has been taken from orbmag_sigij_mesh
        wtk = r0.read_value('kpoint_weights')
        occ = r0.read_value('occupations')
        orbmag_nterms = r0.read_dimvalue('orbmag_nterms')

        self.orbmag_merge_sigij_mesh = np.zeros((orbmag_nterms, nsppol, nkpt, mband, ndir, ndir))

        for iband in range(mband):
            for isppol in range(nsppol):
                for ikpt in range(nkpt):
                    # weight factor for each band and k point
                    trnrm = occ[isppol,ikpt,iband] * wtk[ikpt] / ucvol
                    for iterm in range(orbmag_nterms):
                        # sigij = \sigma_ij the 3x3 shielding tensor term for each sppol, kpt, and band
                        # additional ucvol factor converts to induced dipole moment (was dipole moment density,
                        # that is, magnetization)
                        self.orbmag_merge_sigij_mesh[iterm,isppol,ikpt,iband,0:ndir,0:ndir] = \
                            ucvol * trnrm * self.orbmag_merge_mesh[iterm,0:ndir,0:ndir,isppol,ikpt,iband]

    #def __str__(self) -> str:
    #    """String representation"""
    #    return self.to_string()

    #def to_string(self, verbose: int = 0) -> str:
    #    lines = []; app = lines.append
    #    return "\n".join(lines)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Activated at the end of the with statement. It automatically closes all the files."""
        self.close()

    def close(self):
        """Close all the files."""
        for orb in self.orb_files:
            orb.close()

    @lazy_property
    def structure(self):
        """Structure object."""
        # Perform consistency check
        structure = self.orb_files[0].structure
        if any(orb_file.structure != structure for orb_file in self.orb_files[1:]):
            raise RuntimeError("ORBMAG.nc files have different structures")
        return structure

    @lazy_property
    def has_timrev(self):
        """True if time-reversal symmetry is used in the BZ sampling."""
        has_timrev = self.orb_files[0].ebands.has_timrev
        if any(orb_file.ebands.has_timrev != has_timrev for orb_file in self.orb_files[1:]):
            raise RuntimeError("ORBMAG.nc files have different values of timrev")

    @print_options_decorator(precision=2, suppress=True)
    def report_eigvals(self, report_type):
        """
        """
        #np.set_printoptions(precision=2)

        terms = ['CC   ','VV1  ','VV2  ','NL   ','LR   ','A0An ']
        # Shape: (orbmag_nterms, nsppol, nkpt, mband, ndir, ndir)
        orbmag_merge_sigij_mesh = self.orbmag_merge_sigij_mesh

        total_sigij = orbmag_merge_sigij_mesh.sum(axis=(0, 1, 2, 3))
        eigenvalues = -1.0E6 * np.real(eigvals(total_sigij))
        isotropic = eigenvalues.sum() / 3.0
        span = eigenvalues.max() - eigenvalues.min()
        skew = 3.0 * (eigenvalues.sum() - eigenvalues.max() - eigenvalues.min() - isotropic) / span

        print('\nShielding tensor eigenvalues, ppm : ', eigenvalues)
        print('Shielding tensor iso, span, skew, ppm : %6.2f %6.2f %6.2f \n' % (isotropic,span,skew))

        if report_type == 'T':
            print('Term totals')
            term_sigij = orbmag_merge_sigij_mesh.sum(axis=(1, 2, 3))
            for iterm in range(np.size(orbmag_merge_sigij_mesh, axis=0)):
                eigenvalues = -1.0E6 * np.real(eigvals(term_sigij[iterm]))
                print(terms[iterm] + ': ', eigenvalues)
            print('\n')

        elif report_type == 'B':
            print('Band totals')
            band_sigij = orbmag_merge_sigij_mesh.sum(axis=(0, 1, 2))
            for iband in range(np.size(orbmag_merge_sigij_mesh, axis=3)):
                eigenvalues = -1.0E6 * np.real(eigvals(band_sigij[iband]))
                print('band ' + str(iband) + ' : ', eigenvalues)
            print('\n')

        elif report_type == 'TB':
            print('Terms in each band')
            tband_sigij = orbmag_merge_sigij_mesh.sum(axis=(1, 2))
            for iband in range(np.size(orbmag_merge_sigij_mesh, axis=3)):
                print('band ' + str(iband) + ' : ')
                for iterm in range(np.size(orbmag_merge_sigij_mesh, axis=0)):
                    eigenvalues = -1.0E6 * np.real(eigvals(tband_sigij[iterm, iband]))
                    print('   ' + terms[iterm] + ': ', eigenvalues)
                print('\n')
            print('\n')

    @lazy_property
    def ngkpt_and_shifts(self) -> tuple:
        """
        Return k-mesh divisions and shifts.
        """
        for ifile, orb in enumerate(self.orb_files):
            ksampling = orb.ebands.kpoints.ksampling
            ngkpt, shifts = ksampling.mpdivs, ksampling.shifts

            if ngkpt is None:
                raise ValueError("Non diagonal k-meshes are not supported!")

            if len(shifts) > 1:
                raise ValueError("Multiple shifts are not supported!")

            if ifile == 0:
                _ngkpt, _shifts = ngkpt, shifts
            else:
                # check that all files have the same value.
                if np.any(ngkpt != _ngkpt) or np.any(shifts != _shifts):
                    raise ValueError(f"ORBMAG files have different values of ngkpt: {ngkpt=} {_ngkpt=} or shifts {shifts=}, {_shifts=}")

        return ngkpt, shifts

    def get_value(self, what: str, spin: int, ikpt: int, band: int) -> float:
        """
        Postprocess orbmag_merge_sigij_mesh to compute the quantity associated
        to `what` for the given (spin, ikpt, band).
        """
        if what == "isotropic":
            vals = self.orbmag_merge_sigij_mesh[:, spin, ikpt, band].sum(axis=0)
            eigenvalues = -1.0E6 * vals
            return eigenvalues.sum() / 3.0
            #span = eigenvalues.max() - eigenvalues.min()
            #skew = 3.0 * (eigenvalues.sum() - eigenvalues.max() - eigenvalues.min() - isotropic) / span
        #elif what == "foobar":

        raise ValueError(f"Invalid {what=}")

    def insert_inbox(self, what: str, spin: int) -> tuple:
        """
        Return data, ngkpt, shifts where data is a (mband, nkx, nky, nkz)) array

        Args:
            spin: Spin index.
            what: Strings defining the quantity to insert in the box
        """
        # Need to know the shape of the k-mesh.
        ngkpt, shifts = self.ngkpt_and_shifts
        orb = self.orb_files[0]
        k_indices = kpoints_indices(orb.kpoints.frac_coords, ngkpt, shifts)
        nx, ny, nz = ngkpt

        # I'd start by weighting each band and kpt by trace(sigij)/3.0, the isotropic part of sigij,
        # both as a term-by-term plot and as a single weighting produced by summing over all 6 terms.
        #self.orbmag_merge_sigij_mesh = np.zeros((orbmag_nterms, nsppol, nkpt, mband, ndir, ndir))

        shape = (self.mband, nx, ny, nz)
        data = np.empty(shape, dtype=float)

        for iband in range(self.mband):
            for ikpt, k_inds in zip(range(self.nkpt), k_indices, strict=True):
                ix, iy, iz = k_inds
                value = self.get_value(what, spin, ikpt, iband)
                data[iband, ix, iy, iz] = value

        return data, ngkpt, shifts

    def get_skw_interpolator(self, what: str, lpratio: int, filter_params: None):
        """
        """
        orb = self.orb_files[0]
        ebands = orb.ebands.kpoints
        kpoints = orb.ebands.kpoints

        # Get symmetries from abinit spacegroup (read from file).
        if (abispg := self.structure.abi_spacegroup) is None:
            abispg = self.structure.spgset_abi_spacegroup(has_timerev=self.has_timrev)
        fm_symrel = [s for (s, afm) in zip(abispg.symrel, abispg.symafm, strict=True) if afm == 1]

        from abipy.core.skw import SkwInterpolator
        cell = (self.structure.lattice.matrix, self.structure.frac_coords, self.structure.atomic_numbers)

        interp_spin = [None for _ in range(self.nsppol)]

        for spin in range(self.nsppol):
            values_kb = np.empty((self.nkpt, self.mband))
            for ikpt in range(self.nkpt):
                for iband in range(self.mband):
                    values_kb[ikpt, iband] = self.get_value(what, spin, ikpt, band)

            skw = SkwInterpolator(lpratio, kpoints.frac_coords, self.eigens[:,:,bstart:bstop], ebands.fermie, ebands.nelect,
                                  cell, fm_symrel, self.has_timrev,
                                  filter_params=filter_params, verbose=self.verbose)
            interp_spin[spin] = skw
            #skw.eval_sk(spin, kpt, der1=None, der2=None) -> np.ndarray:

        return interp_spin

    @lazy_property
    def has_full_bz(self) -> bool:
        """True if the list of k-points cover the full BZ."""
        ngkpt, shifts = self.ngkpt_and_shifts
        return np.product(ngkpt) * len(shifts) == self.nkpt

    def get_bz_interpolator_spin(self, what: str, interp_method: str) -> list[BzRegularGridInterpolator]:
        """
        Build and return a list with nsppol interpolators.

        Args:
            what: Strings defining the quantity to insert in the box.
            interp_method: The method of interpolation. Supported are “linear”, “nearest”,
                “slinear”, “cubic”, “quintic” and “pchip”.
        """
        # This part is tricky as we have to handle different cases:
        #
        # 1) k-mesh defined in terms of ngkpt and one shift.
        #    In this case we can interpolate easily although we have to take into account
        #    if the k-points have been reduced by symmetry or not.
        #    If we have the full BZ, we can use BzRegularGridInterpolator
        #    else we have to resort to StarFunction interpolation
        #
        # 2) k-mesh defined in terms of ngkpt and multiple shifts.
        #    If we have the IBZ, we have to resort to StarFunction interpolation
        #    Handling the full BZ is not yet possible due to an internal limitation of BzRegularGridInterpolator

        interp_spin = [None for _ in range(self.nsppol)]

        if self.has_full_bz:
            for spin in range(self.nsppol):
                data, ngkpt, shifts = self.insert_inbox(what, spin)
                interp_spin[spin] = BzRegularGridInterpolator(self.structure, shifts, data, method=interp_method)
        else:
            raise NotImplementedError("k-points must cover the full BZ.")

        return interp_spin

    @add_fig_kwargs
    def plot_fatbands(self, ebands_kpath,
                      what_list="isotropic",
                      ylims=None, scale=50, marker_color="gold", marker_edgecolor="gray",
                      marker_alpha=0.5, fontsize=12, interp_method="linear",
                      ax_mat=None, **kwargs) -> Figure:
        """
        Plot fatbands ...

        Args:
            ebands_kpath: ElectronBands instance with energies along a k-path
                or path to a netcdf file providing it.
            what_list: string or list of strings defining the quantity to show.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
            scale: Scaling factor for fatbands.
            marker_color: Color for markers
            marker_edgecolor: Color for marker edges.
            marker_edgecolor: Marker transparency.
            fontsize: fontsize for legends and titles
            interp_method: The method of interpolation. Supported are “linear”, “nearest”,
                “slinear”, “cubic”, “quintic” and “pchip”.
            ax_mat: matrix of |matplotlib-Axes| or None if a new figure should be created.
        """
        what_list = list_strings(what_list)

        # Build the grid of plots.
        num_plots, ncols, nrows = len(what_list), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        ebands_kpath = ElectronBands.as_ebands(ebands_kpath)

        for what, ax in zip(what_list, ax_list, strict=True):
            # Get interpolator for `what` quantity.
            interp_spin = self.get_bz_interpolator_spin(what, interp_method)

            abs_max = max((interp.get_max_abs_data() for interp in interp_spin))
            scale *= 1. / abs_max

            ymin, ymax = +np.inf, -np.inf
            x, y, s = [], [], []

            for spin in range(self.nsppol):
                for ik, kpoint in enumerate(ebands_kpath.kpoints):
                    enes_n = ebands_kpath.eigens[spin, ik]
                    for e, a2 in zip(enes_n, interp_spin[spin].eval_kpoint(kpoint), strict=True):
                        x.append(ik); y.append(e); s.append(scale * abs(a2))
                        ymin, ymax = min(ymin, e), max(ymax, e)

            # Plot electron bands with markers.
            points = Marker(x, y, s, color=marker_color, edgecolors=marker_edgecolor,
                            alpha=marker_alpha, label=what)

            ebands_kpath.plot(ax=ax, points=points, show=False, linewidth=1.0)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        e0 = self.orb_files[0].ebands.fermie
        if ylims is None:
            # Automatic ylims.
            span = ymax - ymin
            ymin -= 0.1 * span
            ymax += 0.1 * span
            ylims = [ymin - e0, ymax - e0]

        for ax in ax_list:
            set_axlims(ax, ylims, "y")

        return fig


class OrbmagFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands):
    """
    Interface to the ORBMAG.nc file.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram::
    """

    @classmethod
    def from_file(cls, filepath: str) -> OrbmagFile:
        """Initialize the object from a netcdf_ file"""
        return cls(filepath)

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.r = ElectronsReader(filepath)

    @lazy_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def params(self) -> dict:
        """dict with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    def close(self) -> None:
        """Called at the end of the ``with`` context manager."""
        return self.r.close()

    def __str__(self) -> str:
        """String representation"""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=True, title="Electronic Bands"))
        app("")
        #app(marquee("Fatbands Info", mark="="))
        #app("prtdos: %d, prtdosm: %d, mbesslang: %d, pawprtdos: %d, usepaw: %d" % (
        #    self.prtdos, self.prtdosm, self.mbesslang, self.pawprtdos, self.usepaw))
        app("nsppol: %d, nkpt: %d, mband: %d" % (self.nsppol, self.nkpt, self.mband))
        app("")

        if verbose > 1:
            app("")
            app(self.hdr.to_str(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

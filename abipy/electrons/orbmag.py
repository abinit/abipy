# coding: utf-8
"""Post-processing tools to analyze orbital magnetism."""
from __future__ import annotations

import numpy as np

from numpy.linalg import inv, det, eig, eigvals, norm
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

import matplotlib.patches as mpatches

def filter_sigma(sigma, atol):
    """from raw sigma_ij output, convert to ppm shielding and remove small parts"""
    filtered_sigma = -1.0E6 * sigma
    for row in range(3):
        for col in range(3):
            if abs(filtered_sigma[row,col]) < atol:
                filtered_sigma[row,col] = 0.0
    return filtered_sigma


def make_U_from_sigma(gprimd, sigma, atol=0.01):
    """
    """
    fsigma = filter_sigma(sigma, atol)
    sigma_eigs, sigma_evecs = eig(fsigma)
    Dmat = np.array([[norm(gprimd[0]), 0, 0], [0, norm(gprimd[1]), 0], [0, 0, norm(gprimd[2])]])
    intmat=np.matmul(np.transpose(sigma_evecs), gprimd)

    Cmat = np.zeros((3,3))
    for col in range(3):
        for row in range(3):
            Cmat[row,col] = intmat[row,col] / Dmat[col,col]

    Umat= np.matmul(np.transpose(Cmat), np.matmul(fsigma,Cmat))
    Umat = Umat/np.max(Umat)
    return Umat


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

    def __init__(self, filepaths: list, verbose=0):
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
        self.verbose = verbose

        inds = np.array([o.target_atom for o in self.orb_files], dtype=int)
        if not np.all(inds == inds[0]):
            raise RuntimeError(f"ORBMAG.nc files have dipoles on different atoms: {inds}")
        self.target_atom = inds[0]

        # check that the three dipole vectors form a right-handed set
        target_atom_list = [o.target_atom for o in self.orb_files]
        vec_a, vec_b, vec_c = [o.nucdipmom_atom for o in self.orb_files]
        if np.dot(np.cross(vec_a, vec_b), vec_c) < 0.0:
            raise RuntimeError(f"nuclear dipoles do not form a right handed set:\n{vec_a=}\n{vec_b=}\n{vec_c=} ")

        # This piece of code is taken from merge_orbmag_mesh. The main difference
        # is that here ncroots[0] is replaced by the reader instance of the first OrbmagFile.
        r0 = self.orb_files[0].r

        def _read_dim(dimname):
            """Read dimension dimname from files and check for equality."""
            dims = np.array([o.r.read_dimvalue(dimname) for o in self.orb_files], dtype=int)
            if not np.all(dims == dims[0]):
                raise RuntimeError(f"For {dimname=}, files have different dimensions {dims=}")
            return int(dims[0])

        def _read_value(varname):
            """Read variable varname from files and check if they are almost equal."""
            import numpy.testing as nptu
            vals_list = np.array([o.r.read_value(varname) for o in self.orb_files])
            for vals in vals_list[1:]:
                nptu.assert_almost_equal(vals, vals_list[0]) # , decimal, err_msg, verbose)
            return vals_list[0].copy()

        self.mband = mband = _read_dim('mband')
        self.nkpt = nkpt = _read_dim('nkpt')
        self.nsppol = nsppol = _read_dim('nsppol')
        self.ndir = ndir = _read_dim('ndir')
        orbmag_nterms = _read_dim('orbmag_nterms')

        self.rprimd = rprimd = _read_value('primitive_vectors')
        self.gprimd = gprimd = inv(rprimd)
        self.ucvol = ucvol = det(rprimd)

        # merging here means combine the 3 files: each delivers a 3-vector (ndir = 3),
        # output is a 3x3 matrix (ndir x ndir) for each term, sppol, kpt, band
        self.orbmag_merge_mesh = np.zeros((orbmag_nterms, ndir, ndir, nsppol, nkpt, mband))

        # here ncroots have been replaced by a list of reader instances.
        readers = [orb.r for orb in self.orb_files]

        for idir, r in enumerate(readers):
            orbmag_mesh = r.read_value('orbmag_mesh')
            for iterm in range(orbmag_nterms):
                for isppol in range(nsppol):
                    for ikpt in range(nkpt):
                        for iband in range(mband):
                            # convert terms to Cart coords; formulae differ depending on term. First
                            # four were in k space, remaining two already in real space
                            if iterm < 4:
                                omtmp = ucvol * np.matmul(gprimd, orbmag_mesh[iterm,0:ndir,isppol,ikpt,iband])
                            else:
                                omtmp = np.matmul(rprimd, orbmag_mesh[iterm,0:ndir,isppol,ikpt,iband])

                            self.orbmag_merge_mesh[iterm,idir,0:ndir,isppol,ikpt,iband] = omtmp

        # This piece of code has been taken from orbmag_sigij_mesh
        wtk = _read_value('kpoint_weights')
        occ = _read_value('occupations')

        self.orbmag_merge_sigij_mesh = np.zeros((orbmag_nterms, nsppol, nkpt, mband, ndir, ndir))

        for isppol in range(nsppol):
            for ikpt in range(nkpt):
                for iband in range(mband):
                    # weight factor for each band and k point
                    trnrm = occ[isppol,ikpt,iband] * wtk[ikpt] / ucvol
                    for iterm in range(orbmag_nterms):
                        # sigij = \sigma_ij the 3x3 shielding tensor term for each sppol, kpt, and band
                        # additional ucvol factor converts to induced dipole moment (was dipole moment density,
                        # that is, magnetization)
                        self.orbmag_merge_sigij_mesh[iterm,isppol,ikpt,iband,0:ndir,0:ndir] = \
                            ucvol * trnrm * self.orbmag_merge_mesh[iterm,0:ndir,0:ndir,isppol,ikpt,iband]

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level verbose"""
        lines = []; app = lines.append
        for orb in self.orb_files:
            app(orb.to_string(verbose=verbose))
        return "\n".join(lines)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Activated at the end of the with statement. It automatically closes all the files."""
        self.close()

    def close(self) -> None:
        """Close all the files."""
        for orb in self.orb_files:
            orb.close()

    @lazy_property
    def natom(self)  -> int:
        """Number of atoms in the unit cell."""
        return len(self.structure)

    @lazy_property
    def structure(self) -> Structure:
        """Structure object."""
        # Perform consistency check
        structure = self.orb_files[0].structure
        if any(orb_file.structure != structure for orb_file in self.orb_files[1:]):
            raise RuntimeError("ORBMAG.nc files have different structures")
        return structure

    @lazy_property
    def has_nucdipmom(self) -> np.ndarray:
        """
        """
        has_nucdipmom = np.zeros(self.natom, dtype=int)
        for orb in self.orb_files:
            # nctkarr_t("nucdipmom", "dp", "ndir, natom")
            nucdipmom = orb.r.read_value('nucdipmom')
            for iat in range(self.natom):
                if np.any(np.abs(nucdipmom[iat]) > 1e-8):
                    has_nucdipmom[iat] += 1

        if np.count_nonzero(has_nucdipmom) != 1:
            raise RuntimeError(f"Only one site should have nucdipmom activated! {has_nucdipmom=}")

        return has_nucdipmom.astype(bool)

    #@lazy_property
    #def has_timrev(self) -> bool:
    #    """True if time-reversal symmetry is used in the BZ sampling."""
    #    has_timrev = self.orb_files[0].ebands.has_timrev
    #    if any(orb_file.ebands.has_timrev != has_timrev for orb_file in self.orb_files[1:]):
    #        raise RuntimeError("ORBMAG.nc files have different values of timrev")

    def report_eigvals(self, report_type, precision=2) -> None:
        """
        FIXME

        Args:
            report_type:
            precision:
        """
        with np.printoptions(precision=precision):

            terms = ['CC   ','VV1  ','VV2  ','NL   ','LR   ','A0An ']
            # Shape: (orbmag_nterms, nsppol, nkpt, mband, ndir, ndir)

            total_sigij = self.orbmag_merge_sigij_mesh.sum(axis=(0, 1, 2, 3))
            omlamb = self.get_omlamb()
            total_sigij = total_sigij + omlamb
            # note use of -1.0E6 here, we convert from dipole moment (from abinit) to shielding in ppm
            eigenvalues = -1.0E6 * np.real(eigvals(total_sigij))
            isotropic = eigenvalues.sum() / 3.0
            span = eigenvalues.max() - eigenvalues.min()
            skew = 3.0 * (eigenvalues.sum() - eigenvalues.max() - eigenvalues.min() - isotropic) / span

            if report_type=='S':
                print('\nShielding tensor eigenvalues, ppm : ', eigenvalues)
                print('Shielding tensor iso, span, skew, ppm : %6.2f %6.2f %6.2f \n' % (isotropic, span, skew))

            if report_type == 'T':
                print('\nShielding tensor eigenvalues, ppm : ',eigenvalues)
                print('Shielding tensor iso, span, skew, ppm : %6.2f %6.2f %6.2f \n'%(isotropic, span, skew))
                print('Term totals')
                term_sigij = self.orbmag_merge_sigij_mesh.sum(axis=(1, 2, 3))
                for iterm in range(np.size(self.orbmag_merge_sigij_mesh, axis=0)):
                    eigenvalues = -1.0E6 * np.real(eigvals(term_sigij[iterm]))
                    print(terms[iterm] + ': ', eigenvalues)
                print('Lamb  : ', -1.0E6 * np.real(eigvals(omlamb)))
                print('\n')

            elif report_type == 'B':
                print('Band totals')
                band_sigij = self.orbmag_merge_sigij_mesh.sum(axis=(0, 1, 2))
                for iband in range(np.size(self.orbmag_merge_sigij_mesh, axis=3)):
                    eigenvalues = -1.0E6 * np.real(eigvals(band_sigij[iband]))
                    print('band ' + str(iband) + ' : ', eigenvalues)
                print('\n')

            elif report_type == 'TB':
                print('Terms in each band')
                tband_sigij = self.orbmag_merge_sigij_mesh.sum(axis=(1, 2))
                for iband in range(np.size(self.orbmag_merge_sigij_mesh, axis=3)):
                    print('band ' + str(iband) + ' : ')
                    for iterm in range(np.size(self.orbmag_merge_sigij_mesh, axis=0)):
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

    def get_omlamb(self) -> np.ndarray:
        """
        Compute and return ...
        """
        omlamb = np.zeros((3, 3))

        for idir, orb in enumerate(self.orb_files):
            typat = orb.r.read_value("atom_species")
            lambsig = orb.r.read_value('lambsig')
            nucdipmom = orb.r.read_value('nucdipmom')
            for iat in range(self.natom):
                itypat = typat[iat]
                omlamb[idir] += lambsig[itypat-1] * nucdipmom[iat]

        # negative sign to convert to dipole moment, like other quantities
        return -omlamb

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
            what: Strings defining the quantity to insert in the box
            spin: Spin index.
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
                data[iband, ix, iy, iz] = self.get_value(what, spin, ikpt, iband)

        return data, ngkpt, shifts

    @lazy_property
    def has_full_bz(self) -> bool:
        """True if the list of k-points covers the full BZ."""
        ngkpt, shifts = self.ngkpt_and_shifts
        return np.prod(ngkpt) * len(shifts) == self.nkpt

    def get_cif_string(self, symprec=None) -> str:
        """
        Return string with structure and anisotropic U tensor in CIF format.

        print a very minimalisic cif file for use in VESTA,
        where the thermal ellipsoid parameters have been hijacked to
        to show the shielding tensor instead.

        Args:
            symprec (float): If not none, finds the symmetry of the structure
                and writes the cif with symmetry information. Passes symprec
                to the SpacegroupAnalyzer
        """
        # Get string with structure in CIF format.
        # Don't use symprec because it changes the order of the sites
        # and we must be consistent with site_labels when writing aniso_U terms!
        from pymatgen.io.cif import CifWriter
        cif = CifWriter(self.structure, symprec=symprec)

        # @JOE: Note different order of U_ij terms wrt omnc.py
        aniso_u = """loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_23
_atom_site_aniso_U_13
_atom_site_aniso_U_12""".splitlines()

        # TODO: Compute ucif in reduced coords
        # NB: In JOE's script, ucif is called Uaniso
        #raise NotImplementedError("Ucif should be computed.")

        omlamb = self.get_omlamb()
        total_sigij = self.orbmag_merge_sigij_mesh.sum(axis=(0, 1, 2, 3))
        total_sigij = total_sigij + omlamb
        Uaniso = make_U_from_sigma(self.gprimd, total_sigij)

        # Add matrix elements. Use 0 based index
        for iat, site in enumerate(self.structure):
            site_label = "%s%d" % (site.specie.symbol, iat)
            if self.has_nucdipmom[iat]:
                m = Uaniso
                aniso_u.append("%s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" %
                        (site_label, m[0, 0], m[1, 1], m[2, 2], m[1, 2], m[0, 2], m[0, 1]))
            else:
                aniso_u.append("%s %s" % (site_label, " . . . . . . "))

        return str(cif) + "\n".join(aniso_u)

    #def vesta_open(self, temp=300): # pragma: no cover
    #    """
    #    Visualize termal displacement ellipsoids at temperature `temp` (Kelvin) with Vesta_ application.
    #    """
    #    filepath = self.write_cif_file(filepath=None, temp=temp)
    #    cprint("Writing structure + Debye-Waller tensor in CIF format for T = %s (K) to file: %s" % (temp, filepath), "green")
    #    cprint("In the Vesta GUI, select: Properties -> Atoms -> Show as displament ellipsoids.", "green")
    #    from abipy.iotools import Visualizer
    #    visu = Visualizer.from_name("vesta")

    #    return visu(filepath)()

    def write_magres(self, gsr_file=None, filepath="tmp.magres") -> None:
        """
        """
        #def output_magres(ncfiles,gsrroot=None,omroots=None,orbmag_merge_sigij_mesh=None,omlamb=None):

        def get_sig_array(total_sigij, omroots):
            natoms=len(omroots[0].dimensions['number_of_atoms'])
            sigij_array=np.zeros((natoms,3,3))
            for iat in range(natoms):
                 if len(np.nonzero(omroots[0].variables['nucdipmom'][iat])[0]) > 0:
                     sigij_array[iat]=total_sigij
            return sigij_array

        def get_indices_and_labels(ncroot):
            natoms = len(ncroot.dimensions['number_of_atoms'])
            ntypat = len(ncroot.dimensions['ntypat'])
            type_count = np.zeros(ntypat, dtype=int)
            indices = np.zeros(natoms, dtype=int)
            labels = np.empty(natoms, dtype=np.dtype('U2'))
            from netCDF4 import chartostring
            for iat in range(natoms):
                itypat = ncroot.variables['atom_species'][iat]
                type_count[itypat-1] = type_count[itypat-1] + 1
                indices[iat] = type_count[itypat-1]
                labels[iat] = chartostring(ncroot.variables['atom_species_names'][itypat-1])

            return indices, labels

        if not gsrroot and not omroots:
            print("output_magres received no netcdf files, nothing to do.")
            return

        if omroots and not np.any(self.orbmag_merge_sigij_mesh):
            print("magres output for shielding requested but orbmag_merge not provided.")
            return

        if omroots and not np.any(omlamb):
            print("magres output for shielding requested but Lamb shielding not provided.")
            return

        from ase.spacegroup import Spacegroup
        import ase.io.magres as magres

        ase_atoms = self.structure.to_ase_atoms()
        ase_atoms.info.update({'spacegroup': Spacegroup(1, setting=1)})
        if gsrroot:
            indices, labels = get_indices_and_labels(gsrroot[0])
        else:
            indices, labels = get_indices_and_labels(omroots[0])

        ase_atoms.new_array('indices', indices)
        ase_atoms.new_array('labels', labels)

        #if 'magres_units' not in ase_atoms.info.keys():
        ase_atoms.info.update({'magres_units':{'ms':'ppm'}})
        #else:
        #    ase_atoms.info['magres_units'].update({'ms':'ppm'})

        omlamb = self.get_omlamb()
        total_sigij = self.orbmag_merge_sigij_mesh.sum(axis=(0, 1, 2, 3))
        total_sigij = -1.0E6 * (total_sigij + omlamb)
        sig_array = get_sig_array(total_sigij, omroots)
        ase_atoms.new_array('ms',sig_array)

        if gsrroot:
            if 'efg' in gsrroot[0].variables:
                if 'magres_units' not in ase_atoms.info.keys():
                    ase_atoms.info.update({'magres_units': {'efg':'au'}})
                else:
                    ase_atoms.info['magres_units'].update({'efg':'au'})
                efg = gsrroot[0].variables['efg'][:]
                ase_atoms.new_array('efg', efg)

        with open(str(filepath), 'wt') as mout:
            magres.write_magres(mout, ase_atoms)

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

        if not self.has_full_bz:
            raise NotImplementedError("k-points must cover the full BZ.")

        for spin in range(self.nsppol):
            data, ngkpt, shifts = self.insert_inbox(what, spin)
            interp_spin[spin] = BzRegularGridInterpolator(self.structure, shifts, data, method=interp_method)

        return interp_spin

    @add_fig_kwargs
    def plot_fatbands(self, ebands_kpath,
                      what_list="isotropic",
                      ylims=None, scale=50, marker_color="gold", marker_edgecolor="gray",
                      marker_alpha=0.5, fontsize=12, interp_method="linear",
                      ax_mat=None, **kwargs) -> Figure:
        """
        Plot fatbands FIXME ...

        Args:
            ebands_kpath: ElectronBands instance with energies along a k-path
                or path to a netcdf file providing it.
            what_list: string or list of strings defining the quantities to plot as fatbands.
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
                    for e, value in zip(enes_n, interp_spin[spin].eval_kpoint(kpoint), strict=True):
                        # note use of -value: the abinit calculation returns the induced magnetic moment,
                        # but the plot should show shielding, which is -moment by convention
                        # this mimics use of -1.0E6 in all the eigval reports below
                        x.append(ik); y.append(e); s.append(scale * (-value))
                        ymin, ymax = min(ymin, e), max(ymax, e)

            # Compute colors based on sign (e.g., red for positive, blue for negative)
            c=["red" if value >= 0 else "blue" for value in s]

            points = Marker(x, y, s,
                            c=c,marker='s',
                            #color=marker_color,edgecolors=marker_edgecolor,
                            alpha=marker_alpha, label=what)

            # Plot electron bands with markers.
            ebands_kpath.plot(ax=ax, points=points, show=False, linewidth=1.0)
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

            shielding_color=mpatches.Patch(color="red",label="Shielding")
            deshielding_color=mpatches.Patch(color="blue",label="Deshielding")
            ax.legend(handles=[shielding_color,deshielding_color])

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
    .. inheritance-diagram:: OrbmagFile
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
    def target_atom(self) -> tuple:
        return self.target_atom_nucdipmom[0]

    @lazy_property
    def nucdipmom_atom(self) -> tuple:
        return self.target_atom_nucdipmom[1]

    @lazy_property
    def target_atom_nucdipmom(self) -> tuple:
        """
        Return the index of the atom with non-zero nuclear magnetic dipole moment and its values.
        """
        # in C: nucdipmom(natom, ndir)
        nucdipmom = self.r.read_value('nucdipmom')
        target_atom, nfound = -1, 0
        for iat in range(len(self.structure)):
            if np.any(nucdipmom[iat] != 0):
                target_atom = iat
                nfound += 1

        if target_atom == -1:
            raise RuntimeError(f"no dipole found in {self.filepath}")

        if nfound != 1:
            raise RuntimeError(f"Found multiple atoms with nuclear magnetic dipole in {self.filepath}")

        return target_atom, nucdipmom[target_atom].copy()

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

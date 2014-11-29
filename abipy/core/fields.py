# coding: utf-8
"""This module contains the class describing densities in real space on uniform 3D meshes."""
from __future__ import print_function, division, unicode_literals

import numpy as np

from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.core.units import bohr_to_angstrom
from abipy.tools import transpose_last3dims
from abipy.iotools import Visualizer, xsf, ETSF_Reader
from abipy.core.mesh3d import Mesh3D
from abipy.core.mixins import Has_Structure
from abipy.tools import transpose_last3dims
from abipy.iotools import Visualizer, xsf, ETSF_Reader


__all__ = [
    "ScalarField",
    "Density",
    #"Potential",
]


class ScalarField(Has_Structure):
    """
    Base class representing a typical scalar field generated electrons (e.g. densities, potentials).
    The field is represented on a homogenous real-space mesh.
    This class provides helper functions to perform common operations such as FFT transforms.
    """
    def __init__(self, nspinor, nsppol, nspden, datar, structure, iorder="c"):
        """
        Args:
            nspinor:
                Number of spinorial components.
            nsppol:
                Number of spins.
            nspden:
                Number of spin density components.
            datar:
                numpy array with the scalar field in real space. shape [nspden, nx, ny, nz]
                Note that here, unlike, abinit we store the spin-components instead of the
                sum over spins. For spin-polarized calculation, for example, datar[0] contains
                the spin-up term whereas datar[1] stores the spin-up contribution.
            structure:
                `Structure` object describing the crystalline structure.
            iorder:
                Order of the array. "c" for C ordering, "f" for Fortran ordering.
        """
        self.nspinor, self.nsppol, self.nspden = nspinor, nsppol, nspden
        self._structure = structure

        iorder = iorder.lower()
        assert iorder in ["f", "c"]

        if iorder == "f": 
            # (z,x,y) --> (x,y,z)
            datar = transpose_last3dims(datar)

        # Init Mesh3D
        mesh_shape = datar.shape[-3:]
        self._mesh = Mesh3D(mesh_shape, structure.lattice_vectors())

        # Make sure we have the correct shape.
        self._datar = np.reshape(datar, (nspden,) + self.mesh.shape)

    def __len__(self):
        return len(self.datar)

    def __str__(self):
        return self.to_string()

    @property
    def structure(self):
        return self._structure

    def to_string(self, prtvol=0):
        """String representation"""
        lines = ["%s: nspinor = %i, nsppol = %i, nspden = %i" %
                 (self.__class__.__name__, self.nspinor, self.nsppol, self.nspden)]
        app = lines.append
        app(self.mesh.to_string(prtvol))
        if prtvol > 0: app(str(self.structure))

        return "\n".join(lines)

    @property
    def datar(self):
        """`ndarrray` with data in real space."""
        return self._datar

    @lazy_property
    def datag(self):
        """`ndarrray` with data in reciprocal space."""
        # FFT R --> G.
        return self.mesh.fft_r2g(self.datar)

    @property
    def mesh(self):
        """`Mesh3D`"""
        return self._mesh

    @property
    def shape(self):
        """Shape of the array."""
        assert self.datar.shape == self.datag.shape
        return self.datar.shape

    @property
    def nx(self):
        """Number of points along x."""
        return self.mesh.nx

    @property
    def ny(self):
        """Number of points along y."""
        return self.mesh.ny

    @property
    def nz(self):
        """Number of points along z."""
        return self.mesh.nz

    @property
    def is_collinear(self):
        """True if collinear i.e. nspinor==1."""
        return self.nspinor == 1

    #@property
    #def datar_xyz(self):
    #    """
    #    Returns a copy of the real space data with shape [:, nx, ny, nz]. 
    #    Mainly used for post-processing.
    #    """
    #    return self.mesh.reshape(self.datar).copy()

    #@property
    #def datag_xyz(self):
    #    """
    #    Returns a copy of the reciprocal space data with shape [:, nx, ny, nz]. 
    #    Mainly used for post-processing.
    #    """
    #    return self.mesh.reshape(self.datag).copy()

    @staticmethod
    def _check_space(space):
        space = space.lower()
        if space not in ("r", "g"):
            raise ValueError("Wrong space %s" % space)
        return space

    def mean(self, space="r"):
        """Returns the average of the array elements."""
        if "r" == self._check_space(space):
            return self.datar.mean(axis=0)
        else:
            return self.datag.mean(axis=0)

    def std(self, space="r"):
        """Returns the standard deviation."""
        if "r" == self._check_space(space):
            return self.datar.std(axis=0)
        else:
            return self.datag.std(axis=0)

    #def braket_waves(self, bra_wave, ket_wave):
    #    """
    #    Compute the matrix element of <bra_wave|datar|ket_wave> in real space
    #    """

    #    if bra_wave.mesh != self.mesh:
    #       bra_ur = bra_wave.fft_ug(self.mesh)
    #    else:
    #       bra_ur = bra_wave.ur

    #    if ket_wave.mesh != self.mesh:
    #       ket_ur = ket_wave.fft_ug(self.mesh)
    #    else:
    #       ket_ur = ket_wave.ur

    #    if self.nspinor == 1:
    #        assert bra_wave.spin == ket_wave.spin
    #       datar_spin = self.datar[bra_wave.spin]
    #       return self.mesh.integrate(bra_ur.conj() * datar_spin * ket_ur)
    #    else:
    #        raise NotImplemented("nspinor != 1 not implmenented")

    #def map_coordinates(self, rcoords, order=3, frac_coords=True)
    #    """
    #    Interpolate the real space data
    #
    #    Args:
    #        coordinates: array_like
    #           The coordinates at which input is evaluated.
    #        order: int, optional
    #            The order of the spline interpolation, default is 3. The order has to be in the range 0-5.
    #    Returns: 
    #       ndarray with the interpolated results.
    #    """
    #    from scipy.ndimage.interpolation import map_coordinates
    #    # Compute the fractional coordinates at which datar is interpolated.
    #    rcoords = np.asarray(rcoords)
    #    if not frac_coords:
    #       rcoords = self.structure.to_frac_coords(rcoords, in_cell=True)
    #    # Wrap in the unit cell.
    #    rcoords %= 1
    #    coordinates = [rcoords[0], rcoords[1], rcoords[2]]

    #    Interpolate the real part.
    #    interp_data = []
    #    for indata in self.datar_xyz:
    #        assert not np.iscomple(indata)
    #        interp_data.append(map_coordinates(indata.real, coordinates, order=order))

    #    return np.array(interp_data)

    #def fourier_interp(self, new_mesh):
        #intp_datar = self.mesh.fourier_interp(self.datar, new_mesh, inspace="r")
        #return self.__class__(self.nspinor, self.nsppol, self.nspden, self.structure, intp_datar)

    def export(self, filename, visu=None):
        """
        Export the real space data on file filename. 

        Args:
            filename:
                String specifying the file path and the file format.
                The format is defined by the file extension. filename="prefix.xsf", for example, 
                will produce a file in XSF format. An *empty* prefix, e.g. ".xsf" makes the code use a temporary file.
            visu:
               `Visualizer` subclass. By default, this method returns the first available
                visualizer that supports the given file format. If visu is not None, an
                instance of visu is returned. See :class:`Visualizer` for the list of 
                applications and formats supported.

        Returns:
            Instance of :class:`Visualizer`
        """
        if "." not in filename:
            raise ValueError(" Cannot detect file extension in filename: %s " % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]: # filename == ".ext" ==> Create temporary file.
            import tempfile
            filename = tempfile.mkstemp(suffix="." + ext, text=True)[1]

        with open(filename, mode="w") as fh:
            if ext == "xsf":
                # xcrysden
                xsf.xsf_write_structure(fh, self.structure)
                xsf.xsf_write_data(fh, self.structure, self.datar, add_replicas=True)
            else:
                raise NotImplementedError("extension %s is not supported." % ext)

        if visu is None:
            return Visualizer.from_file(filename)
        else:
            return visu(filename)

    def visualize(self, visu_name):
        """
        Visualize data with visualizer.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        visu = Visualizer.from_name(visu_name)

        # Try to export data to one of the formats supported by the visualizer
        # Use a temporary file (note "." + ext)
        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                return self.export(ext, visu=visu)
            except visu.Error:
                pass
        else:
            raise visu.Error("Don't know how to export data for visualizer %s" % visu_name)

    #def get_line(self, line, space="r"):
    #    x, y, z = self.mesh.line_inds(line)
    #    space = self._check_space(space)
    #    if space == "r":
    #       line = self.datar_xyz[:, x, y, z]
    #    elif space == "g":
    #       line = self.datag_xyz[:, x, y, z]
    #    # Return a 2D array.
    #    new_shape = lines.shape[0] + tuple(s for s in shape[-3:] is s)
    #    return np.reshape(line, new_shape)

    #def get_plane(self, plane, h, space="r"):
    #    x, y, z = self.mesh.plane_inds(plane, h=h)
    #    space = self._check_space(space)
    #    if space == "r":
    #       plane = self.datar_xyz[:, x, y, z]
    #    elif space == "g":
    #       plane = self.datag_xyz[:, x, y, z]
    #    # Return a 3D array.
    #    new_shape = lines.shape[0] + tuple(s for s in shape[-3:] is s)
    #    return np.reshape(plane, new_shape)


class Density(ScalarField):
    """
    Electronic density
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netCDF file."""
        with DensityReader(filepath) as r:
            return r.read_density(cls=cls)

    #def __init__(self, nspinor, nsppol, nspden, rhor, structure, iorder="c"):
    #    """
    #    Args:
    #        nspinor:
    #            Number of spinorial components.
    #        nsppol:
    #            Number of spins.
    #        nspden:
    #            Number of spin density components.
    #        datar:
    #            numpy array with the field in real space.
    #        structure:
    #            pymatgen structure
    #        iorder:
    #            Order of the array. "c" for C ordering, "f" for Fortran ordering.
    #    """
    #    super(Density, self).__init__(nspinor, nsppol, nspden, rhor, structure, iorder=iorder)

    def get_nelect(self, spin=None):
        """
        Returns the number of electrons with given spin.

        If spin is None, the total number of electrons is computed.
        """
        if self.is_collinear:
            nelect = self.mesh.integrate(self.datar)
            return np.sum(nelect) if spin is None else nelect[spin]
        else:
            return self.mesh.integrate(self.datar[0])

    @lazy_property
    def total_rhor(self):
        """numpy array with the total density in real space on the FFT mesh"""
        if self.is_collinear:
            if self.nsppol == 1:
                if self.nspden == 2: raise NotImplementedError()
                return self.datar[0]
            elif self.nsppol == 2:
                #tot_rhor = np.sum(self.datar, axis=0)
                return self.datar[0] + self.datar[1]
            else:
                raise ValueError("You should not be here")

        # Non collinear case.
        raise NotImplementedError

    @lazy_property
    def total_rhog(self):
        """numpy array with the total density in G-space."""
        # FFT R --> G.
        return self.mesh.fft_r2g(self.total_rhor)

    @lazy_property
    def magnetization_field(self):
        """
        :return: numpy array with the magnetization field in real space on the FFT mesh:

            #. 0 if spin-unpolarized calculation
            #. spin_up - spin_down if collinear spin-polarized
            #. numpy array with (mx, my, mz) components if non-collinear magnetism
        """
        if self.is_collinear:
            if self.nsppol == 1 and self.nspden == 1:
                # zero magnetization by definition.
                return self.mesh.zeros()
            else:
                # spin_up - spin_down.
                return self.datar[0] - self.datar[1]
        else:
            # mx, my, mz
            return self.datar[1:]

    @lazy_property
    def magnetization(self):
        """
        Magnetization field integrated over the unit cell.
        Scalar if collinear, vector with mx, my, mz components if non-collinear.
        """
        return self.mesh.integrate(self.magnetization_field)

    @lazy_property
    def nelect_updown(self):
        if not self.is_collinear: return None, None

        if self.nsppol == 1:
            if self.nspden == 2: raise NotImplementedError()
            nup = ndown = self.mesh.integrate(self.datar[0]/2)
        else:
            nup = self.mesh.integrate(self.datar[0])
            ndown = self.mesh.integrate(self.datar[1])

        return nup, ndown

    @lazy_property
    def zeta(self):
        """Magnetization(r) / total_density(r)"""
        fact = np.where(self.tot_rhor > 1e-16, 1/self.tot_rhor, 0.0)
        return self.magnetization * fact

    def vhartree(self):
        """
        Solve the Poisson's equation in reciprocal space.

        returns:
            (vhr, vhg) Hartree potential in real, reciprocal space.
        """
        raise NotImplementedError("")
        # Compute |G| for each G in the mesh and treat G=0.
        gvecs = self.mesh.gvecs
        gwork = self.mesh.zeros().ravel()
        gnorm = self.structure.gnorm(gvec)

        for idx, gg in enumerate(gvecs):
            #gnorm = self.structure.gnorm(gg)
            gnorm = 1.0  # self.structure.gnorm(gg)

            #gg = np.atleast_2d(gg)
            #mv = np.dot(self.structure.gmet, gg.T)
            #norm2 = 2*np.pi * np.dot(gg, mv)
            #gnorm = np.sqrt(norm2)

            #print gg, gnorm
            if idx != 0:
                gwork[idx] = 4*np.pi/gnorm
            else:
                gwork[idx] = 0.0

        new_shape = self.mesh.ndivs
        gwork = np.reshape(gwork, new_shape)
        #gwork = self.mesh.reshape(gwork)

        # FFT to obtain vh in real space.
        vhg = self.total_rhog * gwork
        vhr = self.mesh.fft_g2r(vhg, fg_ishifted=False)

        return vhr, vhg

    #@lazy_property
    #def kinden(self):
        #"""Compute the kinetic energy density in real- and reciprocal-space."""
        #return kindr, kindgg

    #def vxc(self, xc_type=None):
        #"""Compute the exchange-correlation potential in real- and reciprocal-space."""
        #return vxcr, vxcg


class DensityReader(ETSF_Reader):
    """This object reads density data from a netcdf file."""

    def read_den_dims(self):
        """Returns an `AttrDict` dictionary with the dimensions characterizing the density"""
        return AttrDict(
            cplex_den=self.read_dimvalue("real_or_complex_density"),
            nspinor=self.read_dimvalue("number_of_spinor_components"),
            nsppol=self.read_dimvalue("number_of_spins"),
            #nspden=self.read_dimvalue("number_of_spin_density_components"),
            nspden=self.read_dimvalue("number_of_components"),
            nfft1=self.read_dimvalue("number_of_grid_points_vector1"),
            nfft2=self.read_dimvalue("number_of_grid_points_vector2"),
            nfft3=self.read_dimvalue("number_of_grid_points_vector3"),
        )

    def read_density(self, cls=Density):
        """Factory function that builds and returns a `Density` object."""
        structure = self.read_structure()
        dims = self.read_den_dims()

        # Abinit conventions:
        # rhor(nfft,nspden) = electron density in r space
        # (if spin polarized, array contains total density in first half and spin-up density in second half)
        # (for non-collinear magnetism, first element: total density, 3 next ones: mx,my,mz in units of hbar/2)
        rhor = self.read_value("density")

        if dims.nspden == 1:
            pass
        elif dims.nspden == 2:
            # Store rho_up, rho_down instead of rho_total, rho_up
            total = rhor[0].copy()
            rhor[0] = rhor[1]
            rhor[1] = total - rhor[1]
        elif dims.nspden == 4:
            raise NotImplementedError("nspden == 4 not coded")
        else:
            raise RuntimeError("You should not be here")

        # use iorder="f" to transpose the last 3 dimensions since ETSF
        # stores data in Fortran order while abipy uses C-ordering.
        if dims.cplex_den == 1:
            # Get rid of fake last dimensions (cplex).
            rhor = np.reshape(rhor, (dims.nspden, dims.nfft1, dims.nfft2, dims.nfft3))

            # Structure uses Angstrom. Abinit uses bohr.
            rhor /= (bohr_to_angstrom ** 3)
            return cls(dims.nspinor, dims.nsppol, dims.nspden, rhor, structure, iorder="f")

        else:
            raise NotImplementedError("cplex_den %s not coded" % dims.cplex_den)

"""This module contains the class describing densities in real space on uniform 3D meshes."""
from __future__ import print_function, division

import numpy as np

from abipy.tools import transpose_last3dims, AttrDict
from abipy.iotools import Visualizer, xsf, ETSF_Reader
from abipy.core.constants import bohr_to_angstrom
from abipy.core.mesh3d import Mesh3D
from abipy.core.structure import Structure


__all__ = [
    "ScalarField",
    "Density",
    #"Potential",
]


class ScalarField(object):
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
                numpy array with the scalar field in real space. shape [..., nx, ny, nz]
            structure:
                `Structure` object describing the crystalline structure.
            iorder:
                Order of the array. "c" for C ordering, "f" for Fortran ordering.
        """
        self.nspinor, self.nsppol, self.nspden = nspinor, nsppol, nspden
        self.structure = structure

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

        # FFT R --> G.
        self._datag = self.mesh.fft_r2g(self.datar)

    def __len__(self):
        return len(self.datar)

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        """String representation"""
        s  = "%s: nspinor = %i, nsppol = %i, nspden = %i" % (
            self.__class__.__name__, self.nspinor, self.nsppol, self.nspden)
        s += "  " + self.mesh.tostring(prtvol)

        if prtvol > 0:
            s += "  " + str(self.structure)

        return s

    @property
    def datar(self):
        """`ndarrray` with data in real space."""
        return self._datar

    @property
    def datag(self):
        """`ndarrray` with data in reciprocal space."""
        return self._datag

    @property
    def mesh(self):
        """`Mesh3D`"""
        return self._mesh

    @property
    def shape(self):
        """Shape of the array."""
        shape_r, shape_g = self.datar.shape, self.datag.shape
        assert shape_r == shape_g
        return shape_r

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

    #    assert self.nspinor == 1
    #    assert bra_wave.spin == ket_wave.spin

    #    datar_spin = self.datar[bra_wave.spin]
    #    return self.mesh.integrate(bra_ur.conj() * datar_spin * ket_ur)

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
        # Insert self.datag in the FFT box of new mesh.
        #intp_datag = np.empty_like(self.datag)
        #intp_datag = new_mesh.transferg_from(self.mesh, self.datag)
        # FFT transform G --> R.
        #intp_datar = new_mesh.fft_g2r(intp_datag)
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
            filename = tempfile.mkstemp(suffix="."+ext, text=True)[1]

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
            raise visu.Error("Don't know how to export data for visualizer %s" % visualizer)

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
    Electron density
    """
    # TODO 
    #  "exchange_functional",      # "exchange_functional"
    #  "valence_charges",          # "valence_charges"
    def __init__(self, nspinor, nsppol, nspden, rhor, structure, iorder="c"):
        """
        Args:
            nspinor:
                Number of spinorial components.
            nsppol:
                Number of spins.
            nspden:
                Number of spin density components.
            datar:
                numpy array with the field in real space.
            structure:
                pymatgen structure
            iorder:
                Order of the array. "c" for C ordering, "f" for Fortran ordering.
        """
        super(Density, self).__init__(nspinor, nsppol, nspden, rhor, structure, iorder=iorder)

    @classmethod
    def from_file(cls, filepath):
        """
        Read density from an external netCDF file.

        Args:
            filepath:
                string or file object.
        """
        with DensityReader(filepath) as r:
            structure = r.read_structure()
            dims = r.read_dendims()
            rhor = r.read_rhor()

        # use iorder="f" to transpose the last 3 dimensions since ETSF
        # stores data in Fortran order while abipy uses C-ordering.
        if dims.cplex_den == 1:

            # Get rid of fake last dimensions (cplex).
            rhor = np.reshape(rhor, (dims.nspden, dims.nfft1, dims.nfft2, dims.nfft3))

            # Structure uses Angstrom. Abinit uses bohr.
            rhor = rhor / bohr_to_angstrom ** 3

            return Density(dims.nspinor, dims.nsppol, dims.nspden, rhor, structure, iorder="f")

        else:
            raise NotImplementedError("cplex_den %s not coded" % dims.cplex_den)

    def get_nelect(self, spin=None):
        """
        Returns the number of electrons with given spin.

        If spin is None, the total number of electrons is computed.
        """
        if not self.is_collinear:
            raise NotImplementedError("Non collinear not implemented")

        nelect = self.mesh.integrate(self.datar)

        return np.sum(nelect) if spin is None else nelect[spin]

    #def get_magnetization(self)

    def get_rhor_tot(self):
        """Returns the total density in real space."""
        if self.nsppol == 2:
            raise NotImplementedError("check whether ETSF-IO uses up-down storage mode")

        if self.is_collinear:
            rhor_tot = np.sum(self.datar, axis=0)
            if self.nspden == 2 and self.nsppol == 1:
                raise NotImplementedError

        else:
            raise NotImplementedError

        return rhor_tot

    def get_rhog_tot(self):
        """Returns the total density in G-space."""
        if self.nsppol == 2:
            raise NotImplementedError("check whether ETSF-IO uses up-down storage mode")

        if self.is_collinear:
            rhog_tot = np.sum(self.datag, axis=0)
            if self.nspden == 2 and self.nsppol == 1: raise NotImplementedError

        else:
            raise NotImplementedError

        return rhog_tot

    def get_vh(self):
        """
        Solve the Poisson's equation in reciprocal space.

        returns:
            (vhr, vhg) Hartree potential in real, reciprocal space.
        """
        raise NotImplementedError("")
        # Compute total density in G-space.
        rhog_tot = self.get_rhog_tot()
        #print rhog_tot

        # Compute |G| for each G in the mesh and treat G=0.
        gvec  = self.mesh.get_gvecs()

        gwork = self.mesh.zeros().ravel()

        gnorm = self.structure.gnorm(gvec)

        for idx, gg in enumerate(gvec):
            #gnorm = self.structure.gnorm(gg)
            gnorm = 1.0 #self.structure.gnorm(gg)

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
        vhg = rhog_tot * gwork

        vhr = self.mesh.fft_g2r(vhg, fg_ishifted=False)
        return vhr, vhg

    #def get_vxc(self, xc_type=None):
        #"""Compute the exchange-correlation potential in real- and reciprocal-space."""
        #return vxcr, vxcg

    #def get_kinden(self):
        #"""Compute the kinetic energy density in real- and reciprocal-space."""
        #return kindr, kindgg


class DensityReader(ETSF_Reader):
    """This object reads density data from a netcdf file."""
    def read_dendims(self):
        """Returns an `AttrDict` dictionary with the basic dimensions."""
        return AttrDict(
            cplex_den=self.read_dimvalue("real_or_complex_density"),
            nspinor=self.read_dimvalue("number_of_spinor_components"),
            nsppol=self.read_dimvalue("number_of_spins"),
            nspden=self.read_dimvalue("number_of_components"),
            nfft1=self.read_dimvalue("number_of_grid_points_vector1"),
            nfft2=self.read_dimvalue("number_of_grid_points_vector2"),
            nfft3=self.read_dimvalue("number_of_grid_points_vector3"),
        )

    def read_rhor(self):
        """Return the density in real space."""
        return self.read_value("density")



# Global variables
#: Flags denoting the potential type. Note bit order.
#VNONE = 1
#VLOC  = 2
#VH    = 4
#VX    = 8
#VC    = 16
#
##: Mapping flag --> potential name.
#_vnames =  {
#VNONE: "None potential",
#VLOC : "Local potential",
#VH   : "Hartree potential",
#VX   : "Exchange potential",
#VC   : "Correlation potential",
#}
#
##: List of allowed potential types.
#VTYPES = _vnames.keys()
#
#
#class PotentialInfo(object):
#
#    def __init__(self, vtype, **kwargs):
#        self.vtype = vtype
#        self._dict = kwargs
#
#    def __str__(self):
#        s = self.vname() + "\n"
#        for k, v in self._dict.items:
#            s += str(k) + " = " + str(v) + "\n"
#        return s
#
#    def __add__(self, other):
#        new_vtype = self.vtype + other.vtype
#
#        new_dict = self._dict.copy()
#
#        for k in other._dict:
#            if k not in self._dict:
#                new_dict[k] = other._dict[k]
#            else:
#                # Append values if key is already present.
#                lst = list()
#                lst.append(new_dict[k])
#                lst.append(other._dict[k])
#                new_dict[k] = lst
#
#        return PotentialInfo(vtype, kwargs=new_dict)
#
#    @property
#    def vname(self):
#        """Return the name of the potential."""
#        s = bin(self.vtype)[2:][::-1]
#        flags = [ 2**idx for idx, c in enumerate(s) if c == "1" ]
#
#        plus = ""
#        if len(flags) > 1: plus = "+"
#        return plus.join( [_vnames[flag] for flag in flags] )
#
#
#class Potential(ScalarField):
#    """This module contains the class describing local potentials in real space on uniform 3D meshes."""
#    #
#    #: Attributes read from the netCDF file.
#    _slots = [
#      "cplex_pot",                # "real_or_complex_potential"
#      "nsppol",                   # "number_of_spins"
#      "nspinor",                  # "number_of_spinor_components"
#      "nspden",                   # "number_of_components"
#      "nfft1",                    # "number_of_grid_points_vector1"
#      "nfft2",                    # "number_of_grid_points_vector2"
#      "nfft3",                    # "number_of_grid_points_vector3"
#      "nelect",                   # "number_of_electrons"
#      "valence_charges",          # "valence_charges"
#    ]
#                                                                    
#    _potopts = [
#      "exchange_potential",              #"vx",
#      "correlation_potential",           #"vc",
#      "exchange_correlation_potential",  #"vxc",
#      "exchange_functional",             #"exchange_functional"
#      "correlation_functional",          #"correlation_functional"
#    ]
#
#    @classmethod
#    def from_file(cls, filepath):
#        """
#        Read density from a the netCDF file.
#
#        Args:
#            filepath:
#                string with the file path
#        returns:
#            :class:`Potential`
#        """
#        raise NotImplementedError("Potential must be rewritten from scrath")
#        # Read all the keywords present in fname so that we know the potential type.
#        #dim_names, var_names = ncread_keys(path)
#
#        nfound = 0
#        for pot in cls._potopts:
#            if pot in var_names:
#                nfound += 1
#                v = Record()
#                pot_vars = ["data"]
#                map["data"] = pot
#                pot_vars += _potopts[pot]
#
#                #missing = ncread_varsdims(v, path, pot_vars, map_names=map)
#                #
#                # Handle the error
#                if missing:
#                    for miss in missing:
#                        print("potential variables %s are missing!" % str(miss[0]))
#                    raise ValueError("Fatal Eror")
#
#        if nfound != 1:  # No potential or more than one.
#            raise ValueError("Number of potential found in file is %s" % nfound )
#
#        # Build potential descriptor.
#        #vinfo = PotentialInfo(vtype, **kwargs)
#
#        structure = Structure.from_file(path)
#
#        rec = Record()
#        #
#        # Handle the error
#        if missing:
#            for miss in missing:
#                print("internal name= " + miss[0] + ", ETSF name= " + miss[1] + " is missing!")
#
#        # use iorder="f" to transpose the last 3 dimensions since ETSF
#        # stores data in Fortran order while abipy uses C-ordering.
#
#        if rec.cplex_pot == 1:
#            # Get rid of fake last dimensions (cplex).
#            vr = np.reshape(v.data, (rec.nspden, rec.nfft3, rec.nfft2, rec.nfft1))
#            return Potential(rec.nspinor, rec.nsppol, rec.nspden, vtype, vr, structure, iorder="f")
#
#        else:
#            raise NotImplementedError("cplex_den = " + str(rec.cplex_den) + "not coded")
#
#
#    def __init__(self, nspinor, nsppol, nspden, vtype, vr, structure, iorder="c"):
#        """
#        Args:
#            nspinor:
#                Number of spinorial components.
#            nsppol:
#                Number of spins.
#            nspden:
#                Number of spin density components.
#            vtype:
#                Flag defining the potential type.
#            vr:
#                numpy array with the potential in real space.
#            structure:
#                pymatgen structure
#            iorder:
#                Order of the array. "c" for C ordering, "f" for Fortran ordering.
#        """
#        super(Potential, self).__init__(nspinor, nsppol, nspden, vr, structure, iorder=iorder)
#
#        if vtype not in VTYPES:
#            raise ValueError("Unknow vtype: " + str(vtype))
#        self.vtype = vtype
#
#    @property
#    def vname(self):
#        """The name of the potential."""
#        s = bin(self.vtype)[2:][::-1]
#        flags = [2**idx for idx, c in enumerate(s) if c == "1"]
#
#        plus = ""
#        if len(flags) > 1: plus = "+"
#        return plus.join([_vnames[flag] for flag in flags])
#
#    def tostring(self, prtvol=0):
#        s = self.vname + "\n"
#        s += super(Potential, self).tostring(self, prtvol)
#        return s
#
#    #def make_vector_field(self):
#    #  """Return vector field."""
#    #  if self.iscollinear: return None


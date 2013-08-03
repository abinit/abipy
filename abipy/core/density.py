"""This module contains the class describing densities in real space on uniform 3D meshes."""
from __future__ import division, print_function

import numpy as np

from .dftscalarfield import DFTScalarField
from .constants import Bohr_Ang
from abipy.iotools import abipy2etsfio, ETSF_Reader

__all__ = [
    "Density",
]


class Density(DFTScalarField):
    """
    Electron density
    """
    #: Attributes read from the netCDF file.
    #_slots = [
    #  "nelect",                   # "number_of_electrons"
    #  "exchange_functional",      # "exchange_functional"
    #  "valence_charges",          # "valence_charges"
    #]

    @classmethod
    def from_file(cls, filename):
        """
        Read density from an external netCDF file.

        Args:
            filename:
                string or file object.
        """
        with ETSF_Reader(filename) as r:

            structure = r.read_structure()
            cplex_den = r.read_dimvalue("real_or_complex_density")
            nspinor = r.read_dimvalue("number_of_spinor_components")
            nsppol = r.read_dimvalue("number_of_spins")
            nspden = r.read_dimvalue("number_of_components")

            nfft1 = r.read_dimvalue("number_of_grid_points_vector1")
            nfft2 = r.read_dimvalue("number_of_grid_points_vector2")
            nfft3 = r.read_dimvalue("number_of_grid_points_vector3")

            rhor = r.read_value("density")

        # use iorder="f" to transpose the last 3 dimensions since ETSF
        # stores data in Fortran order while abipy uses C-ordering.
        if cplex_den == 1:

            # Get rid of fake last dimensions (cplex).
            rhor = np.reshape(rhor, (nspden, nfft1, nfft2, nfft3))

            # Fortran to C, avoid the view.
            #cview = np.transpose(rec.rhor, axes = [0,3,2,1])
            #rec.rhor = np.ascontiguousarray(cview)
            rhor = rhor / Bohr_Ang**3
            return Density(nspinor, nsppol, nspden, rhor, structure, iorder="f")

        else:
            raise NotImplementedError("cplex_den %s not coded" % cplex_den)

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

    def get_nelect(self, spin=None):
        """
        Returns the number of electrons with given spin.

        If spin is None, the total number of electrons is computed.
        """
        if not self.iscollinear:
            raise NotImplementedError("Non collinear not implemented")

        nelect = self.mesh.integrate(self.datar)

        if spin is None:
            return np.sum(nelect)
        else:
            return nelect[spin]

    #def get_magnetization(self)

    def get_rhor_tot(self):
        """Returns the total density in real space."""
        if self.nsppol == 2:
            raise NotImplementedError("check whether ETSF-IO uses up-down storage mode")

        if self.iscollinear:
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

        if self.iscollinear:
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
        gvec  = self.mesh.get_gvec()

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

    #def get_vxc(self, spin=None, xc_type=None):
        #return vxcr, vxcg

    #def get_kinden(self, spin=None):
        #return kindr, kindgg


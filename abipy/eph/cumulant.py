# coding: utf-8
"""
This module contains objects for postprocessing e-ph calculations
using the results stored in the out_EPH_CUMULANT.nc file.

"""
import tempfile
import pickle
import os
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu
from scipy.fft import fft, fftfreq,fftshift,ifft


from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from tabulate import tabulate
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import Kpoint, KpointList, Kpath, IrredZone, has_timrev_from_kptopt, find_points_along_path
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles)
from abipy.tools import duck
from abipy.tools.numtools import gaussian
from abipy.electrons.ebands import ElectronBands, ElectronDos, RobotWithEbands, ElectronBandsPlotter, ElectronDosPlotter
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader
from abipy.eph.sigeph import SigEPhFile,EphSelfEnergy, SigmaPhReader, QpTempState

__all__ = [
    "CumulantQpTempState",
    "CumulantEPhFile",
    "CumulantPhReader",
    "CumulantSelfEnergy",
]






class CumulantQpTempState(QpTempState):
    """
    Quasi-particle result for given (spin, kpoint, band).

    .. Attributes:

        spin: Spin index (C convention, i.e >= 0)
        kpoint: |Kpoint| object.
        band: band index. Global index in the band structure. (C convention, i.e >= 0).
        tmesh: Temperature mesh in Kelvin.
        e0: Initial KS energy.
        qpe: Quasiparticle energy (complex) computed with the linearized perturbative approach (Z factor).
        ze0: Renormalization factor Z computed at e = e0.
        fan0: Fan term (complex) evaluated at e_KS
        dw: Debye-Waller (static, real)
        qpe_oms: Quasiparticle energy (real) in the on-the-mass-shell approximation:
            qpe_oms = e0 + Sigma(e0)

    .. note::

        Energies are in eV.
    """


    def set_energies(self, wmesh, vals_wr, vals_e0ks, ntemp, nwr):
        # Determination of the QP energies
        for it in range(ntemp):
            # Checking where Re Sigma(w) crosses with w => e^QP = e^KS ( set to 0.0 ) + Sigma( e^QP )
            idx = np.argwhere(np.diff(np.sign(vals_wr[it,nwr//4:3*nwr//4].real - wmesh[nwr//4:3*nwr//4] + self.e0))).flatten()
            self.qpe[it] = self.e0 + vals_wr[it,idx]
            # On-the-mass-shell QP energy = e^QP = e^KS ( set to 0.0 ) + Sigma(e^KS)
            self.qpe_oms[it] = self.e0 + vals_e0ks[it]



class CumulantEPhFile(SigEPhFile):
    """
    This file contains the Green's Function using the cumulant expansion from a self-energy eph calculation, the |ElectronBands| on the k-mesh.
    Provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        SigmaEphFile = SigEPhFile("out_SIGEPH.nc")
        with CumulantEPhFile("out_EPH_CUMULANT.nc", SigmaEphFile) as ncfile:
            print(ncfile)
            ncfile.ebands.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: CumulantEPhFile
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf file."""
        return cls(filepath)


    def __init__(self, filepath):
        self._filepath = filepath
        self.reader = r = CumulantPhReader(filepath)
        self.nkcalc = r.nkcalc
        self.ntemp = r.ntemp
        self.nqbz = r.nqbz
        self.nqibz = r.nqibz
        self.ngqpt = r.ngqpt


        self.symsigma = -1


    def __str__(self):
        """String representation."""
        return self.to_string()


    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        ##app("")
        ##app(self.structure.to_string(verbose=verbose, title="Structure"))
        #app("")
        #app(self.ebands.to_string(with_structure=False, verbose=verbose, title="KS Electron Bands"))
        app("")
        # SigmaEPh section.
        app(marquee("CumulantEPh calculation", mark="="))
        app("Number of k-points in Sigma_{nk}: %d" % (self.nkcalc))
        # These variables have added recently
        sigma_ngkpt = self.reader.read_value("sigma_ngkpt", default=None)
        sigma_erange = self.reader.read_value("sigma_erange", default=None)
        #dvdb_add_lr = self.reader.read_value("dvdb_add_lr", default=None)
        #app("sigma_ngkpt: %s, sigma_erange: %s" % (sigma_ngkpt, sigma_erange))
        app("Max bstart: %d, min bstop: %d" % (self.reader.max_bstart, self.reader.min_bstop))
        app("Initial ab-initio q-mesh:\n\tngqpt: %s, with nqibz: %s" % (str(self.ngqpt), self.nqibz))
        #eph_ngqpt_fine = self.reader.read_value("eph_ngqpt_fine")
        #if np.all(eph_ngqpt_fine == 0): eph_ngqpt_fine = self.ngqpt
        #app("q-mesh for self-energy integration (eph_ngqpt_fine): %s" % (str(eph_ngqpt_fine)))
        #app("k-mesh for electrons:")
        #app("\t" + self.ebands.kpoints.ksampling.to_string(verbose=verbose))
        #app("Number of bands included in e-ph self-energy sum: %d" % (self.nbsum))
        #app("zcut: %.5f (Ha), %.3f (eV)" % (self.zcut, self.zcut * abu.Ha_eV))
        app("Number of temperatures: %d, from %.1f to %.1f (K)" % (self.ntemp, self.tmesh[0], self.tmesh[-1]))
        #app("symsigma: %s" % (self.symsigma))
        #app("Has Eliashberg function: %s" % (self.has_eliashberg_function))
        #app("Has Spectral function: %s" % (self.has_spectral_function))

        # Build Table with direct gaps. Only the results for the first and the last T are shown if not verbose.
        #if verbose:
        #    it_list = list(range(self.ntemp))
        #else:
        #    it_list = [0, -1] if self.ntemp != 1 else [0]
        #app("\nPrinting QP results for %d temperatures. Use --verbose to print all results." % len(it_list))

        # QP corrections
        #for it in it_list:
        #    app("\nKS, QP (Z factor) and on-the-mass-shell (OTMS) direct gaps in eV for T = %.1f K:" % self.tmesh[it])
        #    data = []
        #    for ikc, kpoint in enumerate(self.sigma_kpoints):
        #        for spin in range(self.nsppol):
        #            ks_gap = self.ks_dirgaps[spin, ikc]
        #            qp_gap = self.qp_dirgaps_t[spin, ikc, it]
        #            oms_gap = self.qp_dirgaps_otms_t[spin, ikc, it]
        #            data.append([spin, repr(kpoint), ks_gap, qp_gap, qp_gap - ks_gap, oms_gap, oms_gap - ks_gap])
        #    app(str(tabulate(data,
        #        headers=["Spin", "k-point", "KS_gap", "QPZ0_gap", "QPZ0 - KS", "OTMS_gap", "OTMS - KS"],
        #        floatfmt=".3f")))
        #    app("")
        #else:
        # Print info on Lifetimes?

        #if verbose > 1:
        #    app("K-points and bands included in self-energy corrections:")
        #    for spin in range(self.nsppol):
        #        for ikc, kpoint in enumerate(self.sigma_kpoints):
        #            post = "ikc: %d" % (ikc if self.nsppol == 1 else "ikc: %d, spin: %d" % (ikc, spin))
        #            app("\t%s: bstart: %d, bstop: %d, %s" % (
        #                repr(kpoint), self.bstart_sk[spin, ikc], self.bstop_sk[spin, ikc], post))

        return "\n".join(lines)

    def get_cumulant_skb(self, spin, kpoint, band):
        """"Return e-ph self-energy for the given (spin, kpoint, band)."""
        return self.reader.read_cumulant_skb(spin, kpoint, band)




class CumulantPhReader(SigmaPhReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: CumulantPhReader
    """

    def __init__(self, path ):
        super().__init__(path)

        # Check if the cumulant function exists
        # It is an optional output, mostly for debuging
        # If it exists, then G(t) and time_mesh also exists
        variables = self.read_varnames()
        if "ct_vals" in variables:
            self.ct_vals_exists = True
        else:
            self.ct_vals_exists = False


    def read_cumulant_skb(self, spin, kpoint, band):
        """
        Returns the e-ph self energy for the given (spin, k-point, band).

        Args:
            spin: Spin index
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band: band index.

        Return: :class:`EphSelfEnergy` object.
        """
        if self.nwr == 0:
            raise ValueError("%s does not contain spectral function data." % self.path)

        spin, ikc, ib, kpoint = self.get_sigma_skb_kpoint(spin, kpoint, band)

        # Abinit fortran (Ha units)
        # real(dp) :: wrmesh_b(nwr, max_nbcalc, nkcalc, nsppol)
        # Frequency mesh along the real axis (Ha units) used for the different bands
        wmesh = self.read_variable("wrmesh_b")[spin, ikc, ib, :] * abu.Ha_eV

        # Set time mesh, cumulant function and green's function in time domain as None ( default )
        time_mesh =ct_vals = gt_vals = None
        if self.ct_vals_exists:
            # For debugging proposes read time mesh, cumulant function and green's function in time domain

            # real(dp) :: time_mesh(nwr, ntemp, max_nbcalc, nkcalc, nsppol)
            # Time mesh (ps units)
            time_mesh = self.read_variable("time_mesh")[spin,ikc,ib,:,:]
            
            # complex(dpc) :: ct_vals(complex, nwr, ntemp, max_nbcalc, nkcalc, nsppol)
            # Cumulant function with respect to time
            ct_vals = self.read_variable("ct_vals")[spin,ikc,ib,:,:,:]
            ct_vals = ct_vals[:, :, 0] + 1j* ct_vals[:, :, 1]

            # complex(dpc) :: gt_vals(complex, nwr, ntemp, max_nbcalc, nkcalc, nsppol)
            # Green's function with respect to time
            gt_vals = self.read_variable("gt_vals")[spin,ikc,ib,:,:,:]
            gt_vals = gt_vals[:,:, 0] + 1j*gt_vals[:, :, 1]

        # complex(dpc) :: gw_vals(nwr, ntemp, max_nbcalc, nkcalc, nsppol)
        # Green's function in frequency domain
        gw_vals = self.read_variable("gw_vals")[spin,ikc,ib,:,:,:] / abu.Ha_eV
        gw_vals = gw_vals[:, :, 0] + 1j * gw_vals[:, :, 1]

        # Spectral function obtained from the Gren's function in frequency domain
        spfunccumul_wr = -1.0*gw_vals[:, :].imag/np.pi

        # Read QP data. Note band instead of ib index.
        qp = self.read_qp(spin, ikc, band)

        return CumulantSelfEnergy(wmesh, qp, gw_vals, spfunccumul_wr, time_mesh = time_mesh, ct_vals= ct_vals, gt_vals= gt_vals)


    def read_qp(self, spin, kpoint, band, ignore_imag=False):
        """
        Return :class:`QpTempState` for the given (spin, kpoint, band)
        (NB: band is a global index i.e. unshifted)
        Only real part is returned if ``ignore_imag``.
        """
        spin, ikc, ibc, kpoint = self.get_sigma_skb_kpoint(spin, kpoint, band)

        # Quasi-particle energies set initially to zero
        qpe = np.zeros((self.ntemp),dtype = complex)
        qpe_oms = np.zeros((self.ntemp))

        # Fan and Debbye-Waller self-energies set to zero
        fan0 = np.zeros((self.ntemp),dtype = complex)
        dw = np.zeros((self.ntemp))

        # Renormalization of the Quasi-Particle set to zero
        ze0 = 0.0

        # real(dp) :: ks_enes( max_nbcalc, nkcalc, nsppol )
        # Kohn-Sham eigenvalues
        e0 = self.read_variable("ks_enes")[spin, ikc, ibc] * abu.Ha_eV


        return CumulantQpTempState(spin=spin, kpoint=kpoint, band=band, tmesh=self.tmesh,
                           e0=e0, qpe = qpe, ze0=ze0, fan0=fan0, dw=dw, qpe_oms = qpe_oms)






class CumulantSelfEnergy(EphSelfEnergy):
    
    

    def __init__(self, wmesh, qp, gw_vals, spfunccumul_wr, time_mesh = None, ct_vals = None, gt_vals = None, vals_e0ks=None, dvals_de0ks=None, dw_vals=None,
                 frohl_vals_e0ks=None, frohl_dvals_de0ks=None, frohl_spfunc_wr=None):
        
        # Set dimensions
        ntemp = len(qp.tmesh)
        nwr = len(wmesh)

        # In case of debugging, seting:
        # Green's function in time domain, cumulant function and time mesh
        self.gt_vals = gt_vals
        self.ct_vals = ct_vals
        self.time_mesh = time_mesh

        # Set Debye-Waller to zero, probably is the same as Dyson-Migdal since it is added 
        # separately from the cumulant calculation
        dw_vals = np.zeros((ntemp))

        # Calculation of the self-energy (vals_wr)
        # Calculation of the self-energy evaluated at the KS energy (vals_e0ks)
        # Calculation of the derivative of the self-energy evaluated at the KS energy (dvals_de0ks)
        self.gw = gw_vals
        vals_wr, vals_e0ks, dvals_de0ks = self.calculate_sigma_skb_fromgw(wmesh, qp, gw_vals)

        # Setting spectral function and Gre 
        self.spfunccumul_wr = spfunccumul_wr

        # Calculating QP energies
        qp.set_energies(wmesh, vals_wr, vals_e0ks, ntemp, nwr)

        super().__init__(wmesh, qp, vals_e0ks, dvals_de0ks, dw_vals, vals_wr, spfunccumul_wr,
                   frohl_vals_e0ks=None, frohl_dvals_de0ks=None, frohl_spfunc_wr=None)


    def calculate_sigma_skb_fromgw(self, wmesh, qp, gw_vals):

        # Initial setting
        nwr = len(wmesh)
        ntemp = len(qp.tmesh)
        sigma = np.zeros((ntemp,nwr),dtype=complex)
        e0ks = np.zeros((ntemp))
        de0ks = np.zeros((ntemp))
        
        # Identify the index of the KS energy is located ( probably at nwr//2 )
        idx_e0 = np.argmin( np.abs(wmesh - qp.e0))

        for it in range(ntemp):
            # Use Dyson equations to determine self-energy, only good for energies close to KS
            sigma[it,:] = wmesh[:]- qp.e0 - 1/gw_vals[it,:]
            # Evaluate the self-energy at KS energy
            e0ks[it] = sigma[it,idx_e0].real
            # The derivative of the self-energy evaluated at the KS energy.
            de0ks[it] = np.gradient(sigma[it,:].real, wmesh[1]-wmesh[0])[idx_e0]
        return sigma, e0ks, de0ks

    @classmethod
    def calculate_gw_from_sigeph_skb(cls, sigeph, time_tol = 1e-4):
        
        # Initialization
        wmesh_init = sigeph.wmesh
        sigma = sigeph.vals_wr
        sigma0 = sigeph.vals_e0ks.real
        e0 = sigeph.qp.e0
        nwr = sigeph.nwr
        ntemp = sigeph.ntemp
        qpe = np.zeros((ntemp),dtype = complex)
        qpe_oms = np.zeros((ntemp))
        fan0 = np.zeros((ntemp),dtype = complex)
        dw = np.zeros((ntemp))
        ze0 = 0.0
        qp = CumulantQpTempState(spin=sigeph.spin, kpoint=sigeph.kpoint, band=sigeph.band, tmesh=sigeph.tmesh,
                           e0=sigeph.qp.e0, qpe = qpe, ze0=ze0, fan0=fan0, dw=dw, qpe_oms = qpe_oms)

       
        # Frequency mesh step
        w_step = wmesh_init[1] - wmesh_init[0]
        # Shift frequency mesh to set the KS energy close to 0.0, with principal value offseting the mesh slightly.
        wmesh = wmesh_init - e0 - 0.5 * w_step
        # Create the time mesh
        t = np.linspace(0.0,1,nwr)*2*np.pi/w_step
        # Time mesh step
        t_step = t[1] - t[0]
        gw = np.zeros((ntemp,nwr), dtype=complex)
    
        for itemp in range(ntemp):
            #time_max = np.log(time_tol)/ ( -1.0 * np.abs(sigeph.vals_e0ks[itemp].imag))
            #Give warning if time_max < 2*np.pi/w_step
            
         
            beta = np.abs( sigma[itemp].imag)/np.pi

            # Calculate the part of the cumulant function that can be Fourier Transformed
            c1= fft((beta/(wmesh**2)),nwr)*w_step

            # The odd indexes give negative values
            c1[1::2] = -1.0 * c1[1::2]

            # Other terms of the cumulant function
            c2 = - 1j * t * sigma0
            c3 = -1.0 * np.trapz(beta/wmesh**2,x=wmesh) * np.ones(nwr)

            # Calculation of the cumulant function
            Ct = c1 + c2 + c3

            # Green's Function in time domain calculated with the cumulant function
            Gt = -1j * np.exp(Ct)

            # Green's function in frequency domain
            gw_temp = fftshift(ifft(Gt,nwr))*nwr*t_step
            R = 0.5* Gt[0]*t_step
            gw[itemp,:] = gw_temp - R

        # Spectral Function
        spfunc = -1.0 * gw.imag / np.pi


        return cls( wmesh_init, qp, gw, spfunc, time_mesh = t, ct_vals = c1+c2+c3, gt_vals = Gt )



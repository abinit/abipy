# coding: utf-8

import os
import abc
import numpy as np
import abipy.core.abinit_units as abu

from scipy.interpolate import UnivariateSpline
from monty.collections import dict2namedtuple
from monty.functools import lazy_property
from pymatgen.analysis.eos import EOS
from abipy.core.func1d import Function1D
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.electrons.gsr import GsrFile
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.phonons import PhononBandsPlotter, PhononDos, PhdosFile
from abipy.dfpt.gruneisen import GrunsNcFile


class QHA_App(metaclass=abc.ABCMeta):
    """
    Class for the approximations on QHA analysis.
    Provides some basic methods and plotting utils.
    These can be used to obtain other quantities and plots.
    Does not include electronic entropic contributions for metals.
    """

    def __init__(self, structures, structures_from_phdos, index_list, doses,energies, pressures, eos_name='vinet', pressure=0):
        """
        Args:
            structures: list of structures at different volumes.
            energies: list of SCF energies for the structures in eV.
            eos_name: string indicating the expression used to fit the energies. See pymatgen.analysis.eos.EOS.
            pressure: value of the pressure in GPa that will be considered in the p*V contribution to the energy.
        """
        self.structures = structures
        self.energies = np.array(energies)
        self.eos = EOS(eos_name)
        self.eos_name = eos_name
        self.pressure = pressure
        self.pressures = np.array(pressures)

        self.volumes = np.array([s.volume for s in structures])
        self.iv0 = np.argmin(energies)
        self.lattice_a = np.array([s.lattice.abc[0] for s in structures])
        self.lattice_b = np.array([s.lattice.abc[1] for s in structures])
        self.lattice_c = np.array([s.lattice.abc[2] for s in structures])

        self.angles_alpha = np.array([s.lattice.angles[0] for s in structures])
        self.angles_beta  = np.array([s.lattice.angles[1] for s in structures])
        self.angles_gama  = np.array([s.lattice.angles[2] for s in structures])

        self.doses = doses
        self.structures_from_phdos = np.array(structures_from_phdos)
        self.volumes_from_phdos = np.array([s.volume for s in structures_from_phdos])
        self.energies_pdos=self.energies[index_list]
        self.index_list = index_list
        if (len(self.index_list)==5):
            self.iv0_vib=1
            self.iv1_vib=3
            self.V0_vib=self.volumes_from_phdos[2]
        elif (len(self.index_list)==3):
            self.iv0_vib=0
            self.iv1_vib=2
            self.V0_vib=self.volumes_from_phdos[1]
        else :
            self.iv0_vib=0
            self.iv1_vib=1
            self.V0_vib=0.5*(self.volumes_from_phdos[1]+self.volumes_from_phdos[0])
        if abs(self.volumes_from_phdos[self.iv0_vib]+self.volumes_from_phdos[self.iv1_vib]-2*self.volumes[self.iv0])<1e-3 :
            self.scale_points="S"  # Symmetry
        else:
            self.scale_points="D"  # Displaced

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def fit_tot_energies(self, tstart=0, tstop=1000, num=101,tot_energies="energies" ,volumes="volumes"):
        """
        Performs a fit of the energies as a function of the volume at different temperatures.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns:
            `namedtuple` with the following attributes::

                tot_en: numpy array with shape (nvols, num) with the energies used for the fit
                fits: list of subclasses of pymatgen.analysis.eos.EOSBase, depending on the type of
                    eos chosen. Contains the fit for the energies at the different temperatures.
                min_en: numpy array with the minimum energies for the list of temperatures
                min_vol: numpy array with the minimum volumes for the list of temperatures
                temp: numpy array with the temperatures considered
        """
        # Generate a mesh of temperatures
        tmesh = np.linspace(tstart, tstop, num)

        # List of fit objects, one for each temperature
        fits = [self.eos.fit(volumes, e) for e in tot_energies.T]

        # Extract parameters from the fit objects
        v0 = np.array([fit.v0 for fit in fits])
        e0 = np.array([fit.e0 for fit in fits])
        b0 = np.array([fit.b0 for fit in fits])
        b1 = np.array([fit.b1 for fit in fits])

        # Minimum volumes and energies
        min_volumes = np.array([fit.v0 for fit in fits])
        min_energies = np.array([fit.e0 for fit in fits])

        v = min_volumes
        eta = (v / v0) ** (1.0 / 3.0)

        # Calculate the second derivative of free energy
        F2D = ( b0 / v0 * (-2.0 * (eta - 1) * eta ** -5.0 + (1 - 3.0 / 2.0 * (b1 - 1) * (eta - 1)) * eta ** -4.0)
            * np.exp(-3.0 * (b1 - 1.0) * (v ** (1.0 / 3.0) / v0 ** (1 / 3.0) - 1.0) / 2.0))

        return dict2namedtuple(tot_en=tot_energies, fits=fits, min_en=min_energies, min_vol=min_volumes, temp=tmesh , F2D=F2D)

    def second_derivative_energy_v(self, vol="vol"):
        """
        Performs a fit of the energies as a function of the volume at different temperatures.

        Args:
            vol: The volume at which to evaluate the second derivative of the energy.

        Returns:
            E2D_V: The second derivative of the energy with respect to volume.
        """
        # Initialize variables for temperature range (not used in the current function)
        tstart = 0
        tstop = 0
        num = 1

        tot_en = self.energies[np.newaxis, :].T
        fits = [self.eos.fit(self.volumes, e) for e in tot_en.T]

        # Extract parameters from the fit objects
        v0 = np.array([fit.v0 for fit in fits])
        e0 = np.array([fit.e0 for fit in fits])
        b0 = np.array([fit.b0 for fit in fits])
        b1 = np.array([fit.b1 for fit in fits])

        v = vol
        eta = (v / v0) ** (1.0 / 3.0)
        E2D_V = ( b0 / v0 * (-2.0 * (eta - 1) * eta ** -5.0 + (1 - 3.0 / 2.0 * (b1 - 1) * (eta - 1)) * eta ** -4.0) *
            np.exp(-3.0 * (b1 - 1.0) * (v ** (1.0 / 3.0) / v0 ** (1.0 / 3.0) - 1.0) / 2.0))

        return E2D_V
#=============================================================================================
    def vol_E2Vib1(self, tstart=0, tstop=1000, num=101):
        """
        Compute the volume as a function of temperature using the E2Vib1 method with the Vinet equation of state.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: Number of samples to generate. Default is 101.

        Returns:
            vol: The calculated volumes as a function of temperature.
            fits: The list of fit objects for the energies as a function of volume.
        """
        # Generate temperature mesh
        tmesh = np.linspace(tstart, tstop, num)

        # Get phonon free energies
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)

        vol = np.zeros(num)
        dfe_dV1 = np.zeros(num)
        volumes0 = self.volumes_from_phdos
        volumes = self.volumes
        iv0 = self.iv0_vib
        iv1 = self.iv1_vib
        dV = volumes0[iv1] - volumes0[iv0]
        V0 = volumes[self.iv0]
        E2D = self.second_derivative_energy_v(V0)

        # Calculate derivative of free energy with respect to volume and updated volumes
        for i, e in enumerate(ph_energies.T):
            dfe_dV1[i] = (e[iv1] - e[iv0]) / dV
            vol[i] = V0 - dfe_dV1[i] / E2D

        # Calculate total energies
        tot_en = self.energies[self.iv0] + 0.5*(volumes[np.newaxis, :].T - V0)**2*E2D +(volumes[np.newaxis, :].T - V0)*dfe_dV1

        # Fit the energies as a function of volume
        fits = [self.eos.fit(volumes, e) for e in tot_en.T]

        return vol, fits

#=============================================================================================
    def vol_Einf_Vib1(self, tstart=0, tstop=1000, num=101):
        """
        Compute the volume as a function of temperature using the EinfVib1 method with the Vinet equation of state.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: Number of samples to generate. Default is 101.

        Returns:
            vol: The calculated volumes as a function of temperature.
            fits: The list of fit objects for the energies as a function of volume.
        """
        # Generate temperature mesh
        tmesh = np.linspace(tstart, tstop, num)

        # Get phonon free energies
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)

        vol = np.zeros(num)
        dfe_dV1 = np.zeros(num)
        volumes0 = self.volumes_from_phdos
        volumes = self.volumes
        iv0 = self.iv0_vib
        iv1 = self.iv1_vib
        V0 = self.V0_vib

        dV = volumes0[iv1] - volumes0[iv0]

        # Calculate derivative of free energy with respect to volume
        for i, e in enumerate(ph_energies.T):
            dfe_dV1[i] = (e[iv1] - e[iv0]) / dV

        # Calculate total energies
        tot_en = self.energies[np.newaxis, :].T + (volumes[np.newaxis, :].T - V0) * dfe_dV1

        # Fit the energies as a function of volume
        fits = [self.eos.fit(volumes, e) for e in tot_en.T]

        # Extract minimum volumes from the fit objects
        vol = np.array([fit.v0 for fit in fits])

        return vol, fits

#=============================================================================================
    def vol_Einf_Vib2(self, tstart=0, tstop=1000, num=101):
        """
        Compute the volume as a function of temperature using the EinfVib2 method with the Vinet equation of state.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: Number of samples to generate. Default is 101.

        Returns:
            vol: The calculated volumes as a function of temperature.
            fits: The list of fit objects for the energies as a function of volume.
        """
        # Generate temperature mesh
        tmesh = np.linspace(tstart, tstop, num)

        # Get phonon free energies
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)

        vol = np.zeros(num)
        dfe_dV1 = np.zeros(num)
        dfe_dV2 = np.zeros(num)
        fe_V0 = np.zeros(num)

        # Determine index for volume calculations
        iv0 = 2 if len(self.index_list) == 5 else 1

        dV = self.volumes_from_phdos[iv0] - self.volumes_from_phdos[iv0 - 1]

        # Compute derivatives of free energy with respect to volume
        for i, e in enumerate(ph_energies.T):
            dfe_dV1[i] = (e[iv0 + 1] - e[iv0 - 1]) / (self.volumes_from_phdos[iv0 + 1] - self.volumes_from_phdos[iv0 - 1])
            dfe_dV2[i] = (e[iv0 + 1] - 2.0 * e[iv0] + e[iv0 - 1]) / (self.volumes_from_phdos[iv0 + 1] - self.volumes_from_phdos[iv0])**2
            fe_V0[i] = e[iv0]

        # Reference volume
        V0 = self.volumes_from_phdos[iv0]

        # Calculate total energies
        tot_en = ( self.energies[np.newaxis, :].T +fe_V0
              +(self.volumes[np.newaxis, :].T - V0) * dfe_dV1 + 0.5 * (self.volumes[np.newaxis, :].T - V0)**2 * dfe_dV2)

        # Fit the energies as a function of volume
        fits = [self.eos.fit(self.volumes, e) for e in tot_en.T]

        # Extract minimum volumes from the fit objects
        vol = np.array([fit.v0 for fit in fits])

        return vol, fits

#================================================================================================================

    def vol_Einf_Vib4(self, tstart=0, tstop=1000, num=101):
        """
        Compute the volume as a function of temperature using the EinfVib4 method with the Vinet equation of state.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: Number of samples to generate. Default is 101.

        Returns:
            vol: The calculated volumes as a function of temperature.
            fits: The list of fit objects for the energies as a function of volume.
        """

        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        energies = self.energies
        volumes0= self.volumes_from_phdos
        volumes= self.volumes
        vol = np.zeros( num)
        dfe_dV1= np.zeros( num)
        dfe_dV2= np.zeros( num)
        dfe_dV3= np.zeros( num)
        dfe_dV4= np.zeros( num)
        fe_V0= np.zeros( num)
        iv0=2
        dV=volumes0[2]-volumes0[1]

        for i,e in enumerate(ph_energies.T):
            dfe_dV1[i]=(-e[iv0+2]+ 8*e[iv0+1]-8*e[iv0-1]+e[iv0-2])/(12*dV)
            dfe_dV2[i]=(-e[iv0+2]+16*e[iv0+1]-30*e[iv0]+16*e[iv0-1]-e[iv0-2])/(12*dV**2)
            dfe_dV3[i]=(e[iv0+2]-2*e[iv0+1]+2*e[iv0-1]-e[iv0-2])/(2*dV**3)
            dfe_dV4[i]=(e[iv0+2]-4*e[iv0+1]+6*e[iv0]-4*e[iv0-1]+e[iv0-2])/(dV**4)

            fe_V0[i]= e[iv0]
        V0=volumes0[iv0]

        tot_en = (( volumes[np.newaxis, :].T -V0) * dfe_dV1  + 0.5* ( volumes[np.newaxis, :].T -V0)**2*(dfe_dV2)
                 +( volumes[np.newaxis, :].T -V0)**3 *dfe_dV3/6.0  + ( volumes[np.newaxis, :].T -V0)**4*(dfe_dV4/24.0)
                 + fe_V0[:] +energies[np.newaxis, :].T  )

        fits = [self.eos.fit(volumes, e) for e in tot_en.T]
        vol = np.array([fit.v0 for fit in fits])

        return vol , fits
#*********************************************************************************************
    @property
    def nvols(self):
        """Number of volumes"""
        return len(self.volumes_from_phdos)

    def set_eos(self, eos_name):
        """
        Set the EOS model used for the fit.

        Args:
            eos_name: string indicating the expression used to fit the energies. See pymatgen.analysis.eos.EOS.
        """
        self.eos = EOS(eos_name)
        self.eos_name = eos_name
        if eos_name != "vinet":
            raise RuntimeError("This approximation method is only developed for the Vinet equation of state.")

    @add_fig_kwargs
    def plot_energies(self, tstart=0, tstop=1000, num=1, ax=None, **kwargs):
        """
        Plots the BO energy as a function of volume

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 10.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """

        tmesh = np.linspace(tstart, tstop, num)

        f = self.fit_tot_energies(0, 0, 1 ,self.energies[np.newaxis, :].T,self.volumes)

        ax, fig, plt = get_ax_fig_plt(ax)
        xmin, xmax = np.floor(self.volumes.min() * 0.97), np.ceil(self.volumes.max() * 1.03)
        x = np.linspace(xmin, xmax, 100)

        for fit, e, t in zip(f.fits, f.tot_en.T - self.energies[self.iv0], f.temp):
            ax.scatter(self.volumes, e, label=t, color='b', marker='s', s=10)
            ax.plot(x, fit.func(x) - self.energies[self.iv0], color='b', lw=1)

        ax.plot(f.min_vol, f.min_en - self.energies[self.iv0], color='r', linestyle='dashed' , lw=1, marker='o', ms=5)

        ax.set_xlabel(r'V (${\AA}^3$)')
        ax.set_ylabel('E (eV)')
        return fig

    @add_fig_kwargs
    def plot_vol_vs_t(self, tstart=0, tstop=1000, num=101, ax=None, **kwargs):
        """
        Plot the volume as a function of temperature using various methods.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 101.
            ax: Matplotlib Axes object or None. If None, a new figure will be created.
            **kwargs: Additional keyword arguments to pass to plotting functions.

        Returns:
            fig: The Matplotlib Figure object.
        """
        # Get or create the matplotlib Axes and Figure
        ax, fig, plt = get_ax_fig_plt(ax)

        # Generate temperature mesh
        tmesh = np.linspace(tstart, tstop, num)

        # Get phonon free energies
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)

        # Initialize data storage
        iv0 = self.iv0_vib
        iv1 = self.iv1_vib
        volumes = self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']

        # Method 1: E2Vib1
        if self.scale_points=="S":
            vol, _ = self.vol_E2Vib1(tstart=tstart, tstop=tstop, num=num)
            ax.plot(tmesh, vol, color='b', lw=2, label="E2Vib1")
            data_to_save = np.column_stack((data_to_save, vol))
            columns.append('E2vib1')

        # Method 2: Einf_Vib1
        if len(self.index_list) >= 2:
            vol2, _ = self.vol_Einf_Vib1(tstart=tstart, tstop=tstop, num=num)
            ax.plot(tmesh, vol2, color='gold', lw=2, label=r"$E_\infty$ Vib1")
            data_to_save = np.column_stack((data_to_save, vol2))
            columns.append('Einfvib1')

        # Method 3: Einf_Vib2
        if len(self.index_list) >= 3:
            vol3, _ = self.vol_Einf_Vib2(tstart=tstart, tstop=tstop, num=num)
            ax.plot(tmesh, vol3, color='m', lw=2, label=r"$E_\infty$ Vib2")
            data_to_save = np.column_stack((data_to_save, vol3))
            columns.append('Einfvib2')

        # Method 4: Einf_Vib4 and QHA
        if len(self.index_list) == 5:
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_tot_energies(tstart, tstop, num, tot_en, self.volumes_from_phdos)
            vol4, _ = self.vol_Einf_Vib4(tstart=tstart, tstop=tstop, num=num)
            ax.plot(tmesh, vol4, color='c', lw=2, label=r"$E_\infty$ Vib4")
            ax.plot(tmesh, f0.min_vol, color='k', linestyle='dashed', lw=1.5, label="QHA")
            data_to_save = np.column_stack((data_to_save, vol4, f0.min_vol))
            columns.append('Einfvib4')
            columns.append('QHA')

        # Plot V0
        ax.plot(0, self.volumes[self.iv0], color='g', lw=0, marker='o', ms=10, label="V0")

        # Set labels and limits
        ax.set_xlabel('T (K)')
        ax.set_ylabel(r'V (${\AA}^3$)')
        ax.set_xlim(tstart, tstop)
        ax.grid(True)
        ax.legend()

        return fig

###################################################################################################
    def get_thermal_expansion_coeff(self, tstart=0, tstop=1000, num=101, tref=None):
        """
        Calculates the thermal expansion coefficient as a function of temperature

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: Number of samples to generate. Default is 101.

        Returns:
            Function1D: The thermal expansion coefficient as a function of temperature.
        """
        # Get phonon free energies
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
        f = self.fit_tot_energies(tstart, tstop, num, tot_en, self.volumes_from_phdos)

        if tref != None:
            ph_energies2 = self.get_vib_free_energies(tref, tref, 1)
            tot_en2 = self.energies_pdos[np.newaxis, :].T + ph_energies2
            f0 = self.fit_tot_energies(tref, tref, 1, tot_en2, self.volumes_from_phdos)

        dt = f.temp[1] - f.temp[0]

        # Get thermodynamic properties
        thermo = self.get_thermodynamic_properties(tstart, tstop, num)
        entropy = thermo.entropy.T
        df_t = -entropy

        param = np.zeros((num, 4))
        param2 = np.zeros((num, 3))
        d2f_t_v = np.zeros(num)

        for j in range(num):
            param[j] = np.polyfit(self.volumes_from_phdos, df_t[j], 3)
            param2[j] = np.array([3 * param[j][0], 2 * param[j][1], param[j][2]])
            p = np.poly1d(param2[j])
            d2f_t_v[j] = p(f.min_vol[j])

        F2D = f.F2D
        if tref == None:
            alpha = -1 / f.min_vol * d2f_t_v / F2D
        else:
            alpha = -1 / f0.min_vol * d2f_t_v / F2D

        return Function1D(f.temp, alpha)

    @add_fig_kwargs
    def plot_thermal_expansion_coeff(self, tstart=0, tstop=1000, num=101, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        tmesh = np.linspace(tstart, tstop, num)
        thermo = self.get_thermodynamic_properties(tstart=tstart, tstop=tstop, num=num)
        entropy = thermo.entropy.T #* abu.e_Cb * abu.Avogadro
        df_t = np.zeros((num,self.nvols))
        df_t = - entropy
        volumes=self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        dV=volumes[iv0+1]-volumes[iv0]

        if self.scale_points=="S":
            vol,fits = self.vol_E2Vib1(num=num,tstop=tstop,tstart=tstart)
            E2D = self.second_derivative_energy_v(self.volumes[self.iv0])
            #f = self.fit_tot_energies(0, 0, 1 ,self.energies[np.newaxis, :].T,self.volumes)
            #E2D = f.F2D
            if (tref==None):
                alpha_1 = - 1/vol[:] * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D
            else:
                vol_ref,fits = self.vol_E2Vib1(num=1,tstop=tref,tstart=tref)
                b0 = np.array([fit.b0 for fit in fits])
                print("B (E2vib1)   @ ",tref," K =",b0*160.21766208 ,"(GPa)" )
                alpha_1 = - 1/vol_ref * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D
            ax.plot(tmesh, alpha_1,color='b', lw=2, label="E2Vib1")
            data_to_save = np.column_stack((data_to_save,alpha_1))
            columns.append( 'E2vib1')

        if (len(self.index_list)>=2):
            vol2 ,fits = self.vol_Einf_Vib1(num=num,tstop=tstop,tstart=tstart)
            E2D_V = self.second_derivative_energy_v(vol2)
            if (tref==None):
                alpha_2 = - 1/vol2[:] * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D_V[:]
            else:
                vol2_ref,fits  = self.vol_Einf_Vib1(num=1,tstop=tref,tstart=tref)
                b0 = np.array([fit.b0 for fit in fits])
                print("B (Einfvib1) @ ",tref," K =",b0*160.21766208 ,"(GPa)" )
                alpha_2 = - 1/vol2_ref * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D_V[:]
            ax.plot(tmesh, alpha_2,color='gold', lw=2 ,  label=r"$E_\infty Vib1$")
            data_to_save = np.column_stack((data_to_save,alpha_2))
            columns.append( 'Einfvib1')

        if (len(self.index_list)>=3):
            vol3,fits = self.vol_Einf_Vib2(num=num,tstop=tstop,tstart=tstart)
            E2D_V = self.second_derivative_energy_v(vol3)
            dfe_dV2= np.zeros( num)
            for i,e in enumerate(ph_energies.T):
                dfe_dV2[i]=(e[iv0+2]-2.0*e[iv0+1]+e[iv0])/(dV)**2

            ds_dv = (df_t[:,iv0+2]-df_t[:,iv0])/(2*dV)
            ds_dv = ds_dv+ (df_t[:,iv0+2]-2*df_t[:,iv0+1]+df_t[:,iv0])/dV**2 * (vol3[:]-volumes[iv0+1])
            if (tref==None):
                alpha_3 = - 1/vol3[:] * ds_dv / (E2D_V[:]+dfe_dV2[:])
            else:
                vol3_ref,fits = self.vol_Einf_Vib2(num=1,tstop=tref,tstart=tref)
                b0 = np.array([fit.b0 for fit in fits])
                print("B (Einfvib2) @ ",tref," K =",b0*160.21766208 ,"(GPa)" )
                alpha_3 = - 1/vol3_ref * ds_dv / (E2D_V[:]+dfe_dV2[:])
            ax.plot(tmesh, alpha_3,color='m', lw=2 ,  label=r"$E_\infty Vib2$")
            data_to_save = np.column_stack((data_to_save,alpha_3))
            columns.append( 'Einfvib2')

        if (len(self.index_list)==5):
            alpha_qha  = self.get_thermal_expansion_coeff(tstart, tstop, num, tref)
            vol4,fits = self.vol_Einf_Vib4(num=num,tstop=tstop,tstart=tstart)
            E2D_V = self.second_derivative_energy_v(vol4)

            d2fe_dV2= np.zeros( num)
            d3fe_dV3= np.zeros( num)
            d4fe_dV4= np.zeros( num)
            for i,e in enumerate(ph_energies.T):
                d2fe_dV2[i]=(-e[4]+16*e[3]-30*e[2]+16*e[1]-e[0])/(12*dV**2)
                d3fe_dV3[i]=(e[4]-2*e[3]+2*e[1]-e[0])/(2*dV**3)
                d4fe_dV4[i]=(e[4]-4*e[3]+6*e[2]-4*e[1]+e[0])/(dV**4)

            ds_dv =(-df_t[:,4]+ 8*df_t[:,3]-8*df_t[:,1]+df_t[:,0])/(12*dV)
            ds_dv = ds_dv+ (-df_t[:,4]+16*df_t[:,3]-30*df_t[:,2]+16*df_t[:,1]-df_t[:,0])/(12*dV**2) * (vol4[:]-volumes[2])
            ds_dv = ds_dv+ 1.0/2.0*(df_t[:,4]-2*df_t[:,3]+2*df_t[:,1]-df_t[:,0])/(2*dV**3) * (vol4[:]-volumes[2])**2
            ds_dv = ds_dv+ 1.0/6.0* (df_t[:,4]-4*df_t[:,3]+6*df_t[:,2]-4*df_t[:,1]+df_t[:,0])/(dV**4)* (vol4[:]-volumes[2])**3
            D2F=E2D_V[:]+d2fe_dV2[:]+ (vol4[:]-volumes[2])*d3fe_dV3[:]+0.5*(vol4[:]-volumes[2])**2*d4fe_dV4[:]
            if (tref==None):
                alpha_4 = - 1/vol4[:] * ds_dv / D2F
            else:
                vol4_ref,fits = self.vol_Einf_Vib4(num=1,tstop=tref,tstart=tref)
                b0 = np.array([fit.b0 for fit in fits])
                print("B (Einfvib4) @ ",tref," K =",b0*160.21766208 ,"(GPa)" )
                alpha_4 = - 1/vol4_ref * ds_dv / D2F

            ax.plot(tmesh, alpha_4,color='c',linewidth=2 ,  label=r"$E_\infty Vib4$")
            ax.plot(alpha_qha.mesh, alpha_qha.values, color='k',linestyle='dashed', lw=1.5 ,label="QHA")
            data_to_save = np.column_stack((data_to_save,alpha_4,alpha_qha.values))
            columns.append( 'Einfvib4')
            columns.append( 'QHA')


        ax.set_xlabel(r'T (K)')
        ax.set_ylabel(r'$\alpha$ (K$^{-1}$)')
        ax.grid(True)
        ax.legend()

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
    def get_abc(self, tstart=0, tstop=1000, num=101,volumes="volumes"):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        param = np.zeros((num,4))
        param=np.polyfit(self.volumes, self.lattice_a , 3)
        pa = np.poly1d(param)
        aa_qha=pa(volumes)
        param=np.polyfit(self.volumes, self.lattice_b , 3)
        pb = np.poly1d(param)
        bb_qha=pb(volumes)
        param=np.polyfit(self.volumes, self.lattice_c , 3)
        pc = np.poly1d(param)
        cc_qha=pc(volumes)

        return aa_qha,bb_qha,cc_qha

    def get_angles(self, tstart=0, tstop=1000, num=101,volumes="volumes"):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        param = np.zeros((num,4))
        param=np.polyfit(self.volumes, self.angles_alpha, 3)
        pa = np.poly1d(param)
        gamma=pa(volumes)
        param=np.polyfit(self.volumes, self.angles_beta , 3)
        pb = np.poly1d(param)
        beta=pb(volumes)
        param=np.polyfit(self.volumes, self.angles_gama , 3)
        pc = np.poly1d(param)
        alpha=pc(volumes)

        return alpha,beta,gamma

###################################################################################################
    @add_fig_kwargs
    def plot_thermal_expansion_coeff_abc(self, tstart=0, tstop=1000, num=101, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos
        data_to_save = tmesh[1:-1]
        columns = ['#Tmesh']

        if self.scale_points=="S":
            vol2 ,fits= self.vol_E2Vib1(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol2_tref,fits = self.vol_E2Vib1(num=1,tstop=tref,tstart=tref)

        if (len(self.index_list)==2):
            vol ,fits = self.vol_Einf_Vib1(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref,fits  = self.vol_Einf_Vib1(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol,fits = self.vol_Einf_Vib2(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref,fits = self.vol_Einf_Vib2(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_tot_energies(tstart, tstop, num,tot_en , self.volumes_from_phdos)
            method =r"$ (E_\infty Vib4)$"
            vol,fits = self.vol_Einf_Vib4(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref,fits = self.vol_Einf_Vib4(num=1,tstop=tref,tstart=tref)

        #alpha  = self.get_thermal_expansion_coeff(tstart, tstop, num, tref)
        tmesh = np.linspace(tstart, tstop, num)
        dt= tmesh[1] - tmesh[0]

        aa,bb,cc = self.get_abc(tstart, tstop, num,vol)
        if (tref!=None):
            aa_tref,bb_tref,cc_tref = self.get_abc(tref, tref, 1,vol_tref)

        alpha_a = np.zeros( num-2)
        alpha_b = np.zeros( num-2)
        alpha_c = np.zeros( num-2)
        if (tref==None):
            alpha_a = (aa[2:] - aa[:-2]) / (2 * dt) / aa[1:-1]
            alpha_b = (bb[2:] - bb[:-2]) / (2 * dt) / bb[1:-1]
            alpha_c = (cc[2:] - cc[:-2]) / (2 * dt) / cc[1:-1]
        else:
            alpha_a = (aa[2:] - aa[:-2]) / (2 * dt) / aa_tref
            alpha_b = (bb[2:] - bb[:-2]) / (2 * dt) / bb_tref
            alpha_c = (cc[2:] - cc[:-2]) / (2 * dt) / cc_tref

        ax.plot(tmesh[1:-1] ,alpha_a , color='r', lw=2,label = r"$\alpha_a$"+method, **kwargs)
        ax.plot(tmesh[1:-1] ,alpha_b , color='b', lw=2,label = r"$\alpha_b$"+method)
        ax.plot(tmesh[1:-1] ,alpha_c , color='m', lw=2,label = r"$\alpha_c$"+method)

        method_header=method+"  (alpha_a,alpha_b,alpha_c) |"
        data_to_save = np.column_stack((data_to_save,alpha_a,alpha_b,alpha_c))
        columns.append( method_header)

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            aa2,bb2,cc2 = self.get_abc(tstart, tstop, num,vol2)
            if (tref!=None):
                aa2_tref,bb2_tref,cc2_tref = self.get_abc(tref, tref, 1,vol2_tref)

            alpha2_a = np.zeros( num-2)
            alpha2_b = np.zeros( num-2)
            alpha2_c = np.zeros( num-2)
            if (tref==None):
                alpha2_a = (aa2[2:] - aa2[:-2]) / (2 * dt) / aa2[1:-1]
                alpha2_b = (bb2[2:] - bb2[:-2]) / (2 * dt) / bb2[1:-1]
                alpha2_c = (cc2[2:] - cc2[:-2]) / (2 * dt) / cc2[1:-1]
            else:
                alpha2_a = (aa2[2:] - aa2[:-2]) / (2 * dt) / aa2_tref
                alpha2_b = (bb2[2:] - bb2[:-2]) / (2 * dt) / bb2_tref
                alpha2_c = (cc2[2:] - cc2[:-2]) / (2 * dt) / cc2_tref

            ax.plot(tmesh[1:-1] ,alpha2_a ,  linestyle='dashed' ,  color='r', lw=2 ,label = r"$\alpha_a$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_b ,  linestyle='dashed' ,  color='b', lw=2 ,label = r"$\alpha_b$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_c ,  linestyle='dashed' ,  color='m', lw=2 ,label = r"$\alpha_c$"" (E2vib1)")
            data_to_save = np.column_stack((data_to_save,alpha2_a,alpha2_b,alpha2_c))
            columns.append( 'E2vib1 (alpha_a,alpha_b,alpha_c)   ')

        ax.set_xlabel(r'T (K)')
        ax.set_ylabel(r'$\alpha$ (K$^{-1}$)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
    @add_fig_kwargs
    def plot_thermal_expansion_coeff_angles(self, tstart=0, tstop=1000, num=101, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos
        data_to_save = tmesh[1:-1]
        columns = ['#Tmesh']

        if self.scale_points=="S":
            vol2 ,fits= self.vol_E2Vib1(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol2_tref,fits = self.vol_E2Vib1(num=1,tstop=tref,tstart=tref)

        if (len(self.index_list)==2):
            vol ,fits = self.vol_Einf_Vib1(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref,fits  = self.vol_Einf_Vib1(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol,fits = self.vol_Einf_Vib2(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref,fits = self.vol_Einf_Vib2(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_tot_energies(tstart, tstop, num,tot_en , self.volumes_from_phdos)
            method =r"$ (E_\infty Vib4)$"
            vol,fits = self.vol_Einf_Vib4(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref,fits = self.vol_Einf_Vib4(num=1,tstop=tref,tstart=tref)

        #alpha  = self.get_thermal_expansion_coeff(tstart, tstop, num, tref)
        tmesh = np.linspace(tstart, tstop, num)
        dt= tmesh[1] - tmesh[0]

        alpha,beta,cc = self.get_angles(tstart, tstop, num,vol)
        if (tref!=None):
            alpha_tref,beta_tref,cc_tref = self.get_angles(tref, tref, 1,vol_tref)

        alpha_alpha = np.zeros( num-2)
        alpha_beta = np.zeros( num-2)
        alpha_gamma = np.zeros( num-2)
        if (tref==None):
            alpha_alpha = (alpha[2:] - alpha[:-2]) / (2 * dt) / alpha[1:-1]
            alpha_beta = (beta[2:] - beta[:-2]) / (2 * dt) / beta[1:-1]
            alpha_gamma = (cc[2:] - cc[:-2]) / (2 * dt) / cc[1:-1]
        else:
            alpha_alpha = (alpha[2:] - alpha[:-2]) / (2 * dt) / alpha_tref
            alpha_beta = (beta[2:] - beta[:-2]) / (2 * dt) / beta_tref
            alpha_gamma = (cc[2:] - cc[:-2]) / (2 * dt) / cc_tref

        ax.plot(tmesh[1:-1] ,alpha_alpha , color='r', lw=2,label = r"$\alpha_alpha$"+method, **kwargs)
        ax.plot(tmesh[1:-1] ,alpha_beta , color='b', lw=2,label = r"$\alpha_beta$"+method)
        ax.plot(tmesh[1:-1] ,alpha_gamma , color='m', lw=2,label = r"$\alpha_gamma$"+method)

        method_header=method+"  (alpha_alpha,alpha_beta,alpha_gamma) |"
        data_to_save = np.column_stack((data_to_save,alpha_alpha,alpha_beta,alpha_gamma))
        columns.append( method_header)

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            alpha2,beta2,cc2 = self.get_angles(tstart, tstop, num,vol2)
            if (tref!=None):
                alpha2_tref,beta2_tref,cc2_tref = self.get_angles(tref, tref, 1,vol2_tref)

            alpha2_alpha = np.zeros( num-2)
            alpha2_beta = np.zeros( num-2)
            alpha2_gamma = np.zeros( num-2)
            if (tref==None):
                alpha2_alpha = (alpha2[2:] - alpha2[:-2]) / (2 * dt) / alpha2[1:-1]
                alpha2_beta = (beta2[2:] - beta2[:-2]) / (2 * dt) / beta2[1:-1]
                alpha2_gamma = (cc2[2:] - cc2[:-2]) / (2 * dt) / cc2[1:-1]
            else:
                alpha2_alpha = (alpha2[2:] - alpha2[:-2]) / (2 * dt) / alpha2_tref
                alpha2_beta = (beta2[2:] - beta2[:-2]) / (2 * dt) / beta2_tref
                alpha2_gamma = (cc2[2:] - cc2[:-2]) / (2 * dt) / cc2_tref

            ax.plot(tmesh[1:-1] ,alpha2_alpha ,  linestyle='dashed' ,  color='r', lw=2 ,label = r"$\alpha_alpha$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_beta ,  linestyle='dashed' ,  color='b', lw=2 ,label = r"$\alpha_beta$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_gamma ,  linestyle='dashed' ,  color='m', lw=2 ,label = r"$\alpha_gamma$"" (E2vib1)")
            data_to_save = np.column_stack((data_to_save,alpha2_alpha,alpha2_beta,alpha2_gamma))
            columns.append( 'E2vib1 (alpha_alpha,alpha_beta,alpha_gamma)   ')

        ax.set_xlabel(r'T (K)')
        ax.set_ylabel(r'$\alpha$ (K$^{-1}$)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig

    @add_fig_kwargs
    def plot_abc_vs_t(self, tstart=0, tstop=1000, num=101, lattice=None, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']
        if self.scale_points=="S":
            vol2,fits = self.vol_E2Vib1(num=num,tstop=tstop,tstart=tstart)
            aa2,bb2,cc2 = self.get_abc(tstart, tstop, num,vol2)
            data_to_save = np.column_stack((data_to_save,aa2,bb2,cc2))
            columns.append( 'E2vib1 (a,b,c) |            ')

        if (len(self.index_list)==2):
            vol ,fits = self.vol_Einf_Vib1(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol,fits = self.vol_Einf_Vib2(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_tot_energies(tstart, tstop, num,tot_en , self.volumes_from_phdos)
            method =r"$ (E_\infty Vib4)$"
            vol,fits = self.vol_Einf_Vib4(num=num,tstop=tstop,tstart=tstart)
        aa,bb,cc = self.get_abc(tstart, tstop, num,vol)

        method_header=method+"  (a,b,c) |"
        data_to_save = np.column_stack((data_to_save,aa,bb,cc))
        columns.append( method_header)

        if (lattice==None or lattice=="a"):
            ax.plot(tmesh ,aa , color='r', lw=2,label = r"$a(V(T))$"+method,  **kwargs )
        if (lattice==None or lattice=="b"):
            ax.plot(tmesh ,bb , color='b', lw=2,label = r"$b(V(T))$"+method )
        if (lattice==None or lattice=="c"):
            ax.plot(tmesh ,cc , color='m', lw=2,label = r"$c(V(T))$"+method )

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            if (lattice==None or lattice=="a"):
                ax.plot(tmesh ,aa2 ,  linestyle='dashed' , color='r', lw=2,label = r"$a(V(T))$""E2vib1" )
            if (lattice==None or lattice=="b"):
                ax.plot(tmesh ,bb2 ,  linestyle='dashed' , color='b', lw=2,label = r"$b(V(T))$""E2vib1" )
            if (lattice==None or lattice=="c"):
                ax.plot(tmesh ,cc2 ,  linestyle='dashed' , color='m', lw=2,label = r"$c(V(T))$""E2vib1" )

        ax.set_xlabel(r'T (K)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
    @add_fig_kwargs
    def plot_angles_vs_t(self, tstart=0, tstop=1000, num=101, angle=None, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']
        if self.scale_points=="S":
            vol2,fits = self.vol_E2Vib1(num=num,tstop=tstop,tstart=tstart)
            alpha2,beta2,gamma2 = self.get_angles(tstart, tstop, num,vol2)
            data_to_save = np.column_stack((data_to_save,alpha2,beta2,gamma2))
            columns.append( 'E2vib1 (alpha,beta,gamma) |            ')

        if (len(self.index_list)==2):
            vol ,fits = self.vol_Einf_Vib1(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol,fits = self.vol_Einf_Vib2(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_tot_energies(tstart, tstop, num,tot_en , self.volumes_from_phdos)
            method =r"$ (E_\infty Vib4)$"
            vol,fits = self.vol_Einf_Vib4(num=num,tstop=tstop,tstart=tstart)
        alpha,beta,gamma = self.get_angles(tstart, tstop, num,vol)

        method_header=method+"  (alpha,beta,gamm) |"
        data_to_save = np.column_stack((data_to_save,alpha,beta,gamma))
        columns.append( method_header)

        if (angle==None or angle==1):
            ax.plot(tmesh ,alpha , color='r', lw=2,label = r"$alpha(V(T))$"+method,  **kwargs )
        if (angle==None or angle==2):
            ax.plot(tmesh ,beta , color='b', lw=2,label = r"$beta(V(T))$"+method )
        if (angle==None or angle==3):
            ax.plot(tmesh ,gamma , color='m', lw=2,label = r"$gamma(V(T))$"+method )

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            if (angle==None or angle==1):
                ax.plot(tmesh ,alpha2 ,  linestyle='dashed' , color='r', lw=2,label = r"$alpha(V(T))$""E2vib1" )
            if (angle==None or angle==2):
                ax.plot(tmesh ,beta2 ,  linestyle='dashed' , color='b', lw=2,label = r"$beta(V(T))$""E2vib1" )
            if (angle==None or angle==3):
                ax.plot(tmesh ,gamma2 ,  linestyle='dashed' , color='m', lw=2,label = r"$gamma(V(T))$""E2vib1" )

        ax.set_xlabel(r'T (K)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    def fit_forth(self, tstart=0, tstop=1000, num=1,energy="energy",volumes="volumes"):
        """
        Performs a fit of the energies as a function of the volume at different temperatures.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns:
            `namedtuple` with the following attributes::

                tot_en: numpy array with shape (nvols, num) with the energies used for the fit
                fits: list of subclasses of pymatgen.analysis.eos.EOSBase, depending on the type of
                    eos chosen. Contains the fit for the energies at the different temperatures.
                min_en: numpy array with the minimum energies for the list of temperatures
                min_vol: numpy array with the minimum volumes for the list of temperatures
                temp: numpy array with the temperatures considered
        """
        tmesh = np.linspace(tstart, tstop, num)

        param = np.zeros((num,5))
        param2 = np.zeros((num,4))
        param3 = np.zeros((num,3))
        min_vol = np.zeros((num))
        min_en = np.zeros((num))
        F2D_V = np.zeros((num))
        for j,e in enumerate(energy.T):
            param[j]=np.polyfit(volumes,e , 4)
            param2[j]=np.array([4*param[j][0],3*param[j][1],2*param[j][2],param[j][3]])
            param3[j]=np.array([12*param[j][0],6*param[j][1],2*param[j][2]])
            p = np.poly1d(param[j])
            p2 = np.poly1d(param2[j])
            p3 = np.poly1d(param3[j])
            min_vol[j]=self.volumes[self.iv0]
            vv=self.volumes[self.iv0]
            while p2(min_vol[j])**2 > 1e-16 :
                min_vol[j]=min_vol[j]-p2(min_vol[j])*10
            min_en[j]=p(min_vol[j])
            F2D_V[j]=p3(min_vol[j])

        return dict2namedtuple(min_vol=min_vol, temp=tmesh , min_en=min_en , param=param , F2D_V=F2D_V)#, fits=fits)

    def vol_E2Vib1_forth(self, tstart=0, tstop=1000, num=101):

        volumes0= self.volumes_from_phdos
        iv0=self.iv0_vib
        iv1=self.iv1_vib

        dV=volumes0[iv1]-volumes0[iv0]
        V0=self.volumes[self.iv0]

        energy = self.energies[np.newaxis, :].T
        f=self.fit_forth( tstart=0, tstop=0, num=1 ,energy=energy,volumes=self.volumes)
        param3 = np.zeros((num,3))
        param3=np.array([12*f.param[0][0],6*f.param[0][1],2*f.param[0][2]])
        p3 = np.poly1d(param3)
        E2D = p3(V0)

        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        vol = np.zeros( num)

        for i,e in enumerate(ph_energies.T):
            dfe_dV1=(e[iv1]-e[iv0])/dV
            vol[i]=V0-dfe_dV1*E2D**-1

        return vol

    def vol_EinfVib1_forth(self, tstart=0, tstop=1000, num=101):
        """
        Plot the volume as a function of the temperature.
        Returns: |Vol|
        """

        volumes0= self.volumes_from_phdos
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        V0= self.V0_vib

        dV=volumes0[iv1]-volumes0[iv0]

        energy = self.energies[np.newaxis, :].T
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        vol = np.zeros( num)

        dfe_dV= np.zeros( num)

        for i,e in enumerate(ph_energies.T):
            dfe_dV[i]=(e[iv1]-e[iv0])/dV
        tot_en = self.energies[np.newaxis, :].T + ( self.volumes[np.newaxis, :].T -V0) * dfe_dV

        f=self.fit_forth( tstart, tstop, num ,tot_en,self.volumes)
        vol=f.min_vol

        return vol
    def vol_Einf_Vib2_forth(self, tstart=0, tstop=1000, num=101):
        """
        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns: |Vol|
        """
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        energies = self.energies
        volumes0= self.volumes_from_phdos
        volumes= self.volumes
        vol = np.zeros( num)
        dfe_dV= np.zeros( num)
        d2fe_dV2= np.zeros( num)
        fe_V0= np.zeros( num)
        if (len(self.index_list)==5):
            iv0=2
        else :
            iv0=1
        dV=volumes0[iv0]-volumes0[iv0-1]
        V0=volumes0[iv0]
        for i,e in enumerate(ph_energies.T):
            dfe_dV[i]=(e[iv0+1]-e[iv0-1])/(2*dV)
            d2fe_dV2[i]=(e[iv0+1]-2.0*e[iv0]+e[iv0-1])/(dV)**2
            fe_V0[i]= e[iv0]

        tot_en = self.energies[np.newaxis, :].T + ( self.volumes[np.newaxis, :].T -V0) * dfe_dV
        tot_en  = tot_en + 0.5* (( self.volumes[np.newaxis, :].T -V0))**2*(d2fe_dV2)
        tot_en = tot_en+ fe_V0
        f=self.fit_forth( tstart, tstop, num ,tot_en,self.volumes)
        vol=f.min_vol

        return vol
    def vol_Einf_Vib4_forth(self, tstart=0, tstop=1000, num=101):
        """

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns: |Vol|
        """

        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        energies = self.energies
        volumes0= self.volumes_from_phdos
        volumes= self.volumes
        vol = np.zeros( num)
        dfe_dV1= np.zeros( num)
        dfe_dV2= np.zeros( num)
        dfe_dV3= np.zeros( num)
        dfe_dV4= np.zeros( num)
        fe_V0= np.zeros( num)
        iv0=2
        dV=volumes0[2]-volumes0[1]

        for i,e in enumerate(ph_energies.T):
            dfe_dV1[i]=(-e[iv0+2]+ 8*e[iv0+1]-8*e[iv0-1]+e[iv0-2])/(12*dV)
            dfe_dV2[i]=(-e[iv0+2]+16*e[iv0+1]-30*e[iv0]+16*e[iv0-1]-e[iv0-2])/(12*dV**2)
            dfe_dV3[i]=(e[iv0+2]-2*e[iv0+1]+2*e[iv0-1]-e[iv0-2])/(2*dV**3)
            dfe_dV4[i]=(e[iv0+2]-4*e[iv0+1]+6*e[iv0]-4*e[iv0-1]+e[iv0-2])/(dV**4)

            fe_V0[i]= e[iv0]
        V0=volumes0[iv0]

        tot_en = ( volumes[np.newaxis, :].T -V0) * dfe_dV1  + 0.5* ( volumes[np.newaxis, :].T -V0)**2*(dfe_dV2)
        tot_en =  tot_en+( volumes[np.newaxis, :].T -V0)**3 *dfe_dV3/6  +  ( volumes[np.newaxis, :].T -V0)**4*(dfe_dV4/24)
        tot_en  = tot_en + fe_V0 +energies[np.newaxis, :].T

        f=self.fit_forth( tstart, tstop, num ,tot_en,self.volumes)
        vol=f.min_vol

        return vol
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    @add_fig_kwargs
    def plot_vol_vs_t_4th(self, tstart=0, tstop=1000, num=101, ax=None, **kwargs):
        """
        Plot the volume as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']
        if self.scale_points=="S":
            vol_4th = self.vol_E2Vib1_forth(num=num,tstop=tstop,tstart=tstart)
            ax.plot(tmesh, vol_4th,color='b', lw=2, label="E2Vib1")
            data_to_save = np.column_stack((data_to_save,vol_4th))
            columns.append( 'E2vib1')

        if (len(self.index_list)>=2):
            vol2_4th = self.vol_EinfVib1_forth(num=num,tstop=tstop,tstart=tstart)
            ax.plot(tmesh, vol2_4th,color='gold', lw=2 ,  label=r"$E_\infty Vib1$")
            data_to_save = np.column_stack((data_to_save,vol2_4th))
            columns.append( 'Einfvib1')

        if (len(self.index_list)>=3):
            vol3_4th = self.vol_Einf_Vib2_forth(num=num,tstop=tstop,tstart=tstart)
            ax.plot(tmesh, vol3_4th,color='m', lw=2 , label=r"$E_\infty Vib2$")
            data_to_save = np.column_stack((data_to_save,vol3_4th))
            columns.append( 'Einfvib2')

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_forth( tstart, tstop, num ,tot_en,volumes)
            vol4_4th = self.vol_Einf_Vib4_forth(num=num,tstop=tstop,tstart=tstart)
            ax.plot(tmesh, vol4_4th,color='c', lw=2  ,label=r"$E_\infty Vib4$")
            ax.plot(tmesh, f0.min_vol, color='k',linestyle='dashed', lw=1.5 ,label="QHA")
            data_to_save = np.column_stack((data_to_save,vol4_4th,f0.min_vol))
            columns.append( 'Einfvib4')
            columns.append( 'QHA')

        ax.plot(0, self.volumes[self.iv0], color='g', lw=0, marker='o', ms=10,label="V0")
        ax.set_xlabel('T (K)')
        ax.set_ylabel(r'V (${\AA}^3$)')
        ax.set_xlim(tstart, tstop)
        ax.grid(True)
        ax.legend()

        return fig

    def get_thermal_expansion_coeff_4th(self, tstart=0, tstop=1000, num=101 , tref=None):
        """
        Calculates the thermal expansion coefficient as a function of temperature, using
        finite difference on the fitted values of the volume as a function of temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)

            num: int, optional Number of samples to generate. Default is 100.

        Returns: |Function1D|
        """
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
        f = self.fit_forth( tstart, tstop, num ,tot_en,self.volumes_from_phdos)
        if (tref!=None):
            ph_energies2 = self.get_vib_free_energies(tref, tref, 1)
            tot_en2 = self.energies_pdos[np.newaxis, :].T + ph_energies2
            f0 = self.fit_forth(tref, tref , 1 ,tot_en2 , self.volumes_from_phdos)

        dt = f.temp[1] - f.temp[0]
        thermo = self.get_thermodynamic_properties(tstart=tstart, tstop=tstop, num=num)
        entropy = thermo.entropy.T #* abu.e_Cb * abu.Avogadro
        df_t = np.zeros((num,self.nvols))
        df_t = - entropy
        param = np.zeros((num,4))
        param2 = np.zeros((num,3))
        d2f_t_v = np.zeros(num)
        gamma = np.zeros(num)

        #for j in range (1,num-1):
        for j in range (num):
            param[j]=np.polyfit(self.volumes_from_phdos,df_t[j] , 3)
            param2[j] = np.array([3*param[j][0],2*param[j][1],param[j][2]])

            p = np.poly1d(param2[j])
            d2f_t_v[j]= p(f.min_vol[j])

        F2D = f.F2D_V
        if (tref==None):
            #alpha= - 1/f.min_vol[1:-1] *d2f_t_v[1:-1] / F2D[1:-1]
            alpha= - 1/f.min_vol *d2f_t_v / F2D
        else :
            #alpha= - 1/f0.min_vol * d2f_t_v[1:-1] / F2D[1:-1]
            alpha= - 1/f0.min_vol * d2f_t_v / F2D

        return  alpha
    @add_fig_kwargs
    def plot_thermal_expansion_coeff_4th(self, tstart=0, tstop=1000, num=101, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        tmesh = np.linspace(tstart, tstop, num)
        thermo = self.get_thermodynamic_properties(tstart=tstart, tstop=tstop, num=num)
        entropy = thermo.entropy.T #* abu.e_Cb * abu.Avogadro
        df_t = np.zeros((num,self.nvols))
        df_t = - entropy
        volumes=self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        dV=volumes[iv0+1]-volumes[iv0]
        energy = self.energies[np.newaxis, :].T
        f=self.fit_forth( tstart=0, tstop=0, num=1 ,energy=energy,volumes=self.volumes)
        param3 = np.zeros((num,3))
        param3=np.array([12*f.param[0][0],6*f.param[0][1],2*f.param[0][2]])
        p3 = np.poly1d(param3)

        if self.scale_points=="S":
            vol_4th = self.vol_E2Vib1_forth(num=num,tstop=tstop,tstart=tstart)
            E2D = p3(self.volumes[self.iv0])
            if (tref==None):
                alpha_1 = - 1/vol_4th[:] * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D
            else :
                vol_4th_ref = self.vol_E2Vib1_forth(num=1,tstop=tref,tstart=tref)
                alpha_1 = - 1/vol_4th_ref * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D
            ax.plot(tmesh, alpha_1,color='b', lw=2, label="E2Vib1")
            data_to_save = np.column_stack((data_to_save,alpha_1))
            columns.append( 'E2vib1')

        if (len(self.index_list)>=2):
            vol2_4th = self.vol_EinfVib1_forth(num=num,tstop=tstop,tstart=tstart)
            #E2D_V = self.second_derivative_energy_v(vol2)
            E2D_V = p3(vol2_4th)
            if (tref==None):
                alpha_2 = - 1/vol2_4th[:] * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D_V[:]
            else :
                vol2_4th_ref = self.vol_EinfVib1_forth(num=1,tstop=tref,tstart=tref)
                alpha_2 = - 1/vol2_4th_ref * (df_t[:,iv1]-df_t[:,iv0])/(volumes[iv1]-volumes[iv0]) / E2D_V[:]
            ax.plot(tmesh, alpha_2,color='gold', lw=2 ,  label=r"$E_\infty Vib1$")
            data_to_save = np.column_stack((data_to_save,alpha_2))
            columns.append( 'Einfvib1')

        if (len(self.index_list)>=3):
            vol3_4th = self.vol_Einf_Vib2_forth(num=num,tstop=tstop,tstart=tstart)
            E2D_V = p3(vol3_4th)
            dfe_dV2= np.zeros( num)
            for i,e in enumerate(ph_energies.T):
                dfe_dV2[i]=(e[iv0+2]-2.0*e[iv0+1]+e[iv0])/(dV)**2

            ds_dv = (df_t[:,iv0+2]-df_t[:,iv0])/(2*dV)
            ds_dv = ds_dv+ (df_t[:,iv0+2]-2*df_t[:,iv0+1]+df_t[:,iv0])/dV**2 * (vol3_4th[:]-volumes[iv0+1])
            if (tref==None):
                alpha_3 = - 1/vol3_4th[:] * ds_dv / (E2D_V[:]+dfe_dV2[:])
            else :
                vol3_4th_ref = self.vol_Einf_Vib2_forth(num=1,tstop=tref,tstart=tref)
                alpha_3 = - 1/vol3_4th_ref * ds_dv / (E2D_V[:]+dfe_dV2[:])
            ax.plot(tmesh, alpha_3,color='m', lw=2 ,  label=r"$E_\infty Vib2$")
            data_to_save = np.column_stack((data_to_save,alpha_3))
            columns.append( 'Einfvib2')

        if (len(self.index_list)==5):
            vol4_4th = self.vol_Einf_Vib4_forth(num=num,tstop=tstop,tstart=tstart)
            E2D_V = p3(vol4_4th)

            d2fe_dV2= np.zeros( num)
            d3fe_dV3= np.zeros( num)
            d4fe_dV4= np.zeros( num)
            for i,e in enumerate(ph_energies.T):
                d2fe_dV2[i]=(-e[4]+16*e[3]-30*e[2]+16*e[1]-e[0])/(12*dV**2)
                d3fe_dV3[i]=(e[4]-2*e[3]+2*e[1]-e[0])/(2*dV**3)
                d4fe_dV4[i]=(e[4]-4*e[3]+6*e[2]-4*e[1]+e[0])/(dV**4)

            ds_dv =(-df_t[:,4]+ 8*df_t[:,3]-8*df_t[:,1]+df_t[:,0])/(12*dV)
            ds_dv = ds_dv+ (-df_t[:,4]+16*df_t[:,3]-30*df_t[:,2]+16*df_t[:,1]-df_t[:,0])/(12*dV**2) * (vol4_4th[:]-volumes[2])
            ds_dv = ds_dv+ 1.0/2.0*(df_t[:,4]-2*df_t[:,3]+2*df_t[:,1]-df_t[:,0])/(2*dV**3) * (vol4_4th[:]-volumes[2])**2
            ds_dv = ds_dv+ 1.0/6.0* (df_t[:,4]-4*df_t[:,3]+6*df_t[:,2]-4*df_t[:,1]+df_t[:,0])/(dV**4)* (vol4_4th[:]-volumes[2])**3
            D2F=E2D_V[:]+d2fe_dV2[:]+ (vol4_4th[:]-volumes[2])*d3fe_dV3[:]+0.5*(vol4_4th[:]-volumes[2])**2*d4fe_dV4[:]
            if (tref==None):
                alpha_4 = - 1/vol4_4th[:] * ds_dv / D2F
            else :
                vol4_4th_ref = self.vol_Einf_Vib4_forth(num=1,tstop=tref,tstart=tref)
                alpha_4 = - 1/vol4_4th_ref * ds_dv / D2F


            ax.plot(tmesh, alpha_4,color='c',linewidth=2 ,  label=r"$E_\infty Vib4$")

            alpha_qha  = self.get_thermal_expansion_coeff_4th(tstart, tstop, num, tref)
            ax.plot(tmesh, alpha_qha, color='k',linestyle='dashed', lw=1.5 ,label="QHA")
            data_to_save = np.column_stack((data_to_save,alpha_4,alpha_qha))
            columns.append( 'Einfvib4')
            columns.append( 'QHA')


        ax.set_xlabel(r'T (K)')
        ax.set_ylabel(r'$\alpha$ (K$^{-1}$)')
        ax.grid(True)
        ax.legend()

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
    @add_fig_kwargs
    def plot_abc_vs_t_4th(self, tstart=0, tstop=1000, num=101, lattice=None,tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']
        if self.scale_points=="S":
            vol2 = self.vol_E2Vib1_forth(num=num,tstop=tstop,tstart=tstart)
            aa2,bb2,cc2 = self.get_abc(tstart, tstop, num,vol2)
            data_to_save = np.column_stack((data_to_save,aa2,bb2,cc2))
            columns.append( 'E2vib1 (a,b,c) |            ')

        if (len(self.index_list)==2):
            vol = self.vol_EinfVib1_forth(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol = self.vol_Einf_Vib2_forth(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_forth( tstart, tstop, num ,tot_en,volumes)
            method =r"$ (E_\infty Vib4)$"
            vol = self.vol_Einf_Vib4_forth(num=num,tstop=tstop,tstart=tstart)
        aa,bb,cc = self.get_abc(tstart, tstop, num,vol)

        method_header=method+"  (a,b,c) |"
        data_to_save = np.column_stack((data_to_save,aa,bb,cc))
        columns.append( method_header)

        if (lattice==None or lattice=="a"):
            ax.plot(tmesh ,aa , color='r', lw=2,label = r"$a(V(T))$"+method,  **kwargs )
        if (lattice==None or lattice=="b"):
            ax.plot(tmesh ,bb , color='b', lw=2,label = r"$b(V(T))$"+method )
        if (lattice==None or lattice=="c"):
            ax.plot(tmesh ,cc , color='m', lw=2,label = r"$c(V(T))$"+method )

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            if (lattice==None or lattice=="a"):
                ax.plot(tmesh ,aa2 ,  linestyle='dashed' , color='r', lw=2,label = r"$a(V(T))$""E2vib1" )
            if (lattice==None or lattice=="b"):
                ax.plot(tmesh ,bb2 ,  linestyle='dashed' , color='b', lw=2,label = r"$b(V(T))$""E2vib1" )
            if (lattice==None or lattice=="c"):
                ax.plot(tmesh ,cc2 ,  linestyle='dashed' , color='m', lw=2,label = r"$c(V(T))$""E2vib1" )

        ax.set_xlabel(r'T (K)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
    @add_fig_kwargs
    def plot_angles_vs_t_4th(self, tstart=0, tstop=1000, num=101,angle=None,  tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos

        data_to_save = tmesh
        columns = ['#Tmesh']
        if self.scale_points=="S":
            vol2 = self.vol_E2Vib1_forth(num=num,tstop=tstop,tstart=tstart)
            alpha2,beta2,gamma2 = self.get_angles(tstart, tstop, num,vol2)
            data_to_save = np.column_stack((data_to_save,alpha2,beta2,gamma2))
            columns.append( 'E2vib1 (alpha,beta,gamma) |            ')

        if (len(self.index_list)==2):
            vol = self.vol_EinfVib1_forth(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol = self.vol_Einf_Vib2_forth(num=num,tstop=tstop,tstart=tstart)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_forth( tstart, tstop, num ,tot_en,volumes)
            method =r"$ (E_\infty Vib4)$"
            vol = self.vol_Einf_Vib4_forth(num=num,tstop=tstop,tstart=tstart)
        alpha,beta,gamma = self.get_angles(tstart, tstop, num,vol)

        method_header=method+"  (alpha,beta,gamma) |"
        data_to_save = np.column_stack((data_to_save,alpha,beta,gamma))
        columns.append( method_header)

        if (angle==None or angle==1):
            ax.plot(tmesh ,alpha , color='r', lw=2,label = r"$alpha(V(T))$"+method,  **kwargs )
        if (angle==None or angle==2):
            ax.plot(tmesh ,beta , color='b', lw=2,label = r"$beta(V(T))$"+method )
        if (angle==None or angle==3):
            ax.plot(tmesh ,gamma , color='m', lw=2,label = r"$gamma(V(T))$"+method )

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            if (angle==None or angle==1):
                ax.plot(tmesh ,alpha2 ,  linestyle='dashed' , color='r', lw=2,label = r"$alpha(V(T))$""E2vib1" )
            if (angle==None or angle==2):
                ax.plot(tmesh ,beta2 ,  linestyle='dashed' , color='b', lw=2,label = r"$beta(V(T))$""E2vib1" )
            if (angle==None or angle==3):
                ax.plot(tmesh ,gamma2 ,  linestyle='dashed' , color='m', lw=2,label = r"$gamma(V(T))$""E2vib1" )

        ax.set_xlabel(r'T (K)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
###################################################################################################
    @add_fig_kwargs
    def plot_thermal_expansion_coeff_abc_4th(self, tstart=0, tstop=1000, num=101, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos

        data_to_save = tmesh[1:-1]
        columns = ['#Tmesh']
        if self.scale_points=="S":
            vol2 = self.vol_E2Vib1_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol2_tref = self.vol_E2Vib1_forth(num=1,tstop=tref,tstart=tref)

        if (len(self.index_list)==2):
            vol = self.vol_EinfVib1_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref = self.vol_EinfVib1_forth(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol = self.vol_Einf_Vib2_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref = self.vol_Einf_Vib2_forth(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_forth( tstart, tstop, num ,tot_en,volumes)
            method =r"$ (E_\infty Vib4)$"
            vol = self.vol_Einf_Vib4_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref = self.vol_Einf_Vib4_forth(num=1,tstop=tref,tstart=tref)

        #alpha  = self.get_thermal_expansion_coeff(tstart, tstop, num, tref)
        tmesh = np.linspace(tstart, tstop, num)
        dt= tmesh[1] - tmesh[0]

        aa,bb,cc = self.get_abc(tstart, tstop, num,vol)
        if (tref!=None):
            aa_tref,bb_tref,cc_tref = self.get_abc(tref, tref, 1,vol_tref)

        alpha_a = np.zeros( num-2)
        alpha_b = np.zeros( num-2)
        alpha_c = np.zeros( num-2)
        if (tref==None):
            alpha_a = (aa[2:] - aa[:-2]) / (2 * dt) / aa[1:-1]
            alpha_b = (bb[2:] - bb[:-2]) / (2 * dt) / bb[1:-1]
            alpha_c = (cc[2:] - cc[:-2]) / (2 * dt) / cc[1:-1]
        else:
            alpha_a = (aa[2:] - aa[:-2]) / (2 * dt) / aa_tref
            alpha_b = (bb[2:] - bb[:-2]) / (2 * dt) / bb_tref
            alpha_c = (cc[2:] - cc[:-2]) / (2 * dt) / cc_tref

        ax.plot(tmesh[1:-1] ,alpha_a , color='r', lw=2,label = r"$\alpha_a$"+method, **kwargs)
        ax.plot(tmesh[1:-1] ,alpha_b , color='b', lw=2,label = r"$\alpha_b$"+method)
        ax.plot(tmesh[1:-1] ,alpha_c , color='m', lw=2,label = r"$\alpha_c$"+method)

        method_header=method+"  (alpha_a,alpha_b,alpha_c) |"
        data_to_save = np.column_stack((data_to_save,alpha_a,alpha_b,alpha_c))
        columns.append( method_header)

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            aa2,bb2,cc2 = self.get_abc(tstart, tstop, num,vol2)
            if (tref!=None):
                aa2_tref,bb2_tref,cc2_tref = self.get_abc(tref, tref, 1,vol2_tref)

            alpha2_a = np.zeros( num-2)
            alpha2_b = np.zeros( num-2)
            alpha2_c = np.zeros( num-2)
            if (tref==None):
                alpha2_a = (aa2[2:] - aa2[:-2]) / (2 * dt) / aa2[1:-1]
                alpha2_b = (bb2[2:] - bb2[:-2]) / (2 * dt) / bb2[1:-1]
                alpha2_c = (cc2[2:] - cc2[:-2]) / (2 * dt) / cc2[1:-1]
            else:
                alpha2_a = (aa2[2:] - aa2[:-2]) / (2 * dt) / aa2_tref
                alpha2_b = (bb2[2:] - bb2[:-2]) / (2 * dt) / bb2_tref
                alpha2_c = (cc2[2:] - cc2[:-2]) / (2 * dt) / cc2_tref

            ax.plot(tmesh[1:-1] ,alpha2_a ,  linestyle='dashed' ,  color='r', lw=2 ,label = r"$\alpha_a$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_b ,  linestyle='dashed' ,  color='b', lw=2 ,label = r"$\alpha_b$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_c ,  linestyle='dashed' ,  color='m', lw=2 ,label = r"$\alpha_c$"" (E2vib1)")
            data_to_save = np.column_stack((data_to_save,alpha2_a,alpha2_b,alpha2_c))
            columns.append( 'E2vib1 (alpha_a,alpha_b,alpha_c)   ')

        ax.set_xlabel(r'T (K)')
        ax.set_ylabel(r'$\alpha$ (K$^{-1}$)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
    @add_fig_kwargs
    def plot_thermal_expansion_coeff_angles_4th(self, tstart=0, tstop=1000, num=101, tref=None, ax=None, **kwargs):
        """
        Plots the thermal expansion coefficient as a function of the temperature.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            tref: The reference temperature (in Kelvin) used to compute the thermal expansion coefficient 1/V(tref) * dV(T)/dT.
                  (If tref is not available, it uses 1/V(T) * dV(T)/dT instead.)
            num: int, optional Number of samples to generate. Default is 100.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        iv0=self.iv0_vib
        iv1=self.iv1_vib
        volumes=self.volumes_from_phdos

        data_to_save = tmesh[1:-1]
        columns = ['#Tmesh']
        if self.scale_points=="S":
            vol2 = self.vol_E2Vib1_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol2_tref = self.vol_E2Vib1_forth(num=1,tstop=tref,tstart=tref)

        if (len(self.index_list)==2):
            vol = self.vol_EinfVib1_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref = self.vol_EinfVib1_forth(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib1)$"

        if (len(self.index_list)==3):
            vol = self.vol_Einf_Vib2_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref = self.vol_Einf_Vib2_forth(num=1,tstop=tref,tstart=tref)
            method =r"$ (E_\infty Vib2)$"

        if (len(self.index_list)==5):
            tot_en = self.energies_pdos[np.newaxis, :].T + ph_energies
            f0 = self.fit_forth( tstart, tstop, num ,tot_en,volumes)
            method =r"$ (E_\infty Vib4)$"
            vol = self.vol_Einf_Vib4_forth(num=num,tstop=tstop,tstart=tstart)
            if (tref!=None):
                vol_tref = self.vol_Einf_Vib4_forth(num=1,tstop=tref,tstart=tref)

        tmesh = np.linspace(tstart, tstop, num)
        dt= tmesh[1] - tmesh[0]

        alpha,beta,gamma = self.get_angles(tstart, tstop, num,vol)
        if (tref!=None):
            alpha_tref,beta_tref,gamma_tref = self.get_angles(tref, tref, 1,vol_tref)

        alpha_alpha = np.zeros( num-2)
        alpha_beta = np.zeros( num-2)
        alpha_gamma = np.zeros( num-2)
        if (tref==None):
            alpha_alpha = (alpha[2:] - alpha[:-2]) / (2 * dt) / alpha[1:-1]
            alpha_beta = (beta[2:] - beta[:-2]) / (2 * dt) / beta[1:-1]
            alpha_gamma = (gamma[2:] - gamma[:-2]) / (2 * dt) / gamma[1:-1]
        else:
            alpha_alpha = (alpha[2:] - alpha[:-2]) / (2 * dt) / alpha_tref
            alpha_beta = (beta[2:] - beta[:-2]) / (2 * dt) / beta_tref
            alpha_gamma = (gamma[2:] - gamma[:-2]) / (2 * dt) / gamma_tref

        ax.plot(tmesh[1:-1] ,alpha_alpha , color='r', lw=2,label = r"$\alpha_alpha$"+method, **kwargs)
        ax.plot(tmesh[1:-1] ,alpha_beta , color='b', lw=2,label = r"$\alpha_beta$"+method)
        ax.plot(tmesh[1:-1] ,alpha_gamma , color='m', lw=2,label = r"$\alpha_gamma$"+method)

        method_header=method+"  (alpha_alpha,alpha_beta,alpha_gamma) |"
        data_to_save = np.column_stack((data_to_save,alpha_alpha,alpha_beta,alpha_gamma))
        columns.append( method_header)

        if abs(abs(self.volumes[self.iv0]-volumes[iv0])-abs(volumes[iv1]-self.volumes[self.iv0]))<1e-3 :
            alpha2,beta2,gamma2 = self.get_angles(tstart, tstop, num,vol2)
            if (tref!=None):
                alpha2_tref,beta2_tref,gamma2_tref = self.get_angles(tref, tref, 1,vol2_tref)

            alpha2_alpha = np.zeros( num-2)
            alpha2_beta = np.zeros( num-2)
            alpha2_gamma = np.zeros( num-2)
            if (tref==None):
                alpha2_alpha = (alpha2[2:] - alpha2[:-2]) / (2 * dt) / alpha2[1:-1]
                alpha2_beta = (beta2[2:] - beta2[:-2]) / (2 * dt) / beta2[1:-1]
                alpha2_gamma = (gamma2[2:] - gamma2[:-2]) / (2 * dt) / gamma2[1:-1]
            else:
                alpha2_alpha = (alpha2[2:] - alpha2[:-2]) / (2 * dt) / alpha2_tref
                alpha2_beta = (beta2[2:] - beta2[:-2]) / (2 * dt) / beta2_tref
                alpha2_gamma = (gamma2[2:] - gamma2[:-2]) / (2 * dt) / gamma2_tref

            ax.plot(tmesh[1:-1] ,alpha2_alpha ,  linestyle='dashed' ,  color='r', lw=2 ,label = r"$\alpha_alpha$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_beta ,  linestyle='dashed' ,  color='b', lw=2 ,label = r"$\alpha_beta$"" (E2vib1)")
            ax.plot(tmesh[1:-1] ,alpha2_gamma ,  linestyle='dashed' ,  color='m', lw=2 ,label = r"$\alpha_gamma$"" (E2vib1)")
            data_to_save = np.column_stack((data_to_save,alpha2_alpha,alpha2_beta,alpha2_gamma))
            columns.append( 'E2vib1 (alpha_alpha,alpha_beta,alpha_gamma)   ')

        ax.set_xlabel(r'T (K)')
        ax.set_ylabel(r'$\alpha$ (K$^{-1}$)')
        ax.legend()
        ax.grid(True)

        ax.set_xlim(tstart, tstop)
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

        return fig
#*********************************************************************************************
    @classmethod
    def from_files_app(cls, gsr_paths, phdos_paths):
        """
        Creates an instance of QHA from a list of GSR files and a list of PHDOS.nc files.
        The list should have the same size and the volumes should match.

        Args:
            gsr_paths: list of paths to GSR files.
            phdos_paths: list of paths to PHDOS.nc files.

        Returns: A new instance of QHA
        """
        energies = []
        structures = []
        pressures = []
        for gp in gsr_paths:
            with GsrFile.from_file(gp) as g:
                energies.append(g.energy)
                structures.append(g.structure)
                pressures.append(g.pressure)

        #doses = [PhononDos.as_phdos(dp) for dp in phdos_paths]

        doses = []
        structures_from_phdos = []
        for path in phdos_paths:
            with PhdosFile(path) as p:
                doses.append(p.phdos)
                structures_from_phdos.append(p.structure)

       # cls._check_volumes_id(structures, structures_from_phdos)
        vols1 = [s.volume for s in structures]
        vols2 = [s.volume for s in structures_from_phdos]
        dv=np.zeros((len(vols2)-1))
        for j in range(len(vols2)-1):
            dv[j]=vols2[j+1]-vols2[j]
        tolerance = 1e-3
        if (len(vols2)!=2):
            max_difference = np.max(np.abs(dv - dv[0]))
            if max_difference > tolerance:
                raise RuntimeError("Expecting an equal volume change for structures from PDOS." )

        index_list = [i for v2 in vols2 for i, v1 in enumerate(vols1) if abs(v2 - v1) < 1e-3]
        if len(index_list) != len(vols2):
            raise RuntimeError("Expecting the ground state files for all PDOS files!")
        if len(index_list) not in (2, 3, 5):
            raise RuntimeError("Expecting just 2, 3, or 5 PDOS files in the approximation method.")

        return cls(structures,structures_from_phdos,index_list, doses, energies , pressures)

    def get_vib_free_energies(self, tstart=0, tstop=1000, num=101) -> np.ndarray:
        """
        Generates the vibrational free energy from the phonon DOS.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns: A numpy array of `num` values of the vibrational contribution to the free energy
        """
        f = np.zeros((self.nvols, num))

        for i, dos in enumerate(self.doses):
            f[i] = dos.get_free_energy(tstart, tstop, num).values
        return f

    def get_thermodynamic_properties(self, tstart=0, tstop=1000, num=101):
        """
        Generates all the thermodynamic properties corresponding to all the volumes using the phonon DOS.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 100.

        Returns:
            `namedtuple` with the following attributes for all the volumes:

                tmesh: numpy array with the list of temperatures. Shape (num).
                cv: constant-volume specific heat, in eV/K. Shape (nvols, num).
                free_energy: free energy, in eV. Shape (nvols, num).
                entropy: entropy, in eV/K. Shape (nvols, num).
                zpe: zero point energy in eV. Shape (nvols).
        """
        tmesh = np.linspace(tstart, tstop, num)
        cv = np.zeros((self.nvols, num))
        free_energy = np.zeros((self.nvols, num))
        entropy = np.zeros((self.nvols, num))
        internal_energy = np.zeros((self.nvols, num))
        zpe = np.zeros(self.nvols)

        for i, d in enumerate(self.doses):
            cv[i] = d.get_cv(tstart, tstop, num).values
            free_energy[i] = d.get_free_energy(tstart, tstop, num).values
            entropy[i] = d.get_entropy(tstart, tstop, num).values
            zpe[i] = d.zero_point_energy

        return dict2namedtuple(tmesh=tmesh, cv=cv, free_energy=free_energy, entropy=entropy, zpe=zpe)
#=======================================================================================================
    @classmethod
    def from_files_app_ddb(cls, ddb_paths, phdos_paths):
        """
        Creates an instance of QHA from a list of GSR files and a list of PHDOS.nc files.
        The list should have the same size and the volumes should match.

        Args:
            ddb_paths: list of paths to DDB files.
            phdos_paths: list of paths to PHDOS.nc files.

        Returns: A new instance of QHA
        """
        energies = []
        structures = []
        pressures = []
        for gp in ddb_paths:
            with DdbFile.from_file(gp) as g:
                energies.append(g.total_energy)
                structures.append(g.structure)
                #pressures.append(g.pressure)

        #doses = [PhononDos.as_phdos(dp) for dp in phdos_paths]

        doses = []
        structures_from_phdos = []
        for path in phdos_paths:
            with PhdosFile(path) as p:
                doses.append(p.phdos)
                structures_from_phdos.append(p.structure)

        # cls._check_volumes_id(structures, structures_from_phdos)
        vols1 = [s.volume for s in structures]
        vols2 = [s.volume for s in structures_from_phdos]
        dv=np.zeros((len(vols2)-1))
        for j in range(len(vols2)-1):
            dv[j]=vols2[j+1]-vols2[j]
        tolerance = 1e-3

        if (len(vols2)!=2):
            max_difference = np.max(np.abs(dv - dv[0]))
            if max_difference > tolerance:
                raise RuntimeError("Expecting an equal volume change for structures from PDOS." )

        index_list = [i for v2 in vols2 for i, v1 in enumerate(vols1) if abs(v2 - v1) < 1e-3]
        if len(index_list) != len(vols2):
            raise RuntimeError("Expecting the ground state files for all PDOS files!")
        if len(index_list) not in (2, 3, 5):
            raise RuntimeError("Expecting just 2, 3, or 5 PDOS files in the approximation method.")

        return cls(structures, structures_from_phdos, index_list, doses, energies, pressures)


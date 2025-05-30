# coding: utf-8
from __future__ import annotations

import numpy as np
import abipy.core.abinit_units as abu

from pymatgen.io.phonopy import get_pmg_structure
from abipy.tools.plotting import get_ax_fig_plt,add_fig_kwargs
from abipy.tools.typing import Figure
from abipy.embedding.utils_ifc import clean_structure
from abipy.lumi.utils_lumi import A_hw_help,L_hw_help,plot_emission_spectrum_help
from abipy.embedding.utils_ifc import localization_ratio


class Lineshape:
    """

    Object representing a luminescent lineshape, following a multi-phonon mode model (multiD-CCM).
    For 1D-CCM, use plot_lineshape_1D_zero_temp() function of the abipy/lumi/deltaSCF module.

    For equations, notations and formalism, please refer to:
     https://doi.org/10.1103/PhysRevB.96.125132
     https://doi.org/10.1002/adom.202100649
     https://pubs.acs.org/doi/full/10.1021/acs.chemmater.3c00537

    In the 1D-CCM, the vibronic peaks are the ones from a fictious phonon mode that connects
    the atomic relaxation between the ground state and excited state.
    Within this model, the global shape (fwhm) is well represented if the total Huang-Rhys
    factor is large enough (gaussian shaped spectrum). However, if vibronic
    peaks are present in the experimental spectrum, this model is not able to correctly
    reproduce them as it assumes a fictious phonon mode.

    In the multiD-CCM, the atomic relaxation is projected along the phonon eigenvectors of the system,
    allowing a phonon-projected decomposition of the relaxation.
    Better agreement with experimental vibronic peaks is expected.
    """

    @classmethod
    def from_phonopy_phonons(cls,E_zpl,phonopy_ph,dSCF_structure,use_forces=True,dSCF_displacements=None,dSCF_forces=None,coords_defect_dSCF=None,
                             coords_defect_phonons=None,tol=0.3):
        r"""
        Different levels of approximations for the phonons and force/displacements:
        See discussion in the supplementary informations of https://pubs.acs.org/doi/full/10.1021/acs.chemmater.3c00537, section (1).

        - size_supercell deltaSCF = size_supercell phonons (phonons of the bulk structure or phonons of defect structure).
          Use of the forces or the displacemements is allowed.

        - size_supercell dSCF < size_bulk_supercell phonons (bulk)
          Use of the forces only.

        - size_supercell dSCF < size__defect_supercell phonons (embedding)
          Use of the forces only

        The code first extracts the eigenmodes of the phonopy object.
        Then, it tries performs a structure matching between the phonopy
        structure and the dSCF structure (critical part) in order to put the displacements/forces on the right atoms.

        Args:
            E_zpl: Zero-phonon line energy in eV
            phonopy_ph: Phonopy object containing eigenfrequencies and eigenvectors
            dSCF_structure: Delta_SCF structure
            dSCF_displacements: Dispalcements \Delta R induced by the electronic transition
            dSCF_forces: Dispalcements \Delta F induced by the electronic transition
            coords_defect_dSCF: Main coordinates of the defect in defect structure, if defect complex, can be set to the
                center of mass of the complex
            tol: tolerance in Angstrom applied for the matching between the dSCF structure and phonon structure

        Returns: A lineshape object
        """
        ph_modes = phonopy_ph.get_frequencies_with_eigenvectors(q=[0, 0, 0])
        ph_freq_phonopy, ph_vec_phonopy = ph_modes

        freqs = ph_freq_phonopy * (1 / abu.eV_to_THz)   # THz to eV
        vecs = ph_vec_phonopy.transpose() #

        dSCF_structure = clean_structure(dSCF_structure,coords_defect_dSCF)

        phonon_supercell = get_pmg_structure(phonopy_ph.supercell)
        phonon_supercell = clean_structure(phonon_supercell,coords_defect_phonons)

        if not use_forces:
            forces = None
            displacements = get_displacements_on_phonon_supercell(dSCF_supercell=dSCF_structure,
                                                                phonon_supercell=phonon_supercell,
                                                                displacements_dSCF=dSCF_displacements,
                                                                tol=tol)

        if use_forces:
            displacements = None
            forces = get_forces_on_phonon_supercell(dSCF_supercell=dSCF_structure,
                                                 phonon_supercell=phonon_supercell,
                                                 forces_dSCF=dSCF_forces,
                                                 tol=tol)

        return cls(E_zpl=E_zpl,
                   ph_eigvec=vecs,
                   ph_eigfreq=freqs,
                   structure=phonon_supercell,
                   use_forces=use_forces,
                   forces=forces,
                   displacements=displacements)

    def __init__(self, E_zpl, ph_eigvec, ph_eigfreq, structure,
                 forces,displacements,use_forces):
        """
        Args:
            E_zpl: Zero-phonon line energy in eV
            ph_eigvec: phonon eigenvectors, shape: (3 * N_atoms, 3 * N_atoms)
            ph_eigfreq: phonon eigenfrequencies, shape: (3 * N_atoms), in eV
            structure: Structure object
            forces: Forces acting on the atoms in the ground state, with atomic positions of the relaxed excited state, in eV/Ang
            displacements: Atomic relaxation induced by the electronic transition, in Ang
            use_forces: True in order to use the forces, False to use the displacements
        """
        self.E_zpl = E_zpl
        self.ph_eigvec = ph_eigvec
        self.ph_eigfreq = ph_eigfreq
        self.structure = structure
        self.use_forces = use_forces
        self.forces = forces
        self.displacements = displacements

    def n_modes(self):
        """
        Number of phonon modes
        """
        return len(self.ph_eigfreq)

    def mass_list(self):
        """
        List of masses of the atoms in the structure, in amu unit
        """
        #
        amu_list = np.zeros(len(self.structure))
        for i, atom in enumerate(self.structure.species):
            amu_list[i] = atom.atomic_mass
        return (amu_list)

    def Delta_Q_nu(self):
        """
        List of Delta_Q_nu of the atoms in the structure, in SI unit
        """
        Q_nu = np.zeros(self.n_modes())

        masses = np.repeat(self.mass_list(), 3) * 1.66053892173E-27 # amu to kg
        ph_eigvector = self.ph_eigvec
        ph_eigfreq = self.ph_eigfreq * (abu.eV_s)  # eV to rad/s

        if self.use_forces:
            delta_forces = self.forces.flatten() * ((abu.eV_Ha)*(abu.Ha_J)/(1e-10)) # eV/Angstrom to Newtons (Joules/meter)
            for i in range(self.n_modes()):
                if self.ph_eigfreq[i] < 1e-5: # discard phonon freq with too low freq (avoid division by nearly 0)
                    Q_nu[i] = 0
                else:
                    Q_nu[i] = (1/ph_eigfreq[i]**2) * np.sum( delta_forces/np.sqrt(masses) * np.real(ph_eigvector[i]))
                # equation (7) of Alkauskas, A. (2014) New Journal of Physics, 16(7), 073026.
        else:
            displacements = self.displacements.flatten() * (1e-10) # in meters
            for i in range(self.n_modes()):
                Q_nu[i] = np.sum(np.sqrt(masses) * displacements * np.real(ph_eigvector[i]))
                # equation (6) of Alkauskas, A. (2014) New Journal of Physics, 16(7), 073026.`

        Q_nu[0:3] = 0 # acoustic modes are set to 0
        return Q_nu

    def S_nu(self):
        """
        Partial Huang-Rhys factors
        """
        omega = (abu.eV_s) * self.ph_eigfreq # eV to [rad/s]
        Delta_Q = self.Delta_Q_nu()
        hbar = abu.hbar_eVs*((abu.eV_Ha)*(abu.Ha_J)) # hbar in SI
        S_nu = omega * Delta_Q ** 2 / (2 * hbar)
        return (S_nu)

    def S_tot(self):
        """
        Total Huang-Rhys factor = sum of the S_nu.
        """
        return (np.sum(self.S_nu()))

    def Delta_Q(self, unit="SI"):
        """
        Total Delta_Q, "SI" or "atomic" unit.
        """
        dQ = np.sqrt(np.sum(self.Delta_Q_nu() ** 2))
        if unit == "SI":
            return (dQ)
        if unit == "atomic":
            kgm2_amuAng2 = 6.0221366516752*1e26*1e20
            dQ = dQ*np.sqrt(kgm2_amuAng2)
            return (dQ)

    def p_nu(self):
        """
        Contribution of each mode, eq. (11) of  https://doi.org/10.1002/adom.202100649
        """
        return (self.Delta_Q_nu() / self.Delta_Q()) ** 2

    def eff_freq_multiD(self):
        """
        Effective coupling frequency, eq. (13) of  https://doi.org/10.1002/adom.202100649, in eV
        """
        w = np.sqrt(np.sum(self.p_nu() * self.ph_eigfreq ** 2))
        return (w)

    def S_hbarOmega(self, broadening):
        """
        Return (Energy,spectral decomposition S(hw))

        Args:
            broadening: fwhm of the gaussian broadening in meV
        """
        n_step = 100001
        S = np.zeros(n_step)
        sigma = broadening / (1000 * 2.35482)
        freq_eV = self.ph_eigfreq
        omega = np.linspace(0, 1.1 * max(freq_eV), n_step)
        S_nu = self.S_nu()

        for k in np.arange(self.n_modes()):
            S += S_nu[k] * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((omega - freq_eV[k]) ** 2 / (2 * sigma ** 2)))
        return (omega, S)

    def localization_ratio(self):
        """
        array of the phonon mode localisation, see equation (10) of https://pubs.acs.org/doi/10.1021/acs.chemmater.3c00537
        """
        return localization_ratio(self.ph_eigvec)

    #### Generating function ####

    def A_hw(self,T, lamb=3, w=3, model='multi-D'):
        """
        Lineshape function
        Eq. (2) of https://pubs.acs.org/doi/full/10.1021/acs.chemmater.3c00537
        Returns (Energy in eV, Lineshape function )

        Args:
            T: Temperature in K
            lamb: Lorentzian broadening applied to the vibronic peaks, in meV
            w: Gaussian broadening applied to the vibronic peaks, in meV
            model: 'multi-D' for full phonon decomposition, 'one-D' for 1D-CCM PL spectrum.
        """
        S_nu = self.S_nu()
        omega_nu = self.ph_eigfreq
        eff_freq = self.eff_freq_multiD()
        E_zpl = self.E_zpl
        return A_hw_help(S_nu,omega_nu,eff_freq,E_zpl,T, lamb, w, model='multi-D')

    def L_hw(self, T=0,lamb=3, w=3, model='multi-D'):
        """
        Normalized Luminescence intensity (area under the curve = 1)
        Eq. (1) of https://pubs.acs.org/doi/full/10.1021/acs.chemmater.3c00537
        Returns (Energy in eV, Luminescence intensity)

        Args:
            T: Temperature in K
            lamb: Lorentzian broadening applied to the vibronic peaks, in meV
            w: Gaussian broadening applied to the vibronic peaks, in meV
            model: 'multi-D' for full phonon decomposition, 'one-D' for 1D-CCM PL spectrum.
        """
        E_x, A = self.A_hw(T,lamb,w,model)
        E_x, I = L_hw_help(E_x, A)
        return E_x, I

##### Plot functions ######

    @add_fig_kwargs
    def plot_spectral_function(self, broadening=1, ax=None, with_S_nu=False, with_local_ratio=False, **kwargs) -> Figure:
        """
        Plot the Huang-Rhys spectral function S_hbarOmega

        Args:
            broadening: fwhm of the gaussian broadening in meV
            with_S_nu: True to add stem lines associated to the individuals partial Huang-Rhys factors
            with_local_ratio: True to add stem lines with colored based on the mode localisation.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        S_nu = self.S_nu()
        omega_nu = self.ph_eigfreq
        S_x, S_y = self.S_hbarOmega(broadening=broadening)

        local_ratio = localization_ratio(self.ph_eigvec)

        if with_local_ratio:
            with_S_nu = False # such that the plot is fine even if both with_local_ratio and with_S_nu are set to True.
            # reorder for better plot visualisation
            idx_order = np.argsort(local_ratio)#[::-1]
            local_ratio = local_ratio[idx_order]
            omega_nu = omega_nu[idx_order]
            S_nu = S_nu[idx_order]

            import matplotlib as mpl
            cmap = mpl.colormaps["plasma"]
            #norm = mpl.colors.Normalize(vmin=np.min(local_ratio),vmax=np.max(local_ratio))
            norm = mpl.colors.LogNorm(1,vmax=np.max(local_ratio))
            sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

            log_local_ratio = np.log(local_ratio)
            local_ratio_normalized = (log_local_ratio-min(log_local_ratio))/max((log_local_ratio)-min(log_local_ratio))
            color_list = cmap(local_ratio_normalized)

            ax2 = ax.twinx()
            ax2.scatter(omega_nu,S_nu,c=color_list,norm=norm,alpha=0.6)
            ax2.vlines(x=omega_nu, ymin=0, ymax=S_nu,colors=color_list,linestyles="solid",alpha=0.6,norm=norm)
            cbar = plt.colorbar(sm,ax=ax2,location="top",shrink=0.6,label=r'$\beta_{\nu}$')
            ax2.set_ylabel(r'$S_{\nu}$')

        if with_S_nu:
            ax2 = ax.twinx()
            ax2.vlines(x=omega_nu, ymin=0, ymax=S_nu,linestyles="solid",alpha=0.6)
            ax2.scatter(omega_nu,S_nu,alpha=0.6)
            ax2.set_ylabel(r'$S_{\nu}$')

        ax.plot(S_x,S_y,**kwargs)
        ax.set_xlabel('Phonon energy (eV)')
        ax.set_ylabel(r'$S(\hbar\omega)$  (1/eV)')

        return fig

    @add_fig_kwargs
    def plot_emission_spectrum(self,unit='eV',T=0,lamb=3,w=3,max_to_one=False,ax=None,**kwargs):
        """
        Plot the Luminescence intensity

        Args:
            unit: 'eV', 'cm-1', or 'nm'
            T: Temperature in K
            lamb: Lorentzian broadening applied to the vibronic peaks, in meV
            w: Gaussian broadening applied to the vibronic peaks, in meV
            max_to_one: True if max of the curve is normalized to 1.
        """

        x_eV, y_eV = self.L_hw(T=T,lamb=lamb,w=w)
        return plot_emission_spectrum_help(x_eV,y_eV,unit,max_to_one,ax,**kwargs)


def get_forces_on_phonon_supercell(dSCF_supercell,phonon_supercell,forces_dSCF,tol):
    forces_in_supercell = np.zeros(shape=(len(phonon_supercell), 3))
    mapping = get_matching_dSCF_phonon_spcell(dSCF_supercell,phonon_supercell,tol)
    for i in range(len(mapping)):
        forces_in_supercell[mapping[i]] = forces_dSCF[i]

    return forces_in_supercell


def get_displacements_on_phonon_supercell(dSCF_supercell,phonon_supercell,displacements_dSCF,tol):
    displacements_in_supercell = np.zeros(shape=(len(phonon_supercell), 3))
    mapping = get_matching_dSCF_phonon_spcell(dSCF_supercell,phonon_supercell,tol)
    for i in range(len(mapping)):
        displacements_in_supercell[mapping[i]] = displacements_dSCF[i]

    return displacements_in_supercell


def get_matching_dSCF_phonon_spcell(dSCF_spcell,phonon_spcell,tol):
    dSCF_spcell_cart = dSCF_spcell.cart_coords
    phonon_spcell_cart = phonon_spcell.cart_coords
    # perform the matching
    mapping = []
    for i, site_1 in enumerate(dSCF_spcell):  # subset structure
        for j, site_2 in enumerate(phonon_spcell):  # superset structure
            if max(abs(dSCF_spcell_cart[i] - phonon_spcell_cart[j])) < tol:
                mapping.append(j)

    if len(mapping) == len(dSCF_spcell):
        print("Mapping between delta SCF supercell and phonon supercell succeeded. ",len(mapping),"/",len(dSCF_spcell))
    if len(mapping) != len(dSCF_spcell):
        print("Caution... Mapping between delta SCF supercell and phonon supercell did not succeed. ",len(mapping),"/",len(dSCF_spcell))

    return mapping

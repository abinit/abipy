import numpy as np
from numpy import fft
from scipy import signal
from scipy.integrate import simps

from abipy.tools.plotting import get_ax_fig_plt
from abipy.lumi import phonon_folding


class Lineshape():
    """
    define a lineshape object that can be obtained from a LumiWork + optionaly computation of phonons for a multi-D approach.
    """

    @classmethod
    def from_embedded_method(cls,E_zpl,prim_DDB,phonon_spcell_size,dSCF_structure,forces_dSCF):

        ph_freq,ph_vec,structure_phonons=phonon_folding.get_phonons_mapped_and_spcell_structure(phonon_spcell_size,prim_DDB)

        forces=phonon_folding.get_forces_on_phonon_supercell(dSCF_supercell=dSCF_structure,
                                                             phonon_supercell=structure_phonons,
                                                             forces_dSCF=forces_dSCF)

        return cls(E_zpl=E_zpl,
                   ph_eigvec=ph_vec,
                   ph_eigfreq=ph_freq,
                   structuregs=structure_phonons,
                   structureex=structure_phonons,
                   use_forces=True,
                   forces_gs=forces,
                   forces_ex=None)



    def __init__(self, E_zpl, ph_eigvec, ph_eigfreq, structuregs, structureex,
                 use_forces=False,forces_gs=None,forces_ex=None):
        """
        """
        self.E_zpl = E_zpl
        self.ph_eigvec = ph_eigvec
        self.ph_eigfreq = ph_eigfreq
        self.structuregs = structuregs
        self.structureex = structureex
        self.forces_gs = forces_gs
        self.forces_ex = forces_ex
        self.use_forces=use_forces

    def n_modes(self):
        return len(self.ph_eigfreq)

    def n_atoms(self):
        return len(self.structuregs)

    def mass_list(self):
        # return a list of masses of the atoms in the structure, in amu
        amu_list = np.zeros(len(self.structuregs))
        for i, atom in enumerate(self.structuregs.species):
            amu_list[i] = atom.atomic_mass
        return (amu_list)

    def displacements(self):
        return (self.structureex.cart_coords - self.structuregs.cart_coords)

    def delta_forces(self):
        #return(self.forces_ex-self.forces_gs)
        return self.forces_gs

    def Delta_Q_nu(self):
        # return Delta_Q_nu list
        Q_nu = np.zeros(self.n_modes())

        masses = np.repeat(self.mass_list(), 3) * 1.66053892173E-27 # in kg
        displacements = self.displacements().flatten() * (1e-10) # in meters
        ph_eigvector = self.ph_eigvec
        ph_eigfreq= self.ph_eigfreq * 1.5193e+15  # in rad/s

        if self.use_forces==True:
            delta_forces=self.delta_forces().flatten() * 1.602176565E-9 # to Newtons
            for i in range(self.n_modes()):
                if self.ph_eigfreq[i] < 1e-5 :# discard phonon freq with too low freq (avoid division by nearly 0)
                    Q_nu[i] = 0
                else:
                    Q_nu[i] = (1/ph_eigfreq[i]**2) * np.sum( delta_forces/np.sqrt(masses) * np.real(ph_eigvector[i]) )
                # equation (7) of Alkauskas, A. (2014) New Journal of Physics, 16(7), 073026.
        else:
            for i in range(self.n_modes()):
                Q_nu[i] = np.sum(np.sqrt(masses) * displacements * np.real(ph_eigvector[i]))
                # equation (6) of Alkauskas, A. (2014) New Journal of Physics, 16(7), 073026.`

        Q_nu[0:3]=0 # acoustic modes are set to 0
        return Q_nu


    def S_nu(self):
        # Huang Rhys parameter of each mode
        # transform to SI units

        # eV to [rad/s]
        omega = 1.5193e+15 * self.ph_eigfreq
        # amu^1/2.Angstrom to kg.meters
        Delta_Q = self.Delta_Q_nu() #* np.sqrt(1.66053892173e-27) * (1e-10)
        hbar = 1.054571818e-34
        S_nu = omega * Delta_Q ** 2 / (2 * hbar)
        return (S_nu)


    def Delta_Q(self):
        return (np.sqrt(np.sum(self.Delta_Q_nu() ** 2)))



    def p_nu(self):
        # Contribution of each mode
        return (self.Delta_Q_nu() / self.Delta_Q()) ** 2

    def eff_freq_multiD(self):
        w = np.sqrt(np.sum(self.p_nu() * self.ph_eigfreq ** 2))
        return (w)

    def S_hbarOmega(self, broadening):
        # return (Energy,spectral decomposition S(hw))
        # broadening is the fwhm of the gaussian broadening in meV
        n_step = 100001
        S = np.zeros(n_step)
        sigma = broadening / (1000 * 2.35482)
        freq_eV = self.ph_eigfreq
        omega = np.linspace(0, 1.1 * max(freq_eV), n_step)
        S_nu = self.S_nu()

        for k in np.arange(self.n_modes()):
            S += S_nu[k] * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((omega - freq_eV[k]) ** 2 / (2 * sigma ** 2)))
        return (omega, S)

    def S_hbarOmega_SI(self, broadening):
        # return (Energy,spectral decomposition S(hw))
        # broadening is the fwhm of the gaussian broadening in meV
        n_step = 100001
        S = np.zeros(n_step)
        sigma = (broadening / (1000 * 2.35482)) * 1.5193e+15
        freq_eV = self.ph_eigfreq * 1.5193e+15  # in SI
        omega = np.linspace(0, 1.1 * max(freq_eV), n_step)
        S_nu = self.S_nu()

        for k in np.arange(self.n_modes()):
            S += S_nu[k] * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((omega - freq_eV[k]) ** 2 / (2 * sigma ** 2)))
        return (omega, S)

    #### Generating function ####

    def Bose_einstein(self, T, freq):
        # freq in eV
        k_b = 8.6173303e-5  # in eV/K

        if T == 0:
            n = 0
        else:
            n = 1 / (np.exp(freq / (k_b * T)) - 1)

        return n

    def G_t(self, T, S_nu, omega_nu):
        n_step = 10001
        t = np.linspace(-1e-11, +1e-11, n_step)  # time in the fourier domain

        g_plus = np.zeros(n_step, dtype=complex)
        g_minus = np.zeros(n_step, dtype=complex)
        freq = omega_nu
        freq_SI = freq * 1.5193e+15  # in SI rad/sec

        for i in range(len(S_nu)):
            g_plus += self.Bose_einstein(T, freq[i]) * S_nu[i] * np.exp(+1j * freq_SI[i] * t)
            g_minus += (self.Bose_einstein(T, freq[i]) + 1) * S_nu[i] * np.exp(-1j * freq_SI[i] * t)

        index_0 = int((len(t) - 1) / 2)
        g = g_plus[index_0] + g_minus[index_0]

        G_t = np.exp(g_plus + g_minus - g)

        return (t, G_t)

    def A_hw(self, T, lamb=5, w=1, model='multi-D'):

        # lamb is the homogeneous broadening in meV
        # w is the inhomogeneous broadening in meV
        # model : 'multi-D' or 'one-D'
        if model == 'multi-D':
            t, G_t = self.G_t(T, S_nu=self.S_nu(), omega_nu=self.ph_eigfreq)

        elif model == 'one-D':
            t, G_t = self.G_t(T, S_nu=np.array([np.sum(self.S_nu())]), omega_nu=np.array([self.eff_freq_multiD()]))

        n_step = len(t)
        Lambda = 2.4180e+11 * lamb  # For lorentzian :
        delta_t = t[-1] - t[0]

        fourier = fft.fft(G_t * np.exp(-Lambda * np.abs(t)))
        freq = 2 * np.pi * n_step * fft.fftfreq(G_t.size) / delta_t
        # freq change
        freq_2 = np.zeros(n_step)
        fourier_2 = np.zeros(n_step, dtype=complex)

        freq_2[0:n_step // 2] = freq[n_step // 2 + 1:]
        freq_2[n_step // 2] = freq[0]
        freq_2[n_step // 2 + 1:] = freq[1:n_step // 2 + 1]

        fourier_2[0:n_step // 2] = fourier[n_step // 2 + 1:]
        fourier_2[n_step // 2] = fourier[0]
        fourier_2[n_step // 2 + 1:] = fourier[1:n_step // 2 + 1]

        hbar_eV = 6.582119570e-16  # in eV*s
        En = hbar_eV * freq_2
        E_x = En + self.E_zpl

        sigma = w / (2.35482 * 1000)
        gaussian = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((En) ** 2 / (2 * (sigma) ** 2)))
        A_conv = signal.fftconvolve(np.abs(fourier_2), gaussian, mode='same')

        return (E_x, A_conv)

    def L_hw(self, T=0, lamb=5, w=1, model='multi-D'):
        # Normalization such that area under the curve is one
        E_x, A = self.A_hw(T, lamb, w, model)
        C = 1 / (simps(A * E_x ** 3, E_x))
        I = C * A * E_x ** 3  # intensity prop to energy CUBE
        return (E_x, I)


##### Plot functions ######

    def plot_spectral_function(self,broadening=1,ax=None,with_S_nu=False,color='k',legend=None):

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        S_nu=self.S_nu()
        omega_nu=self.ph_eigfreq
        S_x,S_y=self.S_hbarOmega(broadening=broadening)

        ax.plot(S_x,S_y,color=color,label=legend)
        ax.legend()
        ax.set_xlabel('Phonon energy (eV)')
        ax.set_ylabel(r'$S(\hbar\omega)$  (1/eV)')

        if with_S_nu==True:
            ax2=ax.twinx()
            markerline, stemlines, baseline = ax2.stem(omega_nu,S_nu,markerfmt='o',basefmt="k")
            plt.setp(markerline,'color',color)
            plt.setp(stemlines, 'color', plt.getp(markerline, 'color'))
            plt.setp(stemlines, 'linestyle', 'dotted')
            ax2.set_ylabel(r'$S_{\nu}$')

        return fig


    def plot_emission_spectrum(self,unit='eV',T=0,lamb=3,w=3,ax=None):

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        x_eV,y_eV=self.L_hw(T=T,lamb=lamb,w=w) # in eV

        x_cm=x_eV*8065.73
        y_cm=y_eV/8065.73

        x_nm=1239.84193/x_eV
        y_nm=y_eV*(x_eV**2/1239.84193)


        if unit=='eV':
            ax.plot(x_eV,y_eV)
            ax.set_xlim(2,3.2)
            ax.set_xlabel('Photon energy (eV)')
            ax.set_ylabel(r'$L(\hbar\omega)$  (1/eV)')


        elif unit=='cm-1':
            ax.plot(x_cm,y_cm)
            ax.set_xlim(2*8065.73,3.2*8065.73)
            ax.set_xlabel(r'Photon energy ($cm^{-1}$)')
            ax.set_ylabel(r'$L(\hbar\omega)$  (1/$cm^{-1}$)')

        elif unit=='nm':
            ax.plot(x_nm,y_nm)
            ax.set_xlim(350,700)
            ax.set_xlabel(r'Photon wavelength (nm))')
            ax.set_ylabel(r'Intensity (a.u.)')





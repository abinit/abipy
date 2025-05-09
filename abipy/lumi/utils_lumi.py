# help script containing functions implementing the T-dep generating function approach, to be used in both 1D and multi-D case
# and plotting fcts common to 1D and multi-D case
from __future__ import annotations

import numpy as np
import abipy.core.abinit_units as abu

from numpy import fft
from scipy import signal
try:
    from scipy.integrate import simpson as simps
except ImportError:
    from scipy.integrate import simps
from abipy.tools.plotting import get_ax_fig_plt, add_fig_kwargs #, get_axarray_fig_plt


#### Generating function ####

def Bose_einstein(T, freq):
    """
    Bose Einstein average occupation number

    Args:
        T: Temperature in K
        freq: frequency in eV
    """
    k_b = abu.kb_eVK  # in eV/K

    if T == 0:
        n = 0
    else:
        n = 1 / (np.exp(freq / (k_b * T)) - 1)

    return n


def get_G_t(T, S_nu, omega_nu):
    """
    Generation function
    Eq. (3) of https://pubs.acs.org/doi/full/10.1021/acs.chemmater.3c00537

    Args:
        T: Temperature in K
        S_nu: Parial Huang-Rhys factors
        omega_nu: phonon frequencies
    """
    n_step = 100001
    t = np.linspace(-1e-11, +1e-11, n_step)  # time in the fourier domain

    freq = np.array(omega_nu)
    freq_SI = freq * (abu.eV_s) # in SI rad/sec

    S = np.zeros(n_step, dtype=complex)
    C_plus = np.zeros(n_step, dtype=complex)
    C_minus = np.zeros(n_step, dtype=complex)

    for i in range(len(S_nu)):
        S += S_nu[i] * np.exp(-1j * freq_SI[i] * t)
        C_plus += Bose_einstein(T, freq[i]) * S_nu[i] * np.exp(+1j * freq_SI[i] * t)
        C_minus += Bose_einstein(T, freq[i]) * S_nu[i] * np.exp(-1j * freq_SI[i] * t)

    index_0 = int((len(t) - 1) / 2)
    C_0 = 2*C_plus[index_0]
    S_0 = S[index_0]

    G_t = np.exp(S-S_0+C_plus+C_minus-2*C_0)

    return t, G_t


def A_hw_help(S_nu,omega_nu,eff_freq,E_zpl,T, lamb, w, model='multi-D'):
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
    if model == 'multi-D':
        t, G_t = get_G_t(T, S_nu, omega_nu)

    elif model == 'one-D':
        t, G_t = get_G_t(T, S_nu=np.array([np.sum(S_nu)]), omega_nu=np.array(eff_freq))#np.array([self.eff_freq_multiD()]))

    n_step = len(t)
    Lambda = abu.eV_s*0.001/(np.pi*2) * lamb  # meV to Hz
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

    hbar_eV = abu.hbar_eVs  # in eV*s
    En = hbar_eV * freq_2
    E_x = En + E_zpl

    sigma = w / (2.35482 * 1000)
    gaussian = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((En) ** 2 / (2 * (sigma) ** 2)))
    A_conv = signal.fftconvolve(np.abs(fourier_2), gaussian, mode='same')

    return (E_x, A_conv)


def L_hw_help(E_x, A):
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
    C = 1 / (simps(y=A * E_x ** 3, x=E_x))
    I = C * A * E_x ** 3  # intensity prop to energy CUBE
    return (E_x, I)


#### Plotting functions ####


@add_fig_kwargs
def plot_emission_spectrum_help(x_eV, y_eV, unit, max_to_one, ax, **kwargs):
    """
    Plot the Luminescence intensity, computed from the generating fct approach.

    Args:
        unit: 'eV', 'cm-1', or 'nm'
        T: Temperature in K
        lamb: Lorentzian broadening applied to the vibronic peaks, in meV
        w: Gaussian broadening applied to the vibronic peaks, in meV
    """

    ax, fig, plt = get_ax_fig_plt(ax=ax)

    #x_eV,y_eV=self.L_hw(T=T,lamb=lamb,w=w) # in eV

    x_cm = x_eV*8065.73
    y_cm = y_eV/8065.73

    x_nm = 1239.84193/x_eV
    y_nm = y_eV*(x_eV**2/1239.84193)

    if max_to_one:
        y_eV = y_eV/max(y_eV)
        y_cm = y_cm/max(y_cm)
        y_nm = y_nm/max(y_nm)

    if unit == 'eV':
        ax.plot(x_eV,y_eV,**kwargs)
        ax.set_xlabel('Photon energy (eV)')
        ax.set_ylabel(r'$L(\hbar\omega)$  (1/eV)')

    elif unit == 'cm-1':
        ax.plot(x_cm,y_cm,**kwargs)
        ax.set_xlabel(r'Photon energy ($cm^{-1}$)')
        ax.set_ylabel(r'$L(\hbar\omega)$  (1/$cm^{-1}$)')

    elif unit == 'nm':
        ax.plot(x_nm,y_nm,**kwargs)
        ax.set_xlabel(r'Photon wavelength (nm))')
        ax.set_ylabel(r'Intensity (a.u.)')

    else:
        raise ValueError(f"Invalid {unit=}, must be 'eV', 'cm-1', or 'nm'")

    return fig

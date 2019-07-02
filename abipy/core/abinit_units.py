# coding: utf-8
"""
This module defines constants and conversion factors matching those present in abinit that can be used when
it is important to preserve consistency with the results produced by abinit.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np


# taken from abinit/10_defs/defs_basis.F90
# 1 Bohr, in Angstrom
Bohr_Ang = 0.52917720859
# 1 Angstrom in Bohr
Ang_Bohr = 1.0/Bohr_Ang
# 1 Hartree, in cm^-1
Ha_cmm1 = 219474.6313705
# 1 Hartree, in eV
Ha_eV = 27.21138386
Ha_to_eV = Ha_eV
# 1 eV in Hartree
eV_Ha = 1. / Ha_eV
# 1 eV in Rydberg
eV_Ry = 2 * eV_Ha
# 1 Hartree, in meV
Ha_meV = Ha_eV * 1000
# 1 Hartree, in Kelvin
Ha_K = 315774.65
# 1 eV, in Kelvin
eV_to_K = eV_Ha * Ha_K
# 1 Hartree, in THz
Ha_THz = 6579.683920722
# 1 Hartree, in s
Ha_s = Ha_THz * 1e12 * 2 * np.pi 
# 1 eV, in THz
eV_to_THz = eV_Ha * Ha_THz
# 1 eV, in cm-1
eV_to_cm1 = 8065.5440044136285
# 1 Hartree, in J
Ha_J = 4.35974394e-18
# Planck constant and h/2pi in eV x s
h_eVs = 4.13566766225e-15
hbar_eVs = h_eVs / (2 * np.pi)
# minus the electron charge, in Coulomb
e_Cb = 1.602176487e-19
# Boltzmann constant in eV/K and Ha/K
kb_eVK = 8.617343e-5
kb_HaK = kb_eVK / Ha_eV
# 1 atomic mass unit, in electronic mass
amu_emass = 1.660538782e-27 / 9.10938215e-31
# 1 Ha/Bohr^3, in GPa
HaBohr3_GPa = Ha_eV / Bohr_Ang**3 * e_Cb * 1.0e+21
# 1 eV/A^3 to GPa
eVA3_GPa = 160.21766208
# 1 eV in seconds
eV_s = eV_to_THz*1e12 * 2*np.pi
# conversion factor for velocity between atomic units and SI
velocity_at_to_si = 2.1876912633e6

# per mole
Avogadro = 6.02214179e23
# 1 Ohm.cm in atomic units
Ohmcm = 2 * np.pi * Ha_THz * 10 / 9
eps0 = 1 / (4*np.pi*0.0000001*299792458.0**2)
AmuBohr2_Cm2 = e_Cb * 1.0e20 / (Bohr_Ang * Bohr_Ang)
# Inverse of fine structure constant
InvFineStruct = 137.035999679
# speed of light in SI
Sp_Lt_SI = 2.99792458e8
# speed of light in atomic units
Sp_Lt = Sp_Lt_SI / 2.1876912633e6
# Atomic unit of time, in seconds
Time_Sec = 2.418884326505e-17
# Tesla in a.u.
BField_Tesla = 4.254383e-6
# Debye unit in a.u.
dipole_moment_debye = 0.393430307


def phfactor_ev2units(units):
    """
    Return conversion factor eV --> units for phonons (case-insensitive)
    """
    d = {"ev": 1, "mev": 1000, "ha": eV_Ha,
         "cm-1": eV_to_cm1, 'cm^-1': eV_to_cm1,
         "thz": eV_to_THz,
         }
    try:
        return d[units.lower().strip()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


def phunit_tag(units):
    """
    Return latex string from ``units`` (used for phonons)
    """
    d = {"ev": "(eV)", "mev": "(meV)", "ha": '(Ha)',
         "cm-1": "(cm$^{-1}$)", 'cm^-1': "(cm$^{-1}$)", "thz": '(Thz)',
         }
    try:
        return d[units.lower().strip()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


def wlabel_from_units(units):
    """
    Return latex string for phonon frequencies in ``units``.
    """
    d = {'ev': 'Energy (eV)', 'mev': 'Energy (meV)', 'ha': 'Energy (Ha)',
        'cm-1': r'Frequency (cm$^{-1}$)',
        'cm^-1': r'Frequency (cm$^{-1}$)',
        'thz': r'Frequency (Thz)',
    }
    try:
        return d[units.lower().strip()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


def phdos_label_from_units(units):
    """
    Return latex string for phonon DOS values in ``units``.
    """
    d = {"ev": "(states/eV)", "mev": "(states/meV)", "ha": '(states/Ha)',
         "cm-1": "(states/cm$^{-1}$)", 'cm^-1': "(states/cm$^{-1}$)",
         "thz": '(states/Thz)',
         }
    try:
        return d[units.lower().strip()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


def s2itup(comp):
    """
    Convert string in the form ``xx``, ``xyz`` into tuple of two (three) indices
    that can be used to slice susceptibility tensors (numpy array).

    >>> assert s2itup("yy") == (1, 1)
    >>> assert s2itup("xyz") == (0, 1, 2)
    """
    d = {"x": 0, "y": 1, "z": 2}
    comp = str(comp).strip()
    if len(comp) == 2:
        return d[comp[0]], d[comp[1]]
    elif len(comp) == 3:
        return d[comp[0]], d[comp[1]], d[comp[2]]
    else:
        raise ValueError("Expecting component in the form `xy` or `xyz` but got `%s`" % comp)


def itup2s(t):
    """
    Convert tuple of 2 (3) integers into string in the form ``xx`` (``xyz``).
    Assume C-indexing e.g. 0 --> x

    >>> assert itup2s((0, 1)) == "xy"
    >>> assert itup2s((0, 1, 2)) == "xyz"
    """
    if not isinstance(t, tuple) and len(t) not in (2, 3):
        raise TypeError("Expecting tuple of len 2 or 3, got %s" % str(t))
    d = {0: "x", 1: "y", 2: "z"}
    return "".join(d[i] for i in t)


# Use same tolerances as Abinit Fortran version in m_occ.F90.
_maxFDarg = 500.0

@np.vectorize
def occ_fd(ee, kT, mu):
    r"""
    Fermi-Dirac statistic: 1 / [(exp((e - mu)/ KT) + 1]
    Note that occ_fs in [0, 1] so the spin factor is not included,
    unlike the occupations stored in ebands%occ.

    Args:
        ee: Single particle energy in eV
        kT: Value of K_Boltzmann x T in eV.
        mu: Chemical potential in eV.
    """
    ee_mu = ee - mu

    # 1 kelvin [K] = 3.16680853419133E-06 Hartree
    if kT > 1e-6 * Ha_eV:
        arg = ee_mu / kT
        if arg > _maxFDarg:
            occ_fd = 0.0
        elif arg < - _maxFDarg:
            occ_fd = 1.0
        else:
            occ_fd = 1.0 / (np.exp(arg) + 1.0)
    else:
        # Heaviside
        if ee_mu > 0.0:
            occ_fd = 0.0
        elif ee_mu < 0.0:
            occ_fd = 1.0
        else:
            occ_fd = 0.5

    return occ_fd


@np.vectorize
def occ_be(ee, kT, mu=0.0):
    r"""
    Bose-Einstein statistic: 1 / [(exp((e - mu)/ KT) - 1]

    Args:
      ee: Single particle energy in eV.
      kT: Value of K_Boltzmann x T in eV.
      mu: Chemical potential in eV (usually zero)
    """
    ee_mu = ee - mu

    # 1 kelvin [K] = 3.16680853419133E-06 Hartree
    if kT > 1e-12 * Ha_eV:
        arg = ee_mu / kT
        if arg > 1e-12 and arg < 600.0:
            occ_be = 1.0 / (np.exp(arg) - 1.0)
        else:
            occ_be = 0.0
    else:
        # No condensate for T --> 0
        occ_be = 0.0

    return occ_be

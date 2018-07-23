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
# 1 Hartree, in cm^-1
Ha_cmm1 = 219474.6313705
# 1 Hartree, in eV
Ha_eV = 27.21138386
# 1 eV in Hartree
eV_Ha = 1. / Ha_eV
# 1 Hartree, in meV
Ha_meV = Ha_eV * 1000
# 1 Hartree, in Kelvin
Ha_K = 315774.65
# 1 eV, in Kelvin
eV_to_K = eV_Ha * Ha_K
# 1 Hartree, in THz
Ha_THz = 6579.683920722
# 1 eV, in THz
eV_to_THz = eV_Ha * Ha_THz
# 1 eV, in cm-1
eV_to_cm1 = 8065.5440044136285
# 1 Hartree, in J
Ha_J = 4.35974394e-18
# minus the electron charge, in Coulomb
e_Cb = 1.602176487e-19
# Boltzmann constant in eV/K and Ha/K
kb_eVK = 8.617343e-5
kb_HaK = kb_eVK / Ha_eV
# 1 atomic mass unit, in electronic mass
amu_emass = 1.660538782e-27 / 9.10938215e-31
# 1 Ha/Bohr^3, in GPa
HaBohr3_GPa = Ha_eV / Bohr_Ang**3 * e_Cb * 1.0e+21
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

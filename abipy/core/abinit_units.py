# coding: utf-8
"""
This module defines constants and conversion factors matching those present in abinit that can be used when
it is important to preserve consistency with the results produced by abinit.
"""
from __future__ import print_function, division, unicode_literals, absolute_import
import numpy as np


# taken from abinit/10_defs/defs_basis.F90
# 1 Bohr, in Angstrom
Bohr_Ang=0.52917720859
# 1 Hartree, in cm^-1
Ha_cmm1=219474.6313705
# 1 Hartree, in eV
Ha_eV=27.21138386
# 1 eV in Hartree
eV_Ha=1./Ha_eV
# 1 Hartree, in meV
Ha_meV=Ha_eV*1000
# 1Hartree, in Kelvin
Ha_K=315774.65
# 1 Hartree, in THz
Ha_THz=6579.683920722
# 1 eV, in THz
eV_to_THz=eV_Ha*Ha_THz
#1 Hartree, in J
Ha_J=4.35974394e-18
# minus the electron charge, in Coulomb
e_Cb=1.602176487e-19
# Boltzmann constant in eV/K and Ha/K
kb_eVK=8.617343e-5
kb_HaK = kb_eVK / Ha_eV
# 1 atomic mass unit, in electronic mass
amu_emass=1.660538782e-27/9.10938215e-31
# 1 Ha/Bohr^3, in GPa
HaBohr3_GPa=Ha_eV/Bohr_Ang**3*e_Cb*1.0e+21
# per mole
Avogadro=6.02214179e23
# 1 Ohm.cm in atomic units
Ohmcm=2*np.pi*Ha_THz*10/9
eps0=1/(4*np.pi*0.0000001*299792458.0**2)
AmuBohr2_Cm2=e_Cb*1.0e20/(Bohr_Ang*Bohr_Ang)
# Inverse of fine structure constant
InvFineStruct=137.035999679
# speed of light in SI
Sp_Lt_SI=2.99792458e8
# speed of light in atomic units
Sp_Lt=Sp_Lt_SI/2.1876912633e6
#  Atomic unit of time, in seconds
Time_Sec=2.418884326505e-17
# Tesla in a.u.
BField_Tesla=4.254383e-6
# Debye unit in a.u.
dipole_moment_debye=0.393430307
#!/usr/bin/env python
r"""
Thermodinamic properties
========================

This example shows how to compute and plot thermodinamic properties within
the harmonic approximation using the phonon DOS produced by anaddb.
"""
from __future__ import print_function

from abipy.abilab import abiopen
import abipy.data as abidata

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
ncfile = abiopen(abidata.ref_file("trf2_5.out_PHDOS.nc"))
phdos = ncfile.phdos

# Print crystalline structure and zero-point energy.
print(ncfile.structure)
zpe = phdos.zero_point_energy
print("Zero point energy:", zpe, zpe.to("J"), zpe.to("Ha"))

# Compute free energy from 2 to 300 K (20 points)
# By default, energies are is eV and thermodynamic quantities are given
# on a per-unit-cell basis.
f = phdos.get_free_energy(tstart=2, tstop=300, num=20)
#f.plot()

# Plot U, F, S, Cv as a function of T.
# Use J/mol units, results are divided by formula_units.
phdos.plot_harmonic_thermo(units="Jmol", formula_units=1)

ncfile.close()

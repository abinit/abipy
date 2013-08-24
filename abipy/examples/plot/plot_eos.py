#!/usr/bin/env python
#
# This example shows how to compute the equation of state by 
# fitting the total energy as function of the unit cell volume.
from abipy.abilab import EOS

# Extract volumes and energies from the output files of the calculation.
# Here we use hardcoded values.
volumes = [13.72, 14.83, 16.0, 17.23, 18.52]
energies = [-56.29, -56.41, -56.46, -56.46, -56.42]

# For the list of available models, see EOS.MODELS
eos = EOS(eos_name='murnaghan')

# Note that eos.fit expects lengths in Angstrom, energies are in eV.
# To specify different units use len_units and ene_units 
fit = eos.fit(volumes, energies, len_units="Bohr", ene_units="Ha")

print(fit)
fit.plot()

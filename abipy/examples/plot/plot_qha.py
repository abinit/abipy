#!/usr/bin/env python
r"""
quasi-harmonic approximation analysis
=====================================

This example shows how to plot a phonon band structure
with and without imposing the acoustic sum rule at the Gamma point.
"""
import os
import abipy.data as abidata

from abipy.dfpt.phonons import PhononBands
from abipy.dfpt.qha import QHA

# We use a list of GSR.nc and PHDOS.nc files corresponding to different isotropic strains.
# These files are shipped with AbiPy so that we don't need to run calculations from scratch.
strains = [-4, -2, 0, 2, 4, 6]
dirpath = os.path.join(abidata.dirpath, "refs", "si_qha")
gsr_paths = [os.path.join(dirpath, "mp-149_{:+d}_GSR.nc".format(s)) for s in strains]
dos_paths = [os.path.join(dirpath, "mp-149_{:+d}_PHDOS.nc".format(s)) for s in strains]


# Initialize object from files.
qha = QHA.from_files(gsr_paths, dos_paths)
#qha.set_eos("vinet")

qha.plot_energies(title="Energies as a function of volume for different T")

qha.plot_thermal_expansion_coeff(title="Thermal expansion coefficient as a function of T")

qha.plot_vol_vs_t(title="Volume as a function of T")

# fake temperatures to test the plotting function.
phbs_list = [PhononBands.from_file(os.path.join(dirpath, "mp-149_{:+d}_PHBST.nc".format(s))) for s in
             strains[2:4]]

qha.plot_phbs(phbs_list, temperatures=[10, 20], title="Phonon band structures with a color depending on T")

# Build a Phonopy QHA object.
qha_phonopy = qha.get_phonopy_qha(tstop=500, num=11)
qha_phonopy.run()
qha_phonopy.plot_bulk_modulus_temperature().show()

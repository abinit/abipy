#!/usr/bin/env python
r"""
QHA_2D method
=============

This example shows how to use the GSR.nc and PHDOS.nc files computed with different volumes
to compute thermodynamic properties within the quasi-harmonic approximation.
"""
import os
import abipy.data as abidata

from abipy.dfpt.qha_2D import QHA_2D

strains_a = [ 995,1000, 1005, 1010, 1015  ]
strains_c = [ 995,1000, 1005, 1010, 1015  ]
strains_a1 = [ 1000, 1005, 1010 ]
strains_c1 = [ 1000, 1005, 1010 ]

# Root points to the directory in the git submodule with the output results.
root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_approximation")

#gsr_paths = [[f"scale_{s1}_{s3}/out_GSR.nc" for s3 in strains_c] for s1 in strains_a]
gsr_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_GSR_DDB") for s3 in strains_c] for s1 in strains_a]
phdos_paths = [[os.path.join(root, f"scale_{s1}_{s3}/out_PHDOS.nc") for s3 in strains_c1] for s1 in strains_a1]

qha = QHA_2D.from_files(gsr_paths, phdos_paths, gsr_file="DDB")

#%%
qha.plot_energies()

#%%
qha.plot_free_energies(tstop=500, tstart=0, num=6)

#%%
qha.plot_thermal_expansion(tstop=1000, tstart=0, num=101)

#%%
qha.plot_lattice(tstop=1000, tstart=0, num=101)

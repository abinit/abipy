#!/usr/bin/env python
r"""
Quasi-harmonic approximation with v-ZSISA
=========================================

This example shows how to use the GSR.nc and PHDOS.nc files computed with different volumes
to compute thermodynamic properties within the v-ZSISA method.
"""
import os
import abipy.data as abidata

from abipy.dfpt.vzsisa import Vzsisa

# Root points to the directory in the git submodule with the output results.
root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "Si_v_ZSISA_approximation")

strains = [96, 98, 100, 102, 104, 106]
strains2 = [98, 100, 102, 104, 106] # EinfVib4(D)
#strains2 = [96, 98, 100, 102, 104] # EinfVib4(S)
#strains2 = [100, 102, 104] # EinfVib2(D)

gsr_paths = [os.path.join(root, "scale_{:d}_GSR.nc".format(s)) for s in strains]
ddb_paths = [os.path.join(root, "scale_{:d}_GSR_DDB".format(s)) for s in strains]
phdos_paths = [os.path.join(root, "scale_{:d}_PHDOS.nc".format(s)) for s in strains2]

qha = Vzsisa.from_ddb_phdos_Files(ddb_paths, phdos_paths)
tstart, tstop = 0, 800

#%%
# Plot Energies as a function of volume for different T
qha.plot_energies(tstop=tstop, tstart=tstart, num=11)

#%%
# Plot Volume as a function of T
qha.plot_vol_vs_t(tstop=tstop, tstart=tstart, num=101)

#%%
# Plot Lattice as a function of T
qha.plot_abc_vs_t(tstop=tstop, tstart=tstart, num=101)

#%%
# Plot Lattice as a function of T")
qha.plot_abc_vs_t(tstop=tstop, tstart=tstart, num=101, lattice="b")

#%%
# Plot Volumetric thermal expansion coefficient as a function of T
qha.plot_thermal_expansion_coeff(tstop=tstop, tstart=tstart ,num=101)

#%%
# Plot Thermal expansion coefficient as a function of T
qha.plot_thermal_expansion_coeff_abc(tstop=tstop, tstart=tstart ,num=101)

#%%
# Plot Angles as a function of T
qha.plot_angles_vs_t(tstop=tstop, tstart=tstart, num=101)

#%%
#
# Plot Volume as a function of T. 4th order polinomial
qha.plot_vol_vs_t_4th(tstop=tstop, tstart=tstart, num=101)

#%%
# Plot Lattice as a function of T. 4th order polinomial
qha.plot_abc_vs_t_4th(tstop=tstop, tstart=tstart, num=101, lattice="a")

#%%
# Plot Lattice as a function of T. 4th order polinomial
qha.plot_abc_vs_t_4th(tstop=tstop, tstart=tstart)

#%%
# Plot Volumetric thermal expansion coefficient as a function of T
qha.plot_thermal_expansion_coeff_4th(tref=293)

#%%
# Plot Thermal expansion coefficient as a function of T
qha.plot_thermal_expansion_coeff_abc_4th(tstop=tstop, tstart=tstart ,num=101, tref=293)

#%%
# Plot Angles as a function of T.
qha.plot_angles_vs_t_4th(tstop=tstop, tstart=tstart, num=101, angle=3)


#%%
# Create plotter to plot all the phonon DOS.
phdos_plotter = qha.get_phdos_plotter()
phdos_plotter.combiplot()

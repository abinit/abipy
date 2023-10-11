#!/usr/bin/env python
r"""
Luminescent properties (1D-CCM)
===============================

This example shows how to post-process the results of a LumiWork
following a one-dimensional configuration coordinate model (1D-CCM).
Based on NV- center in diamond (64 atoms supercell)
"""

#%%
# Read the 4 points netcdf file produced by a LumiWork

from abipy.lumi.deltaSCF import DeltaSCF
import abipy.data as abidata

NV_center = DeltaSCF.from_four_points_file([abidata.ref_file("relaxed_gs_out_GSR.nc"),
                                            abidata.ref_file("unrelaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("relaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("unrelaxed_gs_out_GSR.nc")])

# %%
# To draw the configuration coordinates diagram.

NV_center.draw_displaced_parabolas(scale_eff_freq=4);

#%%
# To plot the luminescence lineshape following the one effective phonon mode model at 0K

NV_center.plot_lineshape_1D_zero_temp(energy_range=[0.5,3],max_m=20,phonon_width=0.02,
               with_omega_cube='True',
               normalized='Sum');

#%%
# To get a panda dataframe with the main restults :
NV_center.get_dataframe()

#%%
# To plot the magnitude of the displacements/forces, from the N atom.

NV_center.plot_delta_R_distance(defect_symbol="N");
NV_center.plot_delta_F_distance(defect_symbol="N");

#%%
# To visualise displacements in 3D

NV_center.displacements_visu();

#%%
# To plot the four band structures associated to each point
# nscf_files are obtained in a LumiWork by activating the BS computations

nscf_files=[abidata.ref_file("relaxed_gs_nscf_GSR.nc"),
            abidata.ref_file("unrelaxed_ex_nscf_GSR.nc"),
            abidata.ref_file("relaxed_ex_nscf_GSR.nc"),
            abidata.ref_file("unrelaxed_gs_nscf_GSR.nc")]
NV_center.plot_four_BandStructures(nscf_files, ylims=[-4,4]);

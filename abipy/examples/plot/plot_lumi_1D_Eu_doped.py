#!/usr/bin/env python
r"""
Luminescent properties (1D-CCM)
===============================

This example shows how to post-process the results of a LumiWork
following a one-dimensional configuration coordinate model (1D-CCM).
Based on a Eu-doped phosphor SrAlLi3N4:Eu presenting two non-equivalent sites for Eu
See example/flows/run_lumi_Eu_doped_SLA.py
"""

#%%
# Read the 4 points netcdf file produced by the two LumiWork

from abipy.lumi.deltaSCF import DeltaSCF
import abipy.data as abidata
import pandas as pd

SLA_site_1 = DeltaSCF.from_four_points_file([abidata.ref_file("site_1_relaxed_gs_out_GSR.nc"),
                                            abidata.ref_file("site_1_unrelaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("site_1_relaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("site_1_unrelaxed_gs_out_GSR.nc")])

SLA_site_2 = DeltaSCF.from_four_points_file([abidata.ref_file("site_2_relaxed_gs_out_GSR.nc"),
                                            abidata.ref_file("site_2_unrelaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("site_2_relaxed_ex_out_GSR.nc"),
                                            abidata.ref_file("site_2_unrelaxed_gs_out_GSR.nc")])

#%%
# To plot the luminescence lineshape following the one effective phonon mode model at 0K

SLA_site_1.plot_lineshape_1D_zero_temp(energy_range=[0.5,3],max_m=20,phonon_width=0.02,
               with_omega_cube='True',
               normalized='Area');

#%%
# To get a panda dataframe with the main restults :
dataframes=[]
df_1=SLA_site_1.get_dataframe('Site_1')
df_2=SLA_site_2.get_dataframe('Site_2')

pd.concat([df_1,df_2])


#%%
# To plot the magnitude of the displacements
SLA_site_1.plot_delta_R_distance(defect_symbol="Eu")

# %%
# To draw the configuration coordinates diagram.
SLA_site_1.draw_displaced_parabolas(scale_eff_freq=2.5);
SLA_site_2.draw_displaced_parabolas(scale_eff_freq=2.5);


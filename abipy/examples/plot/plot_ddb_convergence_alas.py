#!/usr/bin/env python
r"""
Phonon convergence of ε∞, BECS, phonons and quadrupoles
=======================================================

This example shows how to use the DdbRobot to analyze
the convergence of eps_inf, BECS, phonons and dynamical quadrupole
with respect to the number of k-points.
"""

import os
import abipy.data as abidata

#%%
# Initialize the DdbRobot from a list of paths to DDB files
# Here we use DDB files shipped with the AbiPy package in data/refs/alas_eps_and_becs_vs_ngkpt
# Each DDB has been computed with a different k-mesh.

paths = ["AlAs_222k_DDB", "AlAs_444k_DDB","AlAs_666k_DDB", "AlAs_888k_DDB"]
paths = [os.path.join(abidata.dirpath, "refs", "alas_eps_and_becs_vs_ngkpt", f) for f in paths]

from abipy.dfpt.ddb import DdbRobot
ddb_robot = DdbRobot.from_files(paths)

#%%
# Call anaddb to get ε∞ for all the DDB files.
# Results and metadata are stored in the
# epsinf_data.df dataframe.

epsinf_data = ddb_robot.anacompare_epsinf()

print("ε∞ dataframe:\n", epsinf_data.df)

epsinf_data.plot_conv("nkpt", abs_conv=0.1)

#%%
# Call anaddb to get BECs for all the DDB files.
# Results and metadata are stored in the becs_data.df dataframe.

becs_data = ddb_robot.anacompare_becs()
print("BECS dataframe:\n", becs_data.df)

becs_data.plot_conv("nkpt", abs_conv=0.1)

#%%
# Call anaddb to get phonons at a single q-point for all the DDB files.
# Results and metadata are stored in the ph_data.ph_df dataframe.

ph_data = ddb_robot.get_phdata_at_qpoint(qpoint=(0, 0, 0))

print(ph_data.ph_df)
ph_data.plot_ph_conv("nkpt", abs_conv=0.1)  # meV units.

#%%
# If the DDB contains dynamical quadrupoles, a similar dataframe with Q^*
# is automatically built and made available in ph_data.dyn_quad_df.
# A negative value of abs_conv is interpreted as relative convergence.

print(ph_data.dyn_quad_df)
ph_data.plot_dyn_quad_conv("nkpt", abs_conv=-0.02)

#!/usr/bin/env python
r"""
Effective masses with finite diff
=================================

This example shows how to use finite differences to compute
effective masses and plot the data.
"""

#%%
# Initialize the EffMassAnalyzer from a GSR file with energies
# computed along segments passing through the k-points of interest.
#
# Clearly the accuracy of the finite difference results strongly depends
# on the spacing between the points along the path.
# The GSR file used in this example is not accurate enough!

import abipy.data as abidata
from abipy.electrons.effmass_analyzer import EffMassAnalyzer

emana = EffMassAnalyzer.from_file(abidata.ref_file("si_nscf_GSR.nc"))
print(emana)

#%%
# Before computing the effective masses, we have to select the k-points.

emana.select_vbm()

# Alternatively, one can use:

#emana.select_cbm()
#emana.select_band_edges()

# or the most flexible API:

#emana.select_kpoint_band((0, 0, 0), band=3, spin=0, degtol_ev=0.1)

#%%
# Compute effective masses with different accuracies and print results in tabular format.
# Useful to understand if results are sensitive to the number of points in the finite difference.
# Note however that numerical differentiation is performed with the same delta_k step.

emana.summarize()

#%%
# Print the results to terminal with:

# extract segment.
segment = emana.segments[0]
#repr(segment); str(segment)
#assert segment.to_string(verbose=2)
df = segment.get_dataframe_with_accuracies(acc_list=(2, 4))
print(df)

#assert len(emana.segments) == 1
#for segment in emana.segments[0]:
#    segment.get_effmass_line(acc=2)

#%%
# Plot electronic dispersion and quadratic approximant based on the
# effective masses computed along each segment.

emana.plot_emass(acc=4)

#emana.plot_all_segments()
#emana.segments[0].plot_emass()

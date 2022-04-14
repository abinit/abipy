#!/usr/bin/env python
r"""
Alas mobility
=============

This example shows how to plot the electronic band structure,
phonon dispersion, electron linewidths and 
the convergence of the electron mobility in the case of AlAs 
(see example e-ph mobility flow).
Only electrons are considered, but the script can easily be
extended to holes (and combination of e/h as well). 
"""
#%%

import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt
from abipy.abilab import abiopen
import glob
from abipy.eph.rta import RtaRobot

root = "../../flows/flow_eph_mob/"

fz = 13 # Fontsize for text, labels,...
mz = 2  # Markersize for dots (linewidths plot)

fig, axes = plt.subplots(2, 2)
plt.subplots_adjust(wspace=0.35, hspace=0.33)

axes[0, 0].text(-0.1,1.05, '(a)', size=fz, transform=axes[0, 0].transAxes)
axes[1, 0].text(-0.1,1.05, '(b)', size=fz, transform=axes[1, 0].transAxes)
axes[0, 1].text(-0.1,1.05, '(c)', size=fz, transform=axes[0, 1].transAxes)
axes[1, 1].text(-0.1,1.05, '(d)', size=fz, transform=axes[1, 1].transAxes)

###
### Upper left : electronic band structure
###

# Open the GSR file containing the BS
with abiopen(root + "w0/t1/outdata/out_GSR.nc") as abifile:
    ebands = abifile.ebands
    formula = abifile.structure.formula
    
# Add the title to the whole figure containing the pretty formula
pretty_formula = "".join([i for i in formula if i != str(1) and i != " "]) # Remove '1's and spaces
fig.suptitle(pretty_formula)

# Plot the bands (keep only the grid on the x-axis)
ebands.plot(ax=axes[0,0], show=False, linewidth=1.5, ylims=(-3, 6), ax_grid=False)

# Set label
axes[0, 0].set_ylabel('Energy (eV)', fontsize=fz)
axes[0, 0].set_xlabel(None)

###
### Lower left : phonon dispersion
###

# Open the DDB file and get phonon dispersion
# Can be used to check the phonon dispersion obtained with the e-ph code
with abiopen(root + "w1/outdata/out_DDB") as abifile:
    phbst, phdos = abifile.anaget_phbst_and_phdos_files(ndivsm=10)

# Plot phonon dispersion (keep only the grid on the x-axis)
phbst.phbands.plot(ax=axes[1,0], show=False, units="mev", linewidth=1.5, ax_grid=False)

# Set label
axes[1, 0].set_ylabel('Frequency (meV)', fontsize=fz)
axes[1, 0].set_xlabel(None)

###
### Upper right : linewidths
###

# We get the linewidths for the densest mesh computed
# We open the SIGEPH file corresponding to this task
with abiopen(root + 'w3/t2/outdata/out_SIGEPH.nc') as abifile:
    # First within SERTA
    qparray = abifile.get_qp_array(mode="ks+lifetimes", rta_type="serta")
    qparray = qparray[np.nonzero(qparray)]
    eigs = qparray.real - np.min(qparray.real) # 0 to CBM
    lws = 2000*qparray.imag # Factor or 2x1000 to get from linewidths to 1/tau in meV
    
    # Then within MRTA
    qparray_mrta = abifile.get_qp_array(mode="ks+lifetimes", rta_type="mrta")
    qparray_mrta = qparray_mrta[np.nonzero(qparray_mrta)]
    eigs_mrta = qparray_mrta.real - np.min(qparray_mrta.real) # 0 to CBM
    lws_mrta = 2000*qparray_mrta.imag
    
# Plot 1/tau
axes[0, 1].plot(eigs, lws, 'ob', markersize=mz, label='SERTA')
axes[0, 1].plot(eigs_mrta, lws_mrta, 'xr', markersize=mz, label='MRTA')

axes[0, 1].set_xticks([0, 0.25], ['0', '0.25'])

# Set the axis labels and legend
axes[0, 1].set_xlabel(r'$\varepsilon_{n\mathbf{k}} - \varepsilon_{\mathrm{CBM}}$ (eV)', labelpad=2, fontsize=fz)
axes[0, 1].set_ylabel(r'$\tau_{n\mathbf{k}}^{-1}$ (meV)', fontsize=fz)
axes[0, 1].legend(loc='best', labelcolor='linecolor')


###
### Lower right : convergence of the mobility
###

# First we find all the RTA.nc files
abifiles = glob.glob(root+'w*/t*/outdata/*RTA.nc')

# We create a robot with these files and plot the mobilities
robot = RtaRobot.from_files(abifiles)
robot.plot_mobility_kconv(ax=axes[1, 1], eh=0, bte=["ibte", "mrta", "serta"], show=False, ax_grid=False)

# Tune the plot
axes[1, 1].set_title(None)
axes[1, 1].legend(fontsize=fz, framealpha=0.5)
axes[1, 1].set_xlabel('$N_k \\times$ $N_k \\times$ $N_k$ $\\mathbf{k}$-point grid', fontsize=fz)
axes[1, 1].set_ylabel(r'$\mu_e$'+' (cm$^2$/(V$\\cdot$s))', fontsize=fz)

# Reactivate the grid on the x-axis for the band structures
axes[0, 0].grid(True, axis='x')
axes[1, 0].grid(True, axis='x')

# We save the figure in pdf format
fig.savefig(pretty_formula + ".pdf", bbox_inches='tight')

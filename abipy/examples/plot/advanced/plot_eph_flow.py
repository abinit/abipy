import numpy as np
import matplotlib.pyplot as plt
from abipy.abilab import abiopen
import os
import glob
from abipy.eph.rta import RtaRobot

root = "../../flows/flow_eph_mob/"

fz = 13
mz = 2

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

# Plot the bands
ebands.plot_ax(axes[0,0], e0="fermie", color="black")

# Decorate the axes with the correct ticks and labels
ticks = [tick for tick, label in enumerate(ebands.kpoints.names) if label != None]
labels = [label for tick, label in enumerate(ebands.kpoints.names) if label != None]
axes[0, 0].set_xticks(ticks, [], size=fz)
for tick in ticks:
    axes[0, 0].plot([tick, tick], [-100, 100], '-k', linewidth=1/2, alpha=0.6)
    
# Set limits of the axis
axes[0, 0].set_xlim([ticks[0], ticks[-1]])
axes[0, 0].set_ylim([-3, 6])

# Set label
axes[0, 0].set_ylabel('Energy (eV)', fontsize=fz)

###
### Lower left : phonon dispersion
###

# Open the DDB file and get phonon dispersion
with abiopen(root + "w1/outdata/out_DDB") as abifile:
    phbst, phdos = abifile.anaget_phbst_and_phdos_files(ndivsm=10)

# Plot phonon dispersion
phbst.phbands.plot_ax(axes[1,0], branch=None, units="mev", color="black")

# Decorate the axes with the correct ticks and labels
ticks = [tick for tick, label in enumerate(phbst.qpoints.names) if label != None]
labels = [label for tick, label in enumerate(phbst.qpoints.names) if label != None]
axes[1, 0].set_xticks(ticks, labels, size=fz-2)
for tick in ticks:
    axes[1, 0].plot([tick, tick], [-100, 100], '-k', linewidth=1/2, alpha=0.6)

# Set limits of the axis
axes[1, 0].set_xlim([ticks[0], ticks[-1]])
axes[1, 0].set_ylim([0, np.max(phbst.phbands.phfreqs)*1000+1])

# Set label
axes[1, 0].set_ylabel('Frequency (meV)', fontsize=fz)

###
### Upper right : linewidths
###

# We get the linewidths for the densest mesh computed
# We open the SIGEPH file corresponding to this task
with abiopen(root + 'w3/t2/outdata/out_SIGEPH.nc') as abifile:
    qplist = abifile.reader.read_allqps()[0]
    eigs = qplist.get_field_itemp(field="qpe", itemp=0).real
    lws = 2000*qplist.get_field_itemp(field="fan0", itemp=0).imag # Factor or 2x1000 to get from linewidths to 1/tau in meV

# The zero of the energy axis is the CBM if there are holes, 
zero = np.min(eigs) # CBM

# Plot 1/tau for electrons only : much easier
axes[0, 1].plot(eigs-zero, lws, 'ok', markersize=mz)
xlabel = r'$\varepsilon_{n\mathbf{k}} - \varepsilon_{\mathrm{CBM}}$ (eV)'
axes[0, 1].set_xticks([0, 0.25], ['0', '0.25'])

# Set the axis labels
axes[0, 1].set_xlabel(xlabel, labelpad=2, fontsize=fz)
axes[0, 1].set_ylabel(r'$\tau_{n\mathbf{k}}^{-1}$ (meV)', fontsize=fz)

###
### Lower right : convergence of the mobility
###

# First we find all the RTA.nc files
abifiles = glob.glob(root+'w*/t*/outdata/*RTA.nc')

# We create a robot with these files and get the k-meshes and mobilities
# Alternatively we could sort the files and get the data ourselves
robot = RtaRobot.from_files(abifiles)
figconv = robot.plot_mobility_kconv(eh=0, bte=["ibte"], show=False)

# xdata : ngkpt, ydata: IBTE mobility
xdata = figconv.axes[0].lines[0].get_xdata()
ydata = figconv.axes[0].lines[0].get_ydata()
   
# We plot the convergence on our axes
axes[1, 1].plot(xdata, ydata, '-ok')
# Set the xticks on the ngkpt values
axes[1, 1].set_xticks(xdata)
# Set the yticks at the minimum and maximum mobilities only
axes[1, 1].set_yticks([ydata[0], ydata[-1]], [int(np.round(ydata[0])), int(np.round(ydata[-1]))])
# Set the labels
axes[1, 1].set_xlabel('$N_k \\times$ $N_k \\times$ $N_k$ $\\mathbf{k}$-point grid', fontsize=fz)
axes[1, 1].set_ylabel(r'$\mu_e^\mathrm{IBTE}$'+' (cm$^2$/(V$\\cdot$s))', fontsize=fz)

# We save the figure in pdf format
fig.savefig(pretty_formula + ".pdf", bbox_inches='tight')

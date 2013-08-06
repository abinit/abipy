#!/usr/bin/env python
#
# This example shows how to visualize the QP results 
# stored in the SIGRES produced by the GW code (sigma run)
from abipy import *

sigma_file = SIGRES_File(get_reference_file("tgw1_9o_DS4_SIGRES.nc"))

# Printout of the QPState results
sigma_file.print_qps()

qp = sigma_file.get_qpcorr(spin=0, kpoint=(0,0,0), band=0)
#print(qp)

qplist_spin = sigma_file.qplist_spin
#print(qplist_spin[0])
#import sys
#sys.exit(1)

qplist_spin[0].plot_qps_vs_e0(title="QPState corrections of Si", exclude_fields="vUme")

#print(sigw)
#sigw.plot()

#spfun = sigma_file.get_spfunc(spin, kpoint, band)
#spfun_intg = spfun.cumintegral()
#spfun.plot()
#spfun_intg.plot()

#qp_bands.plot(title="Silicon band structure", klabels = klabels)

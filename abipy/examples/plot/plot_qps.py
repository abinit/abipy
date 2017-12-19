#!/usr/bin/env python
r"""
SIGRES file (GW)
================

This example shows how to visualize the QP results 
stored in the SIGRES produced by the GW code (sigma run)
"""
import abipy.data as abidata
from abipy.abilab import abiopen

sigres = abiopen(abidata.ref_file("tgw1_9o_DS4_SIGRES.nc"))

# Printout of the QPState results
sigres.print_qps()

sigres.plot_qps_vs_e0(tight_layout=True)

qp = sigres.get_qpcorr(spin=0, kpoint=(0, 0, 0), band=0)
print(qp)

#qplist_spin = sigres.qplist_spin
#qplist_spin[0].plot_qps_vs_e0(title="QPState corrections of Si", exclude_fields="vUme")

sigres.close()

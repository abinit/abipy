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
#sigres.print_qps()

sigres.plot_qps_vs_e0(tight_layout=True, title="QP data vs KS energy")

sigres.plot_qpgaps(title="QP direct gaps")

sigres.plot_qpgaps(plot_qpmks=True, title="QP - KS direct gap")

sigres.plot_qpbands_ibz()

sigres.plot_ksbands_with_qpmarkers(fact=1000)

sigres.close()

#!/usr/bin/env python
r"""
self-consistent GW
==================

This example shows how to visualize the SCGW QP amplitudes in the KS basis set.
"""
import abipy.data as abidata
from abipy.abilab import abiopen

sigres = abiopen(abidata.ref_file("QPSC_SIGRES.nc"))

#print("calctyp",sigres.gwcalctyp)
#sigres.print_qps()

#qp = sigres.get_qpcorr(spin=0, kpoint=(0,0,0), band=0)
#print(qp)

# Visualize <KS_b|QP_b'> matrix at the Gamma point
sigres.plot_eigvec_qp(kpoint=[0, 0, 0])

# Plot <KS_b|QP_b'> components for all k-points in the SIGRES.
sigres.plot_eigvec_qp()

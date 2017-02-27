#!/usr/bin/env python
"""
This example shows how to visualize the SCGW QP amplitudes in the KS basis set.
"""
import abipy.data as abidata
from abipy.abilab import abiopen

sigred = abiopen(abidata.ref_file("QPSC_SIGRES.nc"))

print("calctyp",sigred.gwcalctyp)
sigres.print_qps()

#qp = sigres.get_qpcorr(spin=0, kpoint=(0,0,0), band=0)
#print(qp)

sigres.plot_eigvec_qp(spin=0, kpoint=None, title="<KS_b|QP_b'> components")

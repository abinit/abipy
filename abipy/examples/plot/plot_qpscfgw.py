#!/usr/bin/env python
#
# This example shows how to visualize the SCGW QP amplitudes 
# in the KS basis set.
import abipy.data as abidata
from abipy.abilab import abiopen

sigma_file = abiopen(abidata.ref_file("QPSC_SIGRES.nc"))

print("calctyp",sigma_file.gwcalctyp)
sigma_file.print_qps()

#qp = sigma_file.get_qpcorr(spin=0, kpoint=(0,0,0), band=0)
#print(qp)

sigma_file.plot_eigvec_qp(spin=0, kpoint=None, title="<KS_b|QP_b'> components")


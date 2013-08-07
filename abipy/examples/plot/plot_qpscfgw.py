#!/usr/bin/env python
#
# This example shows how to visualize the SCGW QP amplitudes 
# in the KS basis set.
from abipy import *
import abipy.data as data

sigma_file = SIGRES_File(data.ref_file("QPSC_SIGRES.nc"))

print("calctyp",sigma_file.gwcalctyp)
sigma_file.print_qps()

#qp = sigma_file.get_qpcorr(spin=0, kpoint=(0,0,0), band=0)
#print(qp)

sigma_file.plot_eigvec_qp(spin=0, kpoint=None, title="<KS_b|QP_b'> components")


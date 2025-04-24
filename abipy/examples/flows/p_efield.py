#!/usr/bin/env python
import sys
from abipy import abilab
from abipy.abio.outputs import BerryPhasePolarization

#pol = BerryPhasePolarization.from_abo_file("flow_fd_efield/w0/t0/run.abo")
#print(pol)

flow = abilab.abiopen("flow_fd_efield/__AbinitFlow__.pickle")
work = flow[0]
data = work.get_data()

print(data.get_epsinf_df())
data.plot_polarization(what="total")
sys.exit(0)

data.print_eff_charges()
data.plot_etotal()
data.plot_forces()
data.plot_stresses()


#!/usr/bin/env python
import sys

from abipy import abilab

path = "flow_fd_efield/__AbinitFlow__.pickle"
if len(sys.argv) > 1:
    path = sys.argv[1]

print("Reading from:", path)

flow = abilab.abiopen(path)
work = flow[0]
data = work.get_data()
data = data["clamped_ions"]
print(data)
#data.print_relaxed_coords()

print("Exiting now, no plots"); sys.exit(0)
data.plot_forces_vs_field([1, 0, 0], elements="Al")
data.plot_etotal()
data.plot_polarization(what="total")
data.plot_forces()
data.plot_stresses()

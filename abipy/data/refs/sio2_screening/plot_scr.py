#!/usr/bin/env python
"""
This examples shows how to plot the 
"""
import abipy.data as abidata
from abipy.abilab import abiopen

with abiopen(abidata.ref_file("sio2_SCR.nc")) as ncfile:
    ncfile.plot_emacro(title="foo")
    #ncfile.plot_eelf(title="bar")


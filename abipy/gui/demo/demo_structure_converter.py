#!/usr/bin/env python

import abipy.data as abidata

from abipy.gui.wxapps import wxapp_structure_converter

filename = abidata.ref_file("si_nscf_GSR.nc")

app = wxapp_structure_converter(filename)
app.MainLoop()

#!/usr/bin/env python

import abipy.data as abidata

from abipy.gui.wxapps import wxapp_structure_converter

app = wxapp_structure_converter(abidata.ref_file("si_nscf_GSR.nc"))
app.MainLoop()

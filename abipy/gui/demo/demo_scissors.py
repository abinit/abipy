#!/usr/bin/env python
import abipy
import abipy.data

from abipy.gui.wxapps import wxapp_scissors

app = wxapp_scissors(abipy.data.ref_file("tgw1_9o_DS4_SIGRES.nc"))
app.MainLoop()

#!/usr/bin/env python
import abipy
import abipy.data

from abipy.gui.wxapps import wxapp_sigresviewer

app = wxapp_sigresviewer(abipy.data.ref_file("tgw1_9o_DS4_SIGRES.nc"))
app.MainLoop()

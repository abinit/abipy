#!/usr/bin/env python

import abipy

from abipy.gui.wxapps import wxapp_sigresviewer

app = wxapp_sigresviewer(abipy.get_reference_file("tgw1_9o_DS4_SIGRES.nc"))
app.MainLoop()

#!/usr/bin/env python
import os
import abipy
import abipy.data

from abipy.gui.wxapps import wxapp_wfkviewer

datadir = abipy.data.dirpath

wfk_filename = [f for f in os.listdir(datadir) if "_WFK" in f][0]
wfk_filename = os.path.join(datadir, wfk_filename)

app = wxapp_wfkviewer(wfk_filename)
app.MainLoop()


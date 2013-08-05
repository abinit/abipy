#!/usr/bin/env python
import os
import abipy

import abipy.gui.wxapps as wxapps

datadir = abipy.get_datadir()
filepaths = [os.path.join(datadir, f) for f in os.listdir(datadir)]

app = wxapps.wxapp_comparison(filepaths=filepaths).MainLoop()

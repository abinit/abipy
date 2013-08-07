#!/usr/bin/env python
import os
import abipy
import abipy.data

import abipy.gui.wxapps as wxapps

datadir = abipy.data.dirpath
filepaths = [os.path.join(datadir, f) for f in os.listdir(datadir)]

app = wxapps.wxapp_comparison(filepaths=filepaths).MainLoop()

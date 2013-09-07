#!/usr/bin/env python
import os
import abipy
import abipy.data

import abipy.gui.wxapps as wxapps

dirpath = os.path.join(abipy.data.dirpath, "runs", "data_si_ebands")
wxapps.wxapp_events(dirpath).MainLoop()

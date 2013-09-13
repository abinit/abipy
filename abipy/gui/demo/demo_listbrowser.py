#!/usr/bin/env python
import abipy
import abipy.data

import abipy.gui.wxapps as wxapps

wxapps.wxapp_listbrowser(dirpaths=abipy.data.dirpath).MainLoop()


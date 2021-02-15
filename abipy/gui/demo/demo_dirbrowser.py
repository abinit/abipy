#!/usr/bin/env python
import os
import abipy
import abipy.data

import abipy.gui.wxapps as wxapps

wxapps.wxapp_dirbrowser(abipy.data.dirpath).MainLoop()


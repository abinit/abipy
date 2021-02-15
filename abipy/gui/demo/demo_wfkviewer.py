#!/usr/bin/env python
import os
import abipy.data as abidata

from abipy.gui.wxapps import wxapp_wfkviewer

wfk_filename = abidata.WFK_NCFILES[0]

wxapp_wfkviewer(wfk_filename).MainLoop()

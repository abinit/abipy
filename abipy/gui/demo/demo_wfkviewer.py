#!/usr/bin/env python
import os
import abipy.data as data

from abipy.gui.wxapps import wxapp_wfkviewer

wfk_filename = data.WFK_NCFILES[0]

wxapp_wfkviewer(wfk_filename).MainLoop()

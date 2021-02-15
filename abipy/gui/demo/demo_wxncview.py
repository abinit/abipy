#!/usr/bin/env python
import abipy.data as abidata

from abipy.gui.wxapps import wxapp_ncview

wxapp_ncview(filepaths=abidata.WFK_NCFILES[0]).MainLoop()

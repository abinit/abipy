#!/usr/bin/env python

import os
import wx
import abipy

from abipy.gui.comparison import ComparisonFrame

datadir = abipy.get_datadir()

app = wx.App()
filepaths = [os.path.join(datadir, f) for f in os.listdir(datadir)]
frame = ComparisonFrame(None, filepaths=filepaths)
app.SetTopWindow(frame)
frame.Show()
app.MainLoop()

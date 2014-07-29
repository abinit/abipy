#!/usr/bin/env python
import wx
from abipy.gui.awx.elements_gui import WxPeriodicTable

app = wx.App()
frame = WxPeriodicTable(None)
app.SetTopWindow(frame)
frame.Show()
app.MainLoop()

#!/usr/bin/env python
import wx
from abipy.gui.awx.elements_gui import MainFrame

app = wx.App()
frame = MainFrame(None)
app.SetTopWindow(frame)
frame.Show()
app.MainLoop()

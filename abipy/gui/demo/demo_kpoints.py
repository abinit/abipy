#!/usr/bin/env python
import os
import wx
import abipy

from abipy.gui.kpoints import KpointSymmetriesFrame

datadir = abipy.get_datadir()
wfk_filename = [f for f in os.listdir(datadir) if "_WFK" in f][0]
wfk_filename = os.path.join(datadir, wfk_filename)                   

app = wx.App()

wfk = abipy.abiopen(wfk_filename)

frame = KpointSymmetriesFrame(None, wfk.get_structure(), wfk.kpoints)
frame.Show()

app.MainLoop()

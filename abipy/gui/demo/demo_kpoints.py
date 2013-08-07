#!/usr/bin/env python
import os
import wx
import abipy
import abipy.data

from abipy.gui.kpoints import KpointSymmetriesFrame

datadir = abipy.data.dirpath
wfk_filename = [f for f in os.listdir(datadir) if "_WFK" in f][0]
wfk_filename = os.path.join(datadir, wfk_filename)                   

app = wx.App()

wfk = abipy.abiopen(wfk_filename)

frame = KpointSymmetriesFrame(None, wfk.structure, wfk.kpoints)
frame.Show()

app.MainLoop()

#!/usr/bin/env python
import wx
import abipy
import abipy.data as data

from abipy.gui.kpoints import KpointSymmetriesFrame

wfk_filename = data.WFK_NCFILES[0]
wfk = abipy.abiopen(wfk_filename)

app = wx.App()
frame = KpointSymmetriesFrame(None, wfk.structure, wfk.kpoints)
frame.Show()
app.MainLoop()

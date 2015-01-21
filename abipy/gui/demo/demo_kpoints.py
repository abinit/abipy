#!/usr/bin/env python
import wx

import abipy.data as abidata

from abipy.abilab import abiopen
from abipy.gui.kpoints import KpointsPanel

gsr = abiopen(abidata.ref_file("si_nscf_GSR.nc"))

app = wx.App()
frame = wx.Frame(None, -1)
KpointsPanel(frame, gsr.structure, gsr.kpoints)
app.SetTopWindow(frame)
frame.Show()
app.MainLoop()
gsr.close()

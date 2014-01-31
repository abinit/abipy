#!/usr/bin/env python
import wx
import numpy as np

from abipy.core import Function1D
from abipy.gui.awx.func1dframe import Func1dPlotFrame

app = wx.App()
func1d = Function1D.from_func(np.sin, np.arange(1,100,1))
frame = Func1dPlotFrame(None, func1d)
frame.Show()
app.MainLoop()

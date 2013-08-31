#!/usr/bin/env python
import os
import wx

from abipy.gui.notebooks import ReadOnlyTextNotebookFrame

filenames = [f for f in os.listdir(".") if f.endswith(".py")]

text_list = []
for fname in filenames:
    with open(fname, "r") as fh:
        text_list.append(fh.read())

app = wx.App()
frame = ReadOnlyTextNotebookFrame(None, text_list, page_names=filenames)
frame.Show()
app.MainLoop()

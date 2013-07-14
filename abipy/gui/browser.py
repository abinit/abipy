#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import wx

import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy.iotools.files import NcDumper 
from abipy.electrons import ElectronBands
from abipy.waves import WFK_File
from wx.lib.agw.floatspin import FloatSpin
from wxmplot import PlotApp, PlotFrame

from abipy.gui.popupmenus import popupmenu_from_ext

from abipy.gui.wfkviewer import wfk_viewer

_VIEWERS = {
    "WFK-etsf.nc": wfk_viewer,
}

def run_viewer_for_filepath(filepath):
    ext = filepath.split("_")[-1]
    try:
        return _VIEWERS[ext](filepath)
    except KeyError:
        raise KeyError("No wx viewer s has been registered for the extension %s" % ext)


class NcFileDirCtrl(wx.GenericDirCtrl):
    def __init__(self, *args, **kwargs):

        if "filter" not in kwargs:
            kwargs["filter"] = "Netcdf files (*.nc)|*.nc|All files (*.*)|*.*"
        if "dir" not in kwargs:
            kwargs["dir"] = os.getcwd()

        kwargs["style"] = wx.TR_MULTIPLE
        super(NcFileDirCtrl, self).__init__(*args, **kwargs)

        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.OnItemActivated)
        self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClick)

    def OnItemActivated(self, event):
        path = self.GetFilePath()
        if not path: return
        print("in activated with path %s" % path)
        run_viewer_for_filepath(path)

    def OnRightClick(self, event):
        path = self.GetFilePath()
        if not path: return
        print("in right with path %s" % path)

        # Open the popup menum then destroy it to avoid mem leak.
        popmenu = popupmenu_from_ext(path)
        self.PopupMenu(popmenu, event.GetPoint())
        popmenu.Destroy() 

def wxabi_browser(dirpath):
    if dirpath is None: 
        dirpath = ""
    else:
        dirpath = os.path.abspath(dirpath)

    app = wx.App()

    frame = wx.Frame(None, -1)
    dirctrl = NcFileDirCtrl(frame, -1, dir=dirpath)
    frame.Show()

    wx.LogError("Ciao")
    #wx.LoadFileSelector("hello", "txt")
    #wx.InfoMessageBox(None)

    app.MainLoop()
 

if __name__ == "__main__":
    dirpath = None
    if len(sys.argv) > 1:
        dirpath = sys.argv[1]
    wxabi_browser(dirpath)

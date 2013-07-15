#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import wx

import wx.lib.dialogs as wxdg 
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

#from wx.lib.editor.editor import Editor
from abipy.waves import WFK_File
from abipy.iotools.visualizer import supported_visunames

from pymatgen.io.abinitio import EventParser

#class EventPanel(wx.Panel):
#
#    def __init__(self, parent, filepath, **kwargs):
#        super(EventPanel, self).__init__(parent, -1, **kwargs)
#        self.filepath = os.path.abspath(filepath)
#
#        parser = EventParser()
#        self.events = parser.parse(self.filepath)
#        print(self.events)
#                                                                       
#        self.BuildUi()
#                                                                       
#    def BuildUi(self):
#        pass



class EventPanel(wx.Panel):
    def __init__(self, parent, filepath):
        super(EventPanel, self).__init__(parent, -1)

        self.filepath = os.path.abspath(filepath)

        self.BuildUi()

    def BuildUi(self):

        parser = EventParser()
        self.events = events = parser.parse(self.filepath)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        vbox = wx.BoxSizer(wx.VERTICAL)
        panel1 = wx.Panel(self, -1)
        panel2 = wx.Panel(self, -1)

        self.tree = tree = wx.TreeCtrl(panel1, 1, wx.DefaultPosition, (-1,-1), wx.TR_HIDE_ROOT|wx.TR_HAS_BUTTONS)

        root = self.tree.AddRoot('Events')

        err_tree = tree.AppendItem(root, "%d Errors" % events.num_errors)
        warn_tree = tree.AppendItem(root, "%s Warnings" % events.num_warnings)
        com_tree = tree.AppendItem(root, "%s Comments" % events.num_comments)

        for e in self.events.errors:
            tree.AppendItem(err_tree, "line: " + str(e.lineno), data=wx.TreeItemData(e.message))

        for w in self.events.warnings:
            tree.AppendItem(warn_tree, "line: " + str(w.lineno), data=wx.TreeItemData(w.message))

        for c in self.events.comments:
            tree.AppendItem(com_tree, "line: " + str(c.lineno), data=wx.TreeItemData(c.message))

        tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnSelChanged, id=1)

        self.display = wx.StaticText(panel2, -1, '',(10,10), style=wx.ALIGN_LEFT)

        vbox.Add(self.tree, 1, wx.EXPAND)
        hbox.Add(panel1, 1, wx.EXPAND)
        hbox.Add(panel2, 1, wx.EXPAND)
        panel1.SetSizer(vbox)

        self.SetSizer(hbox)
        self.Centre()

    def OnSelChanged(self, event):
        item =  event.GetItem()
        lineno = self.tree.GetItemText(item)
        proxy = self.tree.GetItemData(item) 

        if proxy is not None:
            data = proxy.GetData()
            self.display.SetLabel(data)

def abo_viewer(filepath):
    """Start up the AbinitOutputViewer application."""

    class AboViewerApp(wx.App):
        def OnInit(self): 
            #frame = AbinitOutputViewer(parent=None, filepath=filepath) 
            frame = wx.Frame(None)
            panel = EventPanel(frame, filepath)
            #Editor(panel, -1, style=wx.SUNKEN_BORDER)
            frame.Show() 
            self.SetTopWindow(frame) 
            return True

    AboViewerApp().MainLoop()


if __name__ == "__main__":
    import sys
    filepath = None
    if len(sys.argv) > 1:
        filepath = sys.argv[1] 
    abo_viewer(filepath)

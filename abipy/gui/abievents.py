#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import wx

import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from collections import OrderedDict
from abipy import abiopen
from abipy.waves import WFK_File
from abipy.iotools.visualizer import supported_visunames

from pymatgen.io.abinitio import EventParser


class AbinitEventsPanel(wx.Panel):
    """
    Panel with a TreeCtrl that allows the user to navigate
    the events (WARNINGS/COMMENTS/ERRORS) reported by ABINIT
    in the main output file or in the log file.
    """
    def __init__(self, parent, filepath):
        super(AbinitEventsPanel, self).__init__(parent, -1)

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

    @property
    def has_events(self):
        return len(self.events) > 0

    def OnSelChanged(self, event):
        item =  event.GetItem()
        lineno = self.tree.GetItemText(item)
        proxy = self.tree.GetItemData(item)

        if proxy is not None:
            data = proxy.GetData()
            self.display.SetLabel(data)


class AbinitEventsFrame(wx.Frame):
    """
    Frame with an EventsPanel
    """
    def __init__(self, parent, filepath):
        self.filepath = os.path.abspath(filepath)
        title = "Abinit Events: %s" % os.path.basename(self.filepath)
        super(AbinitEventsFrame, self).__init__(parent, -1, title=title)

        self.BuildUi()

    def BuildUi(self):
        self.event_panel = AbinitEventsPanel(self, self.filepath)

class AbiOutLogDirCtrl(wx.GenericDirCtrl):

    def __init__(self, *args, **kwargs):

        if "filter" not in kwargs:
            kwargs["filter"] = "All files (*.*)|*.*|about files (*.about)|*.out| ablog files (*.ablog)|*.ablog"
        if "dir" not in kwargs:
            kwargs["dir"] = os.getcwd()
        if "style" not in kwargs:
            kwargs["style"] = wx.TR_MULTIPLE

        super(AbiOutLogDirCtrl, self).__init__(*args, **kwargs)

        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.OnItemActivated)
        #self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClick)

    def OnItemActivated(self, event):
        path = self.GetFilePath()
        if not path: return
        print("in activated with path %s" % path)
        frame = EventFrame(self, path).Show()

    #def OnRightClick(self, event):
    #    path = self.GetFilePath()
    #    if not path: return
    #    print("in right with path %s" % path)


class AbinitEventsNotebookFrame(wx.Frame):
    def __init__(self, parent, filenames, **kwargs):

        if "title" not in kwargs:
            kwargs["title"] = "Abinit Events"

        super(AbinitEventsNotebookFrame, self).__init__(parent, **kwargs)

        # Here we create a panel and a notebook on the panel
        p = wx.Panel(self)

        import wx.lib.agw.flatnotebook as fnb
        #nb = wx.Notebook(p)
        nb = fnb.FlatNotebook(p)

        # Add the pages to the notebook with the label to show on the tab
        # Add only files for which we have events.
        for fname in filenames:
            page = AbinitEventsPanel(nb, fname)
            if page.has_events:
                nb.AddPage(page, text=os.path.basename(fname))

        # Finally, put the notebook in a sizer for the panel to manage the layout
        sizer = wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND)
        p.SetSizer(sizer)

def wxapp_events(root):
    """
    Start up the AbinitOutputViewer application.

    Args:
        root:
            Can be: None, filename, directory name or list of filenames.
            None means that we just open the browser without accessing any file.
            If root is a directory, we locate all the output files
            starting from root and we visualize them in the main Frame.
    """
    if root is None:
        filenames = []

    elif isinstance(root, str):
        root = os.path.abspath(root)
        if os.path.isdir(root):
                #filenames = [os.path.join(root, f) for f in ["t01.out", "t02.out"]]
            filenames = [os.path.join(root, f) for f in os.listdir(root) if f.endswith(".out")]
        else:
            filenames = [root]
    else:
        filenames = root
    print(filenames)

    class AbiEventsViewerApp(wx.App):
        def OnInit(self):
            frame = AbinitEventsNotebookFrame(None, filenames)
            self.SetTopWindow(frame)
            frame.Show()
            return True

    return AbiEventsViewerApp()


class AbinitTimerFrame(wx.Frame):
    """
    Frame with controls to plot the timing data.
    """
    def __init__(self, parent, filepath):
        filepath = os.path.abspath(filepath)
        title = "Abinit Timer: %s" % os.path.basename(filepath)
        super(AbinitTimerFrame, self).__init__(parent, -1, title=title)

        abifile = abiopen(filepath)

        try:
            self.timer_data = abifile.timer_data
        except Exception as exc:
            raise awx.Error(str(exc))

        if not self.timer_data:
            raise awx.Error("%s does not contain a valid ABINIT TIMER section!" % filepath)

        self.BuildUi()

    def BuildUi(self):
        timer_data = self.timer_data

        # Callbacks (bound methods of AbiTimerData).
        self.plot_types = OrderedDict([
            ("pie",          timer_data.show_pie),
            ("stacked_hist", timer_data.show_stacked_hist),
        ])

        from abipy.htc.abitimer import Section
        keys = Section.NUMERIC_FIELDS

        main_sizer = wx.BoxSizer( wx.VERTICAL )

        hsizer = wx.BoxSizer( wx.HORIZONTAL )

        self.plot_label = wx.StaticText( self, wx.ID_ANY, "plot type:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.plot_label.Wrap( -1 )
        hsizer.Add( self.plot_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        self.plot_cbox = wx.ComboBox( self, wx.ID_ANY, "pie", wx.DefaultPosition, wx.DefaultSize, self.plot_types.keys(), 0)
        hsizer.Add( self.plot_cbox, 0, wx.ALL, 5 )

        self.key_label = wx.StaticText( self, wx.ID_ANY, "key:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.key_label.Wrap( -1 )
        hsizer.Add( self.key_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        self.key_cbox = wx.ComboBox( self, wx.ID_ANY, "wall_time", wx.DefaultPosition, wx.DefaultSize, keys, 0 )
        hsizer.Add( self.key_cbox, 0, wx.ALL, 5 )

        main_sizer.Add( hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.plot_button = wx.Button( self, wx.ID_ANY, "Plot", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.Bind(wx.EVT_BUTTON, self.OnClick, self.plot_button)
        main_sizer.Add( self.plot_button, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        self.SetSizer( main_sizer )
        self.Layout()

    def OnClick(self, event):
        callback = self.plot_types[self.plot_cbox.GetValue()]
        kwargs = dict(
            key=str(self.key_cbox.GetValue())
        )
        print(callback, kwargs)
        callback(**kwargs)

if __name__ == "__main__":
    import sys
    root = None
    if len(sys.argv) > 1: root = sys.argv[1:]
    awx_events(root)

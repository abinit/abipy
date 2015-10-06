from __future__ import print_function, division, unicode_literals

import os
import wx
import wx.lib.agw.flatnotebook as fnb
import abipy.gui.awx as awx

from collections import OrderedDict
from monty.string import list_strings, is_string
from pymatgen.io.abinit.events import EventsParser


class AbinitEventsPanel(awx.Panel):
    """
    Panel with a TreeCtrl that allows the user to navigate
    the events (WARNINGS/COMMENTS/ERRORS) reported by ABINIT
    in the main output file or in the log file.
    """
    def __init__(self, parent, filepath, **kwargs):
        super(AbinitEventsPanel, self).__init__(parent, -1, **kwargs)

        self.filepath = os.path.abspath(filepath)

        parser = EventsParser()
        self.events = events = parser.parse(self.filepath)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        vbox = wx.BoxSizer(wx.VERTICAL)
        panel1 = awx.Panel(self, -1)
        panel2 = awx.Panel(self, -1)

        self.tree = tree = wx.TreeCtrl(panel1, 1, wx.DefaultPosition, (-1, -1), wx.TR_HIDE_ROOT | wx.TR_HAS_BUTTONS)

        root = self.tree.AddRoot('Events')

        err_tree = tree.AppendItem(root, "%d Errors" % events.num_errors)
        warn_tree = tree.AppendItem(root, "%d Warnings" % events.num_warnings)
        com_tree = tree.AppendItem(root, "%d Comments" % events.num_comments)

        for e in self.events.errors:
            tree.AppendItem(err_tree, "line: " + str(e.lineno), data=wx.TreeItemData(e.message))

        for w in self.events.warnings:
            tree.AppendItem(warn_tree, "line: " + str(w.lineno), data=wx.TreeItemData(w.message))

        for c in self.events.comments:
            tree.AppendItem(com_tree, "line: " + str(c.lineno), data=wx.TreeItemData(c.message))

        tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnSelChanged, id=1)

        self.display = wx.StaticText(panel2, -1, '', (10, 10), style=wx.ALIGN_LEFT)

        vbox.Add(self.tree, 1, wx.EXPAND)
        main_sizer.Add(panel1, 1, wx.EXPAND)
        main_sizer.Add(panel2, 1, wx.EXPAND)
        panel1.SetSizerAndFit(vbox)

        self.SetSizerAndFit(main_sizer)
        self.Centre()

    @property
    def has_events(self):
        return len(self.events) > 0

    def OnSelChanged(self, event):
        item = event.GetItem()
        proxy = self.tree.GetItemData(item)
        if proxy is None: return
        data = proxy.GetData()
        self.display.SetLabel(data)


class AbinitEventsFrame(awx.Frame):
    """
    Frame with an EventsPanel
    """
    def __init__(self, parent, filepath, **kwargs):
        filepath = os.path.abspath(filepath)

        if "title" not in kwargs:
            kwargs["title"] = "Abinit Events: %s" % os.path.basename(filepath)

        super(AbinitEventsFrame, self).__init__(parent, **kwargs)

        self.event_panel = AbinitEventsPanel(self, filepath)


class AbiOutLogDirCtrl(wx.GenericDirCtrl):

    def __init__(self, *args, **kwargs):
        if "filter" not in kwargs:
            kwargs["filter"] = "All files (*.*)|*.*|ABINIT output files (*.abo)|*.out| ABINIT log files (*.log)|*.log"
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
        EventFrame(self, path).Show()

    #def OnRightClick(self, event):
    #    path = self.GetFilePath()
    #    if not path: return
    #    print("in right with path %s" % path)


class AbinitEventsNotebookFrame(awx.Frame):
    """
    This frame receives a list of filenames and displays the ABINIT events in a notebook.
    """
    def __init__(self, parent, filenames, num_dirs=2, **kwargs):
        """
        Args:
            parent:
                Parent Widget.
            filenames:
                List of filenames.
            num_dirs:
                Maximum number of directories that will be shown in the tab.
        """
        if "title" not in kwargs:
            kwargs["title"] = "Abinit Events"

        super(AbinitEventsNotebookFrame, self).__init__(parent, **kwargs)

        filenames = list_strings(filenames)

        # Remove inexistent files.
        filenames = filter(os.path.exists, filenames)

        if not filenames:
            return

        # Here we create a panel and a notebook on the panel
        panel = awx.Panel(self)

        nb = fnb.FlatNotebook(panel)

        for fname in filenames:
            page = AbinitEventsPanel(nb, fname)
            page_name = fname

            if num_dirs > 0:
                tokens = page_name.split(os.path.sep)
                page_name = os.path.join(*tokens[-num_dirs:])

            # Add only files for which we have events.
            #if page.has_events:

            # Add the pages to the notebook with the name to show on the tab
            nb.AddPage(page, text=page_name)

        # Finally, put the notebook in a sizer for the panel to manage the layout
        sizer = wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND)

        panel.SetSizerAndFit(sizer)


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

    elif is_string(root):
        root = os.path.abspath(root)
        if os.path.isdir(root):
            filenames = [os.path.join(root, f) for f in os.listdir(root) if f.endswith(".abo")]
        else:
            filenames = [root]
    else:
        filenames = root

    class AbiEventsViewerApp(awx.App):
        def OnInit(self):
            frame = AbinitEventsNotebookFrame(None, filenames)
            self.SetTopWindow(frame)
            frame.Show()
            return True

    return AbiEventsViewerApp()

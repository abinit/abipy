from __future__ import print_function, division

import os
import wx

import abipy.gui.awx as awx

from collections import OrderedDict
from pymatgen.io.abinitio import EventParser
from abipy import abiopen
from abipy.htc.abitimer import AbinitTimerSection

def is_string(obj):
    try:
        dummy = obj + " "
        return True

    except TypeError:
        return False


class AbinitEventsPanel(awx.Panel):
    """
    Panel with a TreeCtrl that allows the user to navigate
    the events (WARNINGS/COMMENTS/ERRORS) reported by ABINIT
    in the main output file or in the log file.
    """

    def __init__(self, parent, filepath, **kwargs):
        super(AbinitEventsPanel, self).__init__(parent, -1, **kwargs)

        self.filepath = os.path.abspath(filepath)

        self.BuildUi()

    def BuildUi(self):
        parser = EventParser()
        self.events = events = parser.parse(self.filepath)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
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
        hbox.Add(panel1, 1, wx.EXPAND)
        hbox.Add(panel2, 1, wx.EXPAND)
        panel1.SetSizerAndFit(vbox)

        self.SetSizerAndFit(hbox)
        self.Centre()

    @property
    def has_events(self):
        return len(self.events) > 0

    def OnSelChanged(self, event):
        item = event.GetItem()
        proxy = self.tree.GetItemData(item)
        if proxy is not None:
            data = proxy.GetData()
            self.display.SetLabel(data)


class AbinitEventsFrame(awx.Frame):
    """
    Frame with an EventsPanel
    """
    def __init__(self, parent, filepath, **kwargs):
        filepath = os.path.abspath(filepath)
        title = "Abinit Events: %s" % os.path.basename(filepath)

        super(AbinitEventsFrame, self).__init__(parent, -1, title=title, **kwargs)

        self.event_panel = AbinitEventsPanel(self, filepath)


class AbiOutLogDirCtrl(wx.GenericDirCtrl):

    def __init__(self, *args, **kwargs):
        if "filter" not in kwargs:
            kwargs["filter"] = "All files (*.*)|*.*|ABINIT output files (*.abo)|*.out| ABINIT log files (*.abl)|*.ablog"
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
    def __init__(self, parent, filenames, **kwargs):

        if "title" not in kwargs:
            kwargs["title"] = "Abinit Events"

        super(AbinitEventsNotebookFrame, self).__init__(parent, **kwargs)

        # Here we create a panel and a notebook on the panel
        panel = awx.Panel(self)

        import wx.lib.agw.flatnotebook as fnb

        nb = fnb.FlatNotebook(panel)

        # Add the pages to the notebook with the name to show on the tab
        # Add only files for which we have events.
        for fname in filenames:
            page = AbinitEventsPanel(nb, fname)
            if page.has_events:
                nb.AddPage(page, text=os.path.basename(fname))

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

    class AbiEventsViewerApp(wx.App):
        def OnInit(self):
            frame = AbinitEventsNotebookFrame(None, filenames)
            self.SetTopWindow(frame)
            frame.Show()
            return True

    return AbiEventsViewerApp()


class AbinitTimerFrame(awx.Frame):
    """
    Frame with controls to plot the timing data.
    """
    def __init__(self, parent, filepath, **kwargs):
        filepath = os.path.abspath(filepath)
        title = "Abinit Timer: %s" % os.path.basename(filepath)
        super(AbinitTimerFrame, self).__init__(parent, -1, title=title, **kwargs)

        try:
            abifile = abiopen(filepath)
            self.timer_data = abifile.timer_data

        except Exception as exc:
            raise awx.Error(str(exc))

        if not self.timer_data:
            raise awx.Error("%s does not contain a valid ABINIT TIMER section!" % filepath)

        self.BuildUi()

    def BuildUi(self):

        # Set callbacks (bound methods of AbiTimerData).
        timer_data = self.timer_data
        self.plot_types = OrderedDict([
            ("pie", timer_data.show_pie),
            ("stacked_hist", timer_data.show_stacked_hist),
        ])

        keys = AbinitTimerSection.NUMERIC_FIELDS

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.plot_label = wx.StaticText(self, -1, "plot type:", wx.DefaultPosition, wx.DefaultSize, 0)
        self.plot_label.Wrap(-1)
        hsizer.Add(self.plot_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.plot_cbox = wx.ComboBox(self, -1, "pie", wx.DefaultPosition, wx.DefaultSize, self.plot_types.keys(),
                                     0)
        hsizer.Add(self.plot_cbox, 0, wx.ALL, 5)

        self.key_label = wx.StaticText(self, -1, "key:", wx.DefaultPosition, wx.DefaultSize, 0)
        self.key_label.Wrap(-1)
        hsizer.Add(self.key_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.key_cbox = wx.ComboBox(self, -1, "wall_time", wx.DefaultPosition, wx.DefaultSize, keys, 0)
        hsizer.Add(self.key_cbox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.plot_button = wx.Button(self, -1, "Plot", wx.DefaultPosition, wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.OnClick, self.plot_button)
        main_sizer.Add(self.plot_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

    def OnClick(self, event):
        callback = self.plot_types[self.plot_cbox.GetValue()]
        kwargs = dict(
            key=str(self.key_cbox.GetValue())
        )
        callback(**kwargs)

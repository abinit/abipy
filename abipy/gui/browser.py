#!/usr/bin/env python
from __future__ import print_function, division

import os
import wx
import warnings
import abipy.gui.awx as awx

import wx.lib.mixins.listctrl as listmix
from abipy.tools.text import list_strings, is_string, WildCard
from collections import namedtuple

try:
    from wxmplot import PlotApp, PlotFrame
except ImportError:
    #warnings.warn("Error while importing wxmplot. Some features won't be available")
    pass

from abipy.electrons import ElectronBandsPlotter
from abipy.gui.popupmenus import popupmenu_for_filename


def viewerframe_from_filepath(parent, filepath):
    """
    Factory function that returns the viewer (wx frame) associated to the file.
    None if no viewer has been registered for this filename.
    """
    from abipy.gui.wfkviewer import WfkViewerFrame
    from abipy.gui.sigresviewer import SigresViewerFrame
    from abipy.gui.gsrviewer import GsrViewerFrame
    from abipy.gui.mdfviewer import MdfViewerFrame
    from abipy.gui.editor import MyEditorFrame
    from abipy.xwncview import NcViewerFrame

    VIEWER_FRAMES = {
        "WFK-etsf.nc": WfkViewerFrame,
        "SIGRES.nc": SigresViewerFrame,
        "GSR.nc": GsrViewerFrame,
        "MDF.nc": MdfViewer.Frame,
        ".abi": MyEditorFrame, 
        ".abo": MyEditorFrame, 
        ".log": MyEditorFrame, 
        ".sh": MyEditorFrame,  
        ".err": MyEditorFrame, 
        ".files": MyEditorFrame, 
    }

    ext = filepath.split("_")[-1]
    try:
        return VIEWER_FRAMES[ext](parent, filepath)

    except KeyError:
        root, ext = os.path.splitext(filepath)
        try:
            return VIEWER_FRAMES[ext](parent, filepath)
        except KeyError:
            # No frame registered for the file. 
            # Open NcViewer if we have a netcdf file else None
            if filepath.endswith(".nc"): return NcViewerFrame(parent, filepath)
            return None


class NcFileDirCtrl(wx.GenericDirCtrl):
    def __init__(self, *args, **kwargs):

        if "filter" not in kwargs:
            kwargs["filter"] = "Netcdf files (*.nc)|*.nc|All files (*.*)|*.*"

        if "dir" not in kwargs:
            kwargs["dir"] = os.getcwd()

        if "style" not in kwargs:
            kwargs["style"] = wx.TR_MULTIPLE

        super(NcFileDirCtrl, self).__init__(*args, **kwargs)

        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.OnItemActivated)
        self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClick)

    def OnItemActivated(self, event):
        path = self.GetFilePath()
        if not path: return

        frame = viewerframe_from_filepath(self, path)
        if frame is not None:
            frame.Show()

    def OnRightClick(self, event):
        path = self.GetFilePath()
        if not path: return

        # Open the popup menum then destroy it to avoid mem leak.
        popmenu = popupmenu_for_filename(self, path)
        self.PopupMenu(popmenu, event.GetPoint())
        popmenu.Destroy()


class MyListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin):
    """ Mixin class to resize the last column appropriately."""
    def __init__(self, parent, **kwargs):
        wx.ListCtrl.__init__(self, parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)
        listmix.ListCtrlAutoWidthMixin.__init__(self)


class FileDataObj(namedtuple("FileDataObj", "filename type directory")):
    """The fields of the row in `FileListPanel`"""
    @property
    def abspath(self):
        return os.path.join(self.directory, self.filename)
                                                                         
    @classmethod
    def from_abspath(cls, abspath):
        return cls(
            filename=os.path.basename(abspath),
            type=os.path.splitext(abspath)[-1],
            directory=os.path.dirname(abspath),
        )


class FileListPanel(awx.Panel, listmix.ColumnSorterMixin):
    """
    This panel shows a list of files (strings), supports column sorting
    and provides specific popup menus for the different type of files 
    (file type detection is based on file extensions).
    """
    def __init__(self, parent, filepaths, **kwargs):
        """
        Args:
            parent:
                parent window
            filepaths:
                List of file paths.
        """
        super(FileListPanel, self).__init__(parent, -1, **kwargs)

        if filepaths is not None and is_string(filepaths):
            filepaths = [filepaths]

        self.filepaths = filepaths if filepaths is not None else []
        self.filepaths = map(os.path.abspath, self.filepaths)

        self.BuildUi()

    def BuildUi(self):
        self.file_list = file_list = MyListCtrl(self)

        file_list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        file_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

        self.FileDataObj = FileDataObj

        self.id2filedata = {}
        for (index, colname) in enumerate(self.FileDataObj._fields):
            file_list.InsertColumn(index, colname)

        for filepath in self.filepaths:
            self.AppendFilepath(filepath)

        for (index, colname) in enumerate(self.FileDataObj._fields):
            file_list.SetColumnWidth(index, wx.LIST_AUTOSIZE)

        # Now that the list exists we can init the other base class, see wx/lib/mixins/listctrl.py
        self.itemDataMap = self.id2filedata
        listmix.ColumnSorterMixin.__init__(self, len(self.FileDataObj._fields))
        self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.file_list)

        # Pack
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(file_list, 1, wx.ALL | wx.EXPAND, 5)

        self.SetSizerAndFit(sizer)

    def HasAbsPath(self, abspath):
        """True if abspath is already present."""
        return abspath in [f.abspath for f in self.id2filedata.values()]

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self.file_list

    def AppendFilepath(self, abspath):
        """Add a file to the panel."""
        if self.HasAbsPath(abspath):
            awx.showErrorMessage(self, message="%s is already in the list" % abspath)

        self._AppendFilepath(abspath)

    def _AppendFilepath(self, abspath):
        next = len(self.id2filedata)
        # We use next as entry id because we want to be able 
        # to sort the columns with the mixin ColumnSorterMixin.
        entry_id = next
        entry = self.FileDataObj.from_abspath(abspath)
        self.file_list.Append(entry)
        self.file_list.SetItemData(next, entry_id)
        self.id2filedata[entry_id] = entry

    def OnItemActivated(self, event):
        currentItem = event.m_itemIndex
        fd = self.id2filedata[self.file_list.GetItemData(currentItem)]

        frame = viewerframe_from_filepath(self, fd.abspath)
        if frame is not None:
            frame.Show()

    def OnRightClick(self, event):
        currentItem = event.m_itemIndex
        if currentItem == -1:
            return

        fd = self.id2filedata[self.file_list.GetItemData(currentItem)]
        # Open the popup menu then destroy it to avoid mem leak.
        menu = popupmenu_for_filename(self, fd.abspath)

        if menu is not None:
            self.PopupMenu(menu, event.GetPoint())
            menu.Destroy()

    def OnColClick(self, event):
        event.Skip()



class FileListFrame(awx.Frame):
    def __init__(self, parent, dirpaths=None, filepaths=None, walk=True, wildcard="", **kwargs):
        """
        Args:
            parent:
                parent window
            dirpaths:
                List of directories to scan.
            filepaths
                List of filepaths (absolute paths).
            walk:
                True if we have to browse all files and directories starting from filepaths.
            wildcard
                Regular expressions for selecting files (tokens are separated by |).
        """
        super(FileListFrame, self).__init__(parent, -1, **kwargs)

        if dirpaths is not None:
            dirpaths = map(os.path.abspath, list_strings(dirpaths))
        else:
            dirpaths = []

        if filepaths is not None:
            filepaths = map(os.path.abspath, list_strings(filepaths))
        else:
            filepaths = []

        wildcard = WildCard(wildcard)

        self.all_filepaths = filepaths

        if walk:
            for dirpath in dirpaths:
                for root, dirnames, filenames in os.walk(dirpath):
                    fnames = [os.path.join(root, f) for f in filenames]
                    self.all_filepaths += wildcard.filter(fnames)
        else:
            # Select only the files in dirpaths.
            for dirpath in dirpaths:
                fnames = map(os.path.isfile, os.listdir(dirpath))
                self.all_filepaths += wildcard.filter(fnames)

        self.BuildUi()

    def BuildUi(self):
        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        panel = FileListPanel(self, filepaths=self.all_filepaths)
        main_sizer.Add(panel, 1, wx.EXPAND, 5)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        filter_label = wx.StaticText(self, -1, "Filter:", wx.DefaultPosition, wx.DefaultSize, 0)
        filter_label.Wrap(-1)
        hsizer.Add(filter_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        # Combobox to enter shell patterns.
        wildcard_choices = ["*", "*.nc", "*.abo", "*.log"]
        self.filter_combobox = wx.ComboBox(self, id=-1, value="*", style=wx.TE_PROCESS_ENTER, choices=wildcard_choices)
        self.filter_combobox.SetToolTipString("Shell patterns separated by |")

        self.filter_combobox.Bind(wx.EVT_COMBOBOX, self.OnFilterComboBox)
        self.filter_combobox.Bind(wx.EVT_TEXT_ENTER, self.OnFilterComboBox)

        hsizer.Add(self.filter_combobox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizer(main_sizer)
        self.Layout()

    def OnFilterComboBox(self, event):
        wildcard = WildCard(self.filter_combobox.GetValue())

        select_files = wildcard.filter(self.all_filepaths)
        panel = FileListPanel(self, filepaths=select_files)

        main_sizer = self.main_sizer
        main_sizer.Hide(0)
        main_sizer.Remove(0)
        main_sizer.Insert(0, panel, 1, wx.EXPAND, 5)

        self.Layout()
        #self.Fit()


def wxapp_dirbrowser(dirpath):
    if dirpath is None:
        dirpath = " "
    else:
        dirpath = os.path.abspath(dirpath)

    app = awx.App()
    frame = awx.Frame(None, -1)
    NcFileDirCtrl(frame, -1, dir=dirpath)
    frame.Show()
    return app


def wxapp_listbrowser(dirpaths=None, filepaths=None, wildcard=""):
    app = awx.App()
    frame = FileListFrame(None, dirpaths=dirpaths, filepaths=filepaths, wildcard=wildcard)
    frame.Show()
    return app

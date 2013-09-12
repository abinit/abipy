#!/usr/bin/env python
from __future__ import print_function, division

import os
import wx
import fnmatch
import abipy.gui.awx as awx
import wx.lib.mixins.listctrl as listmix

from collections import namedtuple

try:
    from wxmplot import PlotApp, PlotFrame
except ImportError:
    import warnings
    warnings.warn("Error while importing wxmplot. Some features won't be available")

from abipy.electrons import ElectronBandsPlotter
from abipy.gui.popupmenus import popupmenu_for_filename
from abipy.gui.wfkviewer import WfkViewerFrame
from abipy.gui.sigresviewer import SigresViewerFrame


_VIEWER_FRAMES = {
    "WFK-etsf.nc": WfkViewerFrame,
    "SIGRES.nc": SigresViewerFrame,
}

_FRAME_SIZE = (720, 720)

def is_string(obj):
    """True if obj is a string."""
    try:
        dummy = obj + " "
        return True

    except TypeError:
        return False


def viewerframe_from_filepath(parent, filepath):
    """
    Returns the viewer (wx frame) associated to the file.
    None if no viewer has been registered.
    """
    ext = filepath.split("_")[-1]
    try:
        return _VIEWER_FRAMES[ext](parent, filepath)
    except KeyError:
        #print("No WX viewer has been registered for the extension %s" % ext)
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
    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN)
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

    def __init__(self, parent, dirpaths=None, filepaths=None, wildcard="", **kwargs):
        super(FileListPanel, self).__init__(parent, -1, **kwargs)

        if dirpaths is not None and is_string(dirpaths):
            dirpaths = [dirpaths, ]

        if filepaths is not None and is_string(filepaths):
            filepaths = [filepaths, ]

        self.dirpaths = dirpaths if dirpaths is not None else []
        self.filepaths = filepaths if filepaths is not None else []

        self.dirpaths = map(os.path.abspath, self.dirpaths)
        self.filepaths = map(os.path.abspath, self.filepaths)

        if not wildcard:
            self.wildcards = ["*.*"]
        else:
            if "|" in wildcard:
                self.wildcards = wildcard.split("|")
            else:
                self.wildcards = [wildcard]

        self.BuildUi()

    def BuildUi(self):
        self.file_list = file_list = MyListCtrl(self)

        file_list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        file_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

        self.FileDataObj = FileDataObj

        self.id2filedata = {}
        for (index, colname) in enumerate(self.FileDataObj._fields):
            file_list.InsertColumn(index, colname)

        for dirpath in self.dirpaths:
            self.ScanDirectory(dirpath)

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
        self.SetSizer(sizer)

    def HasAbsPath(self, abspath):
        return abspath in [f.abspath for f in self.id2filedata.values()]

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self.file_list

    def AcceptFilename(self, filename):
        for wcard in self.wildcards:
            if not fnmatch.fnmatch(filename, wcard):
                return False
        return True

    def ScanDirectory(self, dirpath):
        for root, dirnames, filenames in os.walk(dirpath):
            wnames = [f for f in filenames if self.AcceptFilename(f)]
            for apath in wnames:
                self._AppendFilepath(os.path.join(root, apath))

    def AppendFilepath(self, abspath):
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

    #def OnFilePicker(self, event):
    #    new_filepath = self.filepicker.GetPath()
    #    self.AppendFilepath(new_filepath)

    def OnItemActivated(self, event):
        currentItem = event.m_itemIndex
        fd = self.id2filedata[self.file_list.GetItemData(currentItem)]
        self.log("In OnItemActivated with filedata %s" % str(fd))

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
    pass


def wxapp_dirbrowser(dirpath):
    if dirpath is None:
        dirpath = " "
    else:
        dirpath = os.path.abspath(dirpath)

    app = wx.App()
    frame = awx.Frame(None, -1)
    NcFileDirCtrl(frame, -1, dir=dirpath)
    app.SetTopWindow(frame)
    frame.Show()
    return app


def wxapp_listbrowser(dirpaths=None, filepaths=None, wildcard=""):
    app = wx.App()
    frame = FileListFrame(None, -1, size=_FRAME_SIZE)
    FileListPanel(frame, dirpaths=dirpaths, filepaths=filepaths, wildcard=wildcard)
    frame.Show()
    return app

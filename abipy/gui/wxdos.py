#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import wx

import wx.lib.dialogs as wxdg 
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy.abifiles import abiopen
from abipy.electrons import ElectronBands
from abipy.waves import WFK_File
from wx.lib.agw.floatspin import FloatSpin
from wxmplot import PlotApp, PlotFrame

from abipy.gui.popupmenus import popupmenu_from_ext

class DosPlotter(wx.Frame):
    VERSION = "0.1"

    # Help menu items.
    ID_ABOUT = wx.NewId() 

    def __init__(self, parent=None, **kwargs):
        super(DosPlotter, self).__init__(parent, -1, self.codename, size=(800,600)) 

        # Create statusbar
        self.statusbar = self.CreateStatusBar() 

        # Setup menu bar.
        menuBar = wx.MenuBar()

        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_OPEN,  "&Open", help="Open an existing WFK file")
        fileMenu.Append(wx.ID_CLOSE, "&Close", help="Close the WFK file")
        fileMenu.Append(wx.ID_EXIT,  "&Quit", help="Exit the application")
        menuBar.Append(fileMenu, "File")

        filehistory = self.filehistory = wx.FileHistory(8)
        self.config = wx.Config(self.codename, style=wx.CONFIG_USE_LOCAL_FILE)
        filehistory.Load(self.config)
        recent = wx.Menu()
        filehistory.UseMenu(recent)
        filehistory.AddFilesToMenu()
        fileMenu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.on_file_history, id=wx.ID_FILE1, id2=wx.ID_FILE9)

        self.helpMenu = wx.Menu()
        self.helpMenu.Append(self.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.helpMenu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        tsize = (15,15)
        self.toolbar = toolbar = self.CreateToolBar()

        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        #toolbar.AddSimpleTool(self.ID_VISTRUCT, wx.Bitmap(awx.path_img("crystal.png")), "Visualize the crystal structure")
        #toolbar.AddSimpleTool(self.ID_VISWAVE, wx.Bitmap(awx.path_img("wave.png")), "Visualize the selected wavefunction")

        self.toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menuHandlers = [
            (wx.ID_OPEN,        self.onOpen),
            (wx.ID_CLOSE,       self.onClose),
            (wx.ID_EXIT,        self.onExit),
            (self.ID_ABOUT,     self.onAboutBox),
        ]

        for combo in menuHandlers:
            id, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=id)

        # Widgets.
        self.panel = panel = FileListPanel(self)
        self.panel.Show()

    @property
    def codename(self):
        return self.__class__.__name__

    def onOpen(self, event):
        dlg = wx.FileDialog(self, message="Choose a Netcdf file", defaultDir=os.getcwd(), 
            defaultFile="", wildcard="Netcdf files (*.nc)|*nc",
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. 
        # If it is the OK response, process the data.
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.statusbar.PushStatusText("Reading %s" % path)
            self.filehistory.AddFileToHistory(path)
            self.filehistory.Save(self.config)
            self.config.Flush()
            self.read_bands(path)

        dlg.Destroy()

    def on_file_history(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        path = self.filehistory.GetHistoryFile(fileNum)
        # move up the list
        self.filehistory.AddFileToHistory(path)  
        self.read_bands(path)

    def read_bands(self, path):
        try:
            self.bands = ElectronBands.from_ncfile(path)
            ewx.ElectronDosFrame(bands=self.bands, parent=self).Show()
        except Exception:
            awx.showErrorMessage(self)

    def onClose(self, event):
        self.bands = None

    def onExit(self, event):
        self.Destroy()

    def onAboutBox(self, event):
        awx.makeAboutBox(codename=self.codename, version=self.VERSION, 
                      description="", developers="M. Giantomassi")


class FileListPanel(wx.Panel):

    def __init__(self, parent, paths=(), **kwargs):
        super(FileListPanel, self).__init__(parent, -1, **kwargs)

        self.file_list = file_list = wx.ListCtrl(self, id=-1, size=(-1,100), 
                                     style=wx.LC_REPORT | wx.BORDER_SUNKEN
                                     )
        file_list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.onRightClick)

        self.ncfiles_by_id = {}
        file_list.InsertColumn(0, "filename")
        #file_list.InsertColumn(1, "filetype")
        for index, path in enumerate(paths):
            self.append_path_to_filelist(path)

        self.filepicker = wx.FilePickerCtrl(self, id=-1, path=os.getcwd(),
            wildcard="Netcdf files (*.nc)|*nc", style = wx.FLP_OPEN | wx.CHANGE_DIR)

        self.Bind(wx.EVT_FILEPICKER_CHANGED, self.onFilePicker)

        # Pack
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(file_list, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(self.filepicker, 0, wx.ALL | wx.CENTER, 5)
        self.SetSizer(sizer)

    def append_path_to_filelist(self, path):
        next = len(self.ncfiles_by_id)
        try:
            #ncfile = abiopen(path)
            #file_type = ncfile.__class__.__name__
            #file_type = ncfile.filetype
            ncid = id(path)
            self.ncfiles_by_id[ncid] = path
            #entry = [os.path.basename(path), file_type]
            entry = [os.path.basename(path)]
            self.file_list.Append(entry)
            self.file_list.SetItemData(next, ncid)
        except:
            awx.showErrorMessage(self)

    def onFilePicker(self, event):
        self.append_path_to_filelist(self.filepicker.GetPath())

    def onRightClick(self, event):
        currentItem = event.m_itemIndex

        if currentItem != -1:
            ncfile = self.ncfiles_by_id[self.file_list.GetItemData(currentItem)]
            print("Select ncfile: ",ncfile)

            menu = popupmenu_from_ext(ncfile)

            # Open the popup menum then destroy it to avoid mem leak.
            self.PopupMenu(menu, event.GetPoint())
            menu.Destroy() 


def popupmenu_from_ext(filepath):
    """
    Factory function that returns the appropriate popmenu menu
    by testing the file extension.
    """
    menu = BasePopupMenu()
    menu.add_target(ncfilepath)
    return menu


def showNcdumpMessage(parent, filepath):
    caption = "ncdump output for file %s" % filepath
    from abipy import abiopen
    ncfile = abiopen(filepath)
    wxdg.ScrolledMessageDialog(parent, ncfile.ncdump(), caption=caption, style=wx.MAXIMIZE_BOX).Show()


class BasePopupMenu(wx.Menu):
    MENU_TITLES = {
        #("Properties",
        "Ncdump": showNcdumpMessage,
        "DOS": ewx.showElectronDosFrame,
    }

    def __init__(self, *args, **kwargs):
        super(BasePopupMenu, self).__init__()

        menu_title_by_id = {}
        for title in self.MENU_TITLES:
            menu_title_by_id[wx.NewId()] = title

        if not hasattr(self, "menu_title_by_id"):
            self.menu_title_by_id = {}

        self.menu_title_by_id.update(menu_title_by_id)

        for (id, title) in self.menu_title_by_id.items():
            self.Append(id, title)
            # registers menu handlers with EVT_MENU, on the menu.
            wx.EVT_MENU(self, id, self.OnMenuSelection)

    def add_target(self, target):
        self._target = target

    @property
    def target(self):
        try:
            return self._target
        except AttributeError:
            return None

    def OnMenuSelection(self, event):
        menu_title = self.menu_title_by_id[event.GetId()]
        operation = self.MENU_TITLES[menu_title]
        print("Perform operation %s on target %s" % (operation, self.target))
        try:
            #ewx.showElectronDosFrame(parent=None, path=self.target)
            operation(parent=None, filepath=self.target)
        except:
            awx.showErrorMessage(parent=None)


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


class AbipyDirCtrl(wx.GenericDirCtrl):
    def __init__(self, *args, **kwargs):
        #style = wx.TR_MULTIPLE
        super(AbipyDirCtrl, self).__init__(*args, **kwargs)

        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.OnItemActivated)
        self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClick)
        #self.Bind(wx.EVT_TREE_SEL_CHANGED, self.dirBrowser_OnSelectionChanged, tree)

    def OnItemActivated(self, event):
        path = self.GetFilePath()
        if not path:
            return
        print("in activated with path %s" % path)
        run_viewer_for_filepath(path)

    def OnRightClick(self, event):
        path = self.GetFilePath()
        if not path:
            return
        print("in right with path %s" % path)

def dos_plotter():
    app = wx.App()
    win = DosPlotter()
    win.Show(True)
    app.MainLoop()

def abi_navigator():
    app = wx.App()
    win = wx.Frame(None, -1)
    filebrowser = AbipyDirCtrl(win, -1, dir=os.getcwd(), filter="Netcdf files (*.nc)|*.nc|All files (*.*)|*.*", name="Netcdf Navigator")
    #filebrowser = wx.GenericDirCtrl(win, -1, dir=os.getcwd(), filter="Netcdf files (*.nc)|*.nc|All files (*.*)|*.*", name="Netcdf Navigator")
    #filebrowser.SetDefaultPath("/Users/gmatteo/Coding/abipy/abipy/tests/data")
    #filebrowser.SetPath("/Users/gmatteo/Coding/abipy/abipy/tests/data")
    win.Show(True)
    app.MainLoop()

if __name__ == "__main__":
    abi_navigator()

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

from abipy.gui.popupmenus import popupmenu_for_filename

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
            self.bands = ElectronBands.from_file(path)
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


def dos_plotter():
    app = wx.App()
    win = DosPlotter()
    win.Show(True)
    app.MainLoop()


if __name__ == "__main__":
    dos_plotter()

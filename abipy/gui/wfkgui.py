#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import traceback
import wx

from abipy.waves import WFK_File
from abipy.iotools.visualizer import supported_visunames
from abipy.gui.tools import build_aboutbox, build_errormessage, straceback, path_img

class WfkGuiApp(wx.App):

    def OnInit(self): 
        # The code for the splash screen.
        #image = wx.Image(path_img("wabi_logo.png"), wx.BITMAP_TYPE_PNG)
        #bmp = image.ConvertToBitmap() 
        #wx.SplashScreen(bmp, wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT, 1000, None, -1) 
        #wx.Yield()

        frame = WfkGui(parent=None) 
        frame.Show(True) 
        self.SetTopWindow(frame) 
        return True


class WfkGui(wx.Frame):
    VERSION = "0.1"
    # Toolbar items.
    ID_VISTRUCT = wx.NewId() 
    ID_VISWAVE  = wx.NewId() 

    def __init__(self, *args, **kwargs):
        super(WfkGui, self).__init__(None, -1, self.codename, size=(800,600)) 

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
        self.helpMenu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.helpMenu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        tsize = (15,15)
        self.toolbar = toolbar = self.CreateToolBar()

        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(self.ID_VISTRUCT, wx.Bitmap(path_img("crystal.png")), "Visualize the crystal structure")
        toolbar.AddSimpleTool(self.ID_VISWAVE, wx.Bitmap(path_img("wave.png")), "Visualize the selected wavefunction")

        self.toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menuHandlers = [
            (wx.ID_OPEN,        self.onOpen),
            (wx.ID_CLOSE,       self.onClose),
            (wx.ID_EXIT,        self.onExit),
            (wx.ID_ABOUT,       self.onAboutBox),
            (self.ID_VISTRUCT,  self.onVisualizeStructure),
            (self.ID_VISWAVE,   self.onVisualizeWave),
        ]

        for combo in menuHandlers:
            id, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=id)

        self.wfk_file = None 

    @property
    def codename(self):
        return self.__class__.__name__

    def destroy_panel(self):
        if hasattr(self, "panel"):
            self.panel.Destroy()

    def build_panel(self):
        if self.wfk_file is None: return
        self.destroy_panel()

        wfk_file = self.wfk_file

        # Widgets.
        self.panel = panel = wx.Panel(self, -1)
        sizer = wx.FlexGridSizer(rows=4, cols=2, vgap=5, hgap=5)

        spin_label = wx.StaticText(panel, -1, "Spin index:")
        self.spin_cbox = wx.ComboBox(panel, id=-1, choices=map(str, range(wfk_file.nsppol)))
        sizer.AddMany([spin_label, self.spin_cbox])

        band_label = wx.StaticText(panel, -1, "Band index:")
        self.band_cbox = wx.ComboBox(panel, id=-1, choices=map(str, range(wfk_file.mband)))
        sizer.AddMany([band_label, self.band_cbox])

        kpt_label = wx.StaticText(panel, -1, "Kpoint:")

        self.kpt_list = kpt_list = wx.ListView(panel, id=-1)

        kpt_list.InsertColumn(0, 'Index')
        kpt_list.InsertColumn(1, 'Reduced Coords')
        kpt_list.InsertColumn(2, 'Weight')
        for index, kpt in enumerate(wfk_file.kpoints):
            kpt_list.InsertStringItem(index, str(index))
            kpt_list.SetStringItem(index, 1, str(kpt.frac_coords))
            kpt_list.SetStringItem(index, 2, str(kpt.weight))

        sizer.AddMany([kpt_label, self.kpt_list])

        visu_label = wx.StaticText(panel, -1, "Visualizer:")
        self.visualizer = wx.ComboBox(panel, id=-1, choices=supported_visunames())
        sizer.AddMany([visu_label, self.visualizer])

        panel.SetSizer(sizer)
        sizer.Layout()

    def onOpen(self, event):
        dlg = wx.FileDialog(self, message="Choose a WFK file", defaultDir=os.getcwd(), 
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
            self.read_wfkfile(path)

        dlg.Destroy()

    def on_file_history(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        path = self.filehistory.GetHistoryFile(fileNum)
        # move up the list
        self.filehistory.AddFileToHistory(path)  
        self.read_wfkfile(path)

    def read_wfkfile(self, path):
        try:
            self.wfk_file = WFK_File(path)
            self.build_panel()
        except Exception:
            build_errormessage(self, straceback())

    def onClose(self, event):
        self.wfk_file = None
        self.destroy_panel()

    def onExit(self, event):
        self.Destroy()

    def onAboutBox(self, event):
        build_aboutbox(codename=self.codename, version=self.VERSION, 
                       description="", developers="M. Giantomassi")

    def onVisualizeStructure(self, event):
        if self.wfk_file is None: return

        visualizer = self.visualizer.GetValue()

        self.statusbar.PushStatusText("Visualizing crystal structure with %s" % visualizer)
        try:
            visu = self.wfk_file.visualize_structure_with(visualizer)
            visu()
        except:
            build_errormessage(self, straceback())

    def onVisualizeWave(self, event):
        if self.wfk_file is None: return

        spin = int(self.spin_cbox.GetValue())
        band = int(self.band_cbox.GetValue())
        kidx  = self.kpt_list.GetFirstSelected()
        kpoint = self.wfk_file.kpoints[kidx]

        visualizer = self.visualizer.GetValue()

        self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
        try:
            visu = self.wfk_file.export_ur2(".xsf", spin, kpoint, band)
            visu()
        except:
            build_errormessage(self, straceback())


def wfkgui_main():
    """ Start up the WfkGui application."""
    app = WfkGuiApp()
    app.MainLoop()

if __name__ == "__main__":
    wfkgui_main()


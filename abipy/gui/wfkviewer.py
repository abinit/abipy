#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import wx

import wx.lib.dialogs as wxdg 
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy import abiopen, WFK_File
from abipy.iotools.visualizer import supported_visunames

class WfkViewerFrame(wx.Frame):
    VERSION = "0.1"
    # Toolbar items.
    ID_VISTRUCT  = wx.NewId() 
    ID_VISWAVE   = wx.NewId() 
    ID_VISBZ     = wx.NewId() 
    ID_DOS       = wx.NewId() 
    ID_SUMMARY   = wx.NewId()
    ID_NCDUMP    = wx.NewId()
    ID_PLOTBANDS = wx.NewId()

    def __init__(self, parent, filename=None):
        super(WfkViewerFrame, self).__init__(parent, -1, self.codename) 

        # Create statusbar
        self.statusbar = self.CreateStatusBar() 

        # Setup menu bar.
        menuBar = wx.MenuBar()

        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_OPEN, "&Open", help="Open an existing WFK file")
        fileMenu.Append(wx.ID_CLOSE,"&Close", help="Close the WFK file")
        fileMenu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")
        fileMenu.Append(self.ID_NCDUMP, "Ncdump", help="ncdump printout")
        menuBar.Append(fileMenu, "File")

        filehistory = self.filehistory = wx.FileHistory(8)
        self.config = wx.Config(self.codename, style=wx.CONFIG_USE_LOCAL_FILE)
        filehistory.Load(self.config)
        recent = wx.Menu()
        filehistory.UseMenu(recent)
        filehistory.AddFilesToMenu()
        fileMenu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)

        self.helpMenu = wx.Menu()
        self.helpMenu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.helpMenu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        self.toolbar = toolbar = self.CreateToolBar()
        artBmp = wx.ArtProvider.GetBitmap

        tsize = (15,15)
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(self.ID_VISTRUCT, wx.Bitmap(awx.path_img("crystal.png")), "Visualize the crystal structure")
        toolbar.AddSimpleTool(self.ID_VISWAVE, wx.Bitmap(awx.path_img("wave.png")), "Visualize the selected wavefunction")
        toolbar.AddSimpleTool(self.ID_VISBZ, wx.Bitmap(awx.path_img("wave.png")), "Visualize the BZ")
        toolbar.AddSimpleTool(self.ID_DOS, wx.Bitmap(awx.path_img("wave.png")), "Compute DOS")
        toolbar.AddSimpleTool(self.ID_PLOTBANDS, wx.Bitmap(awx.path_img("wave.png")), "Plot bands")

        self.toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menuHandlers = [
            (wx.ID_OPEN,        self.OnOpen),
            (wx.ID_CLOSE,       self.OnClose),
            (wx.ID_EXIT,        self.OnExit),
            (wx.ID_ABOUT,       self.OnAboutBox),
            (self.ID_NCDUMP,    self.OnNcdump),
            #
            (self.ID_VISTRUCT,  self.OnVisualizeStructure),
            (self.ID_VISWAVE,   self.OnVisualizeWave),
            (self.ID_VISBZ,     self.OnVisualizeBZ),
            (self.ID_DOS,       self.OnDOS),
            (self.ID_PLOTBANDS, self.OnPlotBands),
            #(self.ID_SUMMARY,  self.OnSummary),
        ]

        for combo in menuHandlers:
            id, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=id)

        self.wfk = None 
        if filename is not None:
            print("in read",filename)
            self.read_wfkfile(filename)

    @property
    def codename(self):
        return self.__class__.__name__

    #def destroy_panel(self):
    #    if hasattr(self, "panel"):
    #        self.panel.Destroy()

    def build_UI(self):
        wfk = self.wfk
        if wfk is None: return
        #self.destroy_panel()

        self.skb_panel = awx.SpinKpointBandPanel(self, wfk.nsppol, wfk.kpoints, wfk.mband)

        #visu_label = wx.StaticText(panel, -1, "Visualizer:")
        #self.visualizer = wx.ComboBox(panel, id=-1, choices=supported_visunames())

        # Python shell
        #from wx.py.shell import Shell
        #pyshell = Shell(panel, locals={"wfk": self.wfk})
        #vsizer.Add(pyshell)

    def OnOpen(self, event):
        dlg = wx.FileDialog(self, message="Choose a WFK file", defaultDir=os.getcwd(), 
            defaultFile="", wildcard="Netcdf files (*.nc)|*nc",
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. 
        # If it is the OK response, process the data.
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.filehistory.AddFileToHistory(path)
            self.filehistory.Save(self.config)
            self.config.Flush()
            self.read_wfkfile(path)

        dlg.Destroy()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        path = self.filehistory.GetHistoryFile(fileNum)
        # move up the list
        self.filehistory.AddFileToHistory(path)  
        self.read_wfkfile(path)

    def read_wfkfile(self, path):
        self.statusbar.PushStatusText("Reading %s" % path)
        try:
            wfkfile = abiopen(path)
            if not isinstance(wfkfile, WFK_File):
                awx.showErrorMessage(self, message="%s is not a valid WFK File" % path)
                return
                
            self.wfk = wfkfile
            self.build_UI()
            self.statusbar.PushStatusText("WFK file %s loaded" % path)

        except Exception:
            awx.showErrorMessage(self)

    def OnClose(self, event):
        self.wfk= None
        self.destroy_panel()

    def OnExit(self, event):
        self.Destroy()

    def OnAboutBox(self, event):
        awx.makeAboutBox(codename=self.codename, version=self.VERSION, 
                         description="", developers="M. Giantomassi")

    def OnVisualizeStructure(self, event):
        if self.wfk is None: return

        visualizer = self.visualizer.GetValue()
        self.statusbar.PushStatusText("Visualizing crystal structure with %s" % visualizer)
        try:
            visu = self.wfk.visualize_structure_with(visualizer)
            visu()
        except:
            awx.showErrorMessage(self)

    def OnVisualizeBZ(self, event):
        if self.wfk is None: return
        self.wfk.structure.show_bz()

    def OnVisualizeWave(self, event):
        if self.wfk is None: return

        spin, kidx, band = self.skb_panel.get_skb()
        kpoint = self.wfk.kpoints[kidx]

        #visualizer = self.visualizer.GetValue()
        self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
        try:
            visu = self.wfk.export_ur2(".xsf", spin, kpoint, band)
            visu()
        except:
            awx.showErrorMessage(self)

    def OnDOS(self, event):
        if self.wfk is None: return
        ewx.ElectronDosFrame(bands=self.wfk.get_bands(), parent=self).Show()

    def OnPlotBands(self, event):
        if self.wfk is None: return
        self.wfk.get_bands().plot()

    def OnNcdump(self, event):
        if self.wfk is None: return
        caption = "ncdump output for WFK file %s" % self.wfk.filepath
        wxdg.ScrolledMessageDialog(self, self.wfk.ncdump(), caption=caption, style=wx.MAXIMIZE_BOX).Show()


class WfkViewerApp(awx.App):

    def OnInit(self):
        frame = WfkViewerFrame(None, filename=None)
        frame.Show()
        self.SetTopWindow(frame) 
        return True

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        awx.PRINT("%s dropped on app %s" % (filename, self.appname)) 
        frame = WfkViewerFrame(parent=None, filename=filename) 
        frame.Show()
        

def wxapp_wfkviewer(wfk_filename):
    app = WfkViewerApp()
    frame = WfkViewerFrame(None, filename=wfk_filename) 
    frame.Show() 
    app.SetTopWindow(frame) 
    return app


if __name__ == "__main__":
    import sys
    wfk_filename = None
    if len(sys.argv) > 1:
        wfk_filename = sys.argv[1] 
    wxapp_wfkviewer(wfk_filename).MainLoop()

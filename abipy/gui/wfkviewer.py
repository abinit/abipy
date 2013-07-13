#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import wx

import wx.lib.dialogs as wxdg 
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy.waves import WFK_File
from abipy.iotools.visualizer import supported_visunames

class WfkViewer(wx.Frame):
    VERSION = "0.1"
    # Toolbar items.
    ID_VISTRUCT = wx.NewId() 
    ID_VISWAVE  = wx.NewId() 
    ID_DOS      = wx.NewId() 
    ID_SUMMARY  = wx.NewId()
    ID_NCDUMP   = wx.NewId()

    def __init__(self, parent, filepath=None):
        super(WfkViewer, self).__init__(parent, -1, self.codename) 

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
        self.Bind(wx.EVT_MENU_RANGE, self.on_file_history, id=wx.ID_FILE1, id2=wx.ID_FILE9)

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
        toolbar.AddSimpleTool(self.ID_DOS, wx.Bitmap(awx.path_img("wave.png")), "Compute DOS")

        self.toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menuHandlers = [
            (wx.ID_OPEN,        self.onOpen),
            (wx.ID_CLOSE,       self.onClose),
            (wx.ID_EXIT,        self.onExit),
            (wx.ID_ABOUT,       self.onAboutBox),
            (self.ID_VISTRUCT,  self.onVisualizeStructure),
            (self.ID_VISWAVE,   self.onVisualizeWave),
            (self.ID_DOS,       self.onDOS),
            (self.ID_NCDUMP,    self.onNcdump),
            #(self.ID_SUMMARY,  self.onSummary),
        ]

        for combo in menuHandlers:
            id, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=id)

        self.wfk = None 
        if filepath is not None:
            self.read_wfkfile(filepath)

    @property
    def codename(self):
        return self.__class__.__name__

    def destroy_panel(self):
        if hasattr(self, "panel"):
            self.panel.Destroy()

    def build_panel(self):
        wfk = self.wfk
        if wfk is None: return

        self.destroy_panel()

        self.panel = panel = wx.Panel(self, -1)

        band_label = wx.StaticText(panel, -1, "Band:")
        self.band_cbox = wx.ComboBox(panel, id=-1, choices=map(str, range(wfk.mband)))

        spin_label = wx.StaticText(panel, -1, "Spin:")
        self.spin_cbox = wx.ComboBox(panel, id=-1, choices=map(str, range(wfk.nsppol)))

        kpt_label = wx.StaticText(panel, -1, "Kpoint:")
        self.kpt_list = kpt_list = wx.ListCtrl(panel, id=-1, size=(-1,100), 
                                     style=wx.LC_REPORT | wx.BORDER_SUNKEN
                                     )

        visu_label = wx.StaticText(panel, -1, "Visualizer:")
        self.visualizer = wx.ComboBox(panel, id=-1, choices=supported_visunames())

        kpt_list.InsertColumn(0, '#')
        kpt_list.InsertColumn(1, 'Reduced Coords')
        kpt_list.InsertColumn(2, 'Weight')
        for index, kpt in enumerate(wfk.kpoints):
            entry = map(str, [index, kpt.frac_coords, kpt.weight])
            kpt_list.Append(entry)

        # Make the first row by adding the label and field  to the first horizontal sizer 
        vsizer = wx.BoxSizer(wx.VERTICAL) 
        field1_sz = wx.BoxSizer(wx.HORIZONTAL) 
        field2_sz = wx.BoxSizer(wx.HORIZONTAL)

        field1_sz.AddSpacer(50) 
        field1_sz.Add(band_label)
        field1_sz.AddSpacer(5)
        field1_sz.Add(self.band_cbox) 

        field1_sz.Add(spin_label)
        field1_sz.AddSpacer(5)
        field1_sz.Add(self.spin_cbox) 

        field1_sz.Add(visu_label)
        field1_sz.AddSpacer(5)
        field1_sz.Add(self.visualizer) 

        field2_sz.AddSpacer(50)
        field2_sz.Add(kpt_label)
        field2_sz.AddSpacer(5)

        BOTH_SIDES = wx.EXPAND|wx.LEFT|wx.RIGHT 
        field2_sz.Add(self.kpt_list, 0, BOTH_SIDES, 50) 

        vsizer.AddSpacer(50) 
        vsizer.Add(field1_sz)
        vsizer.AddSpacer(15) 
        vsizer.Add(field2_sz) 
        vsizer.AddSpacer(50)

        # Python shell
        from wx.py.shell import Shell
        pyshell = Shell(panel, locals={"wfk": self.wfk})
        vsizer.Add(pyshell)

        panel.SetSizer(vsizer)

        #sizer = wx.FlexGridSizer(rows=4, cols=2, vgap=5, hgap=5)
        #sizer.AddMany([band_label, self.band_cbox])
        #sizer.AddMany([spin_label, self.spin_cbox])
        #sizer.AddMany([kpt_label, self.kpt_list])
        #sizer.AddMany([visu_label, self.visualizer])
        #panel.SetSizer(sizer)
        #sizer.Layout()

    def onOpen(self, event):
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

    def on_file_history(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        path = self.filehistory.GetHistoryFile(fileNum)
        # move up the list
        self.filehistory.AddFileToHistory(path)  
        self.read_wfkfile(path)

    def read_wfkfile(self, path):
        self.statusbar.PushStatusText("Reading %s" % path)
        try:
            self.wfk = WFK_File(path)
            self.build_panel()
            self.statusbar.PushStatusText("WFK file %s loaded" % path)
        except Exception:
            awx.showErrorMessage(self)

    def onClose(self, event):
        self.wfk= None
        self.destroy_panel()

    def onExit(self, event):
        self.Destroy()

    def onAboutBox(self, event):
        awx.makeAboutBox(codename=self.codename, version=self.VERSION, 
                         description="", developers="M. Giantomassi")

    def onVisualizeStructure(self, event):
        if self.wfk is None: return

        visualizer = self.visualizer.GetValue()
        self.statusbar.PushStatusText("Visualizing crystal structure with %s" % visualizer)
        try:
            visu = self.wfk.visualize_structure_with(visualizer)
            visu()
        except:
            awx.showErrorMessage(self)

    def onVisualizeWave(self, event):
        if self.wfk is None: return

        spin = int(self.spin_cbox.GetValue())
        band = int(self.band_cbox.GetValue())
        kidx  = self.kpt_list.GetFirstSelected()
        kpoint = self.wfk.kpoints[kidx]

        visualizer = self.visualizer.GetValue()
        self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
        try:
            visu = self.wfk.export_ur2(".xsf", spin, kpoint, band)
            visu()
        except:
            awx.showErrorMessage(self)

    def onDOS(self, event):
        if self.wfk is None: return
        ewx.ElectronDosFrame(bands=self.wfk.get_bands(), parent=self).Show()

    def onNcdump(self, event):
        if self.wfk is None: return
        caption = "ncdump output for WFK file %s" % self.wfk.filepath
        wxdg.ScrolledMessageDialog(self, self.wfk.ncdump(), caption=caption, style=wx.MAXIMIZE_BOX).Show()

    #def onSummary(self, event):
    #    if self.wfk is None: return
    #    summary = self.wfk.summary


def wfk_viewer(filepath):
    """Start up the WfkViewer application."""

    class WfkViewerApp(wx.App):
        def OnInit(self): 
            # The code for the splash screen.
            #image = wx.Image(awx.path_img("wabi_logo.png"), wx.BITMAP_TYPE_PNG)
            #bmp = image.ConvertToBitmap() 
            #wx.SplashScreen(bmp, wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT, 1000, None, -1) 
            #wx.Yield()
            frame = WfkViewer(parent=None,filepath=filepath) 
            frame.Show(True) 
            self.SetTopWindow(frame) 
            return True

    app = WfkViewerApp()
    app.MainLoop()

if __name__ == "__main__":
    import sys
    filepath = None
    if len(sys.argv) > 1:
        filepath = sys.argv[1] 
    wfk_browser(filepath)


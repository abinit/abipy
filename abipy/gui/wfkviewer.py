from __future__ import print_function, division

import os
import wx

import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy import abiopen, WFK_File
from abipy.iotools.visualizer import supported_visunames


ID_VISTRUCT = wx.NewId()
ID_VISWAVE = wx.NewId()
ID_VISBZ = wx.NewId()
ID_DOS = wx.NewId()
ID_JDOS = wx.NewId()
ID_NCDUMP = wx.NewId()
ID_PLOTBANDS = wx.NewId()
ID_TBOX_VIS = wx.NewId()

class WfkViewerFrame(awx.Frame):
    VERSION = "0.1"

    def __init__(self, parent, filename=None):
        super(WfkViewerFrame, self).__init__(parent, -1, self.codename)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing WFK file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the WFK file")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")
        file_menu.Append(ID_NCDUMP, "Ncdump", help="ncdump printout")
        menuBar.Append(file_menu, "File")

        file_history = self.file_history = wx.FileHistory(8)
        self.config = wx.Config(self.codename, style=wx.CONFIG_USE_LOCAL_FILE)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)

        self.help_menu = wx.Menu()
        self.help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.help_menu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        self.toolbar = toolbar = self.CreateToolBar()

        tsize = (15, 15)
        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(ID_VISTRUCT, wx.Bitmap(awx.path_img("struct.png")), "Visualize the crystal structure")
        toolbar.AddSimpleTool(ID_VISBZ, wx.Bitmap(awx.path_img("bz.png")), "Visualize the BZ")
        toolbar.AddSimpleTool(ID_VISWAVE, wx.Bitmap(awx.path_img("wfk.png")), "Visualize the selected wavefunction")
        toolbar.AddSimpleTool(ID_DOS, wx.Bitmap(awx.path_img("dos.png")), "Compute the DOS")
        toolbar.AddSimpleTool(ID_JDOS, wx.Bitmap(awx.path_img("jdos.png")), "Compute the joint DOS")
        toolbar.AddSimpleTool(ID_PLOTBANDS, wx.Bitmap(awx.path_img("bs.png")), "Plot bands")

        toolbar.AddSeparator()
        self.visualizer_cbox = wx.ComboBox(choices=supported_visunames(), id=ID_TBOX_VIS, 
            name='visualizer', parent=toolbar, value='xcrysden') 
        self.visualizer_cbox.Refresh() 

        toolbar.AddControl(control=self.visualizer_cbox) 

        self.toolbar.Realize()
        self.Centre()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_OPEN, self.OnOpen),
            (wx.ID_CLOSE, self.OnClose),
            (wx.ID_EXIT, self.OnExit),
            (wx.ID_ABOUT, self.OnAboutBox),
            #
            (ID_NCDUMP, self.OnNcdump),
            (ID_VISTRUCT, self.OnVisualizeStructure),
            (ID_VISWAVE, self.OnVisualizeWave),
            (ID_VISBZ, self.OnVisualizeBZ),
            (ID_DOS, self.OnDos),
            (ID_JDOS, self.OnJdos),
            (ID_PLOTBANDS, self.OnPlotBands),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.wfk = None
        if filename is not None:
            self.ReadWfkFile(filename)

    @property
    def codename(self):
        return self.__class__.__name__

    def BuildUi(self):
        wfk = self.wfk
        if wfk is None: return

        splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
        splitter.SetMinimumPaneSize(50)
        parent = splitter 
        #parent = self

        self.skb_panel = awx.SpinKpointBandPanel(parent, wfk.nsppol, wfk.kpoints, wfk.mband)

        # Set the callback for double click on k-point row..
        self.skb_panel.SetOnItemActivated(self._visualize_skb)

        # Add Python shell
        from wx.py.shell import Shell
        from abipy.tools import marquee
        msg = "WFK_File object is accessible via the wfk variable. Use wfk.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        pyshell = Shell(parent, introText=msg, locals={"wfk": self.wfk})
        splitter.SplitHorizontally(self.skb_panel, pyshell)

        #main_sizer = wx.BoxSizer(wx.VERTICAL)
        #main_sizer.Add(self.skb_panel, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)
        #main_sizer.Add(pyshell, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)
        #self.SetSizerAndFit(main_sizer)

    def DestroyPanel(self):
        if hasattr(self, "panel"):
            self.panel.Destroy()

    def OnOpen(self, event):
        dlg = wx.FileDialog(self, message="Choose a WFK file", defaultDir=os.getcwd(),
                            wildcard="Netcdf files (*.nc)|*nc",
                            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
        )

        # Show the dialog and retrieve the user response. 
        # If it is the OK response, process the data.
        if dlg.ShowModal() == wx.ID_OK:
            filepath = dlg.GetPath()
            self.file_history.AddFileToHistory(filepath)
            self.file_history.Save(self.config)
            self.config.Flush()
            self.ReadWfkFile(filepath)

        dlg.Destroy()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        # move up the list
        self.file_history.AddFileToHistory(filepath)
        self.ReadWfkFile(filepath)

    def ReadWfkFile(self, filepath):
        """Read the WFK file and build the UI."""
        self.statusbar.PushStatusText("Reading %s" % filepath)

        try:
            wfkfile = abiopen(filepath)
            if not isinstance(wfkfile, WFK_File):
                awx.showErrorMessage(self, message="%s is not a valid WFK File" % filepath)
                return

            self.wfk = wfkfile
            self.BuildUi()
            self.statusbar.PushStatusText("WFK file %s loaded" % filepath)

        except Exception:
            awx.showErrorMessage(self)

    def OnClose(self, event):
        self.wfk = None
        self.DestroyPanel()

    def OnExit(self, event):
        self.Close()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")

    def GetVisualizer(self):
        """Returns a string with the visualizer selected by the user."""
        return self.visualizer_cbox.GetValue()

    def OnVisualizeStructure(self, event):
        """"Call visualizer to visualize the crystalline structure."""
        if self.wfk is None: return

        visualizer = self.GetVisualizer()
        self.statusbar.PushStatusText("Visualizing crystal structure with %s" % visualizer)

        try:
            visu = self.wfk.visualize_structure_with(visualizer)

            thread = awx.WorkerThread(self, target=visu)
            thread.start()

        except:
            awx.showErrorMessage(self)

    def OnVisualizeBZ(self, event):
        """"Visualize the Brillouin zone with matplotlib."""
        if self.wfk is None: return
        self.wfk.structure.show_bz()

    def OnVisualizeWave(self, event):
        """Visualize :math:`|u(r)|^2`."""
        if self.wfk is None: return
        skb = self.skb_panel.GetSKB()
        self._visualize_skb(*skb)

    def _visualize_skb(self, spin, kpoint, band):
        """Calls the visualizer to visualize the specified wavefunction."""
        # To make the Gui responsive one can use the approach described in 
        # http://wiki.wxpython.org/LongRunningTasks
        #visualizer = self.GetVisualizer()
        self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
        try:
            visu = self.wfk.export_ur2(".xsf", spin, kpoint, band)
            #visu = self.wfk.visualize_ur2_with(spin, kpoint, bands, visualizer)
            #visu()

            thread = awx.WorkerThread(self, target=visu)
            thread.start()

        except:
            awx.showErrorMessage(self)

    def OnPlotBands(self, event):
        """Plot band energies with matplotlib."""
        if self.wfk is None: return
        self.wfk.plot_ebands()

    def OnDos(self, event):
        """Open Frame for the computation of the DOS."""
        if self.wfk is None: return
        ewx.ElectronDosFrame(self, bands=self.wfk.ebands()).Show()

    def OnJdos(self, event):
        """Open Frame for the computation of the JDOS."""
        if self.wfk is None: return
        ewx.ElectronJdosFrame(self, bands=self.wfk.ebands).Show()

    def OnNcdump(self, event):
        """Call ncdump and show results in a dialog."""
        if self.wfk is None: return
        caption = "ncdump output for WFK file %s" % self.wfk.filepath
        wxdg.ScrolledMessageDialog(self, self.wfk.ncdump(), caption=caption, style=wx.MAXIMIZE_BOX).Show()


class WfkViewerApp(awx.App):
    def OnInit(self):
        return True

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        self.log("%s dropped on app %s" % (filename, self.appname))
        WfkViewerFrame(parent=None, filename=filename).Show()


def wxapp_wfkviewer(wfk_filename):
    app = WfkViewerApp()
    frame = WfkViewerFrame(None, filename=wfk_filename)
    frame.Show()
    app.SetTopWindow(frame)
    return app


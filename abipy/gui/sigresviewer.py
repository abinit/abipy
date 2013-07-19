from __future__ import print_function, division

import os
import wx

import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy import abiopen, SIGRES_File
from abipy.iotools.visualizer import supported_visunames
from abipy.gui.scissors import ScissorsBuilderFrame

ID_VISTRUCT = wx.NewId()
ID_VISBZ = wx.NewId()
ID_NCDUMP = wx.NewId()
ID_SCISSORS = wx.NewId()

ID_TBOX_VIS = wx.NewId()

class SigresViewerFrame(wx.Frame):
    VERSION = "0.1"

    def __init__(self, parent, filename=None):
        super(SigresViewerFrame, self).__init__(parent, -1, self.codename)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing SIGRES file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the SIGRES file")
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
        toolbar.AddSimpleTool(ID_VISTRUCT, wx.Bitmap(awx.path_img("crystal.png")), "Visualize the crystal structure")
        toolbar.AddSimpleTool(ID_VISBZ, wx.Bitmap(awx.path_img("wave.png")), "Visualize the BZ")
        toolbar.AddSimpleTool(ID_SCISSORS, wx.Bitmap(awx.path_img("wave.png")), "Build energy-dependent scissors from GW correction.")

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
            (ID_VISBZ, self.OnVisualizeBZ),
            (ID_SCISSORS, self.OnScissors),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.sigres = None
        if filename is not None:
            self.ReadSigresFile(filename)

    @property
    def codename(self):
        return self.__class__.__name__

    def BuildUi(self):
        sigres = self.sigres
        if sigres is None: return

        #splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
        #splitter.SetMinimumPaneSize(50)
        #parent = splitter 

        parent = self
        self.skb_panel = awx.SpinKpointBandPanel(parent, sigres.nsppol, sigres.gwkpoints, sigres.max_gwbstop, 
            bstart=sigres.min_gwbstart)

        # Set the callback for double click on k-point row..
        #self.skb_panel.SetOnItemActivated(self._ShowQP_k)

        # Add Python shell
        #from wx.py.shell import Shell
        #pyshell = Shell(parent, introText="SIGRES File methods are available via the variable sigres ", locals={"sigres": sigres})
        #splitter.SplitHorizontally(self.skb_panel, pyshell)

        #main_sizer = wx.BoxSizer(wx.VERTICAL)
        #main_sizer.Add(self.skb_panel, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)
        #main_sizer.Add(pyshell, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)
        #self.SetSizerAndFit(main_sizer)

    def DestroyPanel(self):
        if hasattr(self, "panel"):
            self.panel.Destroy()

    def OnOpen(self, event):
        dlg = wx.FileDialog(self, message="Choose a SIGRES file", defaultDir=os.getcwd(),
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
            self.ReadSigresFile(filepath)

        dlg.Destroy()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        # move up the list
        self.file_history.AddFileToHistory(filepath)
        self.ReadSigresFile(filepath)

    def ReadSigresFile(self, filepath):
        """Read the SIGRES file and build the UI."""
        self.statusbar.PushStatusText("Reading %s" % filepath)

        try:
            sigres = abiopen(filepath)
            if not isinstance(sigres, SIGRES_File):
                awx.showErrorMessage(self, message="%s is not a valid SIGRES file" % filepath)
                return

            self.sigres = sigres
            self.BuildUi()
            self.statusbar.PushStatusText("SIGRES file %s loaded" % filepath)

        except Exception:
            awx.showErrorMessage(self)

    def OnClose(self, event):
        self.sigres = None
        self.DestroyPanel()

    def OnExit(self, event):
        self.Close()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")

    def OnScissors(self, event):
        """Build the scissors operator."""
        if self.sigres is None: return
        ScissorsBuilderFrame(self, self.sigres.filepath).Show()

    def GetVisualizer(self):
        """Returns a string with the visualizer selected by the user."""
        return self.visualizer_cbox.GetValue()

    def OnVisualizeStructure(self, event):
        """"Call visualizer to visualize the crystalline structure."""
        if self.sigres is None: return

        visualizer = self.GetVisualizer()
        self.statusbar.PushStatusText("Visualizing crystal structure with %s" % visualizer)
        structure = self.sigres.get_structure()
        try:
            visu = structure.visualize(visualizer)
            visu()
        except:
            awx.showErrorMessage(self)

    def OnVisualizeBZ(self, event):
        """"Visualize the Brillouin zone with matplotlib."""
        if self.sigres is None: return
        self.sigres.get_structure().show_bz()

    def _ShowQP_k(self, spin, kpoint, band):
        """Calls the visualizer to visualize the specified wavefunction."""
        qpcorr_skb = self.sigres.get_qpcorr(spin, kpoint, band)
        print(qpcorr_skb)

    def OnNcdump(self, event):
        """Call ncdump and show results in a dialog."""
        if self.sigres is None: return
        caption = "ncdump output for SIGRES file %s" % self.sigres.filepath
        wxdg.ScrolledMessageDialog(self, self.sigres.ncdump(), caption=caption, style=wx.MAXIMIZE_BOX).Show()


class SigresViewerApp(awx.App):
    def OnInit(self):
        return True

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        awx.PRINT("%s dropped on app %s" % (filename, self.appname))
        frame = SigresViewerFrame(parent=None, filename=filename)
        frame.Show()


def wxapp_sigresviewer(sigres_filename):
    app = SigresViewerApp()
    frame = SigresViewerFrame(None, filename=sigres_filename)
    frame.Show()
    app.SetTopWindow(frame)
    return app


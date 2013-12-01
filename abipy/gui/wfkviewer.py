from __future__ import print_function, division

import os
import wx

import abipy.gui.awx as awx

from wx.py.shell import Shell
from abipy.tools import marquee
from abipy.abilab import abiopen
from abipy.gui.mixins import Has_Structure, Has_Ebands, Has_Tools, Has_Netcdf


ID_VISWAVE = wx.NewId()

class WfkViewerFrame(awx.Frame, Has_Structure, Has_Ebands, Has_Tools, Has_Netcdf):
    VERSION = "0.1"

    def __init__(self, parent, filename=None, **kwargs):
        super(WfkViewerFrame, self).__init__(parent, -1, title=self.codename, **kwargs)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing WFK file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the WFK file")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")
        menuBar.Append(file_menu, "File")

        menuBar.Append(self.CreateStructureMenu(), "Structure")
        menuBar.Append(self.CreateEbandsMenu(), "Ebands")
        menuBar.Append(self.CreateToolsMenu(), "Tools")
        menuBar.Append(self.CreateNetcdfMenu(), "Netcdf")

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

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        tsize = (48, 48)
        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(ID_VISWAVE, bitmap("wfk.png"), "Visualize the selected wavefunction")
        #toolbar.AddSeparator()
        self.toolbar.Realize()
        self.Centre()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_OPEN, self.OnOpen),
            (wx.ID_CLOSE, self.OnClose),
            (wx.ID_EXIT, self.OnExit),
            (wx.ID_ABOUT, self.OnAboutBox),
            #
            (ID_VISWAVE, self.OnVisualizeWave),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.wfk = None
        if filename is not None:
            self.ReadWfkFile(filename)

    @property
    def codename(self):
        return "WfkViewer"

    @property
    def structure(self):
        """`Structure` object."""
        return self.wfk.structure

    @property
    def ebands(self):
        """`Electron Bands object."""
        return self.wfk.ebands

    @property
    def nc_filepath(self):
        """String with the absolute path of the netcdf file."""
        return self.wfk.filepath

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
        msg = "WFK_File object is accessible via the wfk variable. Use wfk.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        pyshell = Shell(parent, introText=msg, locals={"wfk": self.wfk})
        splitter.SplitHorizontally(self.skb_panel, pyshell)

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
            #if not isinstance(wfkfile, WFK_File):
            #    awx.showErrorMessage(self, message="%s is not a valid WFK File" % filepath)
            #    return

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
            #visu = self.wfk.visualize(".xsf", spin, kpoint, band, visu_name)
            visu = self.wfk.visualize_ur2(spin, kpoint, band, "vesta")

            thread = awx.WorkerThread(self, target=visu)
            thread.start()

        except:
            awx.showErrorMessage(self)


class WfkViewerApp(awx.App):
    def OnInit(self):
        return True

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        #logger.info("%s dropped on app %s" % (filename, self.appname))
        WfkViewerFrame(parent=None, filename=filename).Show()


def wxapp_wfkviewer(wfk_filename):
    app = WfkViewerApp()
    frame = WfkViewerFrame(None, filename=wfk_filename)
    app.SetTopWindow(frame)
    frame.Show()
    return app


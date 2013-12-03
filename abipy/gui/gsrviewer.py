from __future__ import print_function, division

import os
import wx

import abipy.gui.awx as awx

from wx.py.shell import Shell
from abipy.tools import marquee
from abipy.abilab import abiopen
from abipy.gui.mixins import Has_Structure, Has_Ebands, Has_Tools, Has_Netcdf
from abipy.gui.kpoints import KpointsPanel


class GsrViewerFrame(awx.Frame, Has_Structure, Has_Ebands, Has_Tools, Has_Netcdf):
    VERSION = "0.1"

    def __init__(self, parent, filename=None, **kwargs):
        super(GsrViewerFrame, self).__init__(parent, -1, title=self.codename, **kwargs)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing GSR file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the GSR file")
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
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.gsr = None
        if filename is not None:
            self.ReadGsrFile(filename)

    @property
    def codename(self):
        return "GsrViewer"

    @property
    def structure(self):
        """`Structure` object."""
        return self.gsr.structure

    @property
    def ebands(self):
        """`Electron Bands object."""
        return self.gsr.ebands

    @property
    def nc_filepath(self):
        """String with the absolute path of the netcdf file."""
        return self.gsr.filepath

    def BuildUi(self):
        gsr = self.gsr
        if gsr is None: return

        splitter = wx.SplitterWindow(self, id=-1, style=wx.SP_LIVE_UPDATE)
        splitter.SetMinimumPaneSize(50)
        parent = splitter 
        #parent = self

        #self.kpoint_panel = awx.SpinKpointBandPanel(parent, gsr.nsppol, gsr.kpoints, gsr.mband)
        self.kpoints_panel = KpointsPanel(parent, gsr.structure, gsr.kpoints)

        # Add Python shell
        msg = "GSR_File object is accessible via the gsr variable. Use gsr.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        pyshell = Shell(parent, introText=msg, locals={"gsr": self.gsr})
        splitter.SplitHorizontally(self.kpoints_panel, pyshell)

    def DestroyPanel(self):
        if hasattr(self, "panel"):
            self.panel.Destroy()

    def OnOpen(self, event):
        dlg = wx.FileDialog(self, message="Choose a GSR file", defaultDir=os.getcwd(),
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
            self.ReadGsrFile(filepath)

        dlg.Destroy()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        # move up the list
        self.file_history.AddFileToHistory(filepath)
        self.ReadGsrFile(filepath)

    def ReadGsrFile(self, filepath):
        """Read the GSR file and build the UI."""
        self.statusbar.PushStatusText("Reading %s" % filepath)

        try:
            gsrfile = abiopen(filepath)
            #if not isinstance(wfkfile, GSR_File):
            #    awx.showErrorMessage(self, message="%s is not a valid WFK File" % filepath)
            #    return

            self.gsr = gsrfile
            self.BuildUi()
            self.statusbar.PushStatusText("GSR file %s loaded" % filepath)

        except Exception:
            awx.showErrorMessage(self)

    def OnClose(self, event):
        self.gsr = None
        self.DestroyPanel()

    def OnExit(self, event):
        self.Close()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")


class GsrViewerApp(awx.App):
    def OnInit(self):
        return True

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        #logger.info("%s dropped on app %s" % (filename, self.appname))
        GsrViewerFrame(parent=None, filename=filename).Show()


def wxapp_gsrviewer(gsr_filename):
    """Standalone application."""
    app = GsrViewerApp()
    frame = GsrViewerFrame(None, filename=gsr_filename)
    app.SetTopWindow(frame)
    frame.Show()
    return app


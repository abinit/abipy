from __future__ import print_function, division

import os
import wx

import wx.lib.agw.flatnotebook as fnb
import abipy.gui.awx as awx

from wx.py.shell import Shell
from abipy.tools import marquee, list_strings 
from abipy.abilab import abiopen
from abipy.gui import mixins as mix 
from abipy.gui.kpoints import KpointsPanel


class GsrViewerFrame(awx.Frame, mix.Has_Structure, mix.Has_MultipleEbands, mix.Has_Tools, mix.Has_Netcdf):
    VERSION = "0.1"

    def __init__(self, parent, filepaths=(), **kwargs):
        """
        Args:
            parent:
                parent window.
            filepaths:
                String or list of strings with the path of the netcdf WFK files to open
                Empty tuple if no file should be opened during the initialization of the frame.
        """
        super(GsrViewerFrame, self).__init__(parent, -1, title=self.codename, **kwargs)

        # This combination of options for config seems to work on my Mac.
        self.config = wx.FileConfig(appName=self.codename, localFilename=self.codename + ".ini", 
                                    style=wx.CONFIG_USE_LOCAL_FILE)

        # Build menu, toolbar and status bar.
        self.makeMenu()
        self.makeToolBar()
        self.statusbar = self.CreateStatusBar()

        # Open netcdf files.
        filepaths, exceptions = list_strings(filepaths), []
        filepaths = map(os.path.abspath, filepaths)

        # Create the notebook (each file will have its own tab).
        panel = wx.Panel(self, -1)
        self.notebook = fnb.FlatNotebook(panel, -1, style=fnb.FNB_NAV_BUTTONS_WHEN_NEEDED)
                                                                                           
        for path in filepaths:
            gsr = abiopen(path)
            tab = GsrFileTab(self.notebook, gsr)
            self.notebook.AddPage(tab, os.path.basename(path))
                                                                                           
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.EXPAND, 5)
        panel.SetSizerAndFit(sizer)

    @property
    def codename(self):
        return "GsrViewer"

    @property
    def active_tab(self):
        """Returns the active tab. None if notebook is empty."""
        return self.notebook.GetCurrentPage()

    @property
    def active_gsr(self):
        """The active GSR file i.e. the GSR associated to the active tab."""
        return self.active_tab.gsr

    @property
    def structure(self):
        """`Structure` associated to the active tab."""
        return self.active_gsr.structure

    @property
    def ebands(self):
        """`ElectronBands` associated to the active tab."""
        return self.active_gsr.ebands

    @property
    def ebands_list(self):
        """List of `ElectronBands`."""
        ebands_list = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            ebands_list.append(tab.gsr.ebands)
        return ebands_list

    @property
    def ebands_filepaths(self):
        """
        Return a list with the absolute paths of the files 
        from which the `ElectronBands` have been read.
        """
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.gsr.filepath)
        return paths

    @property
    def nc_filepath(self):
        """String with the absolute path of the netcdf file."""
        return self.active_gsr.filepath

    def makeMenu(self):
        """Creates the main menu."""
        # Menu IDs

        self.menu_bar = menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing GSR file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the GSR file")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")

        file_history = self.file_history = wx.FileHistory(8)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)
        menuBar.Append(file_menu, "File")

        menuBar.Append(self.CreateStructureMenu(), "Structure")
        menuBar.Append(self.CreateEbandsMenu(), "Ebands")
        menuBar.Append(self.CreateToolsMenu(), "Tools")
        menuBar.Append(self.CreateNetcdfMenu(), "Netcdf")

        self.help_menu = wx.Menu()
        self.help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.help_menu, "Help")

        self.SetMenuBar(menuBar)

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_OPEN, self.OnOpen),
            (wx.ID_CLOSE, self.OnClose),
            (wx.ID_EXIT, self.OnExit),
            (wx.ID_ABOUT, self.OnAboutBox),
        ]
                                                            
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

    def makeToolBar(self):
        """Creates the toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR), "Open")
        #toolbar.AddSeparator()

        toolbar.Realize()

    def AddFileToHistory(self, filepath):
        """Add the absolute filepath to the file history."""
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

    def read_file(self, filepath):
        """Open netcdf file, create new tab and save the file in the history."""
        self.statusbar.PushStatusText("Reading %s" % filepath)
        try:
            notebook = self.notebook
            gsr = abiopen(filepath)
            #if not isinstance(gsrfile, GSR_File):
            #    awx.showErrorMessage(self, message="%s is not a valid GSR File" % filepath)
            #    return
            self.statusbar.PushStatusText("GSR file %s loaded" % filepath)
            tab = GsrFileTab(notebook, gsr)
            notebook.AddPage(tab, os.path.basename(filepath))
            # don't know why but this does not work!
            notebook.Refresh()
            notebook.SetSelection(notebook.GetPageCount())
            self.AddFileToHistory(filepath)
        except:
            awx.showErrorMessage(self)

    def OnOpen(self, event):
        """Open FileDialog to allow the user to select a WFK.nc file."""
        # Show the dialog and retrieve the user response.
        # If it is the OK response, process the data.
        dialog = wx.FileDialog(self, message="Choose a GSR file", defaultDir=os.getcwd(),
                               wildcard="GSR Netcdf files (*.nc)|*.nc",
                               style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)
        if dialog.ShowModal() == wx.ID_CANCEL: return 

        filepath = os.path.abspath(dialog.GetPath())
        self.read_file(filepath)

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        self.read_file(filepath)

    def OnClose(self, event):
        """
        Remove the active tab from the notebook and 
        close the corresponding netcdf file, 
        """
        notebook = self.notebook
        if notebook.GetPageCount() == 0: return
        idx = notebook.GetSelection()
        if idx == -1: return None

        # Close the file
        tab = notebook.GetPage(idx)
        #tab.gsr.close()

        # Remove tab.
        notebook.DeletePage(idx)
        notebook.Refresh()
        #notebook.SendSizeEvent()

    def OnExit(self, event):
        """Exits the application."""
        # Close open netcdf files.
        #try:
        #    for index in range(self.notebook.GetPageCount()):
        #        tab = self.notebook.GetPage(index)
        #        try:
        #            tab.wfk.close()
        #        except:
        #            pass
        #finally:
        self.Destroy()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")


class GsrFileTab(wx.Panel):
    """Tab showing information on a single GSR file."""
    def __init__(self, parent, gsr, **kwargs):
        """
        Args:
            parent:
                parent window.
            gsr:
        """
        super(GsrFileTab, self).__init__(parent, -1, **kwargs)
        self.gsr = gsr

        splitter = wx.SplitterWindow(self, id=-1, style=wx.SP_3D)
        splitter.SetSashGravity(0.95)

        self.kpoints_panel = KpointsPanel(splitter, gsr.structure, gsr.kpoints)

        # Add Python shell
        msg = "GSR_File object is accessible via the gsr variable. Use gsr.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        pyshell = Shell(splitter, introText=msg, locals={"gsr": self.gsr})
        splitter.SplitHorizontally(self.kpoints_panel, pyshell)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

    @property
    def statusbar(self):
        return self.viewer_frame.statusbar

    def GetVisualizer(self):
        """Returns a string with the visualizer selected by the user."""
        return self.viewer_frame.GetVisualizer()

    def _visualize_skb(self, spin, kpoint, band):
        """Calls the visualizer to visualize the specified wavefunction."""
        # To make the Gui responsive one can use the approach described in 
        # http://wiki.wxpython.org/LongRunningTasks
        visu_name = self.GetVisualizer()
        if visu_name == "None": return
                                                                                                                       
        self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
        try:
            visu = self.wfk.visualize_ur2(spin, kpoint, band, visu_name=visu_name)
                                                                                                                       
            thread = awx.WorkerThread(self, target=visu)
            thread.start()
                                                                                                                       
        except:
            awx.showErrorMessage(self)

    @property
    def viewer_frame(self):
        """The parent frame `GsrViewerFrame`."""
        try:
            return self._viewer_frame
                                                                                    
        except AttributeError:
            self._viewer_frame = self.getParentWithType(GsrViewerFrame)
            return self._viewer_frame


class GsrViewerApp(awx.App):
    def OnInit(self):
        return True

    #def MacOpenFile(self, filename):
    #    """Called for files droped on dock icon, or opened via finders context menu"""
    #    if filename.endswith(".py"):
    #        return
    #    # Open filename in a new frame.
    #    #logger.info("%s dropped on app %s" % (filename, self.appname))
    #    GsrViewerFrame(parent=None, filename=filename).Show()


def wxapp_gsrviewer(gsr_filepaths):
    """Standalone application."""
    app = GsrViewerApp()
    frame = GsrViewerFrame(None, filepaths=gsr_filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    return app

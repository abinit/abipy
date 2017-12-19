from __future__ import print_function, division, unicode_literals, absolute_import

import os
import wx
import abipy.gui.awx as awx

from wx.py.shell import Shell
from monty.string import marquee
from abipy.abilab import abiopen
from abipy.gui import mixins as mix
from abipy.gui.kpoints import KpointsPanel
from abipy.gui.baseviewer import MultiViewerFrame


class GsrViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_MultipleEbands, mix.Has_Tools, mix.Has_NetcdfFiles):
    VERSION = "0.1"

    HELP_MSG = """Quick help:

 Kpoint list:

     Right-Click:  display popup menu with choices.

Also, these key bindings can be used
(For Mac OSX, replace 'Ctrl' with 'Apple'):

  Ctrl-Q:     quit
"""

    @property
    def codename(self):
        return "GsrViewer"

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
    def nc_filepaths(self):
        """String with the absolute path of the netcdf file."""
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.gsr.filepath)
        return paths

    def makeMenu(self):
        """Creates the main menu."""
        # Base menu
        menu_bar = super(GsrViewerFrame, self).makeMenu()

        # Add Mixin menus.
        menu_bar.Append(self.CreateStructureMenu(), "Structure")
        menu_bar.Append(self.CreateEbandsMenu(), "Ebands")
        menu_bar.Append(self.CreateToolsMenu(), "Tools")
        menu_bar.Append(self.CreateNetcdfMenu(), "Netcdf")

        # Help menu
        help_menu = self.makeHelpMenu()
        menu_bar.Append(help_menu, "Help")

        self.SetMenuBar(menu_bar)

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

    def addFileTab(self, parent, filepath):
        """Read data from filepath and create a new notebook tab."""
        gsr = abiopen(filepath)
        tab = GsrFileTab(self.notebook, gsr)
        self.notebook.AddPage(tab, os.path.basename(filepath))


class GsrFileTab(wx.Panel):
    """Tab showing information on a single GSR file."""
    def __init__(self, parent, gsr, **kwargs):
        """
        Args:
            parent:
                parent window.
            gsr:
                `GsrFile` instance.
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
        appname = self.GetVisualizer()
        if appname == "None": return

        self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
        try:
            visu = self.wfk.visualize_ur2(spin, kpoint, band, appname=appname)
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
    pass


def wxapp_gsrviewer(gsr_filepaths):
    """Standalone application."""
    app = GsrViewerApp()
    frame = GsrViewerFrame(None, filepaths=gsr_filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    return app

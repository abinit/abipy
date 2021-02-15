import os
import wx
import abipy.gui.awx as awx

from wx.py.shell import Shell
from monty.string import marquee
from abipy.abilab import abiopen
from abipy.iotools.visualizer import Visualizer
from abipy.gui import mixins as mix
from abipy.gui.kpoints import SpinKpointBandPanel
from abipy.gui.baseviewer import MultiViewerFrame


class WfkViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_MultipleEbands, mix.Has_Tools, mix.Has_NetcdfFiles):
    VERSION = "0.1"

    HELP_MSG = """Quick help:

 Kpoint list:

     Left-Click:   to visualize u(r)^2 for the selected spin, k-point, band
     Right-Click:  display popup menu with choices.

Also, these key bindings can be used
(For Mac OSX, replace 'Ctrl' with 'Apple'):

  Ctrl-Q:     quit
"""
    @property
    def codename(self):
        return "WfkViewer"

    @property
    def active_wfk(self):
        """The active WFK file i.e. the WFK associated to the active tab."""
        return self.active_tab.wfk

    @property
    def structure(self):
        """`Structure` associated to the active tab."""
        return self.active_wfk.structure

    @property
    def ebands(self):
        """`ElectronBands` associated to the active tab."""
        return self.active_wfk.ebands

    @property
    def ebands_list(self):
        """List of `ElectronBands`."""
        ebands_list = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            ebands_list.append(tab.wfk.ebands)
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
            paths.append(tab.wfk.filepath)
        return paths

    @property
    def nc_filepaths(self):
        """String with the absolute path of the active netcdf file."""
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.wfk.filepath)
        return paths

    def makeMenu(self):
        """Creates the main menu."""
        # Base menu.
        menu_bar = super(WfkViewerFrame, self).makeMenu()

        # Add Mixin menus.
        menu_bar.Append(self.CreateStructureMenu(), "Structure")
        menu_bar.Append(self.CreateEbandsMenu(), "Ebands")
        menu_bar.Append(self.CreateToolsMenu(), "Tools")
        menu_bar.Append(self.CreateNetcdfMenu(), "Netcdf")

        # Help menu
        help_menu = self.makeHelpMenu()
        menu_bar.Append(help_menu, "Help")

        #self.menu_bar = menu_bar
        self.SetMenuBar(menu_bar)

    def makeToolBar(self):
        """Creates the toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR), "Open")
        toolbar.AddSeparator()

        # Combo box with the list of visualizers
        avail_visunames = [visu.name for visu in Visualizer.get_available()]
        value = avail_visunames[0] if avail_visunames else "None"
        self.visualizer_cbox = wx.ComboBox(toolbar, id=-1, name='visualizer', choices=avail_visunames, value=value, style=wx.CB_READONLY)
        toolbar.AddControl(control=self.visualizer_cbox)

        toolbar.Realize()

    def addFileTab(self, parent, filepath):
        wfk = abiopen(filepath)
        tab = WfkFileTab(self.notebook, wfk)
        self.notebook.AddPage(tab, os.path.basename(filepath))

    def GetVisualizer(self):
        """Returns a string with the visualizer selected by the user."""
        return self.visualizer_cbox.GetValue()


class WfkFileTab(awx.Panel):
    """Tab showing information on a single WFK file."""
    def __init__(self, parent, wfk, **kwargs):
        """
        Args:
            parent:
                parent window.
            wfk:
        """
        super(WfkFileTab, self).__init__(parent, -1, **kwargs)
        self.wfk = wfk

        splitter = wx.SplitterWindow(self, id=-1, style=wx.SP_3D)
        splitter.SetSashGravity(0.95)

        self.skb_panel = skb_panel = SpinKpointBandPanel(splitter, wfk.structure, wfk.nsppol, wfk.kpoints, wfk.mband)

        # Set the callback for double click on k-point row..
        self.Bind(skb_panel.MYEVT_SKB_ACTIVATED, self.onVisualizeSKB)

        # Add Python shell
        msg = "WFK_File object is accessible via the wfk variable. Use wfk.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        # FIXME <Error>: CGContextRestoreGState: invalid context 0x0
        pyshell = Shell(splitter, introText=msg,locals={"wfk": self.wfk})
        splitter.SplitHorizontally(self.skb_panel, pyshell)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

    @property
    def statusbar(self):
        return self.viewer_frame.statusbar

    def GetVisualizer(self):
        """Returns a string with the visualizer selected by the user."""
        return self.viewer_frame.GetVisualizer()

    def onVisualizeSKB(self, event):
        """
        Calls the visualizer to visualize the specified wavefunction.
        Use the approach described in http://wiki.wxpython.org/LongRunningTasks
        to make the Gui responsive one can use
        """
        #print("in new" ,event.skb)
        spin, kpoint, band = event.skb

        appname = self.GetVisualizer()
        if appname == "None": return

        self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
        try:
            visu = self.wfk.visualize_ur2(spin, kpoint, band, appname=appname)
            thread = awx.WorkerThread(self, target=visu)
            thread.start()

        except:
            awx.showErrorMessage(self)

    #def _visualize_skb(self, spin, kpoint, band):
    #    """Calls the visualizer to visualize the specified wavefunction."""
    #    # To make the Gui responsive one can use the approach described in
    #    # http://wiki.wxpython.org/LongRunningTasks
    #    appname = self.GetVisualizer()
    #    if appname == "None": return
    #
    #    self.statusbar.PushStatusText("Visualizing wavefunction (spin=%d, kpoint=%s, band=%d)" % (spin, kpoint, band))
    #    try:
    #        visu = self.wfk.visualize_ur2(spin, kpoint, band, appname=appname)
    #        thread = awx.WorkerThread(self, target=visu)
    #        thread.start()
    #
    #    except:
    #        awx.showErrorMessage(self)

    @property
    def viewer_frame(self):
        """The parent frame `WfkViewerFrame`."""
        try:
            return self._viewer_frame

        except AttributeError:
            self._viewer_frame = self.getParentWithType(WfkViewerFrame)
            return self._viewer_frame


class WfkViewerApp(awx.App):
    pass


def wxapp_wfkviewer(wfk_filepaths):
    app = WfkViewerApp()
    frame = WfkViewerFrame(None, filepaths=wfk_filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    return app


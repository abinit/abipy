from __future__ import print_function, division, unicode_literals, absolute_import

import os
import wx
import abc
import abipy.gui.awx as awx

from wx.py.shell import Shell
from monty.string import marquee
from abipy.abilab import abiopen
from abipy.gui import mixins as mix
from abipy.gui.kpoints import KpointsPanel
from abipy.gui.baseviewer import MultiViewerFrame


class PhbstViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_MultiplePhbands, mix.Has_Tools, mix.Has_NetcdfFiles):
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
        return "PhbstViewer"

    @property
    def active_phbst(self):
        """The active PHBST file i.e. the PHBST associated to the active tab."""
        return self.active_tab.phbst

    @property
    def structure(self):
        """`Structure` associated to the active tab."""
        return self.active_phbst.structure

    @property
    def phbands(self):
        """`PhononBands` associated to the active tab."""
        return self.active_phbst.phbands

    @property
    def active_phdos_file(self):
        """The active PHDOS file i.e. the PHDOS associated to the active tab."""
        try:
            return self.active_tab.phdos_file
        except AttributeError:
            return None

    @property
    def phdos(self):
        """PHDOS data for the active tab if it has been added. None otherwise"""
        if self.active_phdos_file:
            return self.active_phdos_file.phdos
        else:
            return None

    @property
    def phbands_list(self):
        """List of `PhononBands`."""
        phbands_list = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            phbands_list.append(tab.phbst.phbands)
        return phbands_list

    @property
    def phbands_filepaths(self):
        """
        Return a list with the absolute paths of the files
        from which the `PhononBands` have been read.
        """
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.phbst.filepath)
        return paths

    @property
    def phdos_list(self):
        """List of `PhononDos`."""
        phdos_list = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            if tab.phdos_file:
                phdos_list.append(tab.phdos_file.phdos)
        return phdos_list

    @property
    def phdos_filepaths(self):
        """
        Return a list with the absolute paths of the files
        from which the `PhononDos` have been read.
        """
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            if tab.phdos_file:
                paths.append(tab.phdos_file.filepath)
        return paths

    @property
    def nc_filepaths(self):
        """String with the absolute path of the netcdf file."""
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.phbst.filepath)
        return paths

    def makeMenu(self):
        """Creates the main menu."""
        # Base menu
        menu_bar = super(PhbstViewerFrame, self).makeMenu()

        # Add Mixin menus.
        menu_bar.Append(self.CreateStructureMenu(), "Structure")
        menu_bar.Append(self.CreatePhbandsMenu(), "Phbands")
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
        phbst = abiopen(filepath)
        tab = PhbstFileTab(self.notebook, phbst)
        self.notebook.AddPage(tab, os.path.basename(filepath))


class PhbstFileTab(wx.Panel):
    """Tab showing information on a single PHBST file."""
    def __init__(self, parent, phbst, **kwargs):
        """
        Args:
            parent:
                parent window.
            phbst:
                `PhbstFile` instance.
        """
        super(PhbstFileTab, self).__init__(parent, -1, **kwargs)
        self.phbst = phbst
        self.phdos_file = None

        splitter = wx.SplitterWindow(self, id=-1, style=wx.SP_3D)
        splitter.SetSashGravity(0.95)

        self.qpoints_panel = KpointsPanel(splitter, phbst.structure, phbst.qpoints)

        # Add Python shell
        msg = "PHBST_File object is accessible via the phbst variable. Use phbst.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        pyshell = Shell(splitter, introText=msg, locals={"phbst": self.phbst})
        splitter.SplitHorizontally(self.qpoints_panel, pyshell)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

    @property
    def statusbar(self):
        return self.viewer_frame.statusbar

    def GetVisualizer(self):
        """Returns a string with the visualizer selected by the user."""
        return self.viewer_frame.GetVisualizer()

    @property
    def viewer_frame(self):
        """The parent frame `PhbstViewerFrame`."""
        try:
            return self._viewer_frame

        except AttributeError:
            self._viewer_frame = self.getParentWithType(PhbstViewerFrame)
            return self._viewer_frame


class PhbstViewerApp(awx.App):
    pass


def wxapp_phbstviewer(phbst_filepaths):
    """Standalone application."""
    app = PhbstViewerApp()
    if not phbst_filepaths:
        phbst_filepaths = ()
    frame = PhbstViewerFrame(None, filepaths=phbst_filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    return app

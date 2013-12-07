from __future__ import print_function, division

import os
import wx

import wx.lib.agw.flatnotebook as fnb
import abipy.gui.awx as awx

from wx.py.shell import Shell
from abipy.abilab import abiopen
from abipy.tools import marquee, list_strings 
from abipy.electrons.bse import MDF_Plotter
from abipy.iotools.visualizer import Visualizer
from abipy.gui.kpoints import KpointsPanel
from abipy.gui import mixins as mix 
from abipy.gui.baseviewer import MultiViewerFrame


#class MdfViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_MultipleEbands, mix.Has_Tools, mix.Has_NetcdfFiles):
class MdfViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_Tools, mix.Has_NetcdfFiles):
    VERSION = "0.1"

    HELP_MSG = """Quick help:

 Qpoint list:

     Right-Click:  display popup menu with choices.

Also, these key bindings can be used
(For Mac OSX, replace 'Ctrl' with 'Apple'):

  Ctrl-Q:     quit
"""
    @property
    def codename(self):
        return "MdfViewer"

    @property
    def active_mdf_file(self):
        """The active MDF file i.e. the MDF associated to the active tab."""
        return self.active_tab.mdf_file

    @property
    def structure(self):
        """`Structure` associated to the active tab."""
        return self.active_mdf_file.structure

    @property
    def ebands(self):
        """`ElectronBands` associated to the active tab."""
        return self.active_mdf_file.ebands

    @property
    def ebands_list(self):
        """List of `ElectronBands`."""
        ebands_list = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            ebands_list.append(tab.mdf_file.ebands)
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
            paths.append(tab.mdf_file.filepath)
        return paths

    @property
    def nc_filepaths(self):
        """String with the absolute paths of the netcdf files."""
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.mdf_file.filepath)
        return paths

    @property
    def mdf_filepaths(self):
        """
        Return a list with the absolute paths of the files 
        from which the `MDF_Files` have been read.
        """
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.mdf_file.filepath)
        return paths

    @property
    def mdf_files_list(self):
        """List of `MDF_Files`."""
        mdf_lists = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            mdf_lists.append(tab.mdf_file)
        return mdf_lists

    def makeMenu(self):
        """Creates the main menu."""
        # Base menu.
        menu_bar = super(MdfViewerFrame, self).makeMenu()

        # Add Mixin menus.
        menu_bar.Append(self.CreateStructureMenu(), "Structure")
        #menu_bar.Append(self.CreateEbandsMenu(), "Ebands")
        menu_bar.Append(self.CreateMdfMenu(), "Mdf")
        menu_bar.Append(self.CreateToolsMenu(), "Tools")
        menu_bar.Append(self.CreateNetcdfMenu(), "Netcdf")

        # Help menu
        help_menu = self.makeHelpMenu()
        menu_bar.Append(help_menu, "Help")
                                           
        self.SetMenuBar(menu_bar)

    def CreateMdfMenu(self):
        # MDF Menu ID's
        self.ID_MDF_PLOT = wx.NewId()
        self.ID_MDF_COMPARE = wx.NewId()
                                                                                      
        menu = wx.Menu()
        menu.Append(self.ID_MDF_PLOT, "Plot MDF", "Plot the macroscopic dielectric function")
        self.Bind(wx.EVT_MENU, self.OnMdfPlot, id=self.ID_MDF_PLOT)

        menu.AppendSeparator()

        menu.Append(self.ID_MDF_COMPARE, "Compare MDF", "Compare multiple macroscopic dielectric functions")
        self.Bind(wx.EVT_MENU, self.OnMdfCompare, id=self.ID_MDF_COMPARE)                                                                                                            

        return menu

    def OnMdfPlot(self, event):
        mdf_file = self.active_mdf_file
        mdf_file.plot_mdfs(cplx_mode="Im", mdf_type="all")

    def OnMdfCompare(self, event):
        plotter = MDF_Plotter()

        for path, mdf_file in zip(self.mdf_filepaths, self.mdf_files_list):
            label = os.path.relpath(path)
            # TODO
            #plotter.add_mdf(label, mdf)
            plotter.add_mdf_from_file(mdf_file.filepath)

        plotter.plot()

    def makeToolBar(self):
        """Creates the toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR), "Open")

        toolbar.AddSeparator()

        # Combo box with the list of MDF types
        mdf_types = ["ALL", "EXC", "RPA", "GWRPA"]
        self.mdftype_cbox = wx.ComboBox(toolbar, id=-1, name='MDF type', choices=mdf_types, value="ALL", style=wx.CB_READONLY) 
        toolbar.AddControl(control=self.mdftype_cbox) 

        # Combo box with the list of complex modes
        cplx_modes = ["Re", "Im"]
        self.cplx_cbox = wx.ComboBox(toolbar, id=-1, name='COMPLEX mode', choices=cplx_modes, value="Im", style=wx.CB_READONLY) 
        toolbar.AddControl(control=self.cplx_cbox) 

        toolbar.Realize()

    def addFileTab(self, parent, filepath):
        mdf_file = abiopen(filepath)
        tab = MdfFileTab(self.notebook, mdf_file)
        self.notebook.AddPage(tab, os.path.basename(filepath))

    def getMdfType(self):
        """Return the sting with the MDF type selected by the user."""
        mdf_type = self.mdftype_cbox.GetStringSelection()
        if not mdf_type: mdf_type = "ALL"
        return mdf_type

    def getCplxMode(self):
        """Return the sting with the complex mode used for plotting the spectra."""
        cplx_mode = self.cplx_cbox.GetStringSelection()
        if not cplx_mode: cplx_mode = "Im"
        return cplx_mode



class MdfQpointsPanel(KpointsPanel):
    """Extend KpointsPanel adding popupmenus"""

    def __init__(self, parent, mdf_file, **kwargs):
        KpointsPanel.__init__(self, parent, mdf_file.structure, mdf_file.qpoints)
        self.mdf_file = mdf_file

    def makePopupMenu(self):
        # menu of the base class
        menu = super(MdfQpointsPanel, self).makePopupMenu()

        # Add other options
        self.ID_POPUP_MDF_QPLOT = wx.NewId()
        menu.Append(self.ID_POPUP_MDF_QPLOT, "Plot spectra(q)")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_POPUP_MDF_QPLOT, self.onQpointPlot),
        ]
                                                            
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)
                                                     
        return menu

    def onQpointPlot(self, event):
        """Plot MDF(q)"""
        qpoint = self.getSelectedKpoint()
        if qpoint is None: return
        #print(qpoint)
        self.mdf_file.plot_mdfs(cplx_mode="Im", mdf_type="all", qpoint=qpoint)


class MdfFileTab(wx.Panel):
    """Tab showing information on a single MDF file."""
    def __init__(self, parent, mdf_file, **kwargs):
        """
        Args:
            parent:
                parent window.
            mdf_file:
        """
        super(MdfFileTab, self).__init__(parent, -1, **kwargs)
        self.mdf_file = mdf_file

        splitter = wx.SplitterWindow(self, id=-1, style=wx.SP_3D)
        splitter.SetSashGravity(0.95)

        #self.qpts_panel = KpointsPanel(splitter, mdf_file.structure, mdf_file.qpoints)
        self.qpts_panel = MdfQpointsPanel(splitter, mdf_file)

        # Add Python shell
        msg = "MDF_object is accessible via the mdf_file variable. Use mdf_file.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        # FIXME <Error>: CGContextRestoreGState: invalid context 0x0
        pyshell = Shell(splitter, introText=msg, locals={"mdf_file": mdf_file})
        splitter.SplitHorizontally(self.qpts_panel, pyshell)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

    @property
    def statusbar(self):
        return self.viewer_frame.statusbar

    @property
    def viewer_frame(self):
        """The parent frame `MdfViewerFrame`."""
        try:
            return self._viewer_frame
                                                                                    
        except AttributeError:
            self._viewer_frame = self.getParentWithType(MdfViewerFrame)
            return self._viewer_frame


class MdfViewerApp(awx.App):
    pass


def wxapp_mdfviewer(mdf_filepaths):
    app = MdfViewerApp()
    frame = MdfViewerFrame(None, filepaths=mdf_filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    return app


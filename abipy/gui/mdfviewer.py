import os
import wx
import wx.lib.agw.flatnotebook as fnb
import abipy.gui.awx as awx

from wx.py.shell import Shell
from monty.string import marquee
from abipy.abilab import abiopen
from abipy.electrons.bse import MdfPlotter
from abipy.iotools.visualizer import Visualizer
from abipy.gui.kpoints import KpointsPanel
from abipy.gui import mixins as mix
from abipy.gui.baseviewer import MultiViewerFrame


# TODO Add ebands to MDF.nc
class MdfViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_MultipleEbands, mix.Has_Tools, mix.Has_NetcdfFiles):
#class MdfViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_Tools, mix.Has_NetcdfFiles):
    VERSION = "0.1"

    HELP_MSG = """Quick help:

 Toolbar:
    Click the button to plot the spectrum (averaged over the q-points)
    Use the combo boxes to select the type of spectrum and the quantity to plot.

 Qpoint list:

     Right-Click:  display popup menu with choices.
     Select:  plot the spectra for this qpoint.

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
        self.ID_MDF_PLOT_AVERAGE = wx.NewId()
        self.ID_MDF_COMPARE = wx.NewId()

        menu = wx.Menu()
        menu.Append(self.ID_MDF_PLOT_AVERAGE, "Plot averaged MDF", "Plot the average of the macroscopic dielectric function")
        self.Bind(wx.EVT_MENU, self.onPlotAveragedMdf, id=self.ID_MDF_PLOT_AVERAGE)

        menu.AppendSeparator()

        menu.Append(self.ID_MDF_COMPARE, "Compare MDF", "Compare multiple macroscopic dielectric functions")
        self.Bind(wx.EVT_MENU, self.OnMdfCompare, id=self.ID_MDF_COMPARE)

        return menu

    def onPlotAveragedMdf(self, event):
        """Plot the average of the macroscopic dielectric function"""
        mdf_type = self.getMdfType()
        cplx_mode = self.getCplxMode()
        #print(mdf_type, cplx_mode)

        mdf_file = self.active_mdf_file
        mdf_file.plot_mdfs(cplx_mode=cplx_mode, mdf_type=mdf_type)

    def OnMdfCompare(self, event):
        """Compare multiple averaged macroscopic dielectric functions"""
        mdf_type = self.getMdfType()
        cplx_mode = self.getCplxMode()

        if mdf_type == "ALL":
            return awx.showErrorMessage(self, "ALL is not supported by Compare. Please use EXC, RPA, GWRPA")

        plotter = MdfPlotter()
        for path, mdf_file in zip(self.mdf_filepaths, self.mdf_files_list):
            label = os.path.relpath(path)
            mdf = mdf_file.get_mdf(mdf_type)
            plotter.add_mdf(label, mdf)

        plotter.plot(cplx_mode, qpoint=None)

    def makeToolBar(self):
        """Creates the toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR), "Open")
        toolbar.AddSeparator()

        # Button to plot the averaged MDD
        # TODO: Change icon.
        toolbar.AddSimpleTool(self.ID_MDF_PLOT_AVERAGE, bitmap("wave.png"), "Plot averaged MDF")

        # Combo box with the list of MDF types
        mdf_types = ["ALL", "EXC", "RPA", "GWRPA"]
        self.mdftype_cbox = wx.ComboBox(toolbar, id=-1, name='MDF type', choices=mdf_types, value="ALL", style=wx.CB_READONLY)
        self.mdftype_cbox.SetToolTipString("Select the type of MDF spectra to plot.")
        toolbar.AddControl(control=self.mdftype_cbox, label="MDF Type:")

        # Combo box with the list of complex modes
        cplx_modes = ["Re-Im", "Re", "Im"]
        self.cplx_cbox = wx.ComboBox(toolbar, id=-1, name='COMPLEX mode', choices=cplx_modes, value="Re-Im", style=wx.CB_READONLY)
        self.cplx_cbox.SetToolTipString("Select the component of the MDF spectra to plot (real or imaginary part).")
        toolbar.AddControl(control=self.cplx_cbox, label="Complex Mode:")

        toolbar.Realize()

    def addFileTab(self, parent, filepath):
        """Open a file and add a new tab to the notebook."""
        mdf_file = abiopen(filepath)
        tab = MdfFileTab(self.notebook, mdf_file)
        self.notebook.AddPage(tab, os.path.basename(filepath))

        # List to events triggered by the popup menu in the qpoints_table.
        qpanel = tab.qpoints_panel
        self.Bind(qpanel.MYEVT_COMPARE_SPECTRAQ, self.onCompareSpectraQ)

    def getMdfType(self):
        """Return the sting with the MDF type selected by the user."""
        mdf_type = self.mdftype_cbox.GetStringSelection()
        if not mdf_type: mdf_type = "ALL"
        return mdf_type

    def getCplxMode(self):
        """Return the sting with the complex mode used for plotting the spectra."""
        cplx_mode = self.cplx_cbox.GetStringSelection()
        if not cplx_mode: cplx_mode = "Re-Im"
        return cplx_mode

    def onCompareSpectraQ(self, event):
        """
        Compare different MDF(q) spectra (provided that we have multiple tabs in the notebook).
        This callback is executed when MdfQpointsPanel fires CompareSpectraQEvent.
        """
        qpoint = event.qpoint
        mdf_type = self.getMdfType()
        cplx_mode = self.getCplxMode()

        if mdf_type == "ALL":
            return awx.showErrorMessage(self, "ALL is not supported by Compare. Please use EXC, RPA, GWRPA")

        plotter = MdfPlotter()

        # Extract the type of MDF we are interested in
        for path, mdf_file in zip(self.mdf_filepaths, self.mdf_files_list):
            label = os.path.relpath(path)
            mdf = mdf_file.get_mdf(mdf_type)
            plotter.add_mdf(label, mdf)

        plotter.plot(cplx_mode, qpoint=qpoint)


class MdfQpointsPanel(KpointsPanel):
    """Extend KpointsPanel adding popupmenus"""

    # Command event used to signal that the we should compare multiple spectra (fixed q, move along files).
    CompareSpectraQEvent, MYEVT_COMPARE_SPECTRAQ = wx.lib.newevent.NewCommandEvent()

    def __init__(self, parent, mdf_file, **kwargs):
        KpointsPanel.__init__(self, parent, mdf_file.structure, mdf_file.qpoints)
        self.parent = parent
        self.mdf_file = mdf_file

        # Connect the events fired by klist_ctrl.
        self.klist_ctrl.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.onQpointActivated)

    def makePopupMenu(self):
        """Build the popup menu."""
        # Menu of the base class
        menu = super(MdfQpointsPanel, self).makePopupMenu()

        # Add other options
        self.ID_POPUP_MDF_QCOMPARE = wx.NewId()
        menu.Append(self.ID_POPUP_MDF_QCOMPARE, "Compare spectra(q)")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_POPUP_MDF_QCOMPARE, self.onCompareSpectraQ),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu

    @property
    def viewer_frame(self):
        """The parent frame `MdfViewerFrame`."""
        try:
            return self._viewer_frame

        except AttributeError:
            self._viewer_frame = self.getParentWithType(MdfViewerFrame)
            return self._viewer_frame

    def onQpointActivated(self, event):
        """Plot MDF(q) for the selected q-point."""
        qpoint = self.getSelectedKpoint()
        if qpoint is None: return

        # Get the options from the ViewerFrame.
        mdf_type = self.viewer_frame.getMdfType()
        cplx_mode = self.viewer_frame.getCplxMode()

        #print(qpoint)
        self.mdf_file.plot_mdfs(cplx_mode=cplx_mode, mdf_type=mdf_type, qpoint=qpoint)

    def onCompareSpectraQ(self, event):
        """
        Get the selected q-point and fire CompareSpectraQEvent.
        The parent frame will detect the event and will compare the
        different MDF spectra (provided that we have multiple tabs in the notebook.
        """
        qpoint = self.getSelectedKpoint()
        if qpoint is None: return
        # Post the event.
        event = self.CompareSpectraQEvent(id=-1, qpoint=qpoint)
        wx.PostEvent(self.parent, event)


class MdfFileTab(awx.Panel):
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

        self.qpoints_panel = MdfQpointsPanel(splitter, mdf_file)

        # Add Python shell
        msg = "MDF_object is accessible via the mdf_file variable. Use mdf_file.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        # FIXME <Error>: CGContextRestoreGState: invalid context 0x0
        pyshell = Shell(splitter, introText=msg, locals={"mdf_file": mdf_file})
        splitter.SplitHorizontally(self.qpoints_panel, pyshell)

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
    """Standalone application."""
    app = MdfViewerApp()
    frame = MdfViewerFrame(None, filepaths=mdf_filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    return app


from __future__ import print_function, division

import os
import wx

import wx.lib.agw.flatnotebook as fnb
import abipy.gui.awx as awx

from wx.py.shell import Shell
from abipy.abilab import abiopen
from abipy.tools import AttrDict, marquee, list_strings 
from abipy.electrons import SIGRES_Plotter
from abipy.gui.scissors import ScissorsBuilderFrame
from abipy.gui.baseviewer import MultiViewerFrame
from abipy.gui import mixins as mix
from abipy.gui.baseviewer import MultiViewerFrame


class SigresViewerFrame(MultiViewerFrame, mix.Has_Structure, mix.Has_MultipleEbands, mix.Has_Tools, mix.Has_NetcdfFiles):
    VERSION = "0.1"

    HELP_MSG = """Quick help:

 Left-Click:   to display X,Y coordinates
 Left-Drag:    to zoom in on plot region
 Right-Click:  display popup menu with choices:
                Zoom out 1 level
                Zoom all the way out
                --------------------
                Configure
                Save Image

Also, these key bindings can be used
(For Mac OSX, replace 'Ctrl' with 'Apple'):

  Ctrl-Q:     quit
"""

    @property
    def codename(self):
        return "SigmaResultsViewer"

    @property
    def active_sigres(self):
        """The active SIGRES file i.e. the SIGRES associated to the active tab."""
        return self.active_tab.sigres

    @property
    def structure(self):
        """`Structure` associated to the active tab."""
        return self.active_sigres.structure

    @property
    def ebands(self):
        return self.active_sigres.ebands

    @property
    def ebands_list(self):
        """List of `ElectronBands`."""
        ebands_list = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            ebands_list.append(tab.sigres.ebands)
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
            paths.append(tab.sigres.filepath)
        return paths

    @property
    def nc_filepaths(self):
        """String with the absolute path of the netcdf file."""
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.sigres.filepath)
        return paths

    @property
    def sigres_filepaths(self):
        """
        Return a list with the absolute paths of the files 
        from which the `SigresFile` have been read.
        """
        paths = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            paths.append(tab.sigres.filepath)
        return paths

    def makeMenu(self):
        """Creates the main menu."""
        # Base menu
        menu_bar = super(SigresViewerFrame, self).makeMenu()

        # Add mixin menus.
        menu_bar.Append(self.CreateStructureMenu(), "Structure")
        menu_bar.Append(self.CreateEbandsMenu(), "Ebands")
        menu_bar.Append(self.CreateQpDataMenu(), "QP Data")
        menu_bar.Append(self.CreateToolsMenu(), "Tools")
        menu_bar.Append(self.CreateNetcdfMenu(), "Netcdf")

        # Help menu
        help_menu = self.makeHelpMenu()
        menu_bar.Append(help_menu, "Help")

        self.SetMenuBar(menu_bar)

    def CreateQpDataMenu(self):
        # Structure Menu ID's
        self.ID_QPGAPS_COMPARE = wx.NewId()
        self.ID_QPENES_COMPARE = wx.NewId()

        menu = wx.Menu()
        menu.Append(self.ID_QPGAPS_COMPARE, "Compare QP-gaps", "Compare QP gaps")
        self.Bind(wx.EVT_MENU, self.OnQpGapsCompare, id=self.ID_QPGAPS_COMPARE)

        menu.Append(self.ID_QPENES_COMPARE, "Compare QP-enes", "Compare QP energies")
        self.Bind(wx.EVT_MENU, self.OnQpEnesCompare, id=self.ID_QPENES_COMPARE)

        return menu

    def OnQpGapsCompare(self, event):
        # Instanciate the plotter and add the filepaths to the plotter.
        plotter = SIGRES_Plotter()
        plotter.add_files(self.sigres_filepaths)

        # Plot the convergence of the QP gaps.
        plotter.plot_qpgaps(title="Convergence of QP gaps", hspan=0.05)

    def OnQpEnesCompare(self, event):
        # Instanciate the plotter and add the filepaths to the plotter.
        plotter = SIGRES_Plotter()
        plotter.add_files(self.sigres_filepaths)
                                                                             
        # Plot the convergence of the QP energies.
        plotter.plot_qpenes(title="Convergence of QP energies", hspan=0.05)

    def makeToolBar(self):
        """Creates the toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        self.ID_SCISSORS = wx.NewId()
        self.ID_PLOTQPSE0 = wx.NewId()
        self.ID_PLOTKSWITHMARKS = wx.NewId()

        artBmp = wx.ArtProvider.GetBitmap
        #toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR), "Open")
        toolbar.AddSimpleTool(self.ID_PLOTQPSE0, bitmap("qpresults.png"), "Plot QPState Results.")
        toolbar.AddSimpleTool(self.ID_PLOTKSWITHMARKS, bitmap("qpmarkers.png"), "Plot KS energies with QPState markers.")
        toolbar.AddSimpleTool(self.ID_SCISSORS, bitmap("qpscissor.png"), "Build energy-dependent scissors from GW correction.")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_PLOTQPSE0, self.OnPlotQpsE0),
            (self.ID_PLOTKSWITHMARKS, self.OnPlotKSwithQPmarkers),
            (self.ID_SCISSORS, self.OnScissors),
        ]
                                                                   
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.toolbar.Realize()

    def addFileTab(self, parent, filepath):
        gsr = abiopen(filepath)
        tab = SigresFileTab(self.notebook, gsr)
        self.notebook.AddPage(tab, os.path.basename(filepath))

    def OnPlotQpsE0(self, event):
        """Plot QPState results as function of the KS energy."""
        if self.active_sigres is None: return
        self.active_sigres.plot_qps_vs_e0()

    def OnPlotKSwithQPmarkers(self, event):
        """Plot KS energies with QPState markers."""
        if self.active_sigres is None: return
        sigres = self.active_sigres

        band_range = (sigres.min_gwbstart, sigres.max_gwbstop)
        try:
            EbandsWithMarkersPlotFrame(self, sigres.ebands, band_range=band_range).Show()
        except:
            awx.showErrorMessage(self)

    def OnScissors(self, event):
        """Build the scissors operator."""
        if self.active_sigres is None: return
        ScissorsBuilderFrame(self, self.active_sigres.filepath).Show()

    def ShowQPTable(self, spin, kpoint, band):
        if self.active_sigres is None: return
        qplist = self.active_sigres.get_qplist(spin, kpoint)
        table = qplist.to_table()

        title = "spin: %d, kpoint: %s" % (spin, kpoint)
        awx.SimpleGridFrame(self, table[1:], col_labels=table[0], title=title).Show()


class SigresFileTab(wx.Panel):
    """Tab showing information on a single SIGRES file."""
    def __init__(self, parent, sigres, **kwargs):
        """
        Args:
            parent:
                parent window.
            sigres:
                `SigresFile` object
        """
        super(SigresFileTab, self).__init__(parent, -1, **kwargs)
        self.sigres = sigres

        splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
        splitter.SetSashGravity(0.95)

        self.skb_panel = awx.SpinKpointBandPanel(splitter, sigres.nsppol, sigres.gwkpoints, sigres.max_gwbstop,
            bstart=sigres.min_gwbstart)

        # Set the callback for double click on k-point row..
        self.skb_panel.SetOnItemActivated(self.ShowQPTable)

        # Add Python shell
        msg = "SIGRES_File object is accessible via the sigres variable. Use sigres.<TAB> to access the list of methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        pyshell = Shell(splitter, introText=msg, locals={"sigres": sigres})
        splitter.SplitHorizontally(self.skb_panel, pyshell)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

    @property
    def statusbar(self):
        return self.viewer_frame.statusbar

    @property
    def viewer_frame(self):
        """The parent frame `WfkViewerFrame`."""
        try:
            return self._viewer_frame
                                                                                    
        except AttributeError:
            self._viewer_frame = self.getParentWithType(SigresViewerFrame)
            return self._viewer_frame

    def ShowQPTable(self, spin, kpoint, band):
        qplist = self.sigres.get_qplist(spin, kpoint)
        table = qplist.to_table()
                                                                                       
        title = "spin: %d, kpoint: %s" % (spin, kpoint)
        awx.SimpleGridFrame(self, table[1:], col_labels=table[0], title=title).Show()


class EbandsWithMarkersPlotFrame(awx.Frame):

    def __init__(self, parent, ebands, **kwargs):
        """
        Args:
            parent:
                Parent window
            ebands:
                `ElectronBands` object.
        """
        self.band_range = kwargs.pop("band_range", None)

        super(EbandsWithMarkersPlotFrame, self).__init__(parent, -1, **kwargs)
        self.SetTitle("Select parameters")

        self.ebands = ebands

        if not ebands.markers:
            raise awx.Error("Found empty markers dictionary in ebands %s" % repr(ebands))

        self.BuildUi()

    def BuildUi(self):

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer1 = wx.BoxSizer( wx.HORIZONTAL )

        marker_label = wx.StaticText(self, -1, "Available Markers:")
        marker_label.Wrap(-1)
        hsizer1.Add(marker_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        marker_choices = list(self.ebands.markers.keys())
        self.marker_choice = wx.Choice( self, -1, wx.DefaultPosition, wx.DefaultSize, marker_choices, 0 )
        self.marker_choice.SetSelection(0)
        hsizer1.Add(self.marker_choice, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        scale_label = wx.StaticText(self, -1, "Scale Factor:")
        scale_label.Wrap(-1)
        hsizer1.Add(scale_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        self.scale_ctrl = wx.SpinCtrlDouble(self, id=-1, value=str(50), min=0.0, max=1.e+10, inc=50)
        self.scale_ctrl.SetToolTipString("Multiplicative factor used to render the markers more visible.")
        hsizer1.Add( self.scale_ctrl, 0, wx.ALL, 5 )

        main_sizer.Add( hsizer1, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )

        hsizer2 = wx.BoxSizer( wx.HORIZONTAL )

        ok_button = wx.Button(self, wx.ID_OK, label='Ok')
        cancel_button = wx.Button(self, wx.ID_CANCEL, label='Cancel')

        ok_button.Bind(wx.EVT_BUTTON, self.OnOkButton)
        cancel_button.Bind(wx.EVT_BUTTON, self.OnCloseButton)

        hsizer2.Add(ok_button, 0, wx.ALL, 5 )
        hsizer2.Add(cancel_button, 0, wx.ALL, 5 )

        main_sizer.Add( hsizer2, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.SetSizerAndFit(main_sizer)

    def OnCloseButton(self, event):
        self.Destroy()

    def GetParams(self):
        return AttrDict(
            qpattr=self.marker_choice.GetStringSelection(),
            fact=float(self.scale_ctrl.GetValue()))

    def OnOkButton(self, event):
        p = self.GetParams()
        with_marker = p.qpattr + ":" + str(p.fact)

        self.ebands.plot(marker=with_marker, band_range=self.band_range)


class SigresViewerApp(awx.App):
    pass


def wxapp_sigresviewer(sigres_filepaths):
    app = SigresViewerApp()
    frame = SigresViewerFrame(None, filepaths=sigres_filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    return app


#def wxapp_filewalker(main):
#    """
#    This decorator is used to decorate main functions producing `AbinitFlows`.
#    It adds the initialization of the logger and an argument parser that allows one to select 
#    the loglevel, the workdir of the flow as well as the YAML file with the paramenters of the `TaskManager`.
#    The main function shall have the signature:
#
#        main(options)
#
#    where options in the container with the command line options generated by `ArgumentParser`.
#
#    Args:
#        main:
#            main function.
#    """
#    from functools import wraps
#    @wraps(main)
#    def wrapper(*args, **kwargs):
#        import argparse
#        parser = argparse.ArgumentParser()
#
#        parser.add_argument('--loglevel', default="ERROR", type=str,
#                            help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
#
#        parser.add_argument("-r", '--recurse', default=False, type=bool, help="Recurse the file system starting from the local directory")
#
#        options = parser.parse_args()
#
#        # loglevel is bound to the string value obtained from the command line argument. 
#        # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
#        numeric_level = getattr(logging, options.loglevel.upper(), None)
#        if not isinstance(numeric_level, int):
#            raise ValueError('Invalid log level: %s' % options.loglevel)
#        logging.basicConfig(level=numeric_level)
#
#        if options.recurse:
#            # Walk the filesytem starting from top.
#            files = []
#        else:
#            files = options.args
#
#        return main(files, **kwargs)
#
#    return wrapper

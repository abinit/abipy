from __future__ import print_function, division

import os
import wx

import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx

from abipy.abilab import abiopen
from abipy.tools import AttrDict
from abipy.gui.scissors import ScissorsBuilderFrame
from abipy.gui.mixins import Has_Structure, Has_Ebands, Has_Tools

ID_SCISSORS = wx.NewId()
ID_PLOTQPSE0 = wx.NewId()
ID_PLOTKSWITHMARKS = wx.NewId()
ID_NCDUMP = wx.NewId()


class SigresViewerFrame(awx.Frame, Has_Structure, Has_Ebands):
    VERSION = "0.1"

    def __init__(self, parent, filename=None, **kwargs):
        super(SigresViewerFrame, self).__init__(parent, id=-1, title=self.codename, **kwargs)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing SIGRES file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the SIGRES file")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")
        file_menu.Append(ID_NCDUMP, "Ncdump", help="ncdump printout")
        menuBar.Append(file_menu, "File")

        menuBar.Append(self.CreateStructureMenu(), "Structure")
        menuBar.Append(self.CreateEbandsMenu(), "Ebands")
        menuBar.Append(self.CreateToolsMenu(), "Tools")

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

        tsize = (48, 48)
        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(ID_PLOTQPSE0, wx.Bitmap(awx.path_img("qpresults.png")), "Plot QPState Results.")
        toolbar.AddSimpleTool(ID_PLOTKSWITHMARKS, wx.Bitmap(awx.path_img("qpmarkers.png")), "Plot KS energies with QPState markers.")
        toolbar.AddSimpleTool(ID_SCISSORS, wx.Bitmap(awx.path_img("qpscissor.png")), "Build energy-dependent scissors from GW correction.")

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
            (ID_PLOTQPSE0, self.OnPlotQpsE0),
            (ID_PLOTKSWITHMARKS, self.OnPlotKSwithQPmarkers),
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

    @property
    def structure(self):
        return self.sigres.structure

    @property
    def ebands(self):
        return self.sigres.ebands

    def BuildUi(self):
        sigres = self.sigres
        if sigres is None: return

        splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
        splitter.SetMinimumPaneSize(50)
        parent = splitter
        #parent = self

        self.skb_panel = awx.SpinKpointBandPanel(parent, sigres.nsppol, sigres.gwkpoints, sigres.max_gwbstop,
            bstart=sigres.min_gwbstart)

        # Set the callback for double click on k-point row..
        self.skb_panel.SetOnItemActivated(self.ShowQPTable)

        # Add Python shell
        from wx.py.shell import Shell
        from abipy.tools import marquee
        msg = "SIGRES_File object is accessible via the sigres variable. Use sigres.<TAB> to access the list of the methods."
        msg = marquee(msg, width=len(msg) + 8, mark="#")
        msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        pyshell = Shell(parent, introText=msg, locals={"sigres": sigres})
        splitter.SplitHorizontally(self.skb_panel, pyshell)

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
            #if not isinstance(sigres, SIGRES_File):
            #    awx.showErrorMessage(self, message="%s is not a valid SIGRES file" % filepath)
            #    return

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

    def OnPlotQpsE0(self, event):
        """Plot QPState results as function of the KS energy."""
        if self.sigres is None: return
        self.sigres.plot_qps_vs_e0()

    def OnPlotKSwithQPmarkers(self, event):
        """Plot KS energies with QPState markers."""
        if self.sigres is None: return

        band_range = (self.sigres.min_gwbstart, self.sigres.max_gwbstop)
        try:
            BandsWithMarkersPlotFrame(self, self.sigres.ebands, band_range=band_range).Show()
        except awx.Error as exc:
            awx.showErrorMessage(self)

    def OnScissors(self, event):
        """Build the scissors operator."""
        if self.sigres is None: return
        ScissorsBuilderFrame(self, self.sigres.filepath).Show()

    def ShowQPTable(self, spin, kpoint, band):
        qplist = self.sigres.get_qplist(spin, kpoint)
        table = qplist.to_table()

        title = "spin: %d, kpoint: %s" % (spin, kpoint)
        frame = awx.SimpleGridFrame(self, table[1:], col_labels=table[0], title=title)
        frame.Show()

    def OnNcdump(self, event):
        """Call ncdump and show results in a dialog."""
        if self.sigres is None: return
        caption = "ncdump output for SIGRES file %s" % self.sigres.filepath
        wxdg.ScrolledMessageDialog(self, self.sigres.ncdump(), caption=caption, style=wx.MAXIMIZE_BOX).Show()



class BandsWithMarkersPlotFrame(awx.Frame):

    def __init__(self, parent, bands, **kwargs):
        self.band_range = kwargs.pop("band_range", None)

        super(BandsWithMarkersPlotFrame, self).__init__(parent, -1, **kwargs)
        self.SetTitle("Select parameters")

        self.bands = bands

        if not bands.markers:
            raise awx.Error("Found empty markers dictionary in bands %s" % repr(bands))

        self.BuildUi()

    def BuildUi(self):

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer1 = wx.BoxSizer( wx.HORIZONTAL )

        marker_label = wx.StaticText(self, -1, "Available Markers:")
        marker_label.Wrap(-1)
        hsizer1.Add(marker_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        marker_choices = list(self.bands.markers.keys())
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
            fact=float(self.scale_ctrl.GetValue()),
        )

    def OnOkButton(self, event):
        p = self.GetParams()
        with_marker = p.qpattr + ":" + str(p.fact)

        self.bands.plot(marker=with_marker, band_range=self.band_range)


class SigresViewerApp(awx.App):

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        #logger.info("%s dropped on app %s" % (filename, self.appname))
        frame = SigresViewerFrame(parent=None, filename=filename)
        frame.Show()


def wxapp_sigresviewer(sigres_filename):
    app = SigresViewerApp()
    frame = SigresViewerFrame(None, filename=sigres_filename)
    frame.Show()
    app.SetTopWindow(frame)
    return app

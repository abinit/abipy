from __future__ import print_function, division

import os
import wx

import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy import abiopen, SIGRES_File
from abipy.tools import AttrDict
from abipy.iotools.visualizer import supported_visunames
from abipy.gui.scissors import ScissorsBuilderFrame

ID_VISTRUCT = wx.NewId()
ID_VISBZ = wx.NewId()
ID_NCDUMP = wx.NewId()
ID_SCISSORS = wx.NewId()
ID_PLOTQPSE0 = wx.NewId()
ID_PLOTKSWITHMARKS = wx.NewId()
ID_TBOX_VIS = wx.NewId()

class SigresViewerFrame(awx.Frame):
    VERSION = "0.1"

    def __init__(self, parent, filename=None):
        super(SigresViewerFrame, self).__init__(parent, id=-1, title=self.codename)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing SIGRES file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the SIGRES file")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")
        file_menu.Append(ID_NCDUMP, "Ncdump", help="ncdump printout")
        menuBar.Append(file_menu, "File")

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

        tsize = (15, 15)
        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(ID_VISTRUCT, wx.Bitmap(awx.path_img("crystal.png")), "Visualize the crystal structure")
        toolbar.AddSimpleTool(ID_VISBZ, wx.Bitmap(awx.path_img("wave.png")), "Visualize the BZ")
        toolbar.AddSimpleTool(ID_PLOTQPSE0, wx.Bitmap(awx.path_img("wave.png")), "Plot QP Results.")
        toolbar.AddSimpleTool(ID_PLOTKSWITHMARKS, wx.Bitmap(awx.path_img("wave.png")), "Plot KS energies with QP markers.")
        toolbar.AddSimpleTool(ID_SCISSORS, wx.Bitmap(awx.path_img("wave.png")), "Build energy-dependent scissors from GW correction.")

        toolbar.AddSeparator()
        self.visualizer_cbox = wx.ComboBox(choices=supported_visunames(), id=ID_TBOX_VIS, 
            name='visualizer', parent=toolbar, value='xcrysden') 
        self.visualizer_cbox.Refresh() 

        toolbar.AddControl(control=self.visualizer_cbox) 

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
            (ID_VISTRUCT, self.OnVisualizeStructure),
            (ID_VISBZ, self.OnVisualizeBZ),
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
    def ks_bands(self):
        if self.sigres is None: 
            return None
        else:
            return self.sigres.ks_bands

    def BuildUi(self):
        sigres = self.sigres
        if sigres is None: return

        #splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
        #splitter.SetMinimumPaneSize(50)
        #parent = splitter 

        parent = self
        self.skb_panel = awx.SpinKpointBandPanel(parent, sigres.nsppol, sigres.gwkpoints, sigres.max_gwbstop, 
            bstart=sigres.min_gwbstart)

        # Set the callback for double click on k-point row..
        self.skb_panel.SetOnItemActivated(self.ShowQPTable)

        # Add Python shell
        #from wx.py.shell import Shell
        #pyshell = Shell(parent, introText="SIGRES File methods are available via the variable sigres ", locals={"sigres": sigres})
        #splitter.SplitHorizontally(self.skb_panel, pyshell)

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
            if not isinstance(sigres, SIGRES_File):
                awx.showErrorMessage(self, message="%s is not a valid SIGRES file" % filepath)
                return

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
        """Plot QP results as function of the KS energy."""
        if self.sigres is None: return
        self.sigres.plot_qps_vs_e0()

    def OnPlotKSwithQPmarkers(self, event):
        """Plot KS energies with QP markers."""
        if self.sigres is None: return
        QPAttrPlotFrame(self, self.sigres).Show()

    def OnScissors(self, event):
        """Build the scissors operator."""
        if self.sigres is None: return
        ScissorsBuilderFrame(self, self.sigres.filepath).Show()

    def GetVisualizer(self):
        """Returns a string with the visualizer selected by the user."""
        return self.visualizer_cbox.GetValue()

    def OnVisualizeStructure(self, event):
        """"Call visualizer to visualize the crystalline structure."""
        if self.sigres is None: return

        visualizer = self.GetVisualizer()
        self.statusbar.PushStatusText("Visualizing crystal structure with %s" % visualizer)
        structure = self.sigres.get_structure()
        try:
            visu = structure.visualize(visualizer)
            visu()
        except:
            awx.showErrorMessage(self)

    def OnVisualizeBZ(self, event):
        """"Visualize the Brillouin zone with matplotlib."""
        if self.sigres is None: return
        self.sigres.get_structure().show_bz()

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



class QPAttrPlotFrame(awx.Frame):

    def __init__(self, parent, sigres, **kwargs):
        super(QPAttrPlotFrame, self).__init__(parent, -1, **kwargs)
        self.SetTitle("Select parameters")

        self.sigres = sigres

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Construct a panel for each spin.
        self.panel = QPAttrChoicePanel(self)

        main_sizer.Add(self.panel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        hbox = wx.BoxSizer(wx.HORIZONTAL)

        ok_button = wx.Button(self, wx.ID_OK, label='Ok')
        close_button = wx.Button(self, wx.ID_CANCEL, label='Cancel')

        ok_button.Bind(wx.EVT_BUTTON, self.OnOkButton)
        close_button.Bind(wx.EVT_BUTTON, self.OnCloseButton)

        hbox.Add(ok_button)
        hbox.Add(close_button, flag=wx.LEFT, border=5)

        main_sizer.Add(hbox, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizerAndFit(main_sizer)

    def OnCloseButton(self, event):
        self.Destroy()

    def OnOkButton(self, event):
        p = self.panel.GetParams()
        self.sigres.plot_ksbands_with_qpmarkers(qpattr=p.qpattr, fact=p.fact)

class QPAttrChoicePanel(awx.Panel):

    def __init__(self, parent, **kwargs):
        super(QPAttrChoicePanel, self).__init__(parent, -1, **kwargs)

        hsz1 =  wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, -1, "QP attribute:", wx.DefaultPosition, wx.DefaultSize, 0)
        label.Wrap(-1)
        hsz1.Add(label, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        from abipy.electrons import QP
        qp_choices = QP.get_fields(exclude=("spin", "kpoint", "band"))
        self.attr_choice = attr_choice = wx.Choice(self, -1, wx.DefaultPosition, wx.DefaultSize, qp_choices, 0)
        attr_choice.SetSelection(0)
        attr_choice.SetToolTipString("Select the quantity to use for the markers.")
        hsz1.Add(attr_choice, 0, wx.ALL, 5)

        hsz2 =  wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, -1, "Scale Factor:")
        self.fact_ctrl = wx.SpinCtrlDouble(self, id=-1, value=str(5), min=0.0, inc=1)
        self.fact_ctrl.SetToolTipString("Multiplicative factor used to render the markers more visible.")
        hsz2.Add(label, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        hsz2.Add(self.fact_ctrl, 0, wx.ALL, 5)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(hsz1)
        main_sizer.Add(hsz2)
        self.SetSizerAndFit(main_sizer)

    def GetParams(self):
        return AttrDict(
            qpattr=self.attr_choice.GetStringSelection(),
            fact=float(self.fact_ctrl.GetValue()),
        )


class SigresViewerApp(awx.App):

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        self.log("%s dropped on app %s" % (filename, self.appname))
        frame = SigresViewerFrame(parent=None, filename=filename)
        frame.Show()


def wxapp_sigresviewer(sigres_filename):
    app = SigresViewerApp()
    frame = SigresViewerFrame(None, filename=sigres_filename)
    frame.Show()
    app.SetTopWindow(frame)
    return app


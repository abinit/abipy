from __future__ import print_function, division

import os
import wx

import abipy.gui.awx as awx

#from abipy import abiopen

ID_SHOW_INPUTS = wx.NewId()
ID_SHOW_OUTPUTS = wx.NewId()
ID_SHOW_LOGS = wx.NewId()
ID_BROWSE = wx.NewId()

class WorkflowViewerFrame(awx.Frame):
    VERSION = "0.1"

    def __init__(self, parent, work):
        super(WorkflowViewerFrame, self).__init__(parent, -1, self.codename)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing WFK file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the WFK file")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")
        menuBar.Append(file_menu, "File")

        #file_history = self.file_history = wx.FileHistory(8)
        #self.config = wx.Config(self.codename, style=wx.CONFIG_USE_LOCAL_FILE)
        #file_history.Load(self.config)
        #recent = wx.Menu()
        #file_history.UseMenu(recent)
        #file_history.AddFilesToMenu()
        #file_menu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        #self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)

        self.help_menu = wx.Menu()
        self.help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.help_menu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        self.toolbar = toolbar = self.CreateToolBar()

        tsize = (15, 15)
        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(ID_SHOW_INPUTS, wx.Bitmap(awx.path_img("struct.png")), "Visualize the input files.")
        toolbar.AddSimpleTool(ID_SHOW_OUTPUTS, wx.Bitmap(awx.path_img("struct.png")), "Visualize the output files.")
        toolbar.AddSimpleTool(ID_SHOW_LOGS, wx.Bitmap(awx.path_img("struct.png")), "Visualize the log files.")
        toolbar.AddSimpleTool(ID_BROWSE, wx.Bitmap(awx.path_img("struct.png")), "Browse files.")

        #toolbar.AddSeparator()
        #self.visualizer_cbox = wx.ComboBox(choices=supported_visunames(), id=ID_TBOX_VIS, 
        #    name='visualizer', parent=toolbar, value='xcrysden') 
        #self.visualizer_cbox.Refresh() 
        #toolbar.AddControl(control=self.visualizer_cbox) 

        self.toolbar.Realize()
        self.Centre()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            #(wx.ID_OPEN, self.OnOpen),
            #(wx.ID_CLOSE, self.OnClose),
            #(wx.ID_EXIT, self.OnExit),
            (wx.ID_ABOUT, self.OnAboutBox),
            #
            (ID_SHOW_INPUTS, self.OnShowInputs),
            (ID_SHOW_OUTPUTS, self.OnShowOutputs),
            (ID_SHOW_LOGS, self.OnShowLogs),
            (ID_BROWSE, self.OnBrowse),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.work = work
        #if filename is not None:
        #    self.ReadWfkFile(filename)

    @property
    def codename(self):
        return self.__class__.__name__

    def BuildUi(self):
        if self.work is None: return
        work = self.work

        #splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
        #splitter.SetMinimumPaneSize(50)
        #parent = splitter 
        #parent = self

        #self.skb_panel = awx.SpinKpointBandPanel(parent, wfk.nsppol, wfk.kpoints, wfk.mband)

        # Set the callback for double click on k-point row..
        #self.skb_panel.SetOnItemActivated(self._visualize_skb)

        # Add Python shell
        #from wx.py.shell import Shell
        #from abipy.tools import marquee
        #msg = "WFK_File object is accessible via the wfk variable. Use wfk.<TAB> to access the list of methods."
        #msg = marquee(msg, width=len(msg) + 8, mark="#")
        #msg = "#"*len(msg) + "\n" + msg + "\n" + "#"*len(msg) + "\n"

        #pyshell = Shell(parent, introText=msg, locals={"work": self.work})
        #splitter.SplitHorizontally(self.skb_panel, pyshell)

        #main_sizer = wx.BoxSizer(wx.VERTICAL)
        #main_sizer.Add(self.skb_panel, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)
        #main_sizer.Add(pyshell, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)
        #self.SetSizerAndFit(main_sizer)

    #def DestroyPanel(self):
    #    if hasattr(self, "panel"):
    #        self.panel.Destroy()

    #def OnOpen(self, event):
    #    dlg = wx.FileDialog(self, message="Choose a __workflow__.pickle file", defaultDir=os.getcwd(),
    #                        wildcard="pickle files (*.pickle)|*.pickle",
    #                        style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
    #    )

    #    # Show the dialog and retrieve the user response. 
    #    # If it is the OK response, process the data.
    #    if dlg.ShowModal() == wx.ID_OK:
    #        filepath = dlg.GetPath()
    #        self.file_history.AddFileToHistory(filepath)
    #        self.file_history.Save(self.config)
    #        self.config.Flush()
    #        self.ReadWfkFile(filepath)

    #    dlg.Destroy()

    #def OnFileHistory(self, event):
    #    fileNum = event.GetId() - wx.ID_FILE1
    #    filepath = self.file_history.GetHistoryFile(fileNum)
    #    # move up the list
    #    self.file_history.AddFileToHistory(filepath)
    #    self.ReadWfkFile(filepath)

    #def ReadWfkFile(self, filepath):
    #    """Read the WFK file and build the UI."""
    #    self.statusbar.PushStatusText("Reading %s" % filepath)

    #    try:
    #        work = abiopen(filepath)
    #        #if not isinstance(wfkfile, WFK_File):
    #        #    awx.showErrorMessage(self, message="%s is not a valid WFK File" % filepath)
    #        #    return

    #        self.work = wfkfile
    #        self.BuildUi()
    #        self.statusbar.PushStatusText("Workflow file %s loaded" % filepath)

    #    except Exception:
    #        awx.showErrorMessage(self)

    #def OnClose(self, event):
    #    self.work = None
    #    self.DestroyPanel()

    #def OnExit(self, event):
    #    self.Close()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")

    def OnShowInputs(self, event):
        if self.work is None: return
        self.work.wxshow_inputs()

    def OnShowOutputs(self, event):
        if self.work is None: return
        self.work.wxshow_outputs()

    def OnShowLogs(self, event):
        if self.work is None: return
        self.work.wxshow_logs()

    def OnBrowse(self, event):
        if self.work is None: return
        self.work.wxbrowse()


class WorkflowViewerApp(awx.App):
    def OnInit(self):
        return True

    def MacOpenFile(self, filename):
        """Called for files droped on dock icon, or opened via finders context menu"""
        if filename.endswith(".py"):
            return
        # Open filename in a new frame.
        self.log("%s dropped on app %s" % (filename, self.appname))
        WorkflowViewerFrame(parent=None, work=filename).Show()


def wxapp_workflow_viewer(workflow_or_picklefile):
    app = WorkflowViewerApp()
    frame = WorkflowViewerFrame(None, work=workflow_or_picklefile)
    frame.Show()
    return app

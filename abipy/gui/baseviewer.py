from __future__ import print_function, division

import os
import wx
import abc

import wx.lib.agw.flatnotebook as fnb
import abipy.gui.awx as awx

#from abipy.tools import marquee, list_strings 
#from abipy.abilab import abiopen


class MultiViewerFrame(awx.Frame):
    __metaclass__ = abc.ABCMeta
    """
    Concrete classes must define the following attributes.

        VERSION = "0.1"
        HELP_MSG = 'Quick help'
    """
    def __init__(self, parent, filepaths=(), **kwargs):

        super(MultiViewerFrame, self).__init__(parent, -1, title=self.codename, **kwargs)

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
            self.read_file(path)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.EXPAND, 5)
        panel.SetSizerAndFit(sizer)

    @abc.abstractproperty
    def codename(self):
        """The name of the viewer."""

    @property
    def active_tab(self):
        """Returns the active tab. None if notebook is empty."""
        return self.notebook.GetCurrentPage()

    def makeMenu(self):
        """
        Method of the base class that provides a base menu.
        that can be extended by the subclass.
        """
        self.menu_bar = menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing file in a new tab")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the file associated to the active tab")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")

        file_history = self.file_history = wx.FileHistory(8)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)
        menuBar.Append(file_menu, "File")

        #menuBar.Append(self.CreateStructureMenu(), "Structure")
        #menuBar.Append(self.CreateEbandsMenu(), "Ebands")
        #menuBar.Append(self.CreateToolsMenu(), "Tools")
        #menuBar.Append(self.CreateNetcdfMenu(), "Netcdf")

        self.help_menu = wx.Menu()

        self.ID_HELP_QUICKREF = wx.NewId()
        self.help_menu.Append(self.ID_HELP_QUICKREF, "Quick Reference ", help="Quick reference for " + self.codename)
        self.help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")

        menuBar.Append(self.help_menu, "Help")

        self.SetMenuBar(menuBar)

        # Associate menu/toolbar items with their handlers.
        self.ID_VISWAVE = wx.NewId()

        menu_handlers = [
            (wx.ID_OPEN, self.OnOpen),
            (wx.ID_CLOSE, self.OnClose),
            (wx.ID_EXIT, self.OnExit),
            (wx.ID_ABOUT, self.OnAboutBox),
            (self.ID_QUICKREF, self.OnQuickRef),
        ]
                                                            
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

    @abc.abstractmethod
    def makeToolBar(self):
        """To be provided by the concrete class."""

    @abc.abstractmethod
    def addFileTab(self, parent, file_obj):
        pass

    def read_file(self, filepath):
        """Open netcdf file, create new tab and save the file in the history."""
        self.statusbar.PushStatusText("Reading %s" % filepath)
        try:
            notebook = self.notebook
            wfk = abiopen(filepath)
            #if not isinstance(wfkfile, WFK_File):
            #    awx.showErrorMessage(self, message="%s is not a valid WFK File" % filepath)
            #    return
            self.statusbar.PushStatusText("File %s loaded" % filepath)
            tab = WfkFileTab(notebook, wfk)
            notebook.AddPage(tab, os.path.basename(filepath))
            # don't know why but this does not work!
            notebook.Refresh()
            notebook.SetSelection(notebook.GetPageCount())
            self.AddFileToHistory(filepath)
        except:
            awx.showErrorMessage(self)

    def OnOpen(self, event):
        """Open FileDialog to allow the user to select a file."""
        # Show the dialog and retrieve the user response.
        # If it is the OK response, process the data.
        dialog = wx.FileDialog(self, message="Choose a WFK file", defaultDir=os.getcwd(),
                               wildcard="WFK Netcdf files (*.nc)|*.nc",
                               style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)
        if dialog.ShowModal() == wx.ID_CANCEL: return 
                                                                                          
        self.read_file(dialog.GetPath())

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
        #tab.wfk.close()
                                                                                          
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

    def AddFileToHistory(self, filepath):
        """Add the absolute filepath to the file history."""
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        self.read_file(filepath)

import os
import wx
import abc
import wx.lib.agw.flatnotebook as fnb
import abipy.gui.awx as awx

from monty.string import list_strings


class MultiViewerFrame(awx.Frame, metaclass=abc.ABCMeta):
    """
    Base class for Viewers that can handle the multiple netcdf files
    of the same type. A `MultiViewerFrame` has a notebook where
    each page provides tools to interact with the data stored in the netcd file
    Concrete classes should provide the method that open a new file and creates
    a new page (specific to the file type) that will be added to the notebook.

    Concrete classes must also define the following (class) attributes.

    .. attributes:

        VERSION:
            String with the version of the visualizer e.g  "0.1"
        HELP_MSG:
            Multi line string with a Quick help. For example a list of shortcuts
            or a brief description of the usage of the GUI.
    """
    def __init__(self, parent, filepaths=(), **kwargs):
        """
        Args:
            parent:
                parent window.
            filepaths:
                String or list of strings with the path of the netcdf files to open
                Empty tuple if no file should be opened during the initialization of the frame.
        """
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

        try:
            style = fnb.FNB_NO_X_BUTTON  | fnb.FNB_NAV_BUTTONS_WHEN_NEEDED
        except AttributeError:
            style = fnb.FNB_NO_X_BUTTON 

        self.notebook = fnb.FlatNotebook(panel, -1, style=style)
                                                                                           
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
        Method of the base class that provides a base menu
        that can be extended by the subclass.
        """
        menu_bar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing file in a new tab")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the file associated to the active tab")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")

        file_history = self.file_history = wx.FileHistory(8)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(-1, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)
        menu_bar.Append(file_menu, "File")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_OPEN, self.OnOpen),
            (wx.ID_CLOSE, self.OnClose),
            (wx.ID_EXIT, self.OnExit),
        ]
                                                            
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu_bar

    def makeHelpMenu(self):
        """
        Method of the base class that provides a base help menu
        that can be extended by the subclass.
        """
        help_menu = wx.Menu()
                                                                                                                 
        self.ID_HELP_QUICKREF = wx.NewId()
        help_menu.Append(self.ID_HELP_QUICKREF, "Quick Reference ", help="Quick reference for " + self.codename)
        help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_HELP_QUICKREF, self.onQuickRef),
            (wx.ID_ABOUT, self.OnAboutBox),
        ]
                                                            
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)
                                                     
        return help_menu

    @abc.abstractmethod
    def makeToolBar(self):
        """To be provided by the concrete class."""

    @abc.abstractmethod
    def addFileTab(self, parent, filepath):
        """
        This method must be provided by the subclass. 
        It receives a string with the file path, opens the file
        and adds a new tab to the notebook.

        Example::

            wfk = abiopen(filepath)
            tab = WfkFileTab(self.notebook, wfk)
            self.notebook.AddPage(tab, os.path.basename(filepath))
        """

    def read_file(self, filepath):
        """Open netcdf file, create new tab and save the file in the history."""
        self.statusbar.PushStatusText("Reading %s" % filepath)
        try:
            self.addFileTab(self, filepath)
            # don't know why but this does not work!
            self.notebook.Refresh()
            self.notebook.SetSelection(self.notebook.GetPageCount())
            self.AddFileToHistory(filepath)
        except:
            awx.showErrorMessage(self)

    def OnOpen(self, event):
        """Open FileDialog to allow the user to select a file."""
        # intercept possible problems when determining the current directory (e.g. it has been deleted)
        try:
            def_dir = os.getcwd()
        except OSError:
            def_dir = os.path.expanduser('~')
        # Show the dialog and retrieve the user response.
        # If it is the OK response, process the data.
        dialog = wx.FileDialog(self, message="Choose a netcdf file", defaultDir=def_dir,
                               wildcard="Netcdf files (*.nc)|*.nc",
                               style=wx.OPEN | wx.CHANGE_DIR)

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

    def onQuickRef(self, event=None):
        """Show a dialog with a brief description of the commands."""
        dialog = wx.MessageDialog(self, self.HELP_MSG, self.codename + " Quick Reference",
                               wx.OK | wx.ICON_INFORMATION)
        dialog.ShowModal()
        dialog.Destroy()

    def AddFileToHistory(self, filepath):
        """Add the absolute filepath to the file history."""
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        self.read_file(filepath)

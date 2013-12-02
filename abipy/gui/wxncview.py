from __future__ import print_function, division

import os
import numpy as np
import wx
import netCDF4

import abipy.gui.awx as awx
import wx.lib.mixins.listctrl as listmix
import wx.lib.agw.flatnotebook as fnb
#import wx.lib.dialogs as wxdg

from wxmplot import PlotFrame, ImageFrame
from abipy.tools.text import list_strings #is_string , 
#from collections import OrderedDict


class NcViewFrame(wx.Frame):

    def __init__(self, parent, filepaths=(), **kwargs):
        """
            parent:
                parent window.
            filepaths:
                String or list of strings with the path of the nc files. 
        """
        if "size" not in kwargs:
            kwargs["size"] = (800, 600)

        super(NcViewFrame, self).__init__(parent, -1, **kwargs)

        self.makeMenu()
        self.makeToolBar()
        self.statusbar = self.CreateStatusBar()

        filepaths = list_strings(filepaths)
        self.datasets = []
        for path in filepaths:
            self.datasets.append(netCDF4.Dataset(path, mode="r"))

        # Create the notebook contents
        panel = wx.Panel(self, -1)
        #self.notebook = fnb.FlatNotebook(panel, -1, style=fnb.FNB_NO_X_BUTTON | fnb.FNB_NAV_BUTTONS_WHEN_NEEDED)
        self.notebook = fnb.FlatNotebook(panel, -1, style=fnb.FNB_NAV_BUTTONS_WHEN_NEEDED)

        for dataset in self.datasets:
            tab = NcFileTab(self.notebook, dataset)
            self.notebook.AddPage(tab, dataset.filepath())

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.EXPAND, 5)
        panel.SetSizerAndFit(sizer)

        self.Bind(wx.EVT_CLOSE, self.OnExit)
        #self.Bind(events.EVT_DIRECTORY_CHANGE, self.onDirectoryChange)

    def makeMenu(self):
        """Creates the main menu."""
        self.menu_bar = wx.MenuBar()

        self.menu_bar.Append(self.createFileMenu(), "&File")
        self.menu_bar.Append(self.createHelpMenu(), "&Help")

        self.SetMenuBar(self.menu_bar)

    def createFileMenu(self):
        """Creates the file menu."""
        file_menu = wx.Menu()

        # File Menu ID's
        self.idOPEN = wx.NewId()
        self.idCLOSE = wx.NewId()
        self.idEXIT = wx.NewId()

        file_menu.Append(self.idOPEN, "&Open", "Open netcdf file")
        file_menu.Append(self.idCLOSE, "&Close", "Close netcdf file")
        file_menu.AppendSeparator()
        file_menu.Append(self.idEXIT, "E&xit", "Exit wxnciew")

        self.Bind(wx.EVT_MENU, self.OnOpen, id=self.idOPEN)
        self.Bind(wx.EVT_MENU, self.OnClose, id=self.idCLOSE)
        self.Bind(wx.EVT_MENU, self.OnExit, id=self.idEXIT)
        return file_menu

    def createHelpMenu(self):
        """Create the help menu."""
        help_menu = wx.Menu()
        #help_menu.Append(self.idABOUT, "&About", "About Ftpcube")
        #self.Bind(wx.EVT_MENU, self.onAbout, id=self.idABOUT)
        return help_menu

    def makeToolBar(self):
        """Creates the toolbar."""
        self.toolbar = self.CreateToolBar()
        self.toolbar.SetToolBitmapSize(wx.Size(32, 32))

        #bitmap = icons.connect.getBitmap()
        #self.toolbar.AddSimpleTool(self.idTOOL_CONNECT, bitmap, shortHelpString=_("Normal Connect"))
        #bitmap = icons.quick_connect.getBitmap()

        #self.dir_entry_box = wx.ComboBox(self.toolbar, self.idTOOL_ENTRYBOX,
        #    size=wx.Size(350, -1), style=wx.CB_DROPDOWN)
        #self.binary_button = wx.RadioButton(self.toolbar, self.idTOOL_BINARY, _("BINARY"))
        #self.ascii_button = wx.RadioButton(self.toolbar, self.idTOOL_ASCII, _("ASCII"))
        #config = utils.getApplicationConfiguration()
        #if config['transfer_type'] == protocol.ProtocolInterface.BINARY:
        #    self.binary_button.SetValue(True)
        #else:
        #    self.ascii_button.SetValue(True)

        #utils.addHorizontalSpaceTool(self.toolbar, 75)
        #self.toolbar.AddControl(self.dir_entry_box)
        #utils.addHorizontalSpaceTool(self.toolbar, 5)
        #self.toolbar.AddControl(self.binary_button)
        #utils.addHorizontalSpaceTool(self.toolbar, 5)
        #self.toolbar.AddControl(self.ascii_button)

        #self.Bind(wx.EVT_TOOL, self.onConnect, id=self.idTOOL_CONNECT)
        #self.Bind(wx.EVT_TOOL, self.onQuickConnect, id=self.idTOOL_QUICKCONNECT)
        #self.Bind(wx.EVT_TOOL, self.onDisconnect, id=self.idTOOL_DISCONNECT)
        #self.Bind(wx.EVT_TOOL, self.onAbort, id=self.idTOOL_ABORT)
        #self.Bind(wx.EVT_TOOL, self.onRefresh, id=self.idTOOL_REFRESH)
        #self.Bind(wx.EVT_TOOL, self.onBookmarks, id=self.idTOOL_BOOKMARK)
        #self.Bind(wx.EVT_TOOL, self.onQuit, id=self.idTOOL_QUIT)
        #self.Bind(wx.EVT_COMBOBOX, self.onComboSelect, id=self.idTOOL_ENTRYBOX)
        #self.Bind(wx.EVT_RADIOBUTTON, self.onBinarySelect, id=self.idTOOL_BINARY)
        #self.Bind(wx.EVT_RADIOBUTTON, self.onAsciiSelect, id=self.idTOOL_ASCII)
        self.toolbar.Realize()

    #def onAbout(self, event):
    #    """Displays the about window box."""
    #    about = aboutwin.AboutWindow(self)
    #    ret = about.ShowModal()

    #    awx.makeAboutBox(codename=self.codename, version=self.VERSION,
    #                     description="", developers="M. Giantomassi")

    def OnOpen(self, event):
        """Open netcdf file."""
        dialog = wx.FileDialog(self, wildcard="*.nc")
        if dialog.ShowModal() == wx.ID_CANCEL: return 

        path = dialog.GetPath()

        #self.notebook.AddPage(text_ctrl, filename, select=True)

    def OnClose(self, event):
        """Close netcdf file, remote the tab from the notebook."""
        notebook = self.notebook
        if notebook.GetPageCount() == 0: return
        idx = notebook.GetSelection()
        if idx == -1: return None
        print("idx", idx)

        #notebook.DeletePage(idx)

        # See
        # http://stackoverflow.com/questions/16854332/how-to-delete-a-notebook-page?rq=1
        #def delPage(self, pageTitle):
        #    for index in range(self.dataNoteBook.GetPageCount()):
        #        if self.dataNoteBook.GetPageText(index) == pageTitle:
        #            self.dataNoteBook.DeletePage(index)
        #            self.dataNoteBook.SendSizeEvent()
        #            break


    def OnExit(self, event):
        """Exits the application."""
        self.Destroy()


class NcFileTab(wx.Panel):
    """Tab for one netcdf file."""
    def __init__(self, parent, dataset, **kwargs):
        super(NcFileTab, self).__init__(parent, -1, **kwargs)
        self.dataset = dataset

        #splitter = wx.SplitterWindow(self, -1, style=wx.SP_3DSASH)
        splitter = wx.SplitterWindow(self, -1, style=wx.SP_LIVE_UPDATE)
        splitter.SetMinimumPaneSize(100)

        self.dims_panel = DimsPanel(splitter, dataset)
        self.vars_panel = VarsPanel(splitter, dataset)
        splitter.SplitVertically(self.dims_panel, self.vars_panel)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)


class DimsPanel(wx.Panel, listmix.ColumnSorterMixin):
    """Panel with the netcdf dimensions."""
    def __init__(self, parent, dataset, **kwargs):
        super(DimsPanel, self).__init__(parent, -1, **kwargs)

        columns = ["name", "value"]
        self.list = wx.ListCtrl(self, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)

        for col, name in enumerate(columns):
            self.list.InsertColumn(col, name)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in columns]

        self.id2entry = {}
        for idx, (name, dim) in enumerate(dataset.dimensions.items()):
            entry = map(str, [name, len(dim)])
            self.list.Append(entry)
            self.list.SetItemData(idx, idx)
            self.id2entry[idx] = entry
                                                                                      
            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))
                                                                                      
        for index, col in enumerate(columns):
            self.list.SetColumnWidth(index, column_widths[index])

        # Now that the list exists we can init the other base class, see wx/lib/mixins/listctrl.py
        self.itemDataMap = self.id2entry
        listmix.ColumnSorterMixin.__init__(self, len(columns))
        self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.list)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.list, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self.list

    def OnColClick(self, event):
        event.Skip()


class VarsPanel(wx.Panel, listmix.ColumnSorterMixin):
    """Panel with the netcdf variables.."""
    def __init__(self, parent, dataset, **kwargs):
        super(VarsPanel, self).__init__(parent, -1, **kwargs)

        self.dataset = dataset

        columns = ["name", "dtype", "shape", "dimensions"]
        self.list = wx.ListCtrl(self, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)

        for col, name in enumerate(columns):
            self.list.InsertColumn(col, name)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in columns]

        self.id2entry = {}
        for idx, (name, var) in enumerate(dataset.variables.items()):
            entry = map(str, [name, var.dtype, var.shape, var.dimensions])
            self.list.Append(entry)
            self.list.SetItemData(idx, idx)
            self.id2entry[idx] = entry
                                                                                      
            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))
                                                                                      
        for index, col in enumerate(columns):
            self.list.SetColumnWidth(index, column_widths[index])

        # Now that the list exists we can init the other base class, see wx/lib/mixins/listctrl.py
        self.itemDataMap = self.id2entry
        listmix.ColumnSorterMixin.__init__(self, len(columns))
        self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.list)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.list, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

        # Connect the events whose callback will be set by the client code.
        self.list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)               

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self.list

    def OnColClick(self, event):
        event.Skip()

    def GetNameVarFromEvent(self, event):
        currentItem = event.m_itemIndex
        entry = self.id2entry[self.list.GetItemData(currentItem)]
        name = entry[0]
        return name, self.dataset.variables[name]

    def OnRightClick(self, event):
        name, var = self.GetNameVarFromEvent(event)
        print(var)

    def OnItemActivated(self, event):
        name, var = self.GetNameVarFromEvent(event)
        #print(var)
        #VarInfo(self, var).Show()
        plot_variable(self, name, var, self.dataset)


def plot_variable(parent, var_name, var, dataset):
    #import matplotlib.pyplot as plt

    # Remove fake dimensions.
    shape, dimensions = [], []
    for num, name in zip(var.shape, var.dimensions):
        if num > 1:
            shape.append(num)
            dimensions.append(name)

    ndim = len(shape)
    data = np.reshape(var[:], shape)

    # TODO complex dtype
    if ndim == 1:
        dim_name = dimensions[0]
        xx = range(len(dataset.dimensions[dim_name]))

        frame = PlotFrame(parent)
        frame.plot(xx, data)
        frame.set_xlabel(dim_name)
        frame.set_ylabel(var_name)
        #frame.scatterplot(xx, yy)
        frame.Show()

    elif ndim == 2:
        dim_namex, dim_namey = dimensions
        xx, yy = range(len(dataset.dimensions[dim_namex])), range(len(dataset.dimensions[dim_namey]))

        frame = ImageFrame(mode='intensity')
        frame.display(data, title=var_name, x=xx, y=yy, xlabel=dim_namex, ylabel=dim_namey)
        frame.Show()

    else:
        raise NotImplementedError()

#class VarInfo(wx.Frame):
#    def __init__(self, parent, var, **kwargs):
#        super(VarInfo, self).__init__(parent, -1, **kwargs)
#
#        columns = ["attribute", "value"]
#        self.list = wx.ListCtrl(self, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)
#                                                                                 
#        for col, name in enumerate(columns):
#            self.list.InsertColumn(col, name)
#
#        # Used to store the Max width in pixels for the data in the column.
#        column_widths = [awx.get_width_height(self, s)[0] for s in columns]
#
#
#        d = dict(
#            ncattrs=var.ncattrs(),
#            dimensions=var.dimensions,
#            size=var.size,
#        )
#
#        for k, v in d.items():
#            entry = map(str, [k ,v])
#            self.list.Append(entry)
#
#            w = [awx.get_width_height(self, s)[0] for s in entry]
#            column_widths = map(max, zip(w, column_widths))
#                                                                                      
#        for index, col in enumerate(columns):
#            self.list.SetColumnWidth(index, column_widths[index])


def wxapp_ncview(filepaths=()):
    """Standalone application."""
    app = wx.App()

    frame = NcViewFrame(None, filepaths=filepaths)
    app.SetTopWindow(frame)
    frame.Show()

    return app


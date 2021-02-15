import os
import collections
import wx
import netCDF4
import numpy as np
import abipy.gui.awx as awx
import wx.lib.mixins.listctrl as listmix
import wx.lib.agw.flatnotebook as fnb
import wx.lib.dialogs as wxdg
import wx.lib.newevent

from monty.string import list_strings


#def ndcontract_var(var):
#    """
#    Remove fake dimensions from a Netcdf Variable.
#
#    Returns:
#        data:
#            ndarray
#        dimensions:
#            list with the name of the dimensions in data.
#    """
#    shape, dimensions = [], []
#    for num, name in zip(var.shape, var.dimensions):
#        if num > 1:
#            shape.append(num)
#            dimensions.append(name)
#
#    # Reshape array.
#    data = np.reshape(var[:], shape)
#
#    return data, dimensions


class AttrDict(dict):
    """
    Allows to access dict keys as obj.foo in addition to the traditional way
    obj['foo']"
    """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

    def copy(self):
        newd = super(AttrDict, self).copy()
        return self.__class__(**newd)


class SimpleTextViewer(awx.Frame):
    """Very simple frame that displays text (string) in read-only mode."""
    def __init__(self, parent, text, **kwargs):
        super(SimpleTextViewer, self).__init__(parent, **kwargs)
        wx.TextCtrl(self, -1, text, style=wx.TE_MULTILINE|wx.TE_LEFT|wx.TE_READONLY)


# Command event used to signal that the value of the selected variable/dimension
# should be plotted as function of the file index.
CompareEvent, EVT_COMPARE = wx.lib.newevent.NewCommandEvent()

def getSelected(wxlist):
    """Gets the selected items from a list object."""
    selected = []

    item = -1
    while True:
        item = wxlist.GetNextItem(item, wx.LIST_NEXT_ALL, wx.LIST_STATE_SELECTED)
        if item == -1: break
        selected.append(item)

    return selected


def getColumnText(wxlist, index, col):
    """Gets the text from the specified column entry in a list."""
    item = wxlist.GetItem(index, col)
    return item.GetText()


class NcViewerFrame(wx.Frame):
    """
    This frame allows the user to inspect the dimensions
    and the variables reported in  a netcdf file.
    """
    def __init__(self, parent, filepaths=(), **kwargs):
        """
        Args:
            parent:
                parent window.
            filepaths:
                String or list of strings with the path of the nc files to open
                Empty tuple if no file should be opened during the initialization of the frame.
        """
        if "size" not in kwargs:
            kwargs["size"] = (1200, 800)

        super(NcViewerFrame, self).__init__(parent, id=-1, **kwargs)

        # This combination of options for config seems to work on my Mac.
        self.config = wx.FileConfig(appName=self.codename, localFilename=self.codename + ".ini",
                                    style=wx.CONFIG_USE_LOCAL_FILE)

        # Build menu, toolbar and status bar.
        self.makeMenu()
        self.makeToolBar()
        self.statusbar = self.CreateStatusBar()

        # Open netcdf files.
        filepaths, datasets = list_strings(filepaths), []
        filepaths = map(os.path.abspath, filepaths)

        for path in filepaths:
            datasets.append(netCDF4.Dataset(path, mode="r"))
            self.AddFileToHistory(path)

        # Create the notebook (each file will have its own tab).
        panel = wx.Panel(self, -1)
        try:
            self.notebook = fnb.FlatNotebook(panel, -1, style=fnb.FNB_NAV_BUTTONS_WHEN_NEEDED)
        except AttributeError:
            self.notebook = fnb.FlatNotebook(panel, -1)

        for path, dataset in zip(filepaths, datasets):
            tab = NcFileTab(self.notebook, dataset)
            self.notebook.AddPage(tab, os.path.basename(path))

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.EXPAND, 5)
        panel.SetSizerAndFit(sizer)

        self.Bind(wx.EVT_CLOSE, self.OnExit)

        # Intercept the command event associated to variable/dimension comparison.
        self.Bind(EVT_COMPARE, self.OnCompare)

    @property
    def codename(self):
        """Name of the application."""
        return "wxncview"

    @property
    def num_tabs(self):
        """Number of tabs in the notebook."""
        return self.notebook.GetPageCount()

    @property
    def datasets(self):
        """List of Netcdf4 `Dataset`"""
        datasets = []
        for page in range(self.notebook.GetPageCount()):
            tab = self.notebook.GetPage(page)
            datasets.append(tab.dataset)
        return datasets

    def makeMenu(self):
        """Creates the main menu."""
        self.menu_bar = wx.MenuBar()

        self.menu_bar.Append(self.createFileMenu(), "&File")
        self.menu_bar.Append(self.createPlotMenu(), "&Plot")
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

        file_history = self.file_history = wx.FileHistory(8)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(-1, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)

        file_menu.AppendSeparator()
        file_menu.Append(self.idEXIT, "E&xit", "Exit wxnciew")

        self.Bind(wx.EVT_MENU, self.OnOpen, id=self.idOPEN)
        self.Bind(wx.EVT_MENU, self.OnClose, id=self.idCLOSE)
        self.Bind(wx.EVT_MENU, self.OnExit, id=self.idEXIT)
        return file_menu

    def createPlotMenu(self):
        """Creates the plot menu."""
        # Plot Menu ID's
        plot_menu = wx.Menu()

        # Options for 2D-plots
        self.twod_menu = twod_menu = wx.Menu()
        for plot_mode in ["line", "scatter"]:
            _id =  wx.NewId()
            twod_menu.AppendRadioItem(_id, plot_mode)
        plot_menu.AppendMenu(-1, '2-D plots', twod_menu)

        # Options for image plots
        self.image_menu = image_menu = wx.Menu()
        for image_mode in ["intensity", "contour", "rgb"]:
            _id =  wx.NewId()
            image_menu.AppendRadioItem(_id, image_mode)
        plot_menu.AppendMenu(-1, 'Image mode', image_menu)

        # Options for complex data.
        self.cplx_menu = cplx_menu = wx.Menu()
        for cplx_mode in ["None", "Abs", "Real", "Imag"]:
            _id =  wx.NewId()
            cplx_menu.AppendRadioItem(_id, cplx_mode)
        plot_menu.AppendMenu(-1, 'Complex', cplx_menu)

        return plot_menu

    def createHelpMenu(self):
        """Create the help menu."""
        # Help Menu ID's
        self.idABOUT = wx.NewId()

        help_menu = wx.Menu()
        help_menu.Append(self.idABOUT, "&About", "About wxncview")
        self.Bind(wx.EVT_MENU, self.OnAbout, id=self.idABOUT)

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
        #self.Bind(wx.EVT_COMBOBOX, self.onComboSelect, id=self.idTOOL_ENTRYBOX)
        self.toolbar.Realize()

    def AddFileToHistory(self, filepath):
        """Add the absolute filepath to the file history."""
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

    def GetPlotOptions(self):
        """Returns a dictionary with the options specified in the PlotMenu"""
        d = AttrDict()

        for radio in self.twod_menu.GetMenuItems():
            if radio.IsChecked(): d["plot_mode"] = radio.GetItemLabel()

        for radio in self.image_menu.GetMenuItems():
            if radio.IsChecked(): d["image_mode"] = radio.GetItemLabel()

        for radio in self.cplx_menu.GetMenuItems():
            if radio.IsChecked(): d["cplx_mode"] = radio.GetItemLabel()

        return d

    def OnAbout(self, event):
        """Displays the about window box."""
        # TODO
        #about = aboutwin.AboutWindow(self)
        #ret = about.ShowModal()

        #awx.makeAboutBox(codename=self.codename, version=self.VERSION,
        #                 description="", developers="M. Giantomassi")

    def read_file(self, filepath):
        """Open a netcdf file, create new tab and save the file in the history."""
        try:
            notebook = self.notebook
            dataset = netCDF4.Dataset(filepath, mode="r")
            tab = NcFileTab(notebook, dataset)
            notebook.AddPage(tab, os.path.basename(filepath))
            # don't know why but this does not work!
            notebook.Refresh()
            notebook.SetSelection(notebook.GetPageCount())
            self.AddFileToHistory(filepath)
        except:
            raise

    def OnOpen(self, event):
        """Open FileDialog to allow the user to select a netcdf file."""
        dialog = wx.FileDialog(self, wildcard="*.nc")
        if dialog.ShowModal() == wx.ID_CANCEL: return

        filepath = os.path.abspath(dialog.GetPath())
        self.read_file(filepath)

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
        tab.dataset.close()

        # Remove tab.
        notebook.DeletePage(idx)
        notebook.Refresh()

    def OnExit(self, event):
        """Close the netcdf files and exit the application."""
        try:
            for index in range(self.notebook.GetPageCount()):
                tab = self.notebook.GetPage(index)
                try:
                    tab.dataset.close()
                except:
                    pass
        finally:
            self.Destroy()

    def OnFileHistory(self, event):
        """Open one of the file listed in the file_history."""
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        self.read_file(filepath)

    def OnCompare(self, event):
        """
        Callback triggered by the `compare` option in the popup menu of the panels
        Plote the netcdf variable/dimension as function of the file index.
        """
        # Get the name of the variable.
        name, obj = event.name, event.obj

        # Define the function to use the get the value of the dimension/variable.
        if isinstance(obj, netCDF4.Variable):
            def extract(dataset, name):
                return dataset.variables[name][:]

        elif isinstance(obj, netCDF4.Dimension):
            def extract(dataset, name):
                return np.array(len(dataset.dimensions[name]))

        else:
            raise ValueError("Don't know how to handle %s" % repr(obj))

        # Extract data. Make sure all the shapes are equal.
        data, shape = [], None
        for dataset in self.datasets:
            v = extract(dataset, name)
            data.append(v)
            if shape is None:
                shape = v.shape
            else:
                assert shape == v.shape

        # Plot data. Two Branches for scalars and arrays.
        opts = self.GetPlotOptions()

        is_scalar = not shape or (len(shape) == 1 and shape[0] == 1)
        #print(shape)

        from wxmplot import PlotFrame
        if is_scalar:
            frame = PlotFrame(parent=self)
            xx = range(len(self.datasets))

            if opts.plot_mode == "line":
                frame.plot(xx, data)
            else:
                frame.scatterplot(xx, data)

            frame.set_xlabel("File index")
            frame.set_ylabel(name)

            frame.Show()

        else:
            # Open new frame to allow the user to specify how to handle the array.
            ArrayComparisonFrame(self, name, data).Show()


class NcFileTab(wx.Panel):
    """Tab showing information on the netcdf file."""
    def __init__(self, parent, dataset, **kwargs):
        """
        Args:
            parent:
                parent window.
            dataset:
                Netcdf4 `Dataset`.
        """
        super(NcFileTab, self).__init__(parent, -1, **kwargs)
        self.dataset = dataset

        splitter = wx.SplitterWindow(self, -1, style=wx.SP_3DSASH)

        self.dims_panel = DimsPanel(splitter, dataset)
        self.vars_panel = VarsPanel(splitter, dataset)
        splitter.SplitVertically(self.dims_panel, self.vars_panel)
        splitter.SetSashGravity(0.1)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)


class BasePanel(wx.Panel, listmix.ColumnSorterMixin):
    """Mixin class providing helper functions."""

    @property
    def viewer_frame(self):
        try:
            return self._viewer_frame
        except AttributeError:
            self._viewer_frame = self.getParentWithType(NcViewerFrame)
            return self._viewer_frame

    def makePopupMenu(self):
        """
        Creates a popup menu containing the list of supported actions for local file objects.
        """
        # Popup IDs
        self.idPOPUP_GETINFO = wx.NewId()
        self.idPOPUP_COMPARE = wx.NewId()

        menu = wx.Menu()
        menu.Append(self.idPOPUP_GETINFO, "Get Info")
        self.Bind(wx.EVT_MENU, self.OnGetInfo, id=self.idPOPUP_GETINFO)

        if self.viewer_frame.num_tabs > 1:
            menu.Append(self.idPOPUP_COMPARE, "Compare")
            self.Bind(wx.EVT_MENU, self.OnCompare, id=self.idPOPUP_COMPARE)

        return menu

    def OnRightClick(self, event):
        """Creates a popup menu at the location of the right click."""
        menu = self.makePopupMenu()
        self.PopupMenu(menu, event.GetPosition())

    def OnColClick(self, event):
        event.Skip()

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self.list

    def GetPlotOptions(self):
        """
        Dictionary with the options for the plots taken from
        the MenuBar of the parent `NcViewerFrame`
        """
        frame = self.ncview_frame
        return frame.GetPlotOptions()

    @property
    def ncview_frame(self):
        """The parent frame `NcViewFrame`."""
        try:
            return self._ncview_frame

        except AttributeError:
            parent = self.GetParent()
            while True:
                if parent is None:
                    raise RuntimeError("Cannot find NcViewerFrame, got None parent!")

                if isinstance(parent, NcViewerFrame):
                    break
                else:
                    parent = parent.GetParent()

            self._ncview_frame = parent
            return self._ncview_frame

    def getParentWithType(self, cls):
        """
        Returns the first parent window of type cls.

        Raises:
            RuntimeError if we have reached the head of the linked list.
        """
        parent = self.GetParent()

        while True:
            if parent is None:
                raise RuntimeError("Cannot find parent with class %s, reached None parent!" % cls)

            if isinstance(parent, cls):
                return parent
            else:
                parent = parent.GetParent()


class DimsPanel(BasePanel):
    """Panel with the netcdf dimensions."""
    def __init__(self, parent, dataset, **kwargs):
        """
        Args:
            parent:
                Parent window.
            dataset:
                Netcdf4 `Dataset`.
        """
        super(DimsPanel, self).__init__(parent, -1, **kwargs)
        self.dataset = dataset

        text = wx.StaticText(self, -1, "Number of dimensions: %d" % len(dataset.dimensions))

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
        sizer.Add(text, 0, wx.EXPAND, 5)
        sizer.Add(self.list, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

        # Connect the events.
        self.list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def OnItemActivated(self, event):
        """NoP"""

    def OnCompare(self, event):
        """Trigger CompareEvent."""
        item = getSelected(self.list)[0]
        name = getColumnText(self.list, item, 0)
        dim = self.dataset.dimensions[name]

        # Create and post the event to trigger the refresh of the GUI
        event = CompareEvent(id=-1, name=name, obj=dim)
        wx.PostEvent(self.viewer_frame, event)

    def OnGetInfo(self, event):
        """Show extra information on the selected variables."""
        selected = getSelected(self.list)

        lines = []
        for item in selected:
            var_name = getColumnText(self.list, item, 0)
            var = self.dataset.dimensions[var_name]
            lines.append(str(var))
            #lines.append(ncvar_info(var))

        wxdg.ScrolledMessageDialog(self, "\n".join(lines), caption="Variable Metadata", style=wx.MAXIMIZE_BOX).Show()


class VarsPanel(BasePanel):
    """
    Panel with the netcdf variables.
    Provides popup menu to inspect/plot the selected variable.
    """
    def __init__(self, parent, dataset, **kwargs):
        """
        Args:
            parent:
                Parent window.
            dataset:
                Netcdf4 `Dataset`.
        """
        super(VarsPanel, self).__init__(parent, -1, **kwargs)

        self.dataset = dataset

        text = wx.StaticText(self, -1, "Number of variables: %d" % len(dataset.variables))

        columns = ["name", "dtype", "shape", "value", "dimensions"]
        self.list = wx.ListCtrl(self, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)

        for col, name in enumerate(columns):
            self.list.InsertColumn(col, name)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in columns]

        self.id2entry = {}
        for idx, (name, var) in enumerate(dataset.variables.items()):
            # Show the value of scalars.
            if not var.shape:
                value = str(var[:]).replace("[", "").replace("]", "").strip()
            else:
                value = "[...]"

            entry = map(str, [name, var.dtype, var.shape, value, var.dimensions])
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
        sizer.Add(text, 0, wx.EXPAND, 5)
        sizer.Add(self.list, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(sizer)

        # Connect the events whose callback will be set by the client code.
        self.list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def makePopupMenu(self):
        """
        Creates a popup menu containing the list of
        supported actions for local file objects.
        """
        # Menu of the base class
        menu = super(VarsPanel, self).makePopupMenu()

        self.idPOPUP_VAR_SHOW = wx.NewId()
        menu.Append(self.idPOPUP_VAR_SHOW, "Show")

        # Extend base class.
        #menu.AppendSeparator()

        self.idPOPUP_VAR_EXPORT = wx.NewId()
        menu.Append(self.idPOPUP_VAR_EXPORT, "Export")

        # Extra Popup IDs for arrays
        #self.idPOPUP_ARR_SLICE = wx.NewId()
        #menu.Append(self.idPOPUP_ARR_SLICE, "Slice array")
        #self.Bind(wx.EVT_MENU, self.OnArraySlice, id=self.idPOPUP_ARR_SLICE)

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.idPOPUP_VAR_SHOW, self.onVarShow),
            (self.idPOPUP_VAR_EXPORT, self.onVarExport),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu

    #def GetNameVarFromEvent(self, event):
    #    currentItem = event.m_itemIndex
    #    entry = self.id2entry[self.list.GetItemData(currentItem)]
    #    name = entry[0]
    #    return name, self.dataset.variables[name]

    def OnGetInfo(self, event):
        """Show extra information on the selected variables."""
        selected = getSelected(self.list)

        lines = []
        for item in selected:
            name = getColumnText(self.list, item, 0)
            var = self.dataset.variables[name]
            lines.append(str(var))
            #lines.append(ncvar_info(var))

        wxdg.ScrolledMessageDialog(self, "\n".join(lines), caption="Variable Metadata", style=wx.MAXIMIZE_BOX).Show()

    def onVarShow(self, event):
        item = getSelected(self.list)[0]
        name = getColumnText(self.list, item, 0)
        var = self.dataset.variables[name]
        text = str(var[:])
        SimpleTextViewer(self, text=text).Show()

    def onVarExport(self, event):
        """Export the selected variable to txt format and save it to file."""
        item = getSelected(self.list)[0]
        name = getColumnText(self.list, item, 0)
        var = self.dataset.variables[name]

        dialog = wx.FileDialog(self, message="Save text file", defaultDir="", defaultFile="",
                               wildcard="", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if dialog.ShowModal() == wx.ID_CANCEL: return

        # save the current contents in the file
        path = dialog.GetPath()

        # Save self to a text file. See :func:`np.savetext` for the description of the variables
        np.savetxt(path, var[:], fmt='%.18e')

    def OnCompare(self, event):
        """Trigger CompareEvent."""
        item = getSelected(self.list)[0]
        name = getColumnText(self.list, item, 0)
        var = self.dataset.variables[name]

        # Create and post the event to trigger the comparison.
        event = CompareEvent(id=-1, name=name, obj=var)
        wx.PostEvent(self.viewer_frame, event)

    def OnItemActivated(self, event):
        """Use `wxmplot` to plot the selected variables."""
        selected = getSelected(self.list)
        for item in selected:
            name = getColumnText(self.list, item, 0)
            var = self.dataset.variables[name]
            self.plot_variable(name, var, self.dataset)

    def plot_variable(self, var_name, var, dataset):
        """
        Use `wxmplot` to plot the selected variables.

        Args:
            var_name:
                Name of the variable
            var:
                Netcdf4 `Variable`.
            dataset:
                Netcdf4 `Dataset`.
        """
        # Remove fake dimensions.
        shape, dimensions = [], []
        for num, name in zip(var.shape, var.dimensions):
            if num > 1:
                shape.append(num)
                dimensions.append(name)

        # Get data to plot.
        data = np.reshape(var[:], shape)
        opts = self.GetPlotOptions()

        cplx_mode = opts.cplx_mode
        if cplx_mode != "None":
            if shape[-1] != 2:
                err_msg = "cplx_mode: %s. Expecting 2 as last dimensions but got %d" % (
                    cplx_mode, shape[-1])
                raise ValueError(err_msg)
            # Convert to complex then change shape and dimensions
            data = data[..., 0] + 1j*data[..., 1]
            shape = shape[:-1]
            dimensions = dimensions[:-1]
            if cplx_mode == "Abs":
                data = np.abs(data)
            elif cplx_mode == "Real":
                data = data.real
            elif cplx_mode == "Imag":
                data = data.imag
            else:
                raise ValueError("Wrong value for cplx_mode %s" % cplx_mode)

        # Plotting a scalar?
        if not shape: return
        ndim = len(shape)

        if ndim == 1:
            # Vector
            dim_name = dimensions[0]
            xx = range(len(dataset.dimensions[dim_name]))

            from wxmplot import PlotFrame
            frame = PlotFrame(parent=self)
            if opts.plot_mode == "line":
                frame.plot(xx, data)
            else:
                frame.scatterplot(xx, data)

            frame.set_xlabel(dim_name)
            frame.set_ylabel(var_name)

            frame.Show()

        elif ndim == 2:
            # Two dimensional array.
            dim_namex, dim_namey = dimensions
            xx, yy = range(len(dataset.dimensions[dim_namex])), range(len(dataset.dimensions[dim_namey]))

            mode = opts.image_mode

            if False:
                # 3d plot
                import matplotlib.pyplot as plt
                from mpl_toolkits.mplot3d import Axes3D
                fig = plt.figure()
                ax = Axes3D(fig)
                X, Y = np.meshgrid(xx, yy, sparse=False, indexing='ij')
                print(X.shape, Y.shape, data.shape)
                ax.plot_surface(X, Y, data) # rstride=8, cstride=8, alpha=0.3)
                plt.show()

            from wxmplot import ImageFrame
            frame = ImageFrame(parent=self)
            frame.display(data, title=var_name, style=mode, x=xx, y=yy, xlabel=dim_namex, ylabel=dim_namey)
            frame.Show()

        else:
            raise NotImplementedError()

    #def OnArraySlice(self, event):
    #    print("Onslice")
    #    item = getSelected(self.list)[0]
    #    name = getColumnText(self.list, item, 0)
    #    var = self.dataset.variables[name]
    #    ArraySlicerFrame(self, name, var).Show()


class ArrayComparisonFrame(wx.Frame):
    """This frame allows the user to specify how to handle array comparison."""

    def __init__(self, parent, name, data, **kwargs):
        super(ArrayComparisonFrame, self).__init__(parent, -1, **kwargs)
        self.name, self.data = name, data
        self.panel = ArrayComparisonPanel(self)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        plot_button = wx.Button(self, -1, "Plot")
        hsizer.Add(plot_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        plot_button.Bind(wx.EVT_BUTTON, self.onPlotData)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)
        main_sizer.Add(self.panel, 1, wx.EXPAND, 5)

        self.SetSizerAndFit(main_sizer)

    def onPlotData(self, event):
        operators = self.panel.GetSelectedOperators()

        from wxmplot import PlotFrame
        frame = PlotFrame(parent=self)

        xx = range(len(self.data))

        for oname, op in operators.items():
            #print(oname)
            values = [op(arr) for arr in self.data]

            frame.oplot(xx, values, label=oname)
            #if opts.plot_mode == "line":
            #    frame.oplot(xx, values)
            #else:
            #    frame.scatterplot(xx, values)

        frame.set_xlabel("File index")
        frame.set_ylabel(self.name)

        frame.Show()


class ArrayComparisonPanel(wx.Panel):
    """A panel with checkboxes used to select the scalar quantities to compare."""
    def __init__(self, parent, **kwargs):
        """
        Args:
            parent:
                Parent window.
        """
        super(ArrayComparisonPanel, self).__init__(parent, -1, **kwargs)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        static_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Array contractions"), wx.VERTICAL)

        operators = [
            ("max", np.max),
            ("min", np.min),
            ("mean", np.mean),
            ("std", np.std),
            #("median", np.median),
        ]

        self.operators =  collections.OrderedDict(operators)
        self.check_boxes = collections.OrderedDict()

        for oname, _ in self.operators.items():
            cbox = wx.CheckBox(self, -1, oname, wx.DefaultPosition, wx.DefaultSize, 0)
            cbox.SetValue(True)
            static_sizer.Add(cbox, 0, wx.ALL | wx.EXPAND, 5)
            self.check_boxes[oname] = cbox

        main_sizer.Add(static_sizer, 1, wx.EXPAND, 5)

        # Add buttons to (select|deselect) all checkboxes.
        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        all_button = wx.Button(self, -1, "Select all")
        all_button.Bind(wx.EVT_BUTTON, self.OnSelectAll)
        hsizer.Add(all_button, 0, wx.ALL, 5)

        deselect_button = wx.Button(self, -1, "Deselect all")
        deselect_button.Bind(wx.EVT_BUTTON, self.OnDeselectAll)
        hsizer.Add(deselect_button, 0, wx.ALL, 5)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

    def OnSelectAll(self, event):
        """Select all the check boxes."""
        for cbox in self.check_boxes.values():
            cbox.SetValue(True)

    def OnDeselectAll(self, event):
        """Deselect all the check boxes."""
        for cbox in self.check_boxes.values():
            cbox.SetValue(False)

    def GetSelectedOperators(self):
        """
        Return the list of operators selected by the user.
        """
        ops = collections.OrderedDict()

        for oname, cbox in self.check_boxes.items():
            if cbox.GetValue():
                ops[oname] = self.operators[oname]

        return ops


#class ArraySlicerFrame(wx.Frame):
#
#    def __init__(self, parent, name, var, **kwargs):
#        super(ArraySlicerFrame, self).__init__(parent, -1, **kwargs)
#        self.name, self.var = name, var
#
#        self.panel_1d = Slice1dPanel(self, name, var)
#        #self.panel_2d = Slice1dPanel(self, name, arr)
#        #self.panel_3d = Slice1dPanel(self, name, arr)
#
#        main_sizer = wx.BoxSizer(wx.VERTICAL)
#        main_sizer.Add(self.panel_1d, 1, wx.EXPAND, 5)
#        self.SetSizer(main_sizer)
#
#
#class Slice1dPanel(wx.Panel):
#    """A panel plotting for 1D slices."""
#    def __init__(self, parent, name, var, **kwargs):
#        """
#        Args:
#            parent:
#                Parent window.
#        """
#        super(Slice1dPanel, self).__init__(parent, -1, **kwargs)
#        self.name, self.var = name, var
#
#        label = wx.StaticText(self, -1, "Axis:")
#        label.Wrap(-1)
#        self.axis_choice = wx.ComboBox(self, -1, choices=var.dimensions)
#
#        hsz1 = wx.BoxSizer(wx.HORIZONTAL)
#        hsz1.Add(label, 0, wx.ALL, 5)
#        hsz1.Add(self.axis_choice, 1, wx.ALL, 5)
#
#        plot_button = wx.Button(self, -1, "Plot")
#        plot_button.Bind(wx.EVT_BUTTON, self.onPlotButton)
#        #hsizer = wx.BoxSizer(wx.HORIZONTAL)
#        #hsizer.Add(plot_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
#
#        main_sizer = wx.BoxSizer(wx.VERTICAL)
#        main_sizer.Add(hsz1, 0.4, wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, 5)
#        main_sizer.Add(plot_button, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)
#        self.SetSizerAndFit(main_sizer)
#
#    def onPlotButton(self, event):
#        var = self.var
#        axis = var.dimensions.index(self.axis_choice.GetValue())
#        indices = var.shape
#        #frame = PlotFrame(parent=self)
#        #xx = range(len(self.var))
#        #frame.oplot(xx, values, label=oname)
#        #frame.set_xlabel("File index")
#        #frame.set_ylabel(self.name)
#        #frame.Show()


def wxapp_ncview(filepaths=()):
    """Standalone application."""
    app = wx.App()
    frame = NcViewerFrame(None, filepaths=filepaths)
    app.SetTopWindow(frame)
    frame.Show()

    return app

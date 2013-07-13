#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import wx

from wx.lib.agw.floatspin import FloatSpin
from wxmplot import PlotApp, PlotFrame

from abipy.abifiles import abiopen
from abipy.electrons import ElectronBands
from abipy.waves import WFK_File
from abipy.gui.tools import build_aboutbox, build_errormessage, straceback, path_img

class WxDosPlotter(wx.Frame):
    VERSION = "0.1"

    # Help menu items.
    ID_ABOUT = wx.NewId() 

    def __init__(self, parent=None, **kwargs):
        super(WxDosPlotter, self).__init__(parent, -1, self.codename, size=(800,600)) 

        # Create statusbar
        self.statusbar = self.CreateStatusBar() 

        # Setup menu bar.
        menuBar = wx.MenuBar()

        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_OPEN,  "&Open", help="Open an existing WFK file")
        fileMenu.Append(wx.ID_CLOSE, "&Close", help="Close the WFK file")
        fileMenu.Append(wx.ID_EXIT,  "&Quit", help="Exit the application")
        menuBar.Append(fileMenu, "File")

        filehistory = self.filehistory = wx.FileHistory(8)
        self.config = wx.Config(self.codename, style=wx.CONFIG_USE_LOCAL_FILE)
        filehistory.Load(self.config)
        recent = wx.Menu()
        filehistory.UseMenu(recent)
        filehistory.AddFilesToMenu()
        fileMenu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.on_file_history, id=wx.ID_FILE1, id2=wx.ID_FILE9)

        self.helpMenu = wx.Menu()
        self.helpMenu.Append(self.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.helpMenu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        tsize = (15,15)
        self.toolbar = toolbar = self.CreateToolBar()

        artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        #toolbar.AddSimpleTool(self.ID_VISTRUCT, wx.Bitmap(path_img("crystal.png")), "Visualize the crystal structure")
        #toolbar.AddSimpleTool(self.ID_VISWAVE, wx.Bitmap(path_img("wave.png")), "Visualize the selected wavefunction")

        self.toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menuHandlers = [
            (wx.ID_OPEN,        self.onOpen),
            (wx.ID_CLOSE,       self.onClose),
            (wx.ID_EXIT,        self.onExit),
            (self.ID_ABOUT,     self.onAboutBox),
        ]

        for combo in menuHandlers:
            id, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=id)

        # Widgets.
        self.panel = panel = FileListPanel(self)
        self.panel.Show()

    @property
    def codename(self):
        return self.__class__.__name__

    def onOpen(self, event):
        dlg = wx.FileDialog(self, message="Choose a Netcdf file", defaultDir=os.getcwd(), 
            defaultFile="", wildcard="Netcdf files (*.nc)|*nc",
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. 
        # If it is the OK response, process the data.
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.statusbar.PushStatusText("Reading %s" % path)
            self.filehistory.AddFileToHistory(path)
            self.filehistory.Save(self.config)
            self.config.Flush()
            self.read_ebands(path)

        dlg.Destroy()

    def on_file_history(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        path = self.filehistory.GetHistoryFile(fileNum)
        # move up the list
        self.filehistory.AddFileToHistory(path)  
        self.read_ebands(path)

    def read_ebands(self, path):
        try:
            self.ebands = ElectronBands.from_ncfile(path)
            EdosFrame(ebands=self.ebands, parent=self).Show()
        except Exception:
            build_errormessage(self, straceback())

    def onClose(self, event):
        self.ebands = None

    def onExit(self, event):
        self.Destroy()

    def onAboutBox(self, event):
        build_aboutbox(codename=self.codename, version=self.VERSION, 
                       description="", developers="M. Giantomassi")


class EdosFrame(wx.Frame):

    def __init__(self, ebands, parent=None, **kwargs):
        super(EdosFrame, self).__init__(parent, id=-1, title="Electron DOS", **kwargs)
        self.ebands = ebands

        self.statusbar = self.CreateStatusBar() 

        # Widgets.
        self.panel = panel = wx.Panel(self, -1)
        sizer = wx.FlexGridSizer(rows=3, cols=2, vgap=5, hgap=5)

        label = wx.StaticText(panel, -1, "Broadening [eV]:")

        self.width = 0.2
        self.width_ctrl = FloatSpin(panel, id=-1, value=self.width, min_val=0.0, increment=0.1, digits=3)
        self.Bind(wx.EVT_SPINCTRL, self.onWidth, self.width_ctrl)

        sizer.AddMany([label, self.width_ctrl])

        label = wx.StaticText(panel, -1, "Mesh step [eV]:")
        self.step = 0.1
        self.step_ctrl = FloatSpin(panel, id=-1, value=self.step, min_val=0.0, increment=0.1, digits=3)
        self.Bind(wx.EVT_SPINCTRL, self.onStep, self.step_ctrl)

        sizer.AddMany([label, self.step_ctrl])

        dos_button = wx.Button(panel, -1, "Compute DOS", (20,100))
        self.Bind(wx.EVT_BUTTON, self.onClick, dos_button)

        self.oplot_checkbox = wx.CheckBox(panel, id=-1, label="Same plot")
        self.oplot_checkbox.SetValue(True)

        sizer.AddMany([dos_button, self.oplot_checkbox])

        panel.SetSizer(sizer)
        sizer.Layout()

    def onWidth(self, event):
        self.width = float(self.width_ctrl.GetValue())

    def onStep(self, event):
        self.step = float(self.step_ctrl.GetValue())

    def onClick(self, event):
        try:
            edos = self.ebands.get_dos(step=self.step, width=self.width)
        except:
            build_errormessage(self, straceback())
            return

        tot_dos, tot_idos = edos.dos_idos()
        label = "$\sigma = %s, step = %s$" % (self.width, self.step)

        if self.has_plotframe and self.oplot_checkbox.GetValue():
            self.plotframe.oplot(tot_dos.mesh, tot_dos.values, label=label, draw_legend=True) 
        else:
            self.plotframe = plotframe = PlotFrame(parent=self)
            plotframe.plot(tot_dos.mesh, tot_dos.values, label=label, draw_legend=True)
            plotframe.Show()

    @property
    def has_plotframe(self):
        return hasattr(self, "plotframe")  

    @property
    def destroy_plotframe(self):
        if self.has_plotframe:
            self.plotframe.Destroy()
            del self.plotframe

class FileListPanel(wx.Panel):

    def __init__(self, parent, paths=(), **kwargs):
        super(FileListPanel, self).__init__(parent, -1, **kwargs)

        self.file_list = file_list = wx.ListCtrl(self, id=-1, size=(-1,100), 
                                     style=wx.LC_REPORT | wx.BORDER_SUNKEN
                                     )
        file_list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.onRightClick)

        self.ncfiles_by_id = {}
        file_list.InsertColumn(0, "filename")
        file_list.InsertColumn(1, "filetype")
        for index, path in enumerate(paths):
            self.append_path_to_filelist(path)

        self.filepicker = wx.FilePickerCtrl(self, id=-1, path=os.getcwd(),
            wildcard="Netcdf files (*.nc)|*nc", style = wx.FLP_OPEN | wx.CHANGE_DIR)

        self.Bind(wx.EVT_FILEPICKER_CHANGED, self.onFilePicker)

        # Pack
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(file_list, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(self.filepicker, 0, wx.ALL | wx.CENTER, 5)
        self.SetSizer(sizer)

    def append_path_to_filelist(self, path):
        next = len(self.ncfiles_by_id)
        try:
            ncfile = abiopen(path)
            file_type = ncfile.__class__.__name__
            #file_type = ncfile.filetype
            ncid = id(ncfile)
            self.ncfiles_by_id[ncid] = ncfile 
            entry = [os.path.basename(path), file_type]
            self.file_list.Append(entry)
            self.file_list.SetItemData(next, ncid)
        except:
            build_errormessage(self, straceback())

    def onFilePicker(self, event):
        self.append_path_to_filelist(self.filepicker.GetPath())

    def onRightClick(self, event):
        currentItem = event.m_itemIndex

        if currentItem != -1:
            ncfile = self.ncfiles_by_id[self.file_list.GetItemData(currentItem)]
            print("Select ncfile: ",ncfile)

            menu = popup_menu_for_ncfile(ncfile)

            # Open the popup menum then destroy it to avoid mem leak.
            self.PopupMenu(menu, event.GetPoint())
            menu.Destroy() 

def popup_menu_for_ncfile(ncfile):
    menu = BasePopupMenu()
    menu.add_target(ncfile)
    return menu

class BasePopupMenu(wx.Menu):
    menu_titles = [ 
        "Properties",
    ]

    menu_title_by_id = {}
    for title in menu_titles:
        menu_title_by_id[wx.NewId()] = title

    def __init__(self, *args, **kwargs):
        super(BasePopupMenu, self).__init__()

        for (id, title) in self.menu_title_by_id.items():
            self.Append(id, title)
            # registers menu handlers with EVT_MENU, on the menu.
            wx.EVT_MENU(self, id, self.MenuSelectionCb)

    def add_target(self, target):
        self._target = target

    @property
    def target(self):
        try:
            return self._target
        except AttributeError:
            return None

    def MenuSelectionCb(self, event):
        # do something
        operation = self.menu_title_by_id[event.GetId()]
        print("Perform operation %s on target %s" % (operation, self.target))


if __name__ == "__main__":
    app = wx.App()
    win = WxDosPlotter()
    win.Show(True)
    app.MainLoop()

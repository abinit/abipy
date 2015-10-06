#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals

import os
import wx
import abipy.gui.awx as awx
import wx.lib.mixins.listctrl as listmix
import wx.lib.dialogs as wxdg

from monty.string import is_string 
from pymatgen.io.abinit.pseudos import PseudoParser


class MyListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin):
    """ Mixin class to resize the last column appropriately."""
    def __init__(self, parent, **kwargs):
        wx.ListCtrl.__init__(self, parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)
        listmix.ListCtrlAutoWidthMixin.__init__(self)


class PseudosViewerFrame(awx.Frame):
    def __init__(self, parent, pseudos, **kwargs):
        super(PseudosViewerFrame, self).__init__(parent, -1, **kwargs)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        panel = PseudoListPanel(self, pseudos)
        main_sizer.Add(panel, 1, wx.EXPAND, 5)

        #hsizer = wx.BoxSizer(wx.HORIZONTAL)

        #filter_label = wx.StaticText(self, -1, "Filter:", wx.DefaultPosition, wx.DefaultSize, 0)
        #filter_label.Wrap(-1)
        #hsizer.Add(filter_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        ## Combobox to enter shell patterns.
        #wildcard_choices = ["*", "*.nc", "*.abo", "*.log"]
        #self.filter_combobox = wx.ComboBox(self, id=-1, value="*", style=wx.TE_PROCESS_ENTER, choices=wildcard_choices)
        #self.filter_combobox.SetToolTipString("Shell patterns separated by |")

        #self.filter_combobox.Bind(wx.EVT_COMBOBOX, self.OnFilterComboBox)
        #self.filter_combobox.Bind(wx.EVT_TEXT_ENTER, self.OnFilterComboBox)

        #hsizer.Add(self.filter_combobox, 0, wx.ALL, 5)

        #main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        #panel.SetSizerAndFit(main_sizer)
        self.SetSizer(main_sizer)
        self.Layout()


class PseudoListPanel(awx.Panel, listmix.ColumnSorterMixin):
    """
    This panel shows a list of files (strings), supports column sorting
    and provides specific popup menus for the different type of files 
    (file type detection is based on file extensions).
    """
    def __init__(self, parent, pseudos, **kwargs):
        """
        Args:
            parent:
                parent window
            pseudos:
                List of pseudopotentials.
        """
        super(PseudoListPanel, self).__init__(parent, -1, **kwargs)
        self.pseudos = pseudos

        self.BuildUi()

    def BuildUi(self):
        self.pseudo_list = pseudo_list = MyListCtrl(self)

        pseudo_list.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        pseudo_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

        self.id2pseudo = {}
        columns = self.columns = ["basename", "symbol", "Z_val", "l_max", "rcore", "type"]

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in columns]

        for (index, colname) in enumerate(columns):
            pseudo_list.InsertColumn(index, colname)

        for pseudo in self.pseudos:
            entry = self.AppendPseudo(pseudo)
            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

        for (index, colname) in enumerate(columns):
            pseudo_list.SetColumnWidth(index, column_widths[index])

        # Now that the list exists we can init the other base class, see wx/lib/mixins/listctrl.py
        self.itemDataMap = {pid: [getattr(pseudo, attr) for attr in self.columns] for pid, pseudo in self.id2pseudo.items()}
        listmix.ColumnSorterMixin.__init__(self, len(columns))
        self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.pseudo_list)

        # Pack
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(pseudo_list, 1, wx.ALL | wx.EXPAND, 5)

        self.SetSizerAndFit(sizer)

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self.pseudo_list

    def AppendPseudo(self, pseudo):
        """Add a pseudo to the panel."""
        next = len(self.id2pseudo)

        # We use next as entry id because we want to be able 
        # to sort the columns with the mixin ColumnSorterMixin.
        entry_id = next
        entry = map(str, [getattr(pseudo, attr) for attr in self.columns])
        self.pseudo_list.Append(entry)
        self.pseudo_list.SetItemData(next, entry_id)
        self.id2pseudo[entry_id] = pseudo

        return entry

    def OnColClick(self, event):
        event.Skip()

    def OnItemActivated(self, event):
        currentItem = event.m_itemIndex
        pseudo = self.id2pseudo[self.pseudo_list.GetItemData(currentItem)]
        print(pseudo)

        wxdg.ScrolledMessageDialog(self, str(pseudo), caption=pseudo.filepath, style=wx.MAXIMIZE_BOX).Show()

    #def OnOpenPseudoFile(self, event):
    #    """Open the pseudopotential file in another frame."""
    #    currentItem = event.m_itemIndex
    #    pseudo = self.id2pseudo[self.pseudo_list.GetItemData(currentItem)]
    #    print(pseudo)
    #    filepath = pseudo.filepath
    #    wxdg.ScrolledMessageDialog(self, str(pseudo), caption=pseudo.filepath, style=wx.MAXIMIZE_BOX).Show()

    #def OnPlotPseudo(self, event):
    #    """Open the pseudopotential file in another frame."""
    #    currentItem = event.m_itemIndex
    #    pseudo = self.id2pseudo[self.pseudo_list.GetItemData(currentItem)]
    #    print(pseudo)

    def OnRightClick(self, event):
        currentItem = event.m_itemIndex
        if currentItem == -1:
            return

        #pseudo = self.id2pseudo[self.pseudo_list.GetItemData(currentItem)]
        #print(pseudo)
        ## Open the popup menu then destroy it to avoid mem leak.
        #menu = popupmenu_for_pseudo(self, pseudo)

        #if menu is not None:
        #    self.PopupMenu(menu, event.GetPoint())
        #    menu.Destroy()


#class PopupMenu(wx.Menu):
#    """
#    Base class for popup menus. `A PopupMenu` has a list of callback functions
#    indexed by the menu title. The signature of the callback function is func(parent, filename) where 
#    filename is the name of the file selected in the Widget and parent is the wx
#    Window that will become the parent of the new frame created by the callback.
#    """
#    MENU_TITLES = OrderedDict([
#    ])
#
#    HANDLED_FILES = []
#
#    def __init__(self):
#        super(PopupMenu, self).__init__()
#        self._make_menu()
#
#    @staticmethod
#    def from_filename(filename):
#        """
#        Static factory function that instanciates the appropriate subclass of `NcFilePopupMenu`
#        Returns None if the extesion of filename is not supported.
#        """
#        # Find the AbinitNcFile subclass associated to files.
#        try:
#            file_class = abifile_subclass_from_filename(filename)
#        except KeyError:
#            return None
#
#        # Check whether a subclass handles this file..
#        # Fallback to a simple PopupMenu if no match.
#        def allsubclasses(cls):
#            """Returns the set of subclasses of cls."""
#            children = [cls]
#            for sc in cls.__subclasses__():
#                if sc.__subclasses__():
#                    for k in sc.__subclasses__():
#                        children.extend(allsubclasses(k))
#                else:
#                    children.append(sc)
#            return set(children)
#            
#        for cls in allsubclasses(PopupMenu):
#            if cls.handle_file_class(file_class):
#                return cls()
#        else:
#            return PopupMenu()
#
#    @classmethod
#    def handle_file_class(cls, file_class):
#        """True if the popupmenu is associated to file_class."""
#        return file_class in cls.HANDLED_FILES
#
#    def _make_menu(self):
#        """Build the menu taking into account the options of the superclasses."""
#        base_classes = list(self.__class__.__bases__) + [self.__class__]
#        base_classes.reverse()
#
#        assert not hasattr(self, "menu_title_by_id")
#        assert not hasattr(self, "menu_titles")
#        self.menu_title_by_id, self.menu_titles = OrderedDict(), OrderedDict()
#
#        for cls in base_classes:
#            try:
#                menus = cls.MENU_TITLES
#            except AttributeError as exc:
#                awx.WARNING("exc ",exc," for cls", cls)
#                continue
#
#            self.menu_titles.update(menus)
#
#            for title in menus:
#                self.menu_title_by_id[wx.NewId()] = title
#
#            # Add sentinel for Menu separator.
#            self.menu_title_by_id["separator_" + str(len(self.menu_titles))] = None
#                                                                  
#        for (id, title) in self.menu_title_by_id.items():
#            if title is None:
#                sepid = int(id.split("_")[-1])
#                if sepid != len(self.menu_titles):
#                    self.AppendSeparator()
#            else:
#                # Register menu handlers with EVT_MENU, on the menu.
#                self.Append(id, title)
#                wx.EVT_MENU(self, id, self.OnMenuSelection)
#
#    def set_parent(self, parent):
#        """Set the parent window."""
#        self._parent = parent
#
#    @property
#    def parent(self):
#        """Returns the parent window."""
#        try:
#            return self._parent
#        except AttributeError:
#            awx.WARNING("Popup menu doesn't have parent")
#            return None
#
#    def set_target(self, target):
#        """Set the target of the callback."""
#        self._target = target
#
#    @property
#    def target(self):
#        """The target of the callback."""
#        try:
#            return self._target
#        except AttributeError:
#            return None
#
#    def _get_callback(self, title):
#        return self.menu_titles[title]
#
#    def OnMenuSelection(self, event):
#        title = self.menu_title_by_id[event.GetId()]
#        callback = self._get_callback(title)
#        #print("Calling callback %s on target %s" % (callback, self.target))
#        try:
#            callback(parent=self.parent, filepath=self.target)
#        except:
#            awx.showErrorMessage(parent=self.parent)
#
#
#class AbinitTextFilePopupMenu(PopupMenu):
#    """
#    """
#    MENU_TITLES = OrderedDict([
#        ("events",     showAbinitEventsFrame),
#        ("properties", showFileStat),
#        ("timer",      showAbinitTimerFrame),
#    ])
#
#    HANDLED_FILES = [AbinitLogFile, AbinitOutputFile]
#
#
#class NcFilePopupMenu(PopupMenu):
#    """
#    Base class for popup menus. `A PopupMenu` has a list of callback functions
#    indexed by the menu title and a list of `AbinitNcFile` associated to it.
#    The signature of the callback function is func(filename, parent) where 
#    filename is the name of the file selected in the Widget and parent is the wx
#    Window that will become the parent of the new frame created by the callback.
#
#    How to subclass PopupMenu:
#
#        1. Define a new class that inherits from NcFilePopupMenu. 
#
#        2. Define the callbacks in the class variable MENU_TITLES.
#           Use OrderedDict to have a fixed ordering of the labels.
#
#        3. Define the class variable HANDLED_FILES with the list of 
#           `AbinitNcFile` subclasses associated to the popupmenu.
#
#        4. Done (most of the work is indeed done in the base class and in 
#           the factory function popupmenu_for_filename.
#    """
#    MENU_TITLES = OrderedDict([
#        ("structure",  showStructure),
#        ("ncdump",     showNcdumpMessage),
#        ("properties", showFileStat),
#    ])
#
#    HANDLED_FILES = []
#
#
#class EbandsPopupMenu(NcFilePopupMenu):
#    """Popup menu for Netcdf files that contain the electron band structure."""
#    MENU_TITLES = OrderedDict([
#        ("ePlot", ewx.showElectronBandsPlot),
#        ("eDos",  ewx.showElectronDosFrame),
#        ("eJdos", ewx.showElectronJdosFrame),
#    ])
#
#    HANDLED_FILES = [WFK_File, GSR_File] 


def wxapp_pseudos(dirpath):
    """Standalone application."""
    parser = PseudoParser()
    pseudos = parser.scan_directory(dirpath)

    app = wx.App()
    frame = PseudosViewerFrame(None, pseudos)
    app.SetTopWindow(frame)
    frame.Show()
    return app

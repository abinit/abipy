from __future__ import print_function, division

import wx

import abipy.gui.awx as awx

from collections import OrderedDict

class KpointSymmetriesFrame(awx.Frame):

    def __init__(self, parent, structure, kpoints, **kwargs):
        """
        Args:
            parent:
                Parent window.
            structure:
                `Structure` object.
            kpoints:
                `KpointList` object. 
        """
        super(KpointSymmetriesFrame, self).__init__(parent, **kwargs)

        self.kpoint_panel = KpointSymmetriesPanel(self, structure, kpoints)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.kpoint_panel, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)
        main_sizer.Layout()
        self.SetSizer(main_sizer)

class KpointSymmetriesPanel(awx.Panel):

    def __init__(self, parent, structure, kpoints, **kwargs):
        """
        Args:
            parent:
                Parent window.
            structure:
                `Structure` object.
            kpoints:
                `KpointList` object. 
        """
        super(KpointSymmetriesPanel, self).__init__(parent, -1, **kwargs)

        self.kpoints_listctrl = awx.KpointsListCtrl(parent, kpoints)
        self.structure = structure

        # Connect the events whose callback will be set by the client code.
        self.kpoints_listctrl.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        #self.kpoints_listctrl.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnItemSelected)
        #self.kpoints_listctrl.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def OnRightClick(self, event):
        popmenu = KpointPopupMenu()
        popmenu.set_parent(self)
        popmenu.set_cb_kwargs(structure=self.structure, kpoint=self.kpoints_listctrl.GetKpoint())

        self.PopupMenu(popmenu, event.GetPoint())
        popmenu.Destroy()


##################
### Callbacks ###
##################

def showStarFrame(parent, *args, **kwargs):
    structure, kpoint = kwargs.pop("structure"), kwargs.pop("kpoint")
    #pos = parent.GetPosition()

    star = kpoint.compute_star(structure.fm_symmops)
    KpointSymmetriesFrame(parent, structure, star, title="Star of %s" % str(kpoint)).Show()


class KpointPopupMenu(wx.Menu):
    """
    Base class for popup menus. `A KpointPopupMenu` has a list of callback functions
    indexed by the menu title. The signature of the callback function is func(parent, *args, **kwargs) 
    parent is the wx Window that will become the parent of the new frame created by the callback.
    """
    MENU_TITLES = OrderedDict([
        ("Show Star", showStarFrame),
        #("Little Group", showLittleGroupFrame),
        #("Irreps", showIrrepsFrame),
    ])

    def __init__(self):
        super(KpointPopupMenu, self).__init__()
        self._make_menu()

    def _make_menu(self):
        """Build the menu taking into account the options of the superclasses."""

        assert not hasattr(self, "menu_title_by_id")
        assert not hasattr(self, "menu_titles")
        self.menu_title_by_id, self.menu_titles = OrderedDict(), OrderedDict()

        for title in self.MENU_TITLES:
            self.menu_title_by_id[wx.NewId()] = title

        for (id, title) in self.menu_title_by_id.items():
            # Register menu handlers with EVT_MENU, on the menu.
            self.Append(id, title)
            wx.EVT_MENU(self, id, self.OnMenuSelection)

    def set_parent(self, parent):
        """Set the parent window."""
        self._parent = parent

    @property
    def parent(self):
        """Returns the parent window."""
        try:
            return self._parent
        except AttributeError:
            awx.WARNING("Popup menu doesn't have parent")
            return None

    def set_cb_args(self, *cb_args):
        """Set the positional arguments of the callback."""
        self._cb_args = cb_args

    @property
    def cb_args(self):
        """Callback positional arguments."""
        try: 
            return self._cb_args
        except AttributeError:
            return []

    def set_cb_kwargs(self, **cb_kwargs):
        """Set the kwyword arguments of the callback."""
        self._cb_kwargs = cb_kwargs
                                                  
    @property
    def cb_kwargs(self):
        """Callback keyword arguments."""
        try: 
            return self._cb_kwargs
        except AttributeError:
            return {}

    def OnMenuSelection(self, event):
        """Call the callback selected in the popupmenu."""
        title = self.menu_title_by_id[event.GetId()]
        callback = self.MENU_TITLES[title]

        try:
            callback(parent=self.parent, *self.cb_args, **self.cb_kwargs)
        except:
            awx.showErrorMessage(parent=self)


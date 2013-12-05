from __future__ import print_function, division

import wx

import abipy.gui.awx as awx
import wx.lib.dialogs as wxdg

from collections import OrderedDict


class KpointsPanel(awx.Panel):

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
        super(KpointsPanel, self).__init__(parent, **kwargs)

        self.kpoints_listctrl = awx.KpointsListCtrl(self, kpoints)
        self.structure = structure

        # Connect the events whose callback will be set by the client code.
        self.kpoints_listctrl.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.kpoints_listctrl, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)
        self.SetSizerAndFit(main_sizer)

    def OnRightClick(self, event):
        popmenu = KpointPopupMenu(self, structure=self.structure, kpoint=self.kpoints_listctrl.GetKpoint())
        self.PopupMenu(popmenu, event.GetPoint())
        popmenu.Destroy()


##################
### Callbacks ###
##################


def showLittleGroup(parent, structure, kpoint):
    ltk = structure.spacegroup.find_little_group(kpoint)
    #wxdg.ScrolledMessageDialog(parent, str(ltk), caption=repr(ltk), style=wx.MAXIMIZE_BOX).Show()
    table, header = ltk.bilbao_character_table()
    awx.SimpleGridFrame(parent, table, title=header, labels_from_table=True).Show()


def showStar(parent, structure, kpoint):
    star = kpoint.compute_star(structure.fm_symmops)
    wxdg.ScrolledMessageDialog(parent, str(star), caption="Star of the kpoint", style=wx.MAXIMIZE_BOX).Show()


class KpointPopupMenu(wx.Menu):
    """
    Base class for popup menus. `A KpointPopupMenu` has a list of callback functions
    indexed by the menu title. The signature of the callback function is func(parent, *args, **kwargs) 
    parent is the wx Window that will become the parent of the new frame created by the callback.
    """
    MENU_TITLES = OrderedDict([
        ("Little Group", showLittleGroup),
        ("Star", showStar),
    ])

    def __init__(self, parent, structure, kpoint):
        """
        Args:
            parent:
                Parent window.
            structure:
                `Structure` object
            kpoint:
                `Kpoint` object of vector with the reduced coordinates.
        """
        super(KpointPopupMenu, self).__init__()
        self._parent, self.structure, self.kpoint = parent, structure, kpoint

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

    @property
    def parent(self):
        """Returns the parent window."""
        return self._parent

    def OnMenuSelection(self, event):
        """Call the callback selected in the popupmenu."""
        title = self.menu_title_by_id[event.GetId()]
        callback = self.MENU_TITLES[title]

        try:
            callback(parent=self.parent, structure=self.structure, kpoint=self.kpoint)
        except:
            awx.showErrorMessage(parent=self.parent)


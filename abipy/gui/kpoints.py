from __future__ import print_function, division

import wx

import abipy.gui.awx as awx
import wx.lib.mixins.listctrl as listmix

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

        self.kpoints_listctrl = KpointsListCtrl(self, kpoints)
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


class KpointsFrame(awx.Frame):
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
        super(KpointsFrame, self).__init__(parent, **kwargs)

        panel = KpointsPanel(self, structure, kpoints)


##################
### Callbacks ###
##################


def showLittleGroup(parent, structure, kpoint):
    ltk = structure.spacegroup.find_little_group(kpoint)
    table, header = ltk.bilbao_character_table()
    awx.SimpleGridFrame(parent, table, title=header, labels_from_table=True).Show()


def showStar(parent, structure, kpoint):
    star = kpoint.compute_star(structure.fm_symmops)
    #wxdg.ScrolledMessageDialog(parent, str(star), caption="Star of the kpoint", style=wx.MAXIMIZE_BOX).Show()
    KpointsFrame(parent, structure, star).Show()


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



class SpinKpointBandPanel(awx.Panel):
    """
    This panel shows information on the k-points and the set of bands, spins. 
    Useful if we want to allow the user to select the set of indices (spin, kpt_idx, band).
    """
    def __init__(self, parent, nsppol, kpoints, mband, bstart=0, **kwargs):
        """
        Args:
            nsppol:
                Number of spins.
            kpoints:
                List of `Kpoint` instances.
            mband:
                Maximum band index.
            bstart:
                First band index.
        """
        super(SpinKpointBandPanel, self).__init__(parent, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs) 

        self.nsppol, self.kpoints, self.mband = nsppol, kpoints, mband
        self.bstart = bstart

        self.BuildUi()

    def BuildUi(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer1 = wx.BoxSizer(wx.HORIZONTAL)

        band_label = wx.StaticText(self, -1, "Band:", wx.DefaultPosition, wx.DefaultSize, 0)
        band_label.Wrap(-1)
        hsizer1.Add(band_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.band_cbox = wx.ComboBox(self, -1, choices=map(str, range(self.bstart, self.mband)))
        hsizer1.Add(self.band_cbox, 0, wx.ALL, 5)

        spin_label = wx.StaticText(self, -1, "Spin:", wx.DefaultPosition, wx.DefaultSize, 0)
        spin_label.Wrap(-1)
        hsizer1.Add(spin_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.spin_cbox = wx.ComboBox(self, -1, choices=map(str, range(self.nsppol)))
        hsizer1.Add(self.spin_cbox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer1, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        hsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        kpoint_label = wx.StaticText(self, -1, "Kpoint:", wx.DefaultPosition, wx.DefaultSize, 0)
        kpoint_label.Wrap(-1)
        hsizer2.Add(kpoint_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.kpoint_listctrl = klist = KpointsListCtrl(self, self.kpoints)

        hsizer2.Add(self.kpoint_listctrl, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)
        main_sizer.Add(hsizer2, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

        # Connect the events whose callback will be set by the client code.
        klist.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        klist.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnItemSelected)
        klist.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def GetSKB(self):
        """Returns the tuple (spin, kpoint, band) selected by the user."""
        spin = int(self.spin_cbox.GetValue())
        kpoint = self.kpoint_listctrl.GetKpoint()
        band = int(self.band_cbox.GetValue())
        return spin, kpoint, band

    def SetOnRightClick(self, callback):
        """
        Set the callback when EVT_LIST_ITEM_RIGHT_CLICK is fired.
        The callback expects the tuple (spin, kpoint, band)
        """
        self._on_item_right_click_callback = callback

    def OnRightClick(self, event):
        """Call the callback registered with `SetOnRightClick` (if any)."""
        if hasattr(self, "_on_item_right_click_callback"):
            #print("In OnRightClick with skb %s" % str(skb))
            skb = self.GetSKB()
            self._on_item_right_click_callback(*skb)

    def SetOnItemSelected(self, callback):
        """
        Set the callback when EVT_LIST_ITEM_SELECTED is fired.
        The callback expects the tuple (spin, kpoint, band)
        """
        self._on_item_selected_callback = callback

    def OnItemSelected(self, event):
        """Call the callback registered with `SetOnItemSelected` (if any)."""
        if hasattr(self, "_on_item_selected_callback"):
            #print("In OnItemSelected with skb %s" % str(skb))
            skb = self.GetSKB()
            self._on_item_selected_callback(*skb)

    def SetOnItemActivated(self, callback):
        """
        Set the callback when EVT_LIST_ITEM_ACTIVATED is fired (double click).
        The callback expects the tuple (spin, kpoint, band)
        """
        self._on_item_activated_callback = callback

    def OnItemActivated(self, event):
        """Call the callback registered with `SetOnItemActivated` (if any)."""
        if hasattr(self, "_on_item_activated_callback"):
            skb = self.GetSKB()
            #print("In OnItemActivated with skb %s" % str(skb))
            self._on_item_activated_callback(*skb)


class KpointsListCtrl(wx.ListCtrl, listmix.ColumnSorterMixin):
    """
    ListCtrl with info on the k-points. Support Column sorting.
    """
    def __init__(self, parent, kpoints, **kwargs):
        """
        Args:
            parent:
                Parent window.
            kpoints:
                List of `Kpoint` instances.
        """
        super(KpointsListCtrl, self).__init__(parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)

        self.kpoints = kpoints

        columns = ["#", 'Reduced coordinates', 'Weight', 'Name']

        for (index, col) in enumerate(columns):
            self.InsertColumn(index, col)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in columns]

        # Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
        self.itemDataMap = {}

        for (index, kpt) in enumerate(self.kpoints):
            entry = ["%d\t\t" % index, 
                     "[%.5f,  %.5f,  %.5f]\t\t" % tuple(c for c in kpt.frac_coords), 
                     "%.3f\t\t" % kpt.weight, 
                     "%s" % kpt.name,
                     ]
            self.Append(entry)
            self.itemDataMap[index] = entry
            self.SetItemData(index, index)

            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

        for (index, col) in enumerate(columns):
            self.SetColumnWidth(index, column_widths[index])

        # Now that the list exists we can init the other base class, see wx/lib/mixins/listctrl.py
        listmix.ColumnSorterMixin.__init__(self, len(columns))

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self

    def GetKpoint(self):
        """
        Returns the kpoint selected by the user.
        None if no selection has been done.
        """
        # Get selected index, map to index in kpoints and return the kpoint.
        item = self.GetFirstSelected()
        if item == -1: return None
        index = self.GetItemData(item)
        return self.kpoints[index]

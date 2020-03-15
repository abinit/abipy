"""Collection of widgets that allow the user to interact with list of K-points."""

import wx
import wx.lib.newevent

import abipy.gui.awx as awx
import wx.lib.mixins.listctrl as listmix


class KpointsListCtrl(wx.ListCtrl, listmix.ColumnSorterMixin):
    """
    ListCtrl that allows the user to interact with a list of k-points.
    Support column sorting.
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

    def getSelectedKpoint(self):
        """
        Returns the kpoint selected by the user.
        None if no selection has been done.
        """
        # Get selected index, map to index in kpoints and return the kpoint.
        item = self.GetFirstSelected()
        if item == -1: return None
        index = self.GetItemData(item)
        return self.kpoints[index]


class KpointsPanel(awx.Panel):
    """
    A panel with a list of k-points and a structure.
    Provides popup menus for retrieving information on a single k-point.
    """
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

        self.klist_ctrl = KpointsListCtrl(self, kpoints)
        self.structure = structure

        self.klist_ctrl.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.onRightClick)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.klist_ctrl, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)
        self.SetSizerAndFit(main_sizer)

    def onRightClick(self, event):
        """Generate the popup menu."""
        popup_menu = self.makePopupMenu()
        self.PopupMenu(popup_menu, event.GetPoint())
        popup_menu.Destroy()

    def getSelectedKpoint(self):
        return self.klist_ctrl.getSelectedKpoint()

    def makePopupMenu(self):
        """
        Build and return a popup menu. Subclasses can extend or replace this base method.
        """
        self.ID_POPUP_LITTLEGROUP = wx.NewId()
        self.ID_POPUP_STAR = wx.NewId()

        menu = wx.Menu()
        menu.Append(self.ID_POPUP_LITTLEGROUP, "Little group")
        menu.Append(self.ID_POPUP_STAR, "Show star")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_POPUP_LITTLEGROUP, self.onLittleGroup),
            (self.ID_POPUP_STAR, self.onShowStar),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu

    def onLittleGroup(self, event):
        """Show a table with the character table of the little group."""
        kpoint = self.klist_ctrl.getSelectedKpoint()
        if kpoint is None: return
        ltk = self.structure.abi_spacegroup.find_little_group(kpoint)
        table, header = ltk.bilbao_character_table()
        awx.SimpleGridFrame(self, table, title=header, labels_from_table=True).Show()

    def onShowStar(self, event):
        """Show a new `KpointsFrame` with the list of points in the star of the selected k-point."""
        kpoint = self.klist_ctrl.getSelectedKpoint()
        if kpoint is None: return
        star = kpoint.compute_star(self.structure.abi_spacegroup.fm_symmops)
        KpointsFrame(self, self.structure, star, title="Star of point: " + str(star.base_point)).Show()


class KpointsFrame(awx.Frame):
    """A frame with a `KpointsPanel`."""
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

        self.panel = KpointsPanel(self, structure, kpoints)


class SpinKpointBandPanel(awx.Panel):
    """
    This panel shows information on the k-points and the set of bands, spins.
    Useful if we want to allow the user to select the set of indices (spin, kpt_idx, band).
    """
    # Command event used to signal that the spin-kpoint-band has been selected.
    SkbActivatedEvent, MYEVT_SKB_ACTIVATED = wx.lib.newevent.NewCommandEvent()

    def __init__(self, parent, structure, nsppol, kpoints, mband, bstart=0, **kwargs):
        """
        Args:
            structure:
                Structure object.
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

        self.parent = parent
        self.nsppol, self.mband = nsppol, mband
        self.bstart = bstart

        # Build UI
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

        self.kpanel = KpointsPanel(self, structure, kpoints)
        klist_ctrl = self.kpanel.klist_ctrl

        hsizer2.Add(self.kpanel, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)
        main_sizer.Add(hsizer2, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

        # Connect the events whose callback will be set by the client code.
        klist_ctrl.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)
        #klist_ctrl.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        #klist_ctrl.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnItemSelected)

    def getSelectedKpoint(self):
        return self.kpanel.getSelectedKpoint()

    def GetSKB(self):
        """Returns the tuple (spin, kpoint, band) selected by the user."""
        spin = int(self.spin_cbox.GetValue())
        kpoint = self.getSelectedKpoint()
        band = int(self.band_cbox.GetValue())

        # Default values if no item is selected.:
        if spin == wx.NOT_FOUND: spin = 0
        if band == wx.NOT_FOUND: band = 0

        return spin, kpoint, band

    def OnItemActivated(self, event):
        """Fire SkbActivatedEvent."""
        #print("Firing SkbActivatedEvent")
        event = self.SkbActivatedEvent(id=-1, skb=self.GetSKB())
        wx.PostEvent(self.parent, event)

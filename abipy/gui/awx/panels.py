from __future__ import print_function, division

import wx

import wx.lib.mixins.listctrl as listmix

from .core import Panel

__all__ = [
    "SpinKpointBandPanel",
]


class SpinKpointBandPanel(Panel):
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
        super(SpinKpointBandPanel, self).__init__(parent, **kwargs)

        self.kpoints = kpoints

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer1 = wx.BoxSizer(wx.HORIZONTAL)

        band_label = wx.StaticText(self, -1, "Band:", wx.DefaultPosition, wx.DefaultSize, 0)
        band_label.Wrap(-1)
        hsizer1.Add(band_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.band_cbox = wx.ComboBox(self, -1, choices=map(str, range(bstart, mband)))
        hsizer1.Add(self.band_cbox, 0, wx.ALL, 5)

        spin_label = wx.StaticText(self, -1, "Spin:", wx.DefaultPosition, wx.DefaultSize, 0)
        spin_label.Wrap(-1)
        hsizer1.Add(spin_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.spin_cbox = wx.ComboBox(self, -1, choices=map(str, range(nsppol)))
        hsizer1.Add(self.spin_cbox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer1, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        hsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        kpoint_label = wx.StaticText(self, -1, "Kpoint:", wx.DefaultPosition, wx.DefaultSize, 0)
        kpoint_label.Wrap(-1)
        hsizer2.Add(kpoint_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.kpoint_listctrl = wx.ListCtrl(self, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN)
        hsizer2.Add(self.kpoint_listctrl, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)

        klist = self.kpoint_listctrl
        klist.InsertColumn(0, '#')
        klist.InsertColumn(1, 'Reduced Coordinates')
        klist.InsertColumn(2, 'Weight')
        for index, kpt in enumerate(kpoints):
            entry = map(str, [index, kpt.frac_coords, kpt.weight])
            klist.Append(entry)

        main_sizer.Add(hsizer2, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizer(main_sizer)
        main_sizer.Fit(self)

        # Connect the events whose callback will be set by the client code.
        klist.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        klist.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnItemSelected)
        klist.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def GetSKB(self):
        """Returns the tuple (spin, kpoint, band) selected by the user."""
        spin = int(self.spin_cbox.GetValue())
        kidx = int(self.kpoint_listctrl.GetFirstSelected())
        kpoint = self.kpoints[kidx]
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
            skb = self.GetSKB()
            self.log("In OnRightClick with skb %s" % str(skb))
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
            skb = self.GetSKB()
            self.log("In OnItemSelected with skb %s" % str(skb))
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
            self.log("In OnItemActivated with skb %s" % str(skb))
            self._on_item_activated_callback(*skb)

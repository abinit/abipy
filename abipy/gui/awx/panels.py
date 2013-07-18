from __future__ import print_function, division

import wx

import wx.lib.mixins.listctrl as listmix

__all__ = [
    "SpinKpointBandPanel",
]


class SpinKpointBandPanel(wx.Panel):
    """
    This panel shows information on the k-points and the set of bands, spins. 
    Useful if we want to allow the user to select the set of indices (spin, kpt_idx, band).
    """

    def __init__(self, parent, nsppol, kpoints, mband, **kwargs):
        """
        Args:
            nsppol:
                Number of spins.
            kpoints:
                List of `Kpoint` instances.
            mband:
                Maximum band index.
        """
        wx.Panel.__init__(self, parent, id=-1, **kwargs)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer1 = wx.BoxSizer(wx.HORIZONTAL)

        self.band_label = wx.StaticText(self, -1, u"Band:", wx.DefaultPosition, wx.DefaultSize, 0)
        self.band_label.Wrap(-1)
        hsizer1.Add(self.band_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.band_cbox = wx.ComboBox(self, -1, choices=map(str, range(mband)))
        hsizer1.Add(self.band_cbox, 0, wx.ALL, 5)

        self.spin_label = wx.StaticText(self, -1, u"Spin:", wx.DefaultPosition, wx.DefaultSize, 0)
        self.spin_label.Wrap(-1)
        hsizer1.Add(self.spin_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.spin_cbox = wx.ComboBox(self, -1, choices=map(str, range(nsppol)))
        hsizer1.Add(self.spin_cbox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer1, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        hsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        self.kpoint_label = wx.StaticText(self, -1, u"Kpoint:", wx.DefaultPosition, wx.DefaultSize, 0)
        self.kpoint_label.Wrap(-1)
        hsizer2.Add(self.kpoint_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        class MyListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin):
            """ Mixin class to resize the last column appropriately."""

            def __init__(self, parent):
                wx.ListCtrl.__init__(self, parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN)
                listmix.ListCtrlAutoWidthMixin.__init__(self)

        self.kpoint_list = MyListCtrl(self)
        hsizer2.Add(self.kpoint_list, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)

        klist = self.kpoint_list
        klist.InsertColumn(0, '#')
        klist.InsertColumn(1, 'Reduced Coords')
        klist.InsertColumn(2, 'Weight')
        for index, kpt in enumerate(kpoints):
            entry = map(str, [index, kpt.frac_coords, kpt.weight])
            klist.Append(entry)

        main_sizer.Add(hsizer2, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizer(main_sizer)
        self.Layout()
        main_sizer.Fit(self)

    def get_skb(self):
        """Returns a tuple with the set of indices (spin, kpoint_idx, band)"""
        spin = int(self.spin_cbox.GetValue())
        kidx = int(self.kpoint_list.GetFirstSelected())
        band = int(self.band_cbox.GetValue())
        return spin, kidx, band



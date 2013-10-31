from __future__ import print_function, division

import os
import wx

import abipy.gui.awx as awx

from collections import OrderedDict
from pymatgen.io.abinitio.abitimer import AbinitTimerSection
from abipy.abilab import abiopen


class AbinitTimerFrame(awx.Frame):
    """
    Frame to control to plot the timing data of a *single* ABINIT run.
    """
    def __init__(self, parent, filepath, **kwargs):
        """
        Args:
            parent:
                parent window
            filepath:
                Abinit output file.
        """
        filepath = os.path.abspath(filepath)
        if "title" not in kwargs:
            kwargs["title"] = "Abinit Timer: %s" % os.path.basename(filepath)

        super(AbinitTimerFrame, self).__init__(parent, **kwargs)

        try:
            abifile = abiopen(filepath)
            self.timer_data = abifile.timer_data

        except Exception as exc:
            raise awx.Error(str(exc))

        if not self.timer_data:
            raise awx.Error("%s does not contain a valid ABINIT TIMER section!" % filepath)

        self.BuildUi()

    def BuildUi(self):

        # Set callbacks (bound methods of AbiTimerData).
        timer_data = self.timer_data
        self.plot_types = OrderedDict([
            ("pie", timer_data.show_pie),
            ("stacked_hist", timer_data.show_stacked_hist),
        ])

        keys = AbinitTimerSection.NUMERIC_FIELDS

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        plot_label = wx.StaticText(self, -1, "plot type:", wx.DefaultPosition, wx.DefaultSize, 0)
        plot_label.Wrap(-1)
        hsizer.Add(plot_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.plot_cbox = wx.ComboBox(self, -1, "pie", wx.DefaultPosition, wx.DefaultSize, self.plot_types.keys(),
                                     0)
        hsizer.Add(self.plot_cbox, 0, wx.ALL, 5)

        key_label = wx.StaticText(self, -1, "key:", wx.DefaultPosition, wx.DefaultSize, 0)
        key_label.Wrap(-1)
        hsizer.Add(key_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.key_cbox = wx.ComboBox(self, -1, "wall_time", wx.DefaultPosition, wx.DefaultSize, keys, 0)
        hsizer.Add(self.key_cbox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        plot_button = wx.Button(self, -1, "Plot", wx.DefaultPosition, wx.DefaultSize, 0)
        plot_button.Bind(wx.EVT_BUTTON, self.OnClick)
        main_sizer.Add(plot_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

    def OnClick(self, event):
        callback = self.plot_types[self.plot_cbox.GetValue()]
        kwargs = dict(
            key=str(self.key_cbox.GetValue())
        )
        callback(**kwargs)


class MultiTimerFrame(awx.Frame):
    """
    Frame to control to plot the timing data of a *single* ABINIT run.
    """
    def __init__(self, parent, timers, **kwargs):
        """
        Args:
            parent:
                parent window
            timers:
        """
        super(MultiTimerFrame, self).__init__(parent, **kwargs)
        self.timers = timers

        self.BuildUi()

    def BuildUi(self):
        # Set callbacks (bound methods of AbinitTimer).
        self.plot_types = OrderedDict([
            ("efficiency", self.timers.show_efficiency),
            ("stacked_hist", self.timers.show_stacked_hist),
            ("pie", self.timers.show_pie),
        ])

        keys = AbinitTimerSection.NUMERIC_FIELDS

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        plot_label = wx.StaticText(self, -1, "plot type:", wx.DefaultPosition, wx.DefaultSize, 0)
        plot_label.Wrap(-1)
        hsizer.Add(plot_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.plot_cbox = wx.ComboBox(self, -1, "stacked_hist", wx.DefaultPosition, wx.DefaultSize, self.plot_types.keys(),
                                     0)
        hsizer.Add(self.plot_cbox, 0, wx.ALL, 5)

        key_label = wx.StaticText(self, -1, "key:", wx.DefaultPosition, wx.DefaultSize, 0)
        key_label.Wrap(-1)
        hsizer.Add(key_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        self.key_cbox = wx.ComboBox(self, -1, "wall_time", wx.DefaultPosition, wx.DefaultSize, keys, 0)
        hsizer.Add(self.key_cbox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        plot_button = wx.Button(self, -1, "Plot", wx.DefaultPosition, wx.DefaultSize, 0)
        plot_button.Bind(wx.EVT_BUTTON, self.OnClick)
        main_sizer.Add(plot_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

    def OnClick(self, event):
        callback = self.plot_types[self.plot_cbox.GetValue()]
        kwargs = dict(
            key=str(self.key_cbox.GetValue())
        )
        callback(**kwargs)

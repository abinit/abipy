from __future__ import print_function, division, unicode_literals, absolute_import

import wx
import pymatgen.core.units as units

from abipy.gui import awx


class ConverterFrame(awx.Frame):
    DEFAULT_UTYPE = "energy"

    def __init__(self, parent, **kwargs):
        super(ConverterFrame, self).__init__(parent, -1, title="Unit Converter")
        self.BuildUi()

    def BuildUi(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsz1 = wx.BoxSizer(wx.HORIZONTAL)

        from_text = wx.StaticText(self, -1, "From:")
        from_text.Wrap(-1)
        hsz1.Add(from_text, 0, wx.ALL, 5)

        self.from_textctrl = wx.TextCtrl(self, -1, "1.0")
        hsz1.Add(self.from_textctrl, 1, wx.ALL, 5)

        from_uchoices = units.ALL_UNITS[self.DEFAULT_UTYPE].keys()
        self.from_unit_choice = wx.ComboBox(self, -1, choices=from_uchoices)
        self.from_unit_choice.SetSelection(0)
        hsz1.Add(self.from_unit_choice, 1, wx.ALL, 5)

        main_sizer.Add(hsz1, 0.4, wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, 5)

        hsz2 = wx.BoxSizer(wx.HORIZONTAL)

        to_text = wx.StaticText(self, -1, "To:")
        to_text.Wrap(-1)
        hsz2.Add(to_text, 0, wx.ALL, 5)

        self.to_textctrl = wx.TextCtrl(self, -1, wx.EmptyString)
        hsz2.Add(self.to_textctrl, 1, wx.ALL, 5)

        to_uchoices = units.ALL_UNITS[self.DEFAULT_UTYPE].keys()
        self.to_unit_choice = wx.ComboBox(self, -1, choices=to_uchoices)
        self.to_unit_choice.SetSelection(0)
        hsz2.Add(self.to_unit_choice, 1, wx.ALL, 5)

        main_sizer.Add(hsz2, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, 5)

        hsz3 = wx.BoxSizer(wx.HORIZONTAL)

        utype_text = wx.StaticText(self, -1, "Unit type:")
        utype_text.Wrap(-1)
        hsz3.Add(utype_text, 0, wx.ALL, 5)

        unit_type_choices = units.ALL_UNITS.keys()
        self.unit_type_choice = wx.Choice(self, -1, wx.DefaultPosition, wx.DefaultSize, unit_type_choices, 0)
        self.unit_type_choice.SetSelection(unit_type_choices.index(self.DEFAULT_UTYPE))
        self.Bind(wx.EVT_CHOICE, self.OnUnitTypeChoiced, self.unit_type_choice)

        hsz3.Add(self.unit_type_choice, 0, wx.ALL, 5)

        convert_button = wx.Button(self, -1, "Convert")
        hsz3.Add(convert_button, 0, wx.ALL, 5)

        main_sizer.Add(hsz3, 1, wx.EXPAND, 5)

        self.Bind(wx.EVT_BUTTON, self.OnConvert, convert_button)

        self.SetSizerAndFit(main_sizer)

    def OnUnitTypeChoiced(self, event):
        new_type = self.unit_type_choice.GetStringSelection()
        new_choices = units.ALL_UNITS[new_type].keys()

        self.from_unit_choice.Clear()
        self.from_unit_choice.AppendItems(new_choices)
        self.from_unit_choice.SetSelection(0)

        self.to_unit_choice.Clear()
        self.to_unit_choice.AppendItems(new_choices)
        self.to_unit_choice.SetSelection(0)

    def OnConvert(self, event):
        try:
            fvalue = float(self.from_textctrl.GetValue())
        except:
            awx.showErrorMessage(self)
            return

        from_unit = self.from_unit_choice.GetValue()
        if from_unit == "foo":
            from abipy.gui.awx.eggs import Puzzle
            return Puzzle(self)

        ufloat = units.FloatWithUnit(fvalue, from_unit)
        to_unit = self.to_unit_choice.GetValue()

        try:
            conversion = ufloat.to(to_unit)
            print(ufloat, conversion)
            self.to_textctrl.SetValue(str(float(conversion)))

        except:
            awx.showErrorMessage(self)


def wxapp_converter():
    """Standalon application."""
    app = awx.App()
    frame = ConverterFrame(None)
    frame.Show()
    app.SetTopWindow(frame)
    return app


if __name__ == "__main__":
    wxapp_converter().MainLoop()

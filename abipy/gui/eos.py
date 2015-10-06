from __future__ import print_function, division, unicode_literals

import wx
import abipy.gui.awx as awx

from pymatgen.io.abinit.eos import EOS


class EosFrame(awx.Frame):
    """
    Frame that allows the user to fit E(V).
    """
    def __init__(self, parent, volumes, energies, vol_unit="ang^3", ene_unit="eV", **kwargs):
        super(EosFrame, self).__init__(parent, id=-1, title="EOS Frame", **kwargs)

        self.volumes, self.vol_unit = volumes, vol_unit
        self.energies, self.ene_unit = energies, ene_unit

        self.BuildUi()

    def BuildUi(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        models = EOS.MODELS.keys()
        self.model_choice = wx.Choice(self, -1, wx.DefaultPosition, wx.DefaultSize, models, 0)
        self.model_choice.SetSelection(0)
        main_sizer.Add(self.model_choice, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, 5)

        fit_button = wx.Button(self, -1, "Plot", wx.DefaultPosition, wx.DefaultSize, 0)
        fit_button.Bind(wx.EVT_BUTTON, self.OnFitButton)
        main_sizer.Add(fit_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

    def OnFitButton(self, event):
        model = self.model_choice.GetStringSelection()

        try:
            eos = EOS(eos_name=model)
            fit = eos.fit(self.volumes, self.energies, vol_unit=self.vol_unit, ene_unit=self.ene_unit)
            print(fit)
            fit.plot()

        except:
            awx.showErrorMessage(self)


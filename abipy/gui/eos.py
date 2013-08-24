from __future__ import print_function, division

import wx
import abipy.gui.awx as awx

from pymatgen.io.abinitio.eos import EOS


class EosFrame(awx.Frame):
    """
    Frame that allows the user to fit E(V).
    """
    def __init__(self, parent, volumes, energies, len_units="Ang", ene_units="eV", **kwargs):
        super(EosFrame, self).__init__(parent, -1, title="EOS Frame", **kwargs)

        self.volumes, self.len_units = volumes, len_units
        self.energies, self.ene_units = energies, ene_units

        self.BuildUi()

    def BuildUi(self):
        panel = wx.Panel(self, id=-1)
        models = EOS.MODELS.keys()

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.model_choice = wx.Choice(panel, -1, wx.DefaultPosition, wx.DefaultSize, models, 0)
        self.model_choice.SetSelection(0)
        main_sizer.Add(self.model_choice, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, 5)

        plot_button = wx.Button(panel, -1, "Plot", wx.DefaultPosition, wx.DefaultSize, 0)
        plot_button.Bind(wx.EVT_BUTTON, self.OnFitButton)
        main_sizer.Add(plot_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

    def OnFitButton(self, event):
        model = self.model_choice.GetStringSelection()

        try:
            eos = EOS(eos_name=model)
            fit = eos.fit(self.volumes, self.energies, len_units=self.len_units, ene_units=self.ene_units)
            print(fit)
            fit.plot()

        except EOS.Error as exc:
            awx.showErrorMessage(self, message=str(exc))


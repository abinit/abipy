from __future__ import print_function, division

import wx

from wx.lib.agw.floatspin import FloatSpin
from wxmplot import PlotApp, PlotFrame

from abipy import abiopen
from abipy.electrons import ElectronBands
from abipy.waves import WFK_File
import abipy.gui.awx as awx 


def showElectronDosFrame(parent, filepath):
    bands = ElectronBands.from_ncfile(filepath)
    ElectronDosFrame(bands=bands, parent=parent, title="File: %s" % filepath).Show()

def showElectronBandsPlot(parent, filepath):
    bands = ElectronBands.from_ncfile(filepath)
    bands.plot(title="File: %s" % filepath)

class ElectronDosFrame(wx.Frame):

    def __init__(self, bands, parent=None, title="Electron DOS"):
        super(ElectronDosFrame, self).__init__(parent, id=-1, title=title)
        self.bands = bands

        self.statusbar = self.CreateStatusBar() 

        # Widgets.
        self.panel = panel = wx.Panel(self, -1)
        sizer = wx.FlexGridSizer(rows=3, cols=2, vgap=5, hgap=5)

        label = wx.StaticText(panel, -1, "Broadening [eV]:")

        self.width = 0.2
        self.width_ctrl = FloatSpin(panel, id=-1, value=self.width, min_val=0.0, increment=0.1, digits=3)
        self.Bind(wx.EVT_SPINCTRL, self.onWidth, self.width_ctrl)

        sizer.AddMany([label, self.width_ctrl])

        label = wx.StaticText(panel, -1, "Mesh step [eV]:")
        self.step = 0.1
        self.step_ctrl = FloatSpin(panel, id=-1, value=self.step, min_val=0.0, increment=0.1, digits=3)
        self.Bind(wx.EVT_SPINCTRL, self.onStep, self.step_ctrl)

        sizer.AddMany([label, self.step_ctrl])

        dos_button = wx.Button(panel, -1, "Compute DOS", (20,100))
        self.Bind(wx.EVT_BUTTON, self.onClick, dos_button)

        self.oplot_checkbox = wx.CheckBox(panel, id=-1, label="Same plot")
        self.oplot_checkbox.SetValue(True)

        sizer.AddMany([dos_button, self.oplot_checkbox])

        panel.SetSizer(sizer)
        sizer.Layout()

    def onWidth(self, event):
        self.width = float(self.width_ctrl.GetValue())

    def onStep(self, event):
        self.step = float(self.step_ctrl.GetValue())
        if self.step == 1234: awx.tetris_game()

    def onClick(self, event):
        try:
            edos = self.bands.get_dos(step=self.step, width=self.width)
        except:
            awx.showErrorMessage(self)
            return

        tot_dos, tot_idos = edos.dos_idos()
        label = "$\sigma = %s, step = %s$" % (self.width, self.step)

        if self.has_plotframe and self.oplot_checkbox.GetValue():
            self.plotframe.oplot(tot_dos.mesh, tot_dos.values, label=label, draw_legend=True) 
        else:
            self.plotframe = plotframe = PlotFrame(parent=self)
            plotframe.plot(tot_dos.mesh, tot_dos.values, label=label, draw_legend=True)
            plotframe.Show()

    @property
    def has_plotframe(self):
        return hasattr(self, "plotframe")  

    @property
    def destroy_plotframe(self):
        if self.has_plotframe:
            self.plotframe.Destroy()
            del self.plotframe

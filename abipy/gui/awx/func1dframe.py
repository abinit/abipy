from __future__ import print_function, division

import wx
import inspect
import warnings

from dialogs import showErrorMessage

__all__ = [
    "Func1dPlotFrame",
]


class Func1dPlotFrame(wx.Frame):
    """
    Simple Frame that allows the user to have access to the method
    of a `Function1D` objects and plot the results with wxmplot.
    """

    def __init__(self, parent, func1d, **kwargs):
        """
        Args:
            parent:
                wx parent window.
            func1d:
                `Function1d` object.
        """
        super(Func1dPlotFrame, self).__init__(parent=parent, **kwargs)
        self.func1d = func1d
        self.BuildUi()

        # Store PlotFrames in this list.
        self._pframes = []

    def BuildUi(self):
        panel = wx.Panel(self, id=-1, size=wx.Size(500, 300))

        bSizer10 = wx.BoxSizer(wx.VERTICAL)

        # Get the methods of the objects.
        obj_methods = ["None"]
        obj_methods += [t[0] for t in inspect.getmembers(self.func1d, predicate=inspect.ismethod)
                        if not t[0].startswith("_")]
        sorted(obj_methods)

        self.method_choice = wx.Choice(panel, -1, wx.DefaultPosition, wx.DefaultSize, obj_methods, 0)
        self.method_choice.SetSelection(0)
        bSizer10.Add(self.method_choice, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, 5)

        fgSizer4 = wx.FlexGridSizer(0, 2, 0, 0)
        fgSizer4.SetFlexibleDirection(wx.BOTH)
        fgSizer4.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

        self.replot_checkbox = wx.CheckBox(panel, id=-1, label="Replot")
        self.replot_checkbox.SetValue(True)
        fgSizer4.Add(self.replot_checkbox, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.plot_button = wx.Button(panel, -1, u"Plot", wx.DefaultPosition, wx.DefaultSize, 0)
        fgSizer4.Add(self.plot_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        bSizer10.Add(fgSizer4, 0, 0, 5)

        self.SetSizerAndFit(bSizer10)

        self.plot_button.Bind(wx.EVT_BUTTON, self.OnClick)

    def OnClick(self, event):
        method = self.method_choice.GetStringSelection()

        try:
            from wxmplot import PlotFrame
        except ImportError:
            #warnings.warn("Error while importing wxmplot. Some features won't be available")
            raise

        plotframe = None
        if self.replot_checkbox.GetValue() and len(self._pframes):
            plotframe = self._pframes[-1]

        try:
            if method == "None":
                g = self.func1d
            else:
                g = getattr(self.func1d, method)()

            if plotframe is None:
                plotframe = PlotFrame(self)
                self._pframes.append(plotframe)
                plotframe.plot(g.mesh, g.values, label=method, draw_legend=True)
                plotframe.Show()
            else:
                plotframe.oplot(g.mesh, g.values, label=method, draw_legend=True)

        except:
            showErrorMessage(self)


"""Collections on controls."""
from __future__ import print_function, division

import wx
import collections

from abipy.tools import AttrDict

__all__ = [
    "LinspaceControl",
    "RowMultiCtrl",
    "TableMultiCtrl",
]


class LinspaceControl(wx.Panel):
    """
    This control merges two `SpinCtrlDouble` and a `SpinCtrl` to allow
    the user to specify a range using the `numpy.linspace` syntax.
    """
    # Default parameters passed to SpinCtrl and SpinCtrlDouble.
    SPIN_DEFAULTS = AttrDict(value=str(50), min=0, max=10000, initial=0)

    SPIN_DOUBLE_DEFAULTS = AttrDict(value=str(0.0), min=0, max=10000, initial=0, inc=1)

    def __init__(self, parent, start=None, stop=None, num=None):
        """
        value (string)  Default value (as text).
        min (float)  Minimal value.
        max (float)  Maximal value.
        initial (float)  Initial value.
        inc (float)  Increment value.
        """
        wx.Panel.__init__(self, parent, id=-1)

        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # start
        text = wx.StaticText(self, -1, "Start:")
        text.Wrap(-1)
        text.SetToolTipString("The starting value of the sequence.")

        p = self.SPIN_DOUBLE_DEFAULTS
        if start is not None:
            p.update(start)

        self.start_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        main_sizer.Add(self.start_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        # stop
        text = wx.StaticText(self, -1, "Stop:")
        text.Wrap(-1)
        text.SetToolTipString("The end value of the sequence")

        p = self.SPIN_DOUBLE_DEFAULTS
        if stop is not None:
            p.update(stop)

        self.stop_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        main_sizer.Add(self.stop_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        # num
        text = wx.StaticText(self, -1, "Num:")
        text.Wrap(-1)
        text.SetToolTipString("Number of samples to generate.")

        p = self.SPIN_DEFAULTS
        if num is not None:
            p.update(num)

        self.num_ctrl = wx.SpinCtrl(self, -1, **p)
        # FIXME:
        # There's a problem since the text entered in the SpinCtrl is processed
        # only when the widget looses focus. I tried the solution discussed at
        # https://groups.google.com/forum/#!topic/wxpython-users/Gud8PI6n-4E
        # but it didn't work on my Mac
        #txtctrl = self.num_ctrl.GetChildren[0]
        #txtctrl.WindowStyle |= wx.TE_PROCESS_ENTER

        main_sizer.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        main_sizer.Add(self.num_ctrl, 1, wx.EXPAND, 5)

        self.SetSizerAndFit(main_sizer)

    def GetLinspace(self):
        """Returns the numpy array built with numpy.linspace."""
        # FIXME Values are not updated if I edit the string in the SpinCtrl
        import numpy as np
        p = dict(start=self.start_ctrl.GetValue(),
                 stop=self.stop_ctrl.GetValue(),
                 num=self.num_ctrl.GetValue())
        #print(p)
        return np.linspace(**p)


class RowMultiCtrl(wx.Panel):
    """
    A panel with control widgets for integer, floats, ... placed on a row.
    """
    # Default parameters passed to SpinCtrl and SpinCtrlDouble.
    #SPIN_DEFAULTS = AttrDict(value=str(50), min=0, max=10000, initial=0)

    #SPIN_DOUBLE_DEFAULTS = AttrDict(value=str(0.0), min=0, max=10000, initial=0, inc=1)

    DEFAULT_WIDTH = 0.2
    DEFAULT_STEP = 0.1

    def __init__(self, parent, ctrls):
        """
        Args:
            ctrls:
                List whose items are in the form (dtype: params)
                where dtype is "f" for floats, "i" for integers.
                and params is a dictionary with the arguments used
                to build the controller. Available keys are listed below.
                Mandatory keys are explictly documented.

                ===========  ===================================
                label        label of the controller (mandatory)
                tooltip      tooltip of the controller
                ===========  ===================================

        Example:
            RowMultiCtrl(parent, [
                ("i", dict(label="I'm an integer", tooltip="hello integer)),
                ("f", dict(label="I'm a float")),
            ])
        """
        super(RowMultiCtrl, self).__init__(parent, -1)

        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.ctrls = collections.OrderedDict()
        for c in ctrls:
            dtype = c[0]
            params = AttrDict(**c[1])
            if not hasattr(params, "label"):
                raise ValueError("label must be specified")

            label = wx.StaticText(self, -1, params.label)
            label.Wrap(-1)

            # Set the tooltip
            tooltip = params.get("tooltip", None)
            if tooltip is not None:
                label.SetToolTipString(tooltip)

            # Create the controller and save it in self.ctrls
            if dtype == "f":
                ctrl = wx.SpinCtrlDouble(
                    self, id=-1, value=str(self.DEFAULT_WIDTH), min=self.DEFAULT_WIDTH/1000, inc=0.1)

            elif dtype == "i":
                ctrl = wx.SpinCtrl(self, id=-1, value="1", min=1)

            else:
                raise ValueError("Wrong dtype %s" % str(dtype))

            self.ctrls[params.label] = ctrl

            main_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
            main_sizer.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.SetSizerAndFit(main_sizer)

    def GetParams(self):
        """Return the parameters selected by the user in a `OrderedDict`"""
        tuples = [(label, ctrl.GetValue()) for label, ctrl in self.ctrls.items()]
        return collections.OrderedDict(tuples)

    def SetParams(self, **kwargs):
        """
        Set the value of the controllers from kwargs
        Returns the difference between the number of controllers
        and the number of keys in kwargs that have been set
        """
        count = 0
        for k, v in kwargs.items():
            if k in self.ctrls:
                count += 1
                self.ctrls[k].SetValue(v)

        return len(self.ctrls) - count


class TableMultiCtrl(wx.Panel):
    """
    A panel with a table of control widgets for integer, floats, ....
    """
    def __init__(self, parent, nrows, ctrls, rowmctrl_class=None):
        """
        Args:
            nrows:
                Number of rows
            ctrls:
                List whose items are in the form (dtype: param). See RowMultiCtrl
            rowmctrl_class:
                Subclass of RowMultiCtrl used to customize the behaviour of the
                widget (optional)
        """
        super(TableMultiCtrl, self).__init__(parent, -1)

        # Allow the user to customize RowMultiCtrl
        rowmctrl_class = RowMultiCtrl if rowmctrl_class is None else rowmctrl_class

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        #hsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.ctrl_list = []
        for row in range(nrows):
            rowmc = rowmctrl_class(self, ctrls)
            main_sizer.Add(rowmc)
            self.ctrl_list.append(rowmc)

        self.SetSizerAndFit(main_sizer)

    def GetParams(self):
        """Return the parameters selected by the user in a list of `AttrDict` dictionary"""
        olist = []
        for row in self.ctrl_list:
            od = {label: ctrl.GetValue() for label, ctrl in row.items()}
            olist.append(od)

        return olist

    def SetParams(self, ilist):
        """
        Set the value of the controllers.
        Returns the difference between the number of controls
        and the number of entries that have been set
        """
        count = 0
        for i, d in enumerate(ilist):
            ctrl = self.ctrl_list[i]
            for k, v in d.items():
                if k in ctrl:
                    count += 1
                    ctrl[k].SetValue(v)

        return len(self.ctrl_list) * len(self.ctrl_list[0]) - count

if __name__ == "__main__":
   app = wx.App()
   frame = wx.Frame(None)
   #panel = LinspaceControl(frame)
   #panel = RowMultiCtrl(frame, [
   #    ("f", dict(label="hello", tooltip="Tooltip for hello")),
   #    ("i", dict(label="integer")),
   #])

   panel = TableMultiCtrl(frame, 3, [
       ("f", dict(label="hello", tooltip="Tooltip for hello")),
       ("i", dict(label="integer")),
   ])

   frame.Show()
   app.MainLoop()



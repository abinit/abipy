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

        return np.linspace(**p)


#class CtrlParams(object)
#    def __init__(self):
#        self._params = []
#    def add(self, label, dtype, **kwargs):
#    @property
#    def opts(self):
#        # Make sure value is a string and create the Ctrl
#        _opts["value"] = str(_opts["value"])


class RowMultiCtrl(wx.Panel):
    """
    A panel with control widgets for integer, floats, ... placed on a row.
    """
    # Default parameters passed to SpinCtrl and SpinCtrlDouble.
    SPIN_DEFAULTS = AttrDict(value="0", min=0, max=10000, initial=1)
    SPIN_DOUBLE_DEFAULTS = AttrDict(value="0.0", min=-10000, max=10000, initial=0.0, inc=0.1)

    def __init__(self, parent, ctrl_params):
        """
        Args:
            ctrl_params:
                List whose items are in the form (label: params)
                where label is the name of the Spin button and
                and params is a dictionary with the arguments used
                to build the controller. Available keys are listed below.
                Note that dtype must be specified.

                ===========  ============================================================
                dtype        "f" for floats, "i" for integers, "cbox" for combo box
                tooltip      tooltip of the controller
                choices      list of possible choices (used if dtype == cbox, mandatory
                style        used if dtype == "cbox"
                ===========  ============================================================

        Example:
            RowMultiCtrl(parent, ctrl_params=[
                ("I'm an integer", dict(dtype="i", value="-1", tooltip="hello integer)),
                ("I'm a float", dict(dtype="f", value=str(1/3.))),
            ])
        """
        super(RowMultiCtrl, self).__init__(parent, -1)

        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.ctrls = collections.OrderedDict()

        # Accepts lists or tuples as well.
        if isinstance(ctrl_params, (list, tuple)):
            ctrl_params = collections.OrderedDict(ctrl_params)

        for label, params in ctrl_params.items():
            params = AttrDict(**params)

            dtype = params.pop("dtype", None)
            if dtype is None:
                raise ValueError("dtype must be specified")

            txt = wx.StaticText(self, -1, label)
            txt.Wrap(-1)

            # Set the tooltip
            tooltip = params.get("tooltip", None)
            if tooltip is not None:
                txt.SetToolTipString(tooltip)

            # Create the controller and save it in self.ctrls
            if dtype == "f":
                # Initialize default values then update them with the values in params.
                opts = self.SPIN_DOUBLE_DEFAULTS.copy()
                for k in opts:
                    if k in params:
                        opts[k] = params[k]

                # Make sure value is a string and create the Ctrl
                opts["value"] = str(opts["value"])
                ctrl = wx.SpinCtrlDouble(self, id=-1, **opts)

            elif dtype == "i":
                # Initialize default values then update them with the values in params.
                opts = self.SPIN_DEFAULTS.copy()
                for k in opts:
                    if k in params:
                        opts[k] = params[k]

                # Make sure value is a string and create the Ctrl
                opts["value"] = str(opts["value"])
                ctrl = wx.SpinCtrl(self, id=-1, **opts)

            elif dtype == "cbox":
                # Combo box
                if not hasattr(params, "choices"):
                    raise ValueError("choices must be specified if dtype == cbox")
                choices = params.choices
                ctrl = wx.ComboBox(self, id=-1, choices=choices, value=choices[0],
                                   style=params.get("style", wx.CB_READONLY))
            else:
                raise ValueError("Wrong dtype %s" % str(dtype))

            self.ctrls[label] = ctrl

            main_sizer.Add(txt, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
            main_sizer.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.SetSizerAndFit(main_sizer)

    def GetParams(self):
        """Return the parameters selected by the user in a `OrderedDict`"""
        tuples = [(label, ctrl.GetValue()) for label, ctrl in self.ctrls.items()]
        return collections.OrderedDict(tuples)

    def SetParams(self, d):
        """
        Set the value of the controllers from a dictionary
        Returns the difference between the number of controllers
        and the number of keys in kwargs that have been set
        """
        for k, v in d.items():
            self.ctrls[k].SetValue(v)


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
            olist.append(row.GetParams())

        return olist

    def SetParams(self, ilist):
        """
        Set the value of the controllers from a list of dictionaries
        """
        assert len(ilist) == len(self.ctrl_list)
        count = 0
        for i, d in enumerate(ilist):
            ctrl = self.ctrl_list[i]
            ctrl.SetParams(d)

if __name__ == "__main__":
   app = wx.App()
   frame = wx.Frame(None)
   #panel = LinspaceControl(frame)

   #panel = RowMultiCtrl(frame, [
   #    ("hello", dict(dtype="f", tooltip="Tooltip for hello", value=1/3.0)),
   #    ("integer", dict(dtype="i")),
   #    ("combo", dict(dtype="cbox", choices=["default", "another"])),
   #])

   panel = TableMultiCtrl(frame, 3, [
       ("hello", dict(dtype="f", tooltip="Tooltip for hello")),
       ("integer", dict(dtype="i", value=-1)),
   ])

   frame.Show()
   app.MainLoop()


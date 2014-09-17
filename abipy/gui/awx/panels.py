"""Collections on controls."""
from __future__ import print_function, division

import wx
import collections
import numpy as np

from monty.collections import AttrDict

__all__ = [
    "LinspaceControl",
    "ArangeControl",
    "IntervalControl",
    "RowMultiCtrl",
    "TableMultiCtrl",
    "ListCtrlFromTable",
    "FoldPanelMgr",
]


class LinspaceControl(wx.Panel):
    """
    This control merges two `SpinCtrlDouble` and a `SpinCtrl` to allow
    the user to specify a range using the `numpy.linspace` syntax.
    """
    # Default parameters passed to SpinCtrl and SpinCtrlDouble.
    SPIN_DEFAULTS = dict(value=str(50), min=0, max=10000, initial=0)

    SPIN_DOUBLE_DEFAULTS = dict(value=str(0.0), min=0, max=10000, initial=0, inc=1)

    def __init__(self, parent, start=None, stop=None, num=1, **kwargs):
        """
        value (string)  Default value (as text).
        min (float)  Minimal value.
        max (float)  Maximal value.
        initial (float)  Initial value.
        inc (float)  Increment value.
        """
        super(LinspaceControl, self).__init__(parent, id=-1, **kwargs)

        text_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        ctrl_opts = text_opts

        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # start
        text = wx.StaticText(self, -1, "Start:")
        text.Wrap(-1)
        text.SetToolTipString("The starting value of the sequence.")

        p = self.SPIN_DOUBLE_DEFAULTS.copy()
        if start is not None:
            p["value"] = str(start)

        self.start_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.start_ctrl, **ctrl_opts)

        # stop
        text = wx.StaticText(self, -1, "Stop:")
        text.Wrap(-1)
        text.SetToolTipString("The end value of the sequence")

        p = self.SPIN_DOUBLE_DEFAULTS.copy()
        if stop is not None:
            p["value"] = str(stop)

        self.stop_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.stop_ctrl, **ctrl_opts)

        # num
        text = wx.StaticText(self, -1, "Num:")
        text.Wrap(-1)
        text.SetToolTipString("Number of samples to generate.")

        p = self.SPIN_DEFAULTS.copy()
        if num is not None:
            p["value"] = str(num)

        self.num_ctrl = wx.SpinCtrl(self, -1, **p)
        # FIXME:
        # There's a problem since the text entered in the SpinCtrl is processed
        # only when the widget looses focus. I tried the solution discussed at
        # https://groups.google.com/forum/#!topic/wxpython-users/Gud8PI6n-4E
        # but it didn't work on my Mac
        #txtctrl = self.num_ctrl.GetChildren[0]
        #txtctrl.WindowStyle |= wx.TE_PROCESS_ENTER

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.num_ctrl, **ctrl_opts)

        self.SetSizerAndFit(main_sizer)

    def getValues(self):
        """Returns the numpy array built with numpy.linspace."""
        # FIXME Values are not updated if I edit the string in the SpinCtrl
        p = dict(start=self.start_ctrl.GetValue(),
                 stop=self.stop_ctrl.GetValue(),
                 num=self.num_ctrl.GetValue())

        return np.linspace(**p)


class ArangeControl(wx.Panel):
    """
    This control merges three `SpinCtrlDouble` controls to allow
    the user to specify a range using the `numpy.arange` syntax.
    """
    # Default parameters passed to SpinCtrlDouble.
    SPIN_DOUBLE_DEFAULTS = dict(value=str(0.0), min=0, max=10000, initial=0, inc=1)

    def __init__(self, parent, start=None, stop=None, step=None, **kwargs):
        """
        """
        super(ArangeControl, self).__init__(parent, id=-1, **kwargs)

        text_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        ctrl_opts = text_opts

        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # start
        text = wx.StaticText(self, -1, "Start:")
        text.Wrap(-1)
        text.SetToolTipString("Start of interval. The interval includes this value.")

        p = self.SPIN_DOUBLE_DEFAULTS.copy()
        if start is not None:
            p["value"] = str(start)

        self.start_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.start_ctrl, **ctrl_opts)

        # stop
        text = wx.StaticText(self, -1, "Stop:")
        text.Wrap(-1)
        text.SetToolTipString("""\
End of interval. The interval does not include this value, except
in some cases where `step` is not an integer and floating point round-off affects the length of `out`.""")

        p = self.SPIN_DOUBLE_DEFAULTS.copy()
        if stop is not None:
            p["value"] = str(stop)

        self.stop_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.stop_ctrl, **ctrl_opts)

        # step
        text = wx.StaticText(self, -1, "Step:")
        text.Wrap(-1)
        text.SetToolTipString("""\
Spacing between values. For any output `out`, this is the distance between two adjacent values, out[i+1] - out[i].""")

        p = self.SPIN_DOUBLE_DEFAULTS.copy()
        if step is not None:
            p["value"] = str(step)

        self.step_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.step_ctrl, **ctrl_opts)

        self.SetSizerAndFit(main_sizer)

    def getValues(self):
        """Returns the numpy array built with numpy.linspace."""
        p = dict(start=self.start_ctrl.GetValue(),
                 stop=self.stop_ctrl.GetValue(),
                 step=self.step_ctrl.GetValue())

        return np.arange(**p)


class IntervalControl(wx.Panel):
    """
    This control merges two `SpinCtrlDouble` controls, a `SpinCtrl` and a combo box
    to allow the user to specify a numerical interval
    """
    # Default parameters passed to SpinCtrl and SpinCtrlDouble.
    SPIN_DEFAULTS = dict(value=str(50), min=0, max=10000, initial=0)

    SPIN_DOUBLE_DEFAULTS = dict(value=str(0.0), min=0, max=10000, initial=0, inc=0.1)

    def __init__(self, parent, start, step, num=5, choices=None, **kwargs):
        """
        Args:
            start: Initial value
            step: Step used to generate the mesh
            num: Number of points
        """
        super(IntervalControl, self).__init__(parent, id=-1, **kwargs)

        text_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        ctrl_opts = text_opts

        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # start
        text = wx.StaticText(self, -1, "Start:")
        text.Wrap(-1)
        text.SetToolTipString("Start of interval. The interval includes this value.")

        p = self.SPIN_DOUBLE_DEFAULTS.copy()
        p["value"] = str(start)

        self.start_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.start_ctrl, **ctrl_opts)

        # step
        text = wx.StaticText(self, -1, "Step:")
        text.Wrap(-1)
        text.SetToolTipString("""\
Spacing between values. For any output `out`, this is the distance between two adjacent values, out[i+1] - out[i].""")

        p = self.SPIN_DOUBLE_DEFAULTS.copy()
        p["value"] = str(step)

        self.step_ctrl = wx.SpinCtrlDouble(self, -1, **p)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.step_ctrl, **ctrl_opts)

        # num
        text = wx.StaticText(self, -1, "Num:")
        text.Wrap(-1)
        text.SetToolTipString("Number of samples to generate.")

        p = self.SPIN_DEFAULTS.copy()
        if num is not None:
            p["value"] = str(num)

        self.num_ctrl = wx.SpinCtrl(self, -1, **p)
        # FIXME:
        # There's a problem since the text entered in the SpinCtrl is processed
        # only when the widget looses focus. I tried the solution discussed at
        # https://groups.google.com/forum/#!topic/wxpython-users/Gud8PI6n-4E
        # but it didn't work on my Mac
        #txtctrl = self.num_ctrl.GetChildren[0]
        #txtctrl.WindowStyle |= wx.TE_PROCESS_ENTER

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.num_ctrl, **ctrl_opts)

        # Combo box to select the interval type.
        text = wx.StaticText(self, -1, "Interval type:")
        text.SetToolTipString("""\
Select the interval type:
centered if values are centered on start,
> if the lower bound of the interval is start, < if the uppper bound is start.""")

        if choices is None: choices = ["centered", ">", "<"]
        self.interval_type = wx.ComboBox(
            self, id=-1, name='Interval type type', choices=choices, value=choices[0], style=wx.CB_READONLY)

        main_sizer.Add(text, **text_opts)
        main_sizer.Add(self.interval_type, **ctrl_opts)

        self.SetSizerAndFit(main_sizer)

    def getValues(self):
        """Returns a numpy array with the interval."""
        start = self.start_ctrl.GetValue()
        step = self.step_ctrl.GetValue()
        num = self.num_ctrl.GetValue()
        int_type = self.interval_type.GetValue()

        assert num > 0
        if num == 1:
            return np.array(start)

        if int_type in [">", "<"]:
            sign = 1 if int_type == ">" else -1
            values = [start + i * sign * step for i in range(num)]
            if int_type == "<":
                values.reverse()

        elif int_type == "centered":
            v_left = [start - i * step for i in range(1, num)]
            v_left.reverse()
            v_right = [start + i * step for i in range(num)]
            values = v_left + v_right

        else:
            raise ValueError("Wrong int_type %s" % int_type)

        #print("values", values)
        return np.array(values)


class RowMultiCtrl(wx.Panel):
    """A panel with control widgets for integer, floats, ... placed on a row."""
    # Default parameters passed to SpinCtrl and SpinCtrlDouble.
    SPIN_DEFAULTS = dict(value="0", min=0, max=10000, initial=1)

    SPIN_DOUBLE_DEFAULTS = dict(value="0.0", min=-10000, max=10000, initial=0.0, inc=0.1)

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

        # Accepts lists or tuples as well.
        if isinstance(ctrl_params, (list, tuple)):
            ctrl_params = collections.OrderedDict(ctrl_params)

        self.ctrl_params = ctrl_params
        self.ctrls = collections.OrderedDict()

        self.main_sizer = main_sizer = wx.BoxSizer(wx.HORIZONTAL)

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

            hsizer = wx.BoxSizer(wx.HORIZONTAL)
            hsizer.Add(txt, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
            hsizer.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

            main_sizer.Add(hsizer)

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
        self.rowmctrl_class = RowMultiCtrl if rowmctrl_class is None else rowmctrl_class

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.ctrl_list, self.ctrls = [], ctrls
        for row in range(nrows):
            rowmc = self.rowmctrl_class(self, ctrls)
            self.main_sizer.Add(rowmc, 1, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL)
            self.ctrl_list.append(rowmc)

        self.SetSizerAndFit(self.main_sizer)

    def __getitem__(self, item):
        return self.ctrl_list.__getitem__(item)

    def __iter__(self):
        return self.ctrl_list.__iter__()

    def __len__(self):
        return len(self.ctrl_list)

    def removeRow(self, index=-1):
        """Remove row with the give index (default: last row)."""
        if index == -1: index = len(self) - 1
        if index >= len(self) or len(self) == 0: return
        #print("removing row with index: %s/%s" % (index, len(self)))

        self.main_sizer.Hide(index)
        self.main_sizer.Remove(index)
        self.ctrl_list.pop(index)
        self.main_sizer.Layout()
                                           
    def appendRow(self):
        """Add a new row."""
        rowmc = self.rowmctrl_class(self, self.ctrls)
        self.main_sizer.Add(rowmc, 1, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL)
        self.ctrl_list.append(rowmc)
        self.main_sizer.Layout()

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
        for i, d in enumerate(ilist):
            ctrl = self.ctrl_list[i]
            ctrl.SetParams(d)


import wx.lib.mixins.listctrl as listmix
from abipy.gui.awx.core import get_width_height


class ListCtrlFromTable(wx.ListCtrl, listmix.ColumnSorterMixin, listmix.ListCtrlAutoWidthMixin):
    """
    This ListCtrl displays a list of strings. Support column sorting.
    """
    def __init__(self, parent, table, **kwargs):
        """
        Args:
            parent:
                Parent window.
            table:
                Table of strings (List of lists). The first item contains the names of the columns.
            kwargs:
                enumrows:
                    If True, the first column will contain the row index (default)
        """
        enumrows = kwargs.pop("enumrows", True)
        super(ListCtrlFromTable, self).__init__(parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)

        # Make sure we have a list of strings
        for i, row in enumerate(table):
            table[i] = map(str, row)

        columns = table[0]
        if enumrows:
            columns = ["#"] + columns

        for (index, col) in enumerate(columns):
            self.InsertColumn(index, col)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [get_width_height(self, s)[0] for s in columns]

        # Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
        self.itemDataMap = {}

        for (index, entry) in enumerate(table[1:]):
            if enumrows:
                entry = [str(index)] + entry
            self.Append(entry)
            self.itemDataMap[index] = entry
            self.SetItemData(index, index)

            w = [get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

        for (index, col) in enumerate(columns):
            self.SetColumnWidth(index, column_widths[index])

        # Now that the list exists we can init the other base class, see wx/lib/mixins/listctrl.py
        listmix.ColumnSorterMixin.__init__(self, len(columns))
        listmix.ListCtrlAutoWidthMixin.__init__(self)

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self


import wx.lib.foldpanelbar as foldpanel

class FoldPanelMgr(foldpanel.FoldPanelBar): 
    """
    Fold panel that manages a collection of Panels

    The FoldPanelBar is a custom container class that allows multiple controls 
    to be grouped together into FoldPanelItem controls that allow them to be expanded 
    or contracted by clicking on its CaptionBar. 
    The FoldPanelBar doesn't work with layouts based on a Sizer and as such its API can get a little cumbersome, 
    because it requires you to add each control one by one and set its layout by using various flags. 
    This recipe shows how to create a custom FoldPanelBar that works with Panel objects. 
    This class will allow for you to modularize your code into Panel classes and then just add them to the FoldPanelBar 
    instead of directly adding everything to the FoldPanelBar itself.
    """ 
    def __init__(self, parent, *args, **kwargs):
        super(FoldPanelMgr, self).__init__(parent, *args, **kwargs)

    def addPanel(self, pclass, title="", collapsed=False): 
        """
        Add a panel to the manager 
        @param pclass: Class constructor (callable) 
        @keyword title: foldpanel title
        @keyword collapsed: start with it collapsed @return: pclass instance
        """ 
        fpitem = self.AddFoldPanel(title, collapsed=collapsed) 
        wnd = pclass(fpitem) 
        best = wnd.GetBestSize() 
        wnd.SetSize(best) 
        self.AddFoldPanelWindow(fpitem, wnd) 

        return wnd


if __name__ == "__main__":
    app = wx.App()
    frame = wx.Frame(None)

    #table = [["ciao", "values"], [1,2], ["3", "4"]]
    #ctrl = ListCtrlFromTable(frame, table)

    #panel = LinspaceControl(frame)
    #panel = ArangeControl(frame, start=10, stop=15, step=1)

    #panel = IntervalControl(frame, start=10, step=1, num=2)

    #panel = RowMultiCtrl(frame, [
    #    ("hello", dict(dtype="f", tooltip="Tooltip for hello", value=1/3.0)),
    #    ("integer", dict(dtype="i")),
    #    ("combo", dict(dtype="cbox", choices=["default", "another"])),
    #])

    panel = TableMultiCtrl(frame, 1, [
        ("hello", dict(dtype="f", tooltip="Tooltip for hello")),
        ("integer", dict(dtype="i", value=-1)),
    ])
    panel.appendRow()
    panel.removeRow()

    frame.Show()
    app.MainLoop()

#!/usr/bin/env python
"""Gui for the oncvpsp norm-conserving pseudopotential generator."""
from __future__ import print_function, division

import os
import copy
import time
import abc
import wx
import awx
import numpy as np
import pymatgen.core.periodic_table as periodic_table

import wx.lib.mixins.listctrl as listmix
from collections import OrderedDict
from monty.dev import get_ncpus
from abipy.tools import AttrDict
from abipy.gui.editor import TextNotebookFrame, SimpleTextViewer
from abipy.gui.oncv_tooltips import oncv_tip
from abipy.gui import mixins as mix
from pseudo_dojo.ppcodes.ppgen import OncvGenerator
from pseudo_dojo.ppcodes.oncvpsp import MultiPseudoGenDataPlotter

# TODO
# Change oncvpsp so that 
#   1) we always write the logarithmic derivatie
#   2) better error handling


def add_size(kwargs, size=(800, 600)):
    """Add size to kwargs if not present."""
    if "size" not in kwargs:
        kwargs["size"] = size

    return kwargs

from awx.elements_gui import WxPeriodicTable, PeriodicPanel, ElementButton

def my_periodic_table(parent):
    """
    A periodic table that allows the user to select the element
    before starting the pseudopotential genearation.
    """
    class MyElementButton(ElementButton):

        def makePopupMenu(self):
            # Get the menu of the super class.
            menu = super(MyElementButton, self).makePopupMenu()

            self.ID_POPUP_ONCVPSP = wx.NewId()
            menu.Append(self.ID_POPUP_ONCVPSP, "Generate NC pseudo with oncvpsp")

            # Associate menu/toolbar items with their handlers.
            menu_handlers = [
                (self.ID_POPUP_ONCVPSP, self.onOncvpsp),
            ]
                                                                
            for combo in menu_handlers:
                mid, handler = combo[:2]
                self.Bind(wx.EVT_MENU, handler, id=mid)

            return menu

        def onOncvpsp(self, event):
            """Open a frame for the initialization of oncvpsp."""
            frame = OncvParamsFrame(self, self.Z)
            frame.Show()

    class MyPeriodicPanel(PeriodicPanel):
        element_button_class = MyElementButton

        def OnSelect(self, event):
            # Get the atomic number Z, open a dialog to get basic configuration parameters from the user. 
            # The dialog will then generate the main Frame for the pseudo generation.
            super(MyPeriodicPanel, self).OnSelect(event)
            z = event.GetId() - 100
            print("select", z)

    class MyPeriodicTable(WxPeriodicTable):
        periodic_panel_class = MyPeriodicPanel

    return MyPeriodicTable(parent)


class OncvParamsFrame(awx.Frame):
    """
    This frame allows the user to specify the most important parameters
    used to generate the pseudopotential once the chemical element has been selected.
    """

    HELP_MSG = """\
Quick help:

    Use this window to select the AE reference configuration and how 
    to separate states into core and valence.
"""

    def __init__(self, parent, z, **kwargs):
        super(OncvParamsFrame, self).__init__(parent, **kwargs)
        self.element = periodic_table.Element.from_Z(z)
        self.buildUI()

    def buildUI(self):
        # Build controller with the dimensions.
        panel = wx.Panel(self, -1)
        self.wxdims = awx.RowMultiCtrl(self, ctrl_params=OrderedDict([
                    ("nc", dict(dtype="i", tooltip="Number of core states")),
                    ("nv", dict(dtype="i", tooltip="Number of valence states")),
                    ("lmax", dict(dtype="i", tooltip="Maximum angular momentum for pseudo"))
               ]))

        # Initialize the quantum numbers of the AE atom with the ground-state configuration.
        # E.g., The electronic structure for Fe is represented as:
        # [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 6), (4, "s", 2)]
        ele_struct = self.element.full_electronic_structure

        ctrls = OrderedDict([
            ("n", dict(dtype="i")),
            ("l", dict(dtype="i")),
            ("f", dict(dtype="f"))])

        self.wxaeconf = awx.TableMultiCtrl(self, nrows=len(ele_struct), ctrls=ctrls)

        for wxrow, (n, lchar, f) in zip(self.wxaeconf, ele_struct):
            row = OrderedDict()
            row["n"], row["l"], row["f"] = n, periodic_table.char2l(lchar), f
            wxrow.SetParams(row)

        add_button = wx.Button(self, -1, "Add row")
        add_button.Bind(wx.EVT_BUTTON, self.onAddButton)
        del_button = wx.Button(self, -1, "Delete row")
        del_button.Bind(wx.EVT_BUTTON, self.onDelButton)
        hsz = wx.BoxSizer(wx.HORIZONTAL)
        hsz.Add(add_button)
        hsz.Add(del_button)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(hsz, 0,flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)
        main_sizer.Add(self.wxdims, 0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)
        main_sizer.Add(self.wxaeconf, 0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        help_button = wx.Button(self, wx.ID_HELP)
        help_button.Bind(wx.EVT_BUTTON, self.onHelp)
        main_sizer.Add(help_button, 0, flag=wx.ALL | wx.ALIGN_RIGHT)

        self.SetSizerAndFit(main_sizer)

    def get_oncv_params(self):
        """Return the basic dimensions used in oncvpsp.""" 
        print(self.wxaeconf.GetParams())
        return AttrDict(
            dims=self.wxdims.GetParams(),
        )

    def show_nist_lda_levels(self):
        # Get the LDA levels of the neutral atom.
        # (useful to decide if semicore states should be included in the valence).
        from pseudo_dojo.refdata.nist import nist_database as nist
        entry = nist.get_neutral_entry(symbol=self.element.symbol)
        frame = awx.Frame(self)
        awx.ListCtrlFromTable(frame, table=entry.to_table())
        frame.Show()

    def onAddButton(self, event):
        """Add a new row."""
        self.get_oncv_params()
                                  
    def onDelButton(self, event):
        """Delete last row."""
        self.show_nist_lda_levels()


class WxOncvFrame(awx.Frame, mix.Has_Tools):
    """The main frame of the GUI"""
    VERSION = "0.1"

    HELP_MSG = """\
This window shows a template input file with the variables
used to generated the pseudopotential. The `optimize` buttons
allows you to scan a set of possible values for the generation of the pseudopotential.
"""

    def __init__(self, parent):
        super(WxOncvFrame, self).__init__(parent, id=-1, title=self.codename)

        # This combination of options for config seems to work on my Mac.
        self.config = wx.FileConfig(appName=self.codename, localFilename=self.codename + ".ini", 
                                    style=wx.CONFIG_USE_LOCAL_FILE)

        # Build menu, toolbar and status bar.
        self.SetMenuBar(self.makeMenu())
        self.statusbar = self.CreateStatusBar()

        self.BuildUI()

    @property
    def codename(self):
        """Name of the application."""
        return "WxOncvGui"

    def BuildUI(self):
        """Build user-interface."""
        oncv_dims = dict(nc=1, nv=2, lmax=1, ncfn=0)

        self.makeToolBar()
        self.Centre()

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # The panel with the input variables.
        self.wxoncv_input = WxOncvInput(self, oncv_dims)
        main_sizer.Add(self.wxoncv_input, flag=wx.EXPAND)

        self.SetSizerAndFit(main_sizer)

    def AddFileToHistory(self, filepath):
        """Add the absolute filepath to the file history."""
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        self.file_history.AddFileToHistory(filepath)
        #newpanel = WxOncvInput.from_file(filepath)

    def makeMenu(self):
        """Creates the main menu."""
        menu_bar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an input file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the Gui")
        #file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")

        file_history = self.file_history = wx.FileHistory(8)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(-1, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)
        menu_bar.Append(file_menu, "File")

        # Add Mixin menus.
        menu_bar.Append(self.CreateToolsMenu(), "Tools")

        help_menu = wx.Menu()
        help_menu.Append(wx.ID_HELP, "Help ", help="Quick help")
        help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menu_bar.Append(help_menu, "Help")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_OPEN, self.OnOpen),
            (wx.ID_CLOSE, self.OnClose),
            #(wx.ID_EXIT, self.OnExit),
            (wx.ID_HELP, self.onHelp),
            (wx.ID_ABOUT, self.OnAbout),
        ]
                                                            
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu_bar

    def makeToolBar(self):
        """Creates the toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        self.toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        self.ID_SHOW_INPUT = wx.NewId()
        self.ID_RUN_INPUT = wx.NewId()

        toolbar.AddSimpleTool(self.ID_SHOW_INPUT, bitmap("in.png"), "Visualize the input file(s) of the workflow.")
        toolbar.AddSimpleTool(self.ID_RUN_INPUT, bitmap("run.png"), "Run the input file.")

        toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_SHOW_INPUT, self.OnShowInput),
            (self.ID_RUN_INPUT, self.onRunInput),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

    def OnShowInput(self, event):
        """Show the input file in a new frame."""
        text = self.wxoncv_input.makeInputString()
        SimpleTextViewer(self, text=text, title="Oncvpsp Input").Show()

    def onRunInput(self, event):
        """Build a new generator from the input file, and add it to the queue."""
        text = self.wxoncv_input.makeInputString()
        try:
            psgen = OncvGenerator(text, calc_type=self.wxoncv_input.calc_type)
        except:
            awx.showErrorMessage(self)
            return

        psgen.start()
        psgen.wait()

        if psgen.status == psgen.S_OK:
            psgen.plot_results()
        else:
            SimpleTextViewer(self, text=psgen.get_stdout()).Show()

        #PseudoGeneratorsFrame(self, [psgen], title="Run Input").Show()

    #def _onOptimize_simple_key(self, key):
    #    """Helper function for simple optimizations."""
    #    #if self._warn_inoptrun(): return

    #    template = self.wxoncv_input.makeInput()
    #    # Build the PseudoGeneratorFrame and show it.
    #    # Note how we select the method to call from key.
    #    psgens = []
    #    method = getattr(template, "optimize_" + key)
    #    for inp in method():
    #        psgen = OncvGenerator(str(inp), calc_type=self.wxoncv_input.calc_type)
    #        psgens.append(psgen)

    #    PseudoGeneratorsFrame(self, psgens, title="%s Optimization" % key).Show()

    #def OnOptimizeVloc(self, event):
    #    """Open a new frame for the optimization of Vlocal."""
    #    self._onOptimize_simple_key("vloc")

    def OnOpen(self, event):
        """Open a file"""
        dialog = wx.FileDialog(self, message="Choose an inputfile", defaultDir=os.getcwd(),
                               style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)

        # Show the dialog and retrieve the user response. 
        # If it is the OK response, process the data.
        if dialog.ShowModal() == wx.ID_CANCEL: return 

        filepath = dialog.GetPath()
        dialog.Destroy()

        # Add to the history.
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

        #newpanel = WxOncvInput.from_file(filepath)

    def OnClose(self, event):
        """ Respond to the "Close" menu command."""
        self.Destroy()

    def OnAbout(self, event):
        return awx.makeAboutBox(
            codename=self.codename,
            version=self.VERSION,
            description="oncvgui is a front-end for the pseudopotential generator oncvpsp",
            developers=["Matteo Giantomassi"],
            website="http://www.mat-simresearch.com/")


class OptimizationFrame(awx.Frame):
    """Base class for optimization frames."""
    __metaclass__ = abc.ABCMeta

    def __init__(self, parent, **kwargs):
        super(OptimizationFrame, self).__init__(parent, **kwargs)

        # All optimization buttons are disabled when we start an optimization.
        self.main_frame.wxoncv_input.enable_all_optimize_buttons(False)
        #self.main_frame.Enable(False)
        #self.Bind(wx.EVT_WINDOW_DESTROY, self.onDestroy)

    @property
    def main_frame(self):
        return self.getParentWithType(WxOncvFrame)

    def onDestroy(self, event):
        """Enable all optimize_buttons before destroying the Frame."""
        self.main_frame.wxoncv_input.enable_all_optimize_buttons(True)
        return super(OptimizationFrame, self).Destroy()

    def onCloseButton(self, event):
        self.onDestroy(event)

    def onOkButton(self, event):
        """
        Get input from user, generate new input files by changing some parameters 
        and open a new frame for running the calculations.
        """
        # Build the PseudoGenerators and open a new frame to run them.
        psgens = []
        for inp in self.build_new_inps():
            try:
                psgen = OncvGenerator(str(inp), calc_type=self.wxoncv_input.calc_type)
                psgens.append(psgen)
            except:
                awx.showErrorMessage(self)
                return
                                                                                   
        PseudoGeneratorsFrame(self, psgens, title=self.opt_type).Show()

    def make_buttons(self, parent=None):
        """
        Build the three buttons (Generate, Cancel, Help), binds them and return the sizer.
        """
        parent = self if parent is None else parent

        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        gen_button = wx.Button(parent, wx.ID_OK, label='Generate')
        gen_button.Bind(wx.EVT_BUTTON, self.onOkButton)

        close_button = wx.Button(parent, wx.ID_CANCEL, label='Cancel')
        close_button.Bind(wx.EVT_BUTTON, self.onCloseButton)

        help_button = wx.Button(parent, wx.ID_HELP)
        help_button.Bind(wx.EVT_BUTTON, self.onHelp)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(gen_button, **add_opts)
        hbox.Add(close_button, **add_opts)
        hbox.Add(help_button, **add_opts)

        return hbox

    @abc.abstractproperty
    def opt_type(self):
        """Human-readable string describing the optimization type."""

    @abc.abstractmethod
    def build_new_inps(self):
        """Returns a list of new inputs."""


class LlocOptimizationFrame(OptimizationFrame):
    """
    This frame allows the user to optimize the parameters for the local part of the pseudopotential
    """
    HELP_MSG = """\
This window allows you to change/optimize the parameters governing the local part"""

    def __init__(self, parent, wxoncv_input, **kwargs):
        """
        Args:
            wxoncv_input:
                Instance of WxOncvInput containing the parameters of the template.
        """
        super(LlocOptimizationFrame, self).__init__(parent, **kwargs)

        # Save reference to the input panel.
        self.wxoncv_input = wxoncv_input

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        #self.fcfact_ctl = awx.IntervalControl(self, start=0.25, num=6, step=0.05, choices=[">", "centered", "<"])
        #check_l.SetToolTipString("Enable/Disable optimization for this l-channel")
        #main_sizer.Add(self.fcfact_ctl, 1, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        buttons_sizer = self.make_buttons()
        main_sizer.Add(buttons_sizer, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        self.SetSizerAndFit(main_sizer)

    @property
    def opt_type(self):
        return "lloc optimization"

    def build_new_inps(self):
        #fcfact_list = self.fcfact_ctl.getValues()

        # Generate new list of inputs.
        base_inp = self.wxoncv_input.makeInput()
        return base_inp.optimize_vloc()


class RhoModelOptimizationFrame(OptimizationFrame):
    """
    This frame allows the user to optimize the model core charge 
    """

    HELP_MSG = """\
This window allows you to change/optimize the parameters governing the model core charge"""

    def __init__(self, parent, wxoncv_input, **kwargs):
        """
        Args:
            wxoncv_input:
                Instance of WxOncvInput containing the parameters of the template.
        """
        super(RhoModelOptimizationFrame, self).__init__(parent, **kwargs)

        # Save reference to the input panel.
        self.wxoncv_input = wxoncv_input

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        self.fcfact_ctl = awx.IntervalControl(self, start=0.25, num=6, step=0.05, choices=[">", "centered", "<"])
        #check_l.SetToolTipString("Enable/Disable optimization for this l-channel")
        main_sizer.Add(self.fcfact_ctl, 1, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        buttons_sizer = self.make_buttons()
        main_sizer.Add(buttons_sizer, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        self.SetSizerAndFit(main_sizer)

    @property
    def opt_type(self):
        return "RhoModel optimization"

    def build_new_inps(self):
        fcfact_list = self.fcfact_ctl.getValues()
        print(fcfact_list)

        # Generate new list of inputs.
        base_inp = self.wxoncv_input.makeInput()
        return base_inp.optimize_modelcore(fcfact_list, add_icmod0=True)


class QcutOptimizationFrame(OptimizationFrame):
    """
    This frame allows the user to select the l-channels and
    the list of values of qcut_l to be analyzed.
    """

    HELP_MSG = """\
This window allows you to change/optimize the value of the qcut parameters for 
the different angular channel. Use the checkboxes to select the l-channel(s) to be
analyzed, and the other controls to specify the list of qc values to test. 
"""

    def __init__(self, parent, wxoncv_input, **kwargs):
        """
        Args:
            wxoncv_input:
                Instance of WxOncvInput containing the parameters of the template.
        """
        super(QcutOptimizationFrame, self).__init__(parent, **kwargs)

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Save reference to the input panel.
        self.wxoncv_input = wxoncv_input
        lmax = wxoncv_input.lmax

        # Build list of controls to allow the user to select the list of
        # qcuts for the different l-channels. One can disable particular l-channels via checkboxes.
        self.checkbox_l = [None] * (lmax + 1)
        self.wxqcut_range_l = [None] * (lmax + 1)

        qcut_l = wxoncv_input.makeInput().qcut_l

        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        for l in range(lmax + 1):
            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Angular Channel L=%d" % l), wx.VERTICAL)
            vsz = wx.BoxSizer(wx.HORIZONTAL)

            self.checkbox_l[l] = check_l = wx.CheckBox(self, -1)
            check_l.SetToolTipString("Enable/Disable optimization for this l-channel")
            check_l.SetValue(True)

            self.wxqcut_range_l[l] = qcrange_l = awx.IntervalControl(self, start=qcut_l[l], num=4, step=0.1)

            vsz.Add(check_l, **add_opts)
            vsz.Add(qcrange_l, **add_opts)

            sbox_sizer.Add(vsz, 1, wx.ALL, 5)
            main_sizer.Add(sbox_sizer, 1, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        buttons_sizer = self.make_buttons()
        main_sizer.Add(buttons_sizer, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        self.SetSizerAndFit(main_sizer)

    @property
    def opt_type(self):
        return "Qcut optimization"

    def build_new_inps(self):
        l_list, qcvals_list = [], []
        for l, (checkbox, wxrange) in enumerate(zip(self.checkbox_l, self.wxqcut_range_l)):
            if checkbox.IsChecked():
                l_list.append(l)
                qcvals_list.append(wxrange.getValues())
        #print("\nl_list:", l_list, "\nqcvals_list", qcvals_list)

        # Generate new list of inputs.
        base_inp = self.wxoncv_input.makeInput()
        new_inps = []
        for l, new_qcuts in zip(l_list, qcvals_list):
            new_inps.extend(base_inp.optimize_qcuts_for_l(l=l, new_qcuts=new_qcuts))

        return new_inps


class RcOptimizationFrame(OptimizationFrame):
    """
    This frame allows the user to select the l-channels and
    the list of values of rc_l to be analyzed.
    """

    HELP_MSG = """\
This window allows you to change/optimize the value of the rc parameters (core radius)
for  the different angular channel."""

    def __init__(self, parent, wxoncv_input, **kwargs):
        """
        Args:
            wxoncv_input:
                Instance of WxOncvInput containing the parameters of the template.
        """
        super(RcOptimizationFrame, self).__init__(parent, **kwargs)

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Save reference to the input panel.
        self.wxoncv_input = wxoncv_input
        lmax = wxoncv_input.lmax

        # Build list of controls to allow the user to select the list of
        # qcut values for the different l-channels.
        # One can disable particular l-channels via checkboxes.
        self.checkbox_l = [None] * (lmax + 1)
        self.wxrc_range_l = [None] * (lmax + 1)

        rc_l = wxoncv_input.makeInput().rc_l

        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        for l in range(lmax + 1):
            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Angular Channel L=%d" % l), wx.VERTICAL)
            vsz = wx.BoxSizer(wx.HORIZONTAL)

            self.checkbox_l[l] = check_l = wx.CheckBox(self, -1)
            check_l.SetToolTipString("Enable/Disable optimization for this l-channel")
            check_l.SetValue(True)

            self.wxrc_range_l[l] = qcrange_l = awx.IntervalControl(self, start=rc_l[l], num=4, step=0.1)

            vsz.Add(check_l, **add_opts)
            vsz.Add(qcrange_l, **add_opts)

            sbox_sizer.Add(vsz, 1, wx.ALL, 5)
            main_sizer.Add(sbox_sizer, 1, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        buttons_sizer = self.make_buttons()
        main_sizer.Add(buttons_sizer, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        self.SetSizerAndFit(main_sizer)

    @property
    def opt_type(self):
        return "Rc optimization"

    def build_new_inps(self):
        l_list, rc_list = [], []
        for l, (checkbox, wxrange) in enumerate(zip(self.checkbox_l, self.wxrc_range_l)):
            if checkbox.IsChecked():
                l_list.append(l)
                rc_list.append(wxrange.getValues())

        # Generate new list of inputs.
        base_inp = self.wxoncv_input.makeInput()
        new_inps = []
        for l, new_rcs in zip(l_list, rc_list):
            new_inps.extend(base_inp.optimize_rcs_for_l(l, new_rcs))

        return new_inps


def empty_field(tag, oncv_dims=None):
    """Returns an empty field."""
    # Get the subclass from tag, initialize data with None and call __init__
    cls = _FIELD_LIST[tag]
    data = cls.make_empty_data(oncv_dims=oncv_dims)

    return cls(tag, data, oncv_dims)


class Field(object):
    # Flags used to define the type of field
    # Subclasses should define the class attribute type
    FTYPE_ROW = 1
    FTYPE_TABLE = 2
    FTYPE_RAGGED = 3

    # TODO: Change convention cbox --> slist
    parser_for_dtype = dict(
        i=int,
        f=float,
        cbox=str,
    )

    def __init__(self, tag, data, oncv_dims):
        """
        Args:

            tag:
                Tag used to identify the field.
            data:
                Accepts Ordered dict or List of ordered dicts with varname=value.

            oncv_dims:
                Dictionary with all the dimensions of the calculation.
        """
        #print("tag", tag, "data", data, "oncv_dims", oncv_dims)
        self.oncv_dims = AttrDict(**oncv_dims)
        self.tag = tag
        self.data = data

    @classmethod
    def make_empty_data(cls, oncv_dims=None):
        """Initialize data and fill it with None."""
        if cls.ftype == cls.FTYPE_ROW:
            data = OrderedDict([(k, None) for k in cls.WXCTRL_PARAMS.keys()])

        elif cls.ftype == cls.FTYPE_TABLE:
            nrows = cls.nrows_from_dims(oncv_dims)
            data = nrows * [None]
            for i in range(nrows):
                data[i] = OrderedDict([(k, None) for k in cls.WXCTRL_PARAMS.keys()])

        else:
            raise NotImplementedError(str(cls.ftype))

        return data

    @property
    def nrows(self):
        """
        Number of rows i.e. the number of lines in the section
        as specified in the input file.
        """
        if self.ftype == self.FTYPE_ROW:
            return 1
        elif self.ftype == self.FTYPE_TABLE:
            return len(self.data)
        else:
            raise NotImplementedError()

    #@property
    #def cols_per_row(self):
    #    if self.type == self.FTYPE_RAGGED:
    #    raise NotImplementedError()

    def __str__(self):
        """Returns a string with the input variables."""
        lines = []
        app = lines.append

        # Some fields have a single row but we need a list in the loop below.
        entries = self.data
        if self.ftype == self.FTYPE_ROW:
            entries = [self.data]

        for i, entry in enumerate(entries):
            if i == 0:
                # Put comment with the name of the variables only once.
                app("# " + " ".join((str(k) for k in entry.keys())))
            app(" ".join(str(v) for v in entry.values()))

        return "\n".join(lines)

    def has_var(self, key):
        """Return True if variable belongs to self."""
        if self.ftype == self.FTYPE_ROW:
            return key in self.data
        else:
            return key in self.data[0]

    def set_var(self, key, value):
        """Set the value of a variable."""
        assert self.has_var(key)
        if self.ftype == self.FTYPE_ROW:
            self.data[key] = value

        elif self.ftype == self.FTYPE_TABLE:
            for r in range(self.nrows):
                self.data[r][key] = value

        else:
            raise NotImplementedError()

    def set_vars(self, ord_vars):
        """
        Set the value of the variables inside a field

        Args:
            ord_vars:
                OrderedDict or list of OrderedDict (depending on field_idx).
        """
        assert len(ord_vars) == len(self.data)

        if self.ftype == self.FTYPE_ROW:
            #assert isinstance(ord_vars, OrderedDict)
            for k, v in ord_vars.items():
                self.data[k] = v

        elif self.ftype == self.FTYPE_TABLE:
            # List of ordered dicts.
            for i, od in enumerate(ord_vars):
                for k, v in od.items():
                    self.data[i][k] = v

        else:
            raise NotImplementedError()

    def set_vars_from_lines(self, lines):
        """The the value of the variables from a list of strings."""
        #print("About to read: ", type(self), "\nlines=\n, "\n".join(lines))

        okeys = self.WXCTRL_PARAMS.keys()
        odtypes = [v["dtype"] for v in self.WXCTRL_PARAMS.values()]
        parsers = [self.parser_for_dtype[ot] for ot in odtypes]
        #print("okeys", okeys, "odtypes", odtypes)

        if self.ftype == self.FTYPE_ROW:
            assert len(lines) == 1
            tokens = lines[0].split()
            #print("row tokens", tokens)

            for key, p, tok in zip(okeys, parsers, tokens):
                self.data[key] = p(tok)

        elif self.ftype == self.FTYPE_TABLE:
            assert len(lines) == self.nrows
            for i in range(self.nrows):
                tokens = lines[i].split()
                #print("table tokens: ", tokens)
                for key, p, tok in zip(okeys, parsers, tokens):
                    self.data[i][key] = p(tok)

        else:
            raise NotImplementedError()

    def get_vars(self):
        return self.data

    @classmethod
    def from_wxctrl(cls, wxctrl, tag, oncv_dims):
        """Build the object from a wxpython widget."""
        # Get the variables from the controller
        ord_vars = wxctrl.GetParams()
        # Build empty field and set its variables.
        new = empty_field(tag, oncv_dims)
        new.set_vars(ord_vars)

        return new

    def make_wxctrl(self, parent, **kwargs):
        """"Build the wx controller associated to this field."""
        if self.ftype == self.FTYPE_ROW:
            return awx.RowMultiCtrl(parent, self._customize_wxctrl(**kwargs))

        elif self.ftype == self.FTYPE_TABLE:
            return awx.TableMultiCtrl(parent, self.nrows, self._customize_wxctrl(**kwargs))

        else:
            # Ragged case e.g. test configurations:
            #dims
            raise NotImplementedError()

    def _customize_wxctrl(self, **kwargs):
        """
        Start with the default parameters for the wx controller
        and override them with those given in kwargs
        """
        # Make a deep copy since WXTRL_PARAMS is mutable.
        ctrl_params = copy.deepcopy(self.WXCTRL_PARAMS)

        for label, params in ctrl_params.items():
            value = kwargs.pop("label", None)
            if value is not None:
                params["value"] = str(value)

        return ctrl_params


class RowField(Field):
    """A field made of a single row."""
    ftype = Field.FTYPE_ROW


class TableField(Field):
    """A field made of multiple rows, all with the same number of columns."""
    ftype = Field.FTYPE_TABLE

    @classmethod
    def nrows_from_dims(cls, oncv_dims):
        """Return the number of rows from a dictionary with the dimensions."""
        raise NotImplementedError("Subclasses should define nrows_from_dims")

    def get_col(self, colname):
        """Returns an array with the values of column colname."""
        col = [None] * len(self.data)
        for i, row in enumerate(self.data):
            col[i] = row[colname]

        return col


class RaggedField(Field):
    """
    A field made of ragged rows, i.e. multiple rows with different number of columns.
    """
    ftype = Field.FTYPE_RAGGED

    @classmethod
    def nrows_from_dims(cls, oncv_dims):
        """Return the number of rows from a dictionary with the dimensions."""
        raise NotImplementedError("Subclasses should define nrows_from_dims")

    @classmethod
    def ncols_of_rows(cls, oncv_dims):
        """Return the number of columns in each row from a dictionary with the dimensions."""
        raise NotImplementedError("Subclasses should define nrows_from_dims")


def add_tooltips(cls):
    """Class decorator that add tooltips to WXCTRL_PARAMS."""
    d = cls.WXCTRL_PARAMS
    for key, params in d.items():
        params["tooltip"] = oncv_tip(key)

    return cls


@add_tooltips
class AtomConfField(RowField):
    name = "ATOMIC CONFIGURATION"

    WXCTRL_PARAMS = OrderedDict([
        ("atsym", dict(dtype="cbox", choices=periodic_table.all_symbols())),
        ("z", dict(dtype="i")),
        ("nc", dict(dtype="i", tooltip="number of core states"),),
        ("nv", dict(dtype="i", tooltip="number of valence states")),
        ("iexc", dict(dtype="i", value="4", default="4", tooltip="xc functional")),
        ("psfile", dict(dtype="cbox", choices=["psp8", "upf"]))])


@add_tooltips
class RefConfField(TableField):
    name = "REFERENCE CONFIGURATION"

    WXCTRL_PARAMS = OrderedDict([
        ("n", dict(dtype="i")),
        ("l", dict(dtype="i")),
        ("f", dict(dtype="f"))])

    #@classmethod
    #def neutral_from_symbol(cls, symbol):
    #    # TODO
    #    element = periodic_table.Element(symbol)
    #    # E.g., The electronic structure for Fe is represented as:
    #    # [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 6), (4, "s", 2)]
    #    #new = empty_field(cls.tag, oncv_dims=dict(nc))

    #    for row, (n, lchar, f) in zip(new, element.full_electronic_structure):
    #        row["n"], row["l"], row["f"] = n, periodic_table.char2l(lchar), f

    #    return new

    @classmethod
    def nrows_from_dims(cls, oncv_dims):
        return oncv_dims["nv"] + oncv_dims["nc"]

    @property
    def nrows(self):
        return self.oncv_dims["nv"] + self.oncv_dims["nc"]


@add_tooltips
class PseudoConfField(TableField):
    name = "PSEUDOPOTENTIAL AND OPTIMIZATION"

    WXCTRL_PARAMS = OrderedDict([
        ("l", dict(dtype="i")),
        ("rc", dict(dtype="f")),
        ("ep", dict(dtype="f")),
        ("ncon", dict(dtype="i", value="4")),
        ("nbas", dict(dtype="i", value="7")),
        ("qcut", dict(dtype="f"))])

    @classmethod
    def nrows_from_dims(cls, oncv_dims):
        return oncv_dims["lmax"] + 1

    @property
    def nrows(self):
        return self.oncv_dims["lmax"] + 1


@add_tooltips
class LmaxField(RowField):
    name = "LMAX"

    WXCTRL_PARAMS = OrderedDict([
        ("lmax", dict(dtype="i"))])


@add_tooltips
class VlocalField(RowField):
    name = "LOCAL POTENTIAL"

    WXCTRL_PARAMS = OrderedDict([
        ("lloc", dict(dtype="i", value="4")),
        ("lpopt", dict(dtype="i", value="5")),
        ("rc5", dict(dtype="f")),
        ("dvloc0", dict(dtype="f"))])


@add_tooltips
class VkbConfsField(TableField):
    name = "VANDERBILT-KLEINMAN-BYLANDER PROJECTORs"

    WXCTRL_PARAMS = OrderedDict([
        ("l", dict(dtype="i")),
        ("nproj", dict(dtype="i", value="2")),
        ("debl", dict(dtype="f"))])

    @classmethod
    def nrows_from_dims(cls, oncv_dims):
        return oncv_dims["lmax"] + 1

    @property
    def nrows(self):
        return self.oncv_dims["lmax"] + 1


@add_tooltips
class ModelCoreField(RowField):
    name = "MODEL CORE CHARGE"

    WXCTRL_PARAMS = OrderedDict([
        ("icmod", dict(dtype="i", value="1")),
        ("fcfact", dict(dtype="f", value="0.25"))])


@add_tooltips
class LogDerField(RowField):
    name = "LOG DERIVATIVE ANALYSIS"

    WXCTRL_PARAMS = OrderedDict([
        ("epsh1", dict(dtype="f", value="-2.0")),
        ("epsh2", dict(dtype="f", value="+2.0")),
        ("depsh", dict(dtype="f", value="0.02"))])


@add_tooltips
class RadGridField(RowField):
    name = "OUTPUT GRID"

    WXCTRL_PARAMS = OrderedDict([
        ("rlmax", dict(dtype="f", value="6.0")),
        ("drl", dict(dtype="f", value="0.01"))])


#@add_tooltips
#class TestConfigsField(RaggedField):
    #    name = "TEST CONFIGURATIONS"
    #    WXCTRL_PARAMS = OrderedDict([
    #        ("ncnf", dict(dtype="i", value="0")),
    #        ("nvcnf", dict(dtype="i", value="0")),
    #        ("n", dict(dtype="i")),
    #        ("l", dict(dtype="i")),
    #        ("f", dict(dtype="f"))])

    #@classmethod
    #def from_oxidation_states(cls, symbol, only_common=True):
    #    """
    #    Initialize the test configurations with the most common oxidation states.

    #    Args:
    #        symbol:
    #            Chemical symbol.:w
    #        only_common:
    #            If False all the known oxidations states are considered, else only
    #            the most common ones.
    #    """
    #    element = periodic_table.Element(symbol)

    #    if only_common:
    #        oxi_states = element.common_oxidation_states
    #    else:
    #        oxi_states = element.oxidation_states

    #    for oxi in oxi_states:
    #        # Get the electronic configuration of atom with Z = Z + oxi
    #        if oxi == 0:
    #            continue
    #        oxiele = periodic_table.Element.from_Z(element.Z + oxi)

    #        # Here we found the valence configuration by comparing
    #        # the full configuration of oxiele and the one of the initial element.


    #    return new

    #@property
    #def nrows(self):
    #    return self.oncv_dims["ncnf"]

    #@classmethod
    #def nrows_from_dims(cls, oncv_dims):
    #    return oncv_dims["ncnf"]

    #@property
    #def nlines_for_row(self, row):


# List with the field in the same order as the one used in the input file.
_FIELD_LIST = [
    AtomConfField,
    RefConfField,
    LmaxField,
    PseudoConfField,
    VlocalField,
    VkbConfsField,
    ModelCoreField,
    LogDerField,
    RadGridField,
    #TestConfigsField,
]

_NFIELDS = len(_FIELD_LIST)


class OncvInput(object):
    """
    This object stores the variables needed for generating a pseudo with oncvsps.
    One can initialize this object either from a prexisting file
    or programmatically from the input provided by the user in a GUI.

    An input consistst of _NFIELDS fields. Each field is either a OrderedDict
    or a list of ordered dicts with the input variables.
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from an external file."""
        # Read input lines: ignore empty lines or line starting with #
        lines = []
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    lines.append(line)

        # Read dimensions
        # nc and nv from the first line.
        tokens = lines[0].split()
        atsym, z, nc, nv = tokens[0:4]
        z, nc, nv = map(int, (z, nc, nv))

        # Read lmax and ncfn
        lmax = int(lines[nc + nv + 1])

        # TODO
        # number of tests and number of rows for the different configurations.
        ncnf = 0

        # Initialize the object
        inp = OncvInput(oncv_dims=dict(atsym=atsym, nc=nc, nv=nv, lmax=lmax, ncnf=ncnf))

        # TODO
        # Fill it
        start = 0
        for field in inp:
            stop = start + field.nrows
            field.set_vars_from_lines(lines[start:stop])
            start = stop

        return inp

    def __init__(self, oncv_dims, fields=None):
        """
        Initialize the object from a dict with the fundamental dimensions.
        If fields is None, we create an empty dict, else we use fields.
        """
        self.dims = AttrDict(**oncv_dims)
        #print("oncv_dims", self.dims)

        if fields is None:
            # Default fields.
            self.fields = _NFIELDS * [None]
            for i in range(_NFIELDS):
                new = empty_field(i, self.dims)
                self.fields[i] = new

        else:
            self.fields = fields

        # 1) ATOM CONFIGURATION
        # atsym, z, nc, nv, iexc   psfile
        # O    8     1   2   3   psp8
        #
        # 2) REFERENCE CONFIGURATION
        # n, l, f  (nc+nv lines)
        # 1    0    2.0
        # 2    0    2.0
        # 2    1    4.0
        #
        # 3) LMAX FIELD
        # lmax
        # 1
        # 4) PSEUDOPOTENTIAL AND OPTIMIZATION
        # l, rc, ep, ncon, nbas, qcut  (lmax+1 lines, l's must be in order)
        # 0    1.60    0.00    4    7    8.00
        # 1    1.60    0.00    4    7    8.00
        #
        # 5) LOCAL POTENTIAL
        # lloc, lpopt, rc(5), dvloc0
        # 4    5    1.4    0.0
        #
        # 6) VANDERBILT-KLEINMAN-BYLANDER PROJECTORs
        # l, nproj, debl  (lmax+1 lines, l's in order)
        # 0    2    1.50
        # 1    2    1.00
        #
        # 7) MODEL CORE CHARGE
        # icmod, fcfact
        # 0    0.0
        #
        # 8) LOG DERIVATIVE ANALYSIS
        # epsh1, epsh2, depsh
        # -2.0  2.0  0.02
        #
        # 9) OUTPUT GRID
        # rlmax, drl
        # 4.0  0.01

        # TODO
        # 10) TEST CONFIGURATIONS
        # ncnf
        # 2
        # nvcnf    (repeated ncnf times)
        # n, l, f  (nvcnf lines, repeated follwing nvcnf's ncnf times)
        # 2
        # 2    0    2.0
        # 2    1    3.0
        #
        # 2
        # 2    0    1.0
        # 2    1    4.0
        #ncnf = 2
        #nvcnf = 2

    @property
    def lmax(self):
        return self.dims.lmax

    def __iter__(self):
        return self.fields.__iter__()

    def __str__(self):
        """Returns a string with the input variables."""
        lines = []
        app = lines.append

        for i, field in enumerate(self):
            # FIXME This breaks the output parser!!!!!
            #pre = "\n#" if i > 0 else ""
            #app(pre + field.name)
            app(str(field))

        s = "\n".join(lines)
        # FIXME needed to bypass problems with tests
        return s + "\n 0"

    def __setitem__(self, key, value):
        ncount = 0
        for f in self.fields:
            if f.has_var(key):
                ncount += 1
                f.set_var(key, value)

        assert ncount == 1

    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    def optimize_vloc(self):
        """Produce a list of new input files by changing the lloc option for vloc."""
        # Test all possible vloc up to lmax
        inps, new = [], self.deepcopy()
        for il in range(self.lmax+1):
            new["lloc"] = il
            inps.append(new.deepcopy())

        # Add option for smooth polynomial
        new["lloc"] = 4
        inps.append(new)

        return inps

    def optimize_modelcore(self, fcfact_list, add_icmod0=True):
        """Produce a list of new input files by changing the icmod option for model core."""
        inps, new = [], self.deepcopy()

        if add_icmod0:
            new["icmod"] = 0
            inps.append(new.deepcopy())
            
        for fcfact in fcfact_list:
            new["icmod"] = 1
            new["fcfact"] = fcfact
            inps.append(new.deepcopy())

        return inps

    @property
    def qcut_l(self):
        """List with the values of qcuts as function of l."""
        i = _FIELD_LIST.index(PseudoConfField)
        return self.fields[i].get_col("qcut")


    @property
    def rc_l(self):
        """List with the values of rc as function of l."""
        i = _FIELD_LIST.index(PseudoConfField)
        return self.fields[i].get_col("rc")

    def optimize_qcuts_for_l(self, l, new_qcuts):
        """
        Returns a list of new input objects in which the qcut parameter for
        the given l has been replaced by the values listed in new_qcuts.

        Args:
            l:
                Angular momentum
            new_qcuts:
                Iterable with the new values of qcut.
                The returned list will have len(new_qcuts) input objects.
        """
        # Find the field with the configuration parameters.
        i = _FIELD_LIST.index(PseudoConfField)

        # Find the row with the given l.
        for irow, row in enumerate(self.fields[i].data):
            if row["l"] == l:
                break
        else:
            raise ValueError("Cannot find l %s in the PseudoConfField" % l)

        # This is the dict we want to change
        inps = []
        for qc in new_qcuts:
            new_inp = self.deepcopy()
            new_inp.fields[i].data[irow]["qcut"] = qc
            inps.append(new_inp)

        return inps

    def optimize_rcs_for_l(self, l, new_rcs):
        """
        Returns a list of new input objects in which the rc parameter for
        the given l has been replaced by the values listed in new_rcs.

        Args:
            l:
                Angular momentum
            new_rcs:
                Iterable with the new values of rcs.
                The returned list will have len(new_rcs) input objects.
        """
        # Find the field with the configuration parameters.
        i = _FIELD_LIST.index(PseudoConfField)

        # Find the row with the given l.
        for irow, row in enumerate(self.fields[i].data):
            if row["l"] == l:
                break
        else:
            raise ValueError("Cannot find l %s in the PseudoConfField" % l)

        # This is the dict we want to change
        inps = []
        for rc in new_rcs:
            new_inp = self.deepcopy()
            new_inp.fields[i].data[irow]["rc"] = rc
            inps.append(new_inp)

        return inps


# TODO Try this
#import wx.lib.foldpanelbar as foldpanel
#class WxOncvInput(foldpanel.FoldPanelBar):
from wx.lib.scrolledpanel import ScrolledPanel

#class WxOncvInput(awx.Panel):
class WxOncvInput(ScrolledPanel):
    """
    Panel with widgets allowing the user to select the input parameters.
    """
    def __init__(self, parent, oncv_dims, **kwargs):
        """
        Args:
            oncv_dims:
                Basic dimensions of the calculation.
        """
        #super(WxOncvInput, self).__init__(parent, id=-1)
        super(WxOncvInput, self).__init__(parent, id=-1, style=wx.VSCROLL)

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        stext = wx.StaticText(self, -1, "Calculation type:")
        choices = ["scalar-relativistic", "fully-relativistic", "non-relativistic"]
        self.calctype_cbox = wx.ComboBox(
            self, id=-1, name='Calculation type', choices=choices, value=choices[0], style=wx.CB_READONLY)

        add_opts = dict(proportion=0, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        hbox0.Add(stext, **add_opts)
        hbox0.Add(self.calctype_cbox)

        self.main_sizer.Add(hbox0, **add_opts)

        # Set the dimensions and build the widgets.
        self.oncv_dims = oncv_dims
        sz, sizer_addopts = wx.BoxSizer(wx.VERTICAL), dict(proportion=0, flag=wx.ALL, border=5)

        # We have nfields sections in the input file.
        # Each field has a widget that returns the variables in a dictionary
        self.wxctrls = _NFIELDS * [None]
        self.sbox_sizers = _NFIELDS * [None]

        # Keep an internal list of buttons so that we can disable them easily.
        self.all_optimize_buttons = []

        # FieldClass: [(label, OptimizationFrame), ....]
        self.fields_with_optimization = {
            PseudoConfField: [
                ("Optimize rc", RcOptimizationFrame),
                ("Optimize qcut", QcutOptimizationFrame),
                ],
            VlocalField: [("Optimize lloc", LlocOptimizationFrame)],
            ModelCoreField: [("Optimize fcfact", RhoModelOptimizationFrame)],
            #VkbConfsField: [("Optimize pseudo", [])],
        }

        for i in range(_NFIELDS):
            f = empty_field(i, oncv_dims)
            wxctrl = f.make_wxctrl(self)
            self.wxctrls[i] = wxctrl

            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, f.name + ":"), wx.VERTICAL)
            sbox_sizer.Add(wxctrl, **sizer_addopts)

            # add optimization button if the field supports it.
            if f.__class__ in self.fields_with_optimization:
                hsz = wx.BoxSizer(wx.HORIZONTAL)
                for label, opt_frame in self.fields_with_optimization[f.__class__]:
                    optimize_button = wx.Button(self, -1, label)
                    optimize_button.Bind(wx.EVT_BUTTON, self.onOptimize)
                    optimize_button.field_class = f.__class__
                    optimize_button.opt_frame = opt_frame

                    self.all_optimize_buttons.append(optimize_button)
                    hsz.Add(optimize_button, 0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=5)

                #sbox_sizer.Add(optimize_button, 0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)
                sbox_sizer.Add(hsz, 0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

            sz.Add(sbox_sizer, **sizer_addopts)
            self.sbox_sizers[i] = sbox_sizer

        self.main_sizer.Add(sz, **add_opts)

        self.fill_from_file("08_O.dat")

        self.SetupScrolling()
        self.SetSizerAndFit(self.main_sizer)

    def enable_all_optimize_buttons(self, enable=True):
        """Enable/disable the optimization buttons."""
        for button in self.all_optimize_buttons:
            button.Enable(enable)

    def onOptimize(self, event):
        button = event.GetEventObject()
        self.enable_all_optimize_buttons(False)
        opt_frame = button.opt_frame
        try:
            opt_frame(self, wxoncv_input=self).Show()
        finally:
            self.enable_all_optimize_buttons(True)

    @property
    def lmax(self):
        return self.oncv_dims["lmax"]

    def fill_from_file(self, filename):
        """Build a panel from an input file."""
        inp = OncvInput.from_file(filename)
        for field, wxctrl in zip(inp, self.wxctrls):
            wxctrl.SetParams(field.data)

    @property
    def calc_type(self):
        """"Return a string with the calculation type."""
        return self.calctype_cbox.GetValue()

    def makeInput(self):
        """Build an instance of OncvInput from the data specified in the Wx controllers."""
        inp = OncvInput(self.oncv_dims)
        for i, field in enumerate(self.wxctrls):
            inp.fields[i].set_vars(field.GetParams())

        return inp

    def makeInputString(self):
        """Return a string with the input passed to the pp generator."""
        return str(self.makeInput())

# Event posted when we start an optimization.
#EVT_OPTIMIZATION_TYPE = wx.NewEventType()
#EVT_OPTIMIZATION = wx.PyEventBinder(EVT_CONSOLE_TYPE, 1)
#
#class OptimizationEvent(wx.PyEvent):
#    """
#    This event is triggered when we start/end the optimization process
#    """
#    def __init__(self, kind, msg):
#        wx.PyEvent.__init__(self)
#        self.SetEventType(EVT_OPTIMIZATION_TYPE)
#        self.kind, self.msg = kind, msg
#
#    #@classmethod
#    #def start_optimization(cls, msg)
#    #@classmethod
#    #def end_optimization(cls, msg)
#    #wx.PostEvent(self.console, event)



class PseudoGeneratorListCtrl(wx.ListCtrl, listmix.ColumnSorterMixin, listmix.ListCtrlAutoWidthMixin):
    """
    ListCtrl that allows the user to interact with a list of pseudogenerators. Supports column sorting 
    """
    # List of columns
    _COLUMNS = ["#", 'status', "max_ecut", "atan_logder_err", "max_psexc_abserr", "herm_err"]

    def __init__(self, parent, psgens=(), **kwargs):
        """
        Args:
            parent:
                Parent window.
            psgens:
                List of `PseudoGenerator` instances.
        """
        super(PseudoGeneratorListCtrl, self).__init__(
            parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)

        self.psgens = psgens if psgens else []

        for index, col in enumerate(self._COLUMNS):
            self.InsertColumn(index, col)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in self._COLUMNS]

        # Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
        self.itemDataMap = {}

        for index, psgen in enumerate(self.psgens):
            entry = self.make_entry(index, psgen)
            self.Append(entry)
            self.SetItemData(index, index)
            self.itemDataMap[index] = entry

            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

        for index, col in enumerate(self._COLUMNS):
            self.SetColumnWidth(index, column_widths[index])

        # Now that the list exists we can init the other base class, see wx/lib/mixins/listctrl.py
        listmix.ColumnSorterMixin.__init__(self, len(self._COLUMNS))
        listmix.ListCtrlAutoWidthMixin.__init__(self)

        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.onItemActivated) 
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.onRightClick)

    @property
    def wxinput(self):
        """Reference to the Wx window with the input file."""
        try:
            return self._wxinput
        except AttributeError:
            parent = self.GetParent()
            cls = WxOncvFrame

            while True:
                if parent is None:
                    raise RuntimeError("Cannot find parent with class %s, reached None parent!" % cls)

                if isinstance(parent, cls):
                    break
                else:
                    parent = parent.GetParent()

            # The input we need is an attribute of WxOncvFrame.
            self._wxinput = parent.wxoncv_input
            return self._wxinput

    @staticmethod
    def make_entry(index, psgen):
        """Returns the entry associated to the generator psgen with the given index."""
        max_ecut = None
        max_atan_logder_l1err = None
        max_psexc_abserr = None
        herm_err = None

        if psgen.results is not None:
            max_ecut = psgen.results.max_ecut
            max_atan_logder_l1err = psgen.results.max_atan_logder_l1err
            max_psexc_abserr = psgen.results.max_psexc_abserr
            herm_err = psgen.results.herm_err

        return [
            "%d\t\t" % index,
            "%s" % psgen.status,
            "%s" % max_ecut,
            "%s" % max_atan_logder_l1err,
            "%s" % max_psexc_abserr,
            "%s" % herm_err,
        ]

    def doRefresh(self):
        """Refresh the panel and redraw it."""
        column_widths = [awx.get_width_height(self, s)[0] for s in self._COLUMNS]

        for index, psgen in enumerate(self.psgens):
            entry = self.make_entry(index, psgen)
            self.SetItemData(index, index)
            self.itemDataMap[index] = entry

            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

        for index, col in enumerate(self._COLUMNS):
            self.SetColumnWidth(index, column_widths[index])

    def add_psgen(self, psgen):
        """Add a PseudoGenerator to the list."""
        index = len(self.psgens)
        entry = self.make_entry(index, psgen)
        self.Append(entry)
        self.SetItemData(index, index)
        self.itemDataMap[index] = entry

        # Add it to the list and update column widths.
        self.psgens.append(psgen)
        self.doRefresh()

    def GetListCtrl(self):
        """Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py"""
        return self

    def getSelectedPseudoGen(self):
        """
        Returns the PseudoGenerators selected by the user.
        None if no selection has been done.
        """
        # Get selected index, map to index in psgens and return the object.
        item = self.GetFirstSelected()
        if item == -1: return None
        index = self.GetItemData(item)
        return self.psgens[index]

    def onItemActivated(self, event):
        """Call psgen.plot_results."""
        psgen = self.getSelectedPseudoGen()
        if psgen is None: return
        psgen.plot_results()

    def onPlotSubMenu(self, event):
        """Called by plot submenu."""
        psgen = self.getSelectedPseudoGen()
        if psgen is None or psgen.plotter is None:
            return

        key = self._id2plotdata[event.GetId()]
        psgen.plotter.plot_key(key)

    def onRightClick(self, event):
        """Generate the popup menu."""
        popup_menu = self.makePopupMenu()
        self.PopupMenu(popup_menu, event.GetPoint())
        popup_menu.Destroy()

    def makePopupMenu(self):
        """
        Build and return a popup menu. Subclasses can extend or replace this base method.
        """
        self.ID_POPUP_STDIN = wx.NewId()
        self.ID_POPUP_STDOUT = wx.NewId()
        self.ID_POPUP_STDERR = wx.NewId()
        self.ID_POPUP_CHANGE_INPUT = wx.NewId()

        menu = wx.Menu()

        # Make sub-menu with the list of supported quantities
        plot_submenu = wx.Menu()

        # TODO: this list could be taken from the class or from the plotter instance.
        all_keys = [
            "radial_wfs",
            "projectors",
            "densities",
            "potentials",
            "atan_logders",
            "ene_vs_ecut",
        ]

        self._id2plotdata = {}

        for aqt in all_keys:
            _id = wx.NewId()
            plot_submenu.Append(_id, aqt)
            self._id2plotdata[_id] = aqt
            self.Bind(wx.EVT_MENU, self.onPlotSubMenu, id=_id)
        menu.AppendMenu(-1, 'Plot', plot_submenu)

        menu.Append(self.ID_POPUP_STDIN, "Show standard input")
        menu.Append(self.ID_POPUP_STDOUT, "Show standard output")
        menu.Append(self.ID_POPUP_STDERR, "Show standard error")
        menu.Append(self.ID_POPUP_CHANGE_INPUT, "Use these variables as new template")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_POPUP_STDIN, self.onShowStdin),
            (self.ID_POPUP_STDOUT, self.onShowStdout),
            (self.ID_POPUP_STDERR, self.onShowStderr),
            (self.ID_POPUP_CHANGE_INPUT, self.onChangeInput),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu

    def onChangeInput(self, event):
        """Change the template input file."""
        psgen = self.getSelectedPseudoGen()
        if psgen is None: return

        # Change the parameters in wxinput using those in psgen.stdin_path
        if not os.path.exists(psgen.stdin_path):
            awx.showErrorMessage(self, "Input file %s does not exist" % psgen.stdin_path)

        self.wxinput.fill_from_file(psgen.stdin_path)

    #def onAddToHistory(self, event):
    #    psgen = self.getSelectedPseudoGen()
    #    if psgen is None: return
    #    # Add it to history.
    #    self.history.append(psgen)

    def _showStdfile(self, event, stdfile):
        """
        Helper function used to access the std files
        of the PseudoGenerator selected by the user.
        """
        psgen = self.getSelectedPseudoGen()
        if psgen is None: return
        call = dict(
            stdin=psgen.get_stdin,
            stdout=psgen.get_stdout,
            stderr=psgen.get_stderr,
        )[stdfile]

        SimpleTextViewer(self, text=call()).Show()

    def onShowStdin(self, event):
        """Open a frame with the input file."""
        self._showStdfile(event, "stdin")

    def onShowStdout(self, event):
        """Open a frame with the output file."""
        self._showStdfile(event, "stdout")

    def onShowStderr(self, event):
        """Open a frame with the stderr file."""
        self._showStdfile(event, "stderr")

    def plot_columns(self, **kwargs):
        """Use matplotlib to plot the values reported in the columns."""
        keys = self._COLUMNS[2:]
        table = OrderedDict([(k, []) for k in keys])

        for index, psgen in enumerate(self.psgens):
            entry = self.make_entry(index, psgen)
            for k, v in zip(keys, entry[2:]):
                #print("k", k, "v", v)
                try:
                    v = float(v)
                except ValueError:
                    # This happens if the run is not completed.
                    v = None

                table[k].append(v)

        # Return immediately if all entries are None i.e.
        # if all calculations are still running.
        count = 0
        for items in table.values():
            if all(item is None for item in items): count += 1
        if count == len(table): return

        import matplotlib.pyplot as plt

        # Build grid of plots.
        fig, ax_list = plt.subplots(nrows=len(table), ncols=1, sharex=False, squeeze=True)

        for key, ax in zip(table.keys(), ax_list):
            ax.grid(True)
            ax.set_title(key)
            row = table[key]
            xs = np.arange(len(row))
            ys = np.array(table[key]).astype(np.double)
            # Use mask to exclude None values from the plot.
            mask = np.isfinite(ys)
            line, = ax.plot(xs[mask], ys[mask], linewidth=2, markersize=10, marker="o")

        plt.show()


class PseudoGeneratorsPanel(awx.Panel):
    """
    A panel with a list control providing info on the status pseudogenerators
    Provides popup menus for interacting with the generators.
    """
    def __init__(self, parent, psgens=(), **kwargs):
        """
        Args:
            parent:
                Parent window.
            psgens:
                List of `PseudoGenerator` objects.
        """
        super(PseudoGeneratorsPanel, self).__init__(parent, **kwargs)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.psgen_list_ctrl = PseudoGeneratorListCtrl(self, psgens)

        main_sizer.Add(self.psgen_list_ctrl, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 5)
        self.SetSizerAndFit(main_sizer)

    @property
    def psgens(self):
        return self.psgen_list_ctrl.psgens

    def plot_columns(self, **kwargs):
        self.psgen_list_ctrl.plot_columns(**kwargs)


class PseudoGeneratorsFrame(awx.Frame):
    """
    This frame contains a list of pseudopotential generators,
    It provides controls to run the calculation, and inspect/plot the results.
    """
    REFRESH_INTERVAL = 120

    HELP_MSG = """\
This window allows you to generate and analyze multiple pseudopotentials.
"""

    def __init__(self, parent, psgens=(), **kwargs):
        """
        Args:
            wxoncv_input:
                Instance of WxOncvInput
            psgens:
                List of `PseudoGenerators`.
        """
        super(PseudoGeneratorsFrame, self).__init__(parent, -1, **add_size(kwargs))

        # Build menu, toolbar and status bar.
        self.makeToolBar()
        self.statusbar = self.CreateStatusBar()
        self.Centre()

        self.panel = panel = wx.Panel(self, -1)
        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.psgens_wxlist = PseudoGeneratorsPanel(panel, psgens)
        main_sizer.Add(self.psgens_wxlist, 1, wx.ALL | wx.EXPAND, 5)

        submit_button = wx.Button(panel, -1, label='Submit')
        submit_button.Bind(wx.EVT_BUTTON, self.OnSubmitButton)

        text = wx.StaticText(panel, -1, "Max nlaunch:")
        text.Wrap(-1)
        text.SetToolTipString("Maximum number of tasks that can be submitted. Use -1 for unlimited launches.")
        self.max_nlaunch = wx.SpinCtrl(panel, -1, value=str(get_ncpus()), min=-1)

        help_button = wx.Button(panel, wx.ID_HELP)
        help_button.Bind(wx.EVT_BUTTON, self.onHelp)
        main_sizer.Add(help_button, 0, flag=wx.ALL | wx.ALIGN_RIGHT)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(submit_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        hsizer.Add(text, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        hsizer.Add(self.max_nlaunch, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        hsizer.Add(help_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)

        main_sizer.Add(hsizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        panel.SetSizerAndFit(main_sizer)

        # Register this event when the GUI is IDLE
        self.last_refresh = time.time()
        self.Bind(wx.EVT_IDLE, self.OnIdle)

    @property
    def psgens(self):
        """List of PseudoGenerators."""
        return self.psgens_wxlist.psgens

    def makeToolBar(self):
        """Create toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        self.ID_SHOW_INPUTS = wx.NewId()
        self.ID_SHOW_OUTPUTS = wx.NewId()
        self.ID_SHOW_ERRS = wx.NewId()
        self.ID_CHECK_STATUS = wx.NewId()
        self.ID_MULTI_PLOTTER = wx.NewId()
        self.ID_PLOT_COLUMNS = wx.NewId()

        toolbar.AddSimpleTool(self.ID_SHOW_INPUTS, bitmap("in.png"), "Visualize the input file(s) of the generators.")
        toolbar.AddSimpleTool(self.ID_SHOW_OUTPUTS, bitmap("out.png"), "Visualize the output file(s) of the generators.")
        toolbar.AddSimpleTool(self.ID_SHOW_ERRS, bitmap("log.png"), "Visualize the errors file(s) of the generators.")
        toolbar.AddSimpleTool(self.ID_MULTI_PLOTTER, bitmap("log.png"), "Multi plotter.")
        toolbar.AddSimpleTool(self.ID_PLOT_COLUMNS, bitmap("log.png"), "Plot columns.")
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(self.ID_CHECK_STATUS, bitmap("refresh.png"), "Check the status of the workflow(s).")

        toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_SHOW_INPUTS, self.OnShowInputs),
            (self.ID_SHOW_OUTPUTS, self.OnShowOutputs),
            (self.ID_SHOW_ERRS, self.OnShowErrors),
            (self.ID_CHECK_STATUS, self.OnCheckStatus),
            (self.ID_MULTI_PLOTTER, self.OnMultiPlotter),
            (self.ID_PLOT_COLUMNS, self.onPlotColumns),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

    def OnSubmitButton(self, event):
        """
        Called when Run button is pressed.
        Run the calculation in a subprocess in non-blocking mode and add it to
        the list containing the generators in executions
        Submit up to max_nlauch tasks (-1 to run'em all)
        """
        max_nlaunch = int(self.max_nlaunch.GetValue())

        nlaunch = 0
        for psgen in self.psgens:
            nlaunch += psgen.start()
            if nlaunch == max_nlaunch:
                break

        self.statusbar.PushStatusText("Submitted %d tasks" % nlaunch)

    def OnCheckStatus(self, event):
        """
        Callback triggered by the checkstatus button.
        Check the status of the `PseudoGenerators` and refresh the panel.
        """
        self.CheckStatusAndRedraw()

    def CheckStatusAndRedraw(self):
        """Check the status of all the workflows and redraw the panel."""
        self.statusbar.PushStatusText("Checking status...")
        start = time.time()

        for psgen in self.psgens:
            psgen.check_status()
        self.statusbar.PushStatusText("Check completed in %.1f [s]" % (time.time() - start))

        # Redraw the panel
        main_sizer = self.main_sizer
        main_sizer.Hide(0)
        main_sizer.Remove(0)
        new_psgen_wxlist = PseudoGeneratorsPanel(self.panel, self.psgens)
        main_sizer.Insert(0, new_psgen_wxlist, 1, wx.EXPAND, 5)
        self.psgen_wxlist = new_psgen_wxlist

        self.panel.Layout()

        # Write number of jobs with given status.
        #message = ", ".join("%s: %s" % (k, v) for (k, v) in counter.items())
        #self.statusbar.PushStatusText(message)

    def OnIdle(self, event):
        """Function executed when the GUI is idle."""
        now = time.time()
        if (now - self.last_refresh) > self.REFRESH_INTERVAL:
            self.CheckStatusAndRedraw()
            self.last_refresh = time.time()

    def OnShowInputs(self, event):
        """Show all input files."""
        TextNotebookFrame(self, text_list=[psgen.get_stdin() for psgen in self.psgens], 
                          page_names=["PSGEN # %d" % i for i in range(len(self.psgens))]).Show()

    def OnShowOutputs(self, event):
        """Show all output files."""
        TextNotebookFrame(self, text_list=[psgen.get_stdout() for psgen in self.psgens],
                          page_names=["PSGEN # %d" % i for i in range(len(self.psgens))]).Show()

    def OnShowErrors(self, event):
        """Show all error files."""
        TextNotebookFrame(self, text_list=[psgen.get_stderr() for psgen in self.psgens], 
                          page_names=["PSGEN # %d" % i for i in range(len(self.psgens))]).Show()

    def onPlotColumns(self, event):
        self.psgens_wxlist.plot_columns()

    def OnMultiPlotter(self, event):
        """Open a dialog that allows the user to plot the results of multiple generators."""
        multi_plotter = MultiPseudoGenDataPlotter()

        # Add psgen only if run is OK.
        for i, psgen in enumerate(self.psgens):
            if psgen.status == psgen.S_OK:
                multi_plotter.add_psgen(label="%d" % i, psgen=psgen)

        # Return immediately if no calculation is OK.
        if not len(multi_plotter):
            return

        keys = list(multi_plotter.keys())

        class MyFrame(awx.FrameWithChoice):
            """Get a string with the quantity to plot and call multi_plotter.plot_key"""
            def onOkButton(self, event):
                multi_plotter.plot_key(key=self.getChoice())

        MyFrame(self, choices=keys, title="MultiPlotter").Show()


class OncvApp(awx.App):

    def OnInit(self):
        # The code for the splash screen.
        #image = wx.Image(path_img("wabi_logo.png"), wx.BITMAP_TYPE_PNG)
        #    bmp = image.ConvertToBitmap()
        #    wx.SplashScreen(bmp, wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT, 1000, None, -1)
        #    wx.Yield()

        frame = WxOncvFrame(None)
        #frame = my_periodic_table(None)
        #frame = OncvParamsFrame(None, z=12)
        frame.Show(True)
        self.SetTopWindow(frame)
        return True


def wxapp_oncvpsp():
    """Standalone application."""
    app = OncvApp()
    return app

if __name__ == "__main__":
    import sys
    onc_inp = OncvInput.from_file("08_O.dat")
    #print(onc_inp)

    #print("optimizing qcut")
    #for inp in onc_inp.optimize_qcuts_for_l(l=0, new_qcuts=[1, 9]):
    #    a = 1
    #    #print("new model\n", inp)

    #for inp in onc_inp.optimize_vloc():
    #    print("new\n", inp)
    #for inp in onc_inp.optimize_modelcore():
    #    print("new model\n", inp)
    #sys.exit(0)

    wxapp_oncvpsp().MainLoop()
    #app = awx.App()
    #frame = QcutOptimizationFrame(None, lmax=1)
    #frame.Show()
    #app.MainLoop()


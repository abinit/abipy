#!/usr/bin/env python
"""Gui for the oncvpsp norm-conserving pseudopotential generator."""
import os
import copy
import time
import shutil
import abc
import sys
import wx
import wx.lib.mixins.listctrl as listmix
import numpy as np

from collections import OrderedDict
from monty.dev import get_ncpus
from monty.collections import AttrDict
from pymatgen.core.periodic_table import Element
from abipy.gui.editor import TextNotebookFrame, SimpleTextViewer
from abipy.gui.oncvtooltips import oncv_tip
from abipy.gui import mixins as mix
from abipy.gui import awx
from abipy.gui.awx.elements_gui import WxPeriodicTable, PeriodicPanel, ElementButton
from abipy.flowtk import Pseudo

try:
    from pseudo_dojo.core.dojoreport import DojoReport
    from pseudo_dojo.refdata.nist import database as nist
    from pseudo_dojo.ppcodes.ppgen import OncvGenerator
    from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, MultiPseudoGenDataPlotter
except ImportError as exc:
    print("Error while trying to import pseudo_dojo modules:\n%s" % str(exc))
    #raise

# TODO
# Change oncvpsp so that
#   1) we always write the logarithmic derivative
#   2) better error handling


_char2l = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
    "g": 4,
    "h": 5,
    "i": 6,
}

def char2l(char):
    return _char2l[char]


def all_symbols():
    return [e.symbol for e in Element]


def add_size(kwargs, size=(800, 600)):
    """Add size to kwargs if not present."""
    if "size" not in kwargs:
        kwargs["size"] = size

    return kwargs



def my_periodic_table(parent):
    """
    A periodic table that allows the user to select the element
    before starting the pseudopotential generation.
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
    to generate the pseudopotential once the chemical element has been selected.
    """

    HELP_MSG = """\
Quick help:
    Use this window to select the AE reference configuration and how
    to separate states into core and valence.
"""

    def __init__(self, parent, z, **kwargs):
        super(OncvParamsFrame, self).__init__(parent, **kwargs)
        self.element = Element.from_Z(z)
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
            row["n"], row["l"], row["f"] = n, char2l(lchar), f
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

    def __init__(self, parent, filepath=None):
        super(WxOncvFrame, self).__init__(parent, id=-1, title=self.codename)

        # This combination of options for config seems to work on my Mac.
        self.config = wx.FileConfig(appName=self.codename, localFilename=self.codename + ".ini",
                                    style=wx.CONFIG_USE_LOCAL_FILE)

        # Build menu, toolbar and status bar.
        self.SetMenuBar(self.makeMenu())
        self.statusbar = self.CreateStatusBar()
        self.Centre()
        self.makeToolBar()
        #self.toolbar.Enable(False)

        self.input_file = None
        if filepath is not None:
            if os.path.exists(filepath):
                self.input_file = os.path.abspath(filepath)
                self.BuildUI(notebook=OncvNotebook.from_file(self, filepath))
            else:
                # Assume symbol
                self.BuildUI(notebook=OncvNotebook.from_symbol(self, filepath))
        else:
            self.BuildUI()

    @property
    def codename(self):
        """Name of the application."""
        return "WxOncvGui"

    def BuildUI(self, notebook=None):
        """Build user-interface."""
        old_selection = None
        if hasattr(self, "main_sizer"):
            # Remove old notebook
            main_sizer = self.main_sizer
            main_sizer.Hide(0)
            main_sizer.Remove(0)

            # Save the active tab so that we can set it afterwards.
            old_selection = self.notebook.GetSelection()
            del self.notebook
        else:
            self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        if notebook is not None:
            self.notebook = notebook
        else:
            # Default is, of course, silicon.
            self.notebook = OncvNotebook.from_symbol(self, "Si")

        main_sizer.Add(self.notebook, flag=wx.EXPAND)
        self.SetSizerAndFit(main_sizer)
        main_sizer.Layout()

        # Reinstate the old selection
        if old_selection is not None:
            self.notebook.SetSelection(old_selection)

    def AddFileToHistory(self, filepath):
        """Add the absolute filepath to the file history."""
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

    def onFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        self.file_history.AddFileToHistory(filepath)

        self.BuildUI(notebook=OncvNotebook.from_file(self, filepath))

    def makeMenu(self):
        """Creates the main menu."""
        menu_bar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an input file")
        file_menu.Append(wx.ID_SAVE, "&Save", help="Save the input file")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the Gui")
        #file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")

        # Make sub-menu with the list of supported visualizers.
        symbol_menu = wx.Menu()
        self._id2symbol = {}

        for symbol in all_symbols():
            _id = wx.NewId()
            symbol_menu.Append(_id, symbol)
            self._id2symbol[_id] = symbol
            self.Bind(wx.EVT_MENU, self.onNewNotebookFromSymbol, id=_id)

        file_menu.AppendMenu(-1, 'Template from element', symbol_menu)

        file_history = self.file_history = wx.FileHistory(8)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(-1, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.onFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)
        menu_bar.Append(file_menu, "File")

        # Add Mixin menus.
        menu_bar.Append(self.CreateToolsMenu(), "Tools")

        help_menu = wx.Menu()
        help_menu.Append(wx.ID_HELP, "Help ", help="Quick help")
        help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menu_bar.Append(help_menu, "Help")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_OPEN, self.onOpen),
            (wx.ID_CLOSE, self.onClose),
            #(wx.ID_EXIT, self.onExit),
            (wx.ID_SAVE, self.onSave),
            (wx.ID_HELP, self.onHelp),
            (wx.ID_ABOUT, self.onAbout),
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

        toolbar.AddSimpleTool(self.ID_SHOW_INPUT, bitmap("in.png"), "Visualize the input file")
        toolbar.AddSimpleTool(self.ID_RUN_INPUT, bitmap("run.png"), "Run the input file.")

        toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_SHOW_INPUT, self.onShowInput),
            (self.ID_RUN_INPUT, self.onRunInput),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

    def onNewNotebookFromSymbol(self, event):
        symbol = self._id2symbol[event.GetId()]
        self.BuildUI(notebook=OncvNotebook.from_symbol(self, symbol))

    def onShowInput(self, event):
        """Show the input file in a new frame."""
        text = self.notebook.makeInputString()
        SimpleTextViewer(self, text=text, title="Oncvpsp Input").Show()

    def onRunInput(self, event):
        """Build a new generator from the input file, and add it to the queue."""
        text = self.notebook.makeInputString()
        try:
            psgen = OncvGenerator(text, calc_type=self.notebook.calc_type)
        except:
            return awx.showErrorMessage(self)

        frame = PseudoGeneratorsFrame(self, [psgen], title="Run Input")
        frame.launch_psgens()
        frame.Show()

    def onOpen(self, event):
        """Open a file"""
        dialog = wx.FileDialog(self, message="Choose an input file", style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_CANCEL: return

        filepath = dialog.GetPath()
        dialog.Destroy()

        # Add to the history.
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

        self.BuildUI(notebook=OncvNotebook.from_file(self, filepath))

    def onSave(self, event):
        """Save a file"""
        dialog = wx.FileDialog(self, message="Save file as...", style=wx.SAVE | wx.OVERWRITE_PROMPT,
                               wildcard="Dat files (*.dat)|*.dat")
        if dialog.ShowModal() == wx.ID_CANCEL: return

        filepath = dialog.GetPath()
        dialog.Destroy()

        # Add to the history.
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

        with open(filepath, "w") as fh:
            fh.write(self.notebook.makeInputString())

    def onClose(self, event):
        """ Respond to the "Close" menu command."""
        self.Destroy()

    def onAbout(self, event):
        return awx.makeAboutBox(
            codename=self.codename,
            version=self.VERSION,
            description="oncvgui is a front-end for the pseudopotential generator oncvpsp",
            developers=["Matteo Giantomassi"],
            website="http://www.mat-simresearch.com/")


class OptimizationFrame(awx.Frame, metaclass=abc.ABCMeta):
    """Base class for optimization frames."""
    def __init__(self, parent, **kwargs):
        super(OptimizationFrame, self).__init__(parent, **kwargs)

        # All optimization buttons are disabled when we start an optimization.
        #self.main_frame.notebook.enable_all_optimize_buttons(False)
        #self.main_frame.Enable(False)
        #self.Bind(wx.EVT_WINDOW_DESTROY, self.onDestroy)

    @property
    def main_frame(self):
        return self.getParentWithType(WxOncvFrame)

    def onDestroy(self, event):
        """Enable all optimize_buttons before destroying the Frame."""
        #self.main_frame.notebook.enable_all_optimize_buttons(True)
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
                psgen = OncvGenerator(str(inp), calc_type=self.notebook.calc_type)
                psgens.append(psgen)
            except:
                return awx.showErrorMessage(self)

        frame = PseudoGeneratorsFrame(self, psgens, title=self.opt_type)
        frame.launch_psgens()
        frame.Show()

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

    def __init__(self, parent, notebook, **kwargs):
        """
        Args:
            notebook:
                `OncvNotebool` containing the parameters of the template.
        """
        super(LlocOptimizationFrame, self).__init__(parent, **kwargs)

        # Save reference to the input panel.
        self.notebook = notebook

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
        # Generate new list of inputs.
        base_inp = self.notebook.makeInput()
        return base_inp.optimize_vloc()


class Rc5OptimizationFrame(OptimizationFrame):
    """
    This frame allows the user to optimize the parameters for the local part of the pseudopotential
    """
    HELP_MSG = """\
This window allows you to change/optimize the rc5 parameter"""

    def __init__(self, parent, notebook, **kwargs):
        """
        Args:
            notebook:
                `OncvNotebool` containing the parameters of the template.
        """
        super(Rc5OptimizationFrame, self).__init__(parent, **kwargs)

        # Save reference to the input panel.
        self.notebook = notebook

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        self.rc5_ctlr = awx.IntervalControl(self, start=0.25, num=6, step=0.05, choices=[">", "centered", "<"])
        main_sizer.Add(self.rc5_ctlr, 1, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        buttons_sizer = self.make_buttons()
        main_sizer.Add(buttons_sizer, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        self.SetSizerAndFit(main_sizer)

    @property
    def opt_type(self):
        return "rc5 optimization"

    def build_new_inps(self):
        # Generate new list of inputs.
        base_inp = self.notebook.makeInput()
        # TODO
        return base_inp.optimize_rc5()


class DeblOptimizationFrame(OptimizationFrame):
    """
    This frame allows the user to optimize the VKB projectors
    """

    HELP_MSG = """\
This window allows you to optimize the parameters used to construct the VKB projectors."""

    def __init__(self, parent, notebook, **kwargs):
        """
        Args:
            notebook:
                Notebook containing the parameters of the template.
        """
        super(DeblOptimizationFrame, self).__init__(parent, **kwargs)

        # Save reference to the input panel.
        self.notebook = notebook
        lmax = notebook.lmax

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        # Build list of controls to allow the user to select the list of
        # debl for the different l-channels. One can disable l-channels via checkboxes.
        self.checkbox_l = [None] * (lmax + 1)
        self.debl_range_l = [None] * (lmax + 1)

        debl_l = notebook.makeInput().debl_l

        for l in range(lmax + 1):
            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Angular Channel L=%d" % l), wx.VERTICAL)
            vsz = wx.BoxSizer(wx.HORIZONTAL)

            self.checkbox_l[l] = check_l = wx.CheckBox(self, -1)
            check_l.SetToolTipString("Enable/Disable optimization for this l-channel")
            check_l.SetValue(True)
            self.debl_range_l[l] = debl_range_l = awx.IntervalControl(self, start=debl_l[l], num=3, step=0.5)

            vsz.Add(check_l, **add_opts)
            vsz.Add(debl_range_l, **add_opts)

            sbox_sizer.Add(vsz, 1, wx.ALL, 5)
            main_sizer.Add(sbox_sizer, 1, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        buttons_sizer = self.make_buttons()
        main_sizer.Add(buttons_sizer, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

        self.SetSizerAndFit(main_sizer)

    @property
    def opt_type(self):
        return "Vkb optimization"

    def build_new_inps(self):
        # Get the list of angular channels activated and the corresponding arrays.
        l_list, deblvals_list = [], []
        for l, (checkbox, wxrange) in enumerate(zip(self.checkbox_l, self.debl_range_l)):
            if checkbox.IsChecked():
                l_list.append(l)
                deblvals_list.append(wxrange.getValues())

        # Generate new list of inputs.
        base_inp = self.notebook.makeInput()
        new_inps = []
        for l, new_debls in zip(l_list, deblvals_list):
            new_inps.extend(base_inp.optimize_debls_for_l(l=l, new_debls=new_debls))

        return new_inps


class FcfactOptimizationFrame(OptimizationFrame):
    """
    This frame allows the user to optimize the model core charge
    """

    HELP_MSG = """\
This window allows you to change/optimize the parameters governing the model core charge"""

    def __init__(self, parent, notebook, **kwargs):
        """
        Args:
            notebook:
                Notebook containing the parameters of the template.
        """
        super(FcfactOptimizationFrame, self).__init__(parent, **kwargs)

        # Save reference to the input panel.
        self.notebook = notebook

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
        return "fcfact optimization"

    def build_new_inps(self):
        fcfact_list = self.fcfact_ctl.getValues()
        print(fcfact_list)

        # Generate new list of inputs.
        base_inp = self.notebook.makeInput()
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

    def __init__(self, parent, notebook, **kwargs):
        """
        Args:
            notebook:
                Notebook containing the parameters of the template.
        """
        super(QcutOptimizationFrame, self).__init__(parent, **kwargs)

        # Save reference to the input panel.
        self.notebook = notebook
        lmax = notebook.lmax

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        # Build list of controls to allow the user to select the list of
        # qcuts for the different l-channels. One can disable l-channels via checkboxes.
        self.checkbox_l = [None] * (lmax + 1)
        self.wxqcut_range_l = [None] * (lmax + 1)

        qcut_l = notebook.makeInput().qcut_l

        for l in range(lmax + 1):
            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Angular Channel L=%d" % l), wx.VERTICAL)
            vsz = wx.BoxSizer(wx.HORIZONTAL)

            self.checkbox_l[l] = check_l = wx.CheckBox(self, -1)
            check_l.SetToolTipString("Enable/Disable optimization for this l-channel")
            check_l.SetValue(True)

            self.wxqcut_range_l[l] = qcrange_l = awx.IntervalControl(self, start=qcut_l[l], num=4, step=0.5)

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
        # Get the list of angular channels activated and the corresponding arrays.
        l_list, qcvals_list = [], []
        for l, (checkbox, wxrange) in enumerate(zip(self.checkbox_l, self.wxqcut_range_l)):
            if checkbox.IsChecked():
                l_list.append(l)
                qcvals_list.append(wxrange.getValues())
        #print("\nl_list:", l_list, "\nqcvals_list", qcvals_list)

        # Generate new list of inputs.
        base_inp = self.notebook.makeInput()
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

    def __init__(self, parent, notebook, **kwargs):
        """
        Args:
            notebook:
                Notebook containing the parameters of the template.
        """
        super(RcOptimizationFrame, self).__init__(parent, **kwargs)

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Save reference to the input panel.
        self.notebook = notebook
        lmax = notebook.lmax

        # Build list of controls to allow the user to select the list of
        # qcut values for the different l-channels.
        # One can disable particular l-channels via checkboxes.
        self.checkbox_l = [None] * (lmax + 1)
        self.wxrc_range_l = [None] * (lmax + 1)

        rc_l = notebook.makeInput().rc_l

        add_opts = dict(flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        for l in range(lmax + 1):
            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Angular Channel L=%d" % l), wx.VERTICAL)
            vsz = wx.BoxSizer(wx.HORIZONTAL)

            self.checkbox_l[l] = check_l = wx.CheckBox(self, -1)
            check_l.SetToolTipString("Enable/Disable optimization for this l-channel")
            check_l.SetValue(True)

            self.wxrc_range_l[l] = qcrange_l = awx.IntervalControl(self, start=rc_l[l], num=4, step=0.05)

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
        base_inp = self.notebook.makeInput()
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
            data = OrderedDict([(k, v.get("value")) for k, v in cls.WXCTRL_PARAMS.items()])

        elif cls.ftype == cls.FTYPE_TABLE:
            nrows = cls.nrows_from_dims(oncv_dims)
            data = nrows * [None]
            for i in range(nrows):
                data[i] = OrderedDict([(k, v.get("value")) for k, v in cls.WXCTRL_PARAMS.items()])

        else:
            raise NotImplementedError(str(cls.ftype))

        return data

    #@property
    #def filepos(self):
    #    for pos, field in enumerate(_FIELD_LIST):
    #        if isinstance(self, field):
    #            return pos
    #    else:
    #        raise ValueError("Cannot find position of class:" % self.__class__.__name__)

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
        #print("About to read: ", type(self), "\nlines=\n", "\n".join(lines))
        okeys = self.WXCTRL_PARAMS.keys()
        odtypes = [v["dtype"] for v in self.WXCTRL_PARAMS.values()]
        parsers = [self.parser_for_dtype[ot] for ot in odtypes]
        #print("okeys", okeys, "odtypes", odtypes)

        if self.ftype == self.FTYPE_ROW:
            assert len(lines) == 1
            tokens = lines[0].split()
            #print("row tokens", tokens)
            #if self.__class__ == VlocalField: tokens[-1] = int(tokens[-1])

            for key, p, tok in zip(okeys, parsers, tokens):
                #print(key)
                try:
                    self.data[key] = p(tok)
                except Exception:
                    print("Exception while trying to convert: key= %s, tok= %s" % (key, tok))
                    raise

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
        ("atsym", dict(dtype="cbox", choices=all_symbols())),
        ("z", dict(dtype="i")),
        ("nc", dict(dtype="i", value=0, tooltip="number of core states"),),
        ("nv", dict(dtype="i", value=0, tooltip="number of valence states")),
        #("iexc", dict(dtype="i", value=4, tooltip="xc functional")),
        # GGA-PBE
        #("iexc", dict(dtype="f", value=4, tooltip="xc functional")),
        #("iexc", dict(dtype="cbox", value="4", choices=["-001013", "4"], tooltip="xc functional")),
        # LDA
        ("iexc", dict(dtype="cbox", value="4", choices=["-001012", "4"], tooltip="xc functional")),
        # PBEsol
        #("iexc", dict(dtype="cbox", value="4", choices=["-116133", "4"], tooltip="xc functional")),
        ("psfile", dict(dtype="cbox", value="psp8", choices=["psp8", "upf", "both"]))])


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
        ("l", dict(dtype="i", value=0)),
        ("rc", dict(dtype="f", value=3.0)),
        ("ep", dict(dtype="f", value=0.0)),
        ("ncon", dict(dtype="i", value=4)),
        ("nbas", dict(dtype="i", value=7)),
        ("qcut", dict(dtype="f", value=6.0))])

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
        ("lmax", dict(dtype="i", value=2))])


@add_tooltips
class VlocalField(RowField):
    name = "LOCAL POTENTIAL"

    WXCTRL_PARAMS = OrderedDict([
        ("lloc", dict(dtype="i", value=4)),
        ("lpopt", dict(dtype="i", value=5)),
        ("rc5", dict(dtype="f", value=3.0)),
        ("dvloc0", dict(dtype="f", value=0))])


@add_tooltips
class VkbConfsField(TableField):
    name = "VANDERBILT-KLEINMAN-BYLANDER PROJECTORs"

    WXCTRL_PARAMS = OrderedDict([
        ("l", dict(dtype="i", value=0)),
        ("nproj", dict(dtype="i", value=2)),
        ("debl", dict(dtype="f", value=1.0))])

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
        ("icmod", dict(dtype="i", value=0)),
        ("fcfact", dict(dtype="f", value=0.25)),
        ("rcfact", dict(dtype="f", value=0.0)),
        ])


@add_tooltips
class LogDerField(RowField):
    name = "LOG DERIVATIVE ANALYSIS"

    WXCTRL_PARAMS = OrderedDict([
        ("epsh1", dict(dtype="f", value=-12.0)),
        ("epsh2", dict(dtype="f", value=+12.0)),
        ("depsh", dict(dtype="f", value=0.02))])


@add_tooltips
class RadGridField(RowField):
    name = "OUTPUT GRID"

    WXCTRL_PARAMS = OrderedDict([
        ("rlmax", dict(dtype="f", value=6.0, step=1.0)),
        ("drl", dict(dtype="f", value=0.01))])


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
        """Initialize the object from an external input file."""
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
        #print("lmax = ",lmax)

        # TODO
        # number of tests and number of rows for the different configurations.
        ncnf = 0

        # Initialize the object
        new = OncvInput(oncv_dims=dict(atsym=atsym, nc=nc, nv=nv, lmax=lmax, ncnf=ncnf))

        # TODO
        # Fill it
        start = 0
        for field in new:
            stop = start + field.nrows
            #print(type(field))
            field.set_vars_from_lines(lines[start:stop])
            start = stop

        return new

    @classmethod
    def from_symbol(cls, symbol):
        """
        Return a tentative input file for generating a pseudo for the given chemical symbol

        .. note:
              Assume default values that might not be optimal.
        """
        nc, nv, lmax = 0, 0, 0
        #atom = nist.get_neutral_entry(symbol=symbol)
        #for state in atom.states:
        #    lmax = max(lmax, state.l)

        # E.g., The electronic structure for Fe is represented as:
        # [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 6), (4, "s", 2)]
        element = Element[symbol]
        for (n, lchar, f) in element.full_electronic_structure:
            nc += 1
            lmax = max(lmax, char2l(lchar))

        # FIXME
        lmax = 1
        lmax = 2
        #lmax = 3

        ncnf = 0
        oncv_dims = dict(atsym=symbol, nc=nc, nv=nv, lmax=lmax, ncnf=ncnf)

        new = cls(oncv_dims)

        field = new.fields[_FIELD_LIST.index(RefConfField)]
        for row, (n, lchar, f) in zip(field.data, element.full_electronic_structure):
            row["n"], row["l"], row["f"] = n, char2l(lchar), f

        return new

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
                self.fields[i] = empty_field(i, self.dims)

            header = self.fields[_FIELD_LIST.index(AtomConfField)]
            header.set_var("atsym", self.dims.atsym)
            header.set_var("z", Element[self.dims.atsym].Z)
            header.set_var("nc", self.dims.nc)
            header.set_var("nv", self.dims.nv)

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
        return s + "\n 0\n"

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
    def debl_l(self):
        """List with the values of debl as function of l."""
        i = _FIELD_LIST.index(VkbConfsField)
        return self.fields[i].get_col("debl")

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

    def optimize_debls_for_l(self, l, new_debls):
        """
        Returns a list of new input objects in which the debls parameter for
        the given l has been replaced by the values listed in new_debls.

        Args:
            l:
                Angular momentum
            new_debls:
                Iterable with the new values of debls.
                The returned list will have len(new_debls) input objects.
        """
        # Find the field with the configuration parameters.
        i = _FIELD_LIST.index(VkbConfsField)

        # Find the row with the given l.
        for irow, row in enumerate(self.fields[i].data):
            if row["l"] == l:
                break
        else:
            raise ValueError("Cannot find l %s in the VkbConfsField" % l)

        # This is the dict we want to change
        inps = []
        for debl in new_debls:
            new_inp = self.deepcopy()
            new_inp.fields[i].data[irow]["debl"] = debl
            inps.append(new_inp)

        return inps


class OncvNotebook(wx.Notebook):

    @classmethod
    def from_file(cls, parent, filename):
        inp = OncvInput.from_file(filename)
        new = cls(parent, inp.dims)

        for field in inp:
            wxctrl = new.wxctrls[field.__class__]
            wxctrl.SetParams(field.data)

        return new

    @classmethod
    def from_symbol(cls, parent, symbol):
        inp = OncvInput.from_symbol(symbol)
        new = cls(parent, inp.dims)

        for field in inp:
            wxctrl = new.wxctrls[field.__class__]
            wxctrl.SetParams(field.data)

        return new

    #def fromInput(cls)

    def __init__(self, parent, oncv_dims):
        super(OncvNotebook, self).__init__(parent)

        # Build tabs
        self.oncv_dims = oncv_dims
        self.ae_tab = AeConfTab(self, oncv_dims)
        self.ps_tab = PsConfTab(self, oncv_dims)
        self.pstests_tab = PsTestsTab(self, oncv_dims)

        # Add tabs
        self.AddPage(self.ae_tab, "AE config")
        self.AddPage(self.ps_tab, "PP config")
        self.AddPage(self.pstests_tab, "Tests")

    @property
    def tabs(self):
        return (self.ae_tab, self.ps_tab, self.pstests_tab)

    @property
    def wxctrls(self):
        d = {}
        for tab in self.tabs:
            d.update(tab.wxctrls)
        return d

    @property
    def calc_type(self):
        return self.ae_tab.calctype_cbox.GetValue()

    @property
    def lmax(self):
        # TODO: property or method?
        return self.ps_tab.wxctrls[LmaxField].GetParams()["lmax"]

    def getElement(self):
        symbol = self.ae_tab.wxctrls[AtomConfField].GetParams()["atsym"]
        return Element[symbol]

    def makeInput(self):
        """Build an OncvInput instance from the values specified in the controllers."""
        inp = OncvInput(self.oncv_dims)

        for cls, wxctrl in self.wxctrls.items():
            i = _FIELD_LIST.index(cls)
            inp.fields[i].set_vars(wxctrl.GetParams())

        return inp

    def makeInputString(self):
        """Return a string with the input passed to the pp generator."""
        return str(self.makeInput())


class AeConfTab(awx.Panel):
    def __init__(self, parent, oncv_dims):
        super(AeConfTab, self).__init__(parent)

        # Set the dimensions and build the widgets.
        self.oncv_dims = oncv_dims
        self.buildUI()

    def buildUI(self):
        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        stext = wx.StaticText(self, -1, "Calculation type:")
        choices = ["scalar-relativistic", "fully-relativistic", "non-relativistic"]
        self.calctype_cbox = wx.ComboBox(
            self, id=-1, name='Calculation type', choices=choices, value=choices[0], style=wx.CB_READONLY)

        add_opts = dict(proportion=0, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
        sizer_addopts = dict(proportion=0, flag=wx.ALL, border=5)

        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        hbox0.Add(stext, **add_opts)
        hbox0.Add(self.calctype_cbox)

        main_sizer.Add(hbox0, **add_opts)

        #(list_of_classes_on_row, show)
        layout = [
            AtomConfField,
            RefConfField,
            RadGridField,
            ]

        # Each field has a widget that returns the variables in a dictionary
        self.wxctrls = {}
        self.sbox_sizers = {}
        #self.sboxes = {}

        for cls in layout:
            i = _FIELD_LIST.index(cls)
            f = empty_field(i, self.oncv_dims)
            wxctrl = f.make_wxctrl(self)
            sbox = wx.StaticBox(self, -1, f.name + ":")
            sbox_sizer = wx.StaticBoxSizer(sbox, wx.VERTICAL)
            #sbox_sizer = wx.BoxSizer(wx.VERTICAL)
            sbox_sizer.Add(wxctrl, **sizer_addopts)
            #sz = wx.FlexGridSizer(wx.VERTICAL)
            #sz.Add(sbox_sizer)

            self.wxctrls[cls] = wxctrl
            self.sbox_sizers[cls] = sbox_sizer
            #self.sboxes[cls] = sz

            main_sizer.Add(sbox_sizer, **add_opts)

        lda_levels_button = wx.Button(self, -1, "LDA levels (NIST)")
        lda_levels_button.Bind(wx.EVT_BUTTON, self.onShowLdaLevels)
        main_sizer.Add(lda_levels_button, **add_opts)

        add_button = wx.Button(self, -1, "Add level")
        add_button.action = "add"
        add_button.Bind(wx.EVT_BUTTON, self.onAddRemoveLevel)
        main_sizer.Add(add_button, **add_opts)

        remove_button = wx.Button(self, -1, "Remove level")
        remove_button.action = "remove"
        remove_button.Bind(wx.EVT_BUTTON, self.onAddRemoveLevel)
        main_sizer.Add(remove_button, **add_opts)

        self.SetSizerAndFit(main_sizer)

    @property
    def atomconf(self):
        return self.wxctrls[AtomConfField].GetParams()

    @property
    def nc(self):
        return self.atomconf["nc"]

    @property
    def nv(self):
        return self.atomconf["nv"]

    @property
    def symbol(self):
        return self.atomconf["atsym"]

    def get_core(self):
        # E.g., The electronic structure for Fe is represented as:
        # [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 6), (4, "s", 2)]
        core = []
        refconf = self.wxctrls[RefConfField]
        for ic in range(self.nc):
            d = refconf[ic].GetParams()
            t = [d[k] for k in ("n", "l", "f")]
            core.append(tuple(t))

        return core

    def get_valence(self):
        valence = []
        refconf = self.wxctrls[RefConfField]
        for iv in range(self.nc, self.nc + self.nv):
            d = refconf[iv].GetParams()
            t = [d[k] for k in ("n", "l", "f")]
            valence.append(tuple(t))

        return valence

    def onAddRemoveLevel(self, event):
        button = event.GetEventObject()
        sbox_sizer = self.sbox_sizers[RefConfField]
        old = self.wxctrls[RefConfField]
        sbox_sizer.Hide(0)
        sbox_sizer.Remove(0)

        if button.action == "add":
            old.appendRow()
        elif button.action == "remove":
            old.removeRow()

        sizer_addopts = dict(proportion=0, flag=wx.ALL, border=5)
        sbox_sizer.Insert(0, old, **sizer_addopts)

        #self.sboxes[RefConfField].Layout()
        sbox_sizer.Show(0)
        sbox_sizer.Layout()
        #self.main_sizer.Layout()
        self.Layout()
        self.Fit()

        #frame = self.GetParent().GetParent()
        #frame.fSizer.Layout()
        #frame.Fit()

    def onShowLdaLevels(self, event):
        # Get the LDA levels of the neutral atom.
        # (useful to decide if semicore states should be included in the valence).
        entry = nist.get_neutral_entry(self.symbol)
        frame = awx.Frame(self, title="LDA levels for neutral %s (NIST database)" % self.symbol)
        awx.ListCtrlFromTable(frame, table=entry.to_table())
        frame.Show()


class PsConfTab(awx.Panel):
    def __init__(self, parent, oncv_dims):
        super(PsConfTab, self).__init__(parent)

        self.notebook = parent

        # Set the dimensions and build the widgets.
        self.oncv_dims = oncv_dims
        self.buildUI()

    def buildUI(self):
        sizer_addopts = dict(proportion=0, flag=wx.ALL, border=5)
        add_opts = dict(proportion=0, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        # list_of_classes
        layout = [
            LmaxField,
            PseudoConfField,
            VkbConfsField,
            VlocalField,
            ModelCoreField,
        ]

        # FieldClass: [(label, OptimizationFrame), ....]
        fields_with_optimization = {
            PseudoConfField: [
                ("Change rc", RcOptimizationFrame),
                ("Change qcut", QcutOptimizationFrame),
                ],
            VlocalField: [
                ("Change lloc", LlocOptimizationFrame),
                ("Change rc5", Rc5OptimizationFrame),
                ],
            ModelCoreField: [("Change fcfact", FcfactOptimizationFrame)],
            VkbConfsField: [("Change debl", DeblOptimizationFrame)],
        }

        # Each field has a widget that returns the variables in a dictionary
        self.wxctrls = {}

        # Keep an internal list of buttons so that we can disable them easily.
        self.all_optimize_buttons = []
        self.sboxes = {}

        for cls in layout:
            i = _FIELD_LIST.index(cls)
            f = empty_field(i, self.oncv_dims)
            wxctrl = f.make_wxctrl(self)

            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, f.name + ":"), wx.VERTICAL)
            #sbox_sizer = wx.BoxSizer(wx.VERTICAL)
            sbox_sizer.Add(wxctrl, **sizer_addopts)

            self.sboxes[cls] = sbox_sizer
            self.wxctrls[cls] = wxctrl

            # add optimization button if the field supports it.
            if f.__class__ in fields_with_optimization:
                hsz = wx.BoxSizer(wx.HORIZONTAL)
                for label, opt_frame in fields_with_optimization[f.__class__]:
                    optimize_button = wx.Button(self, -1, label)
                    optimize_button.Bind(wx.EVT_BUTTON, self.onOptimize)
                    optimize_button.field_class = f.__class__
                    optimize_button.opt_frame = opt_frame

                    self.all_optimize_buttons.append(optimize_button)
                    hsz.Add(optimize_button, 0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=5)

                sbox_sizer.Add(hsz, 0, flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL)

            main_sizer.Add(sbox_sizer, **add_opts)

        #self.old_lmax = self.lmax_ctrl.GetParams()["lmax"]
        #self.wxctrls[LmaxField].Bind(

        add_button = wx.Button(self, -1, "Add row")
        add_button.Bind(wx.EVT_BUTTON, self.onLmaxChanged)
        main_sizer.Add(add_button)

        self.SetSizerAndFit(main_sizer)

    @property
    def lmax_ctrl(self):
        return self.wxctrls[LmaxField]

    def onLmaxChanged(self, event):
        #self.wxctrls[PseudoConfField].appendRow()
        #self.wxctrls[VkbConfsField].appendRow()

        #main_sizer = self.main_sizer
        #main_sizer.Hide(1)
        #main_sizer.Remove(1)
        #main_sizer.Insert(1, self.wxctrls[PseudoConfField], 0, wx.EXPAND | wx.ALL, 5)

        #for sbox in self.sboxes.values():
        #    sbox.Layout()
        self.main_sizer.Layout()

    def onOptimize(self, event):
        button = event.GetEventObject()
        #self.enable_all_optimize_buttons(False)
        opt_frame = button.opt_frame

        try:
            opt_frame(self, notebook=self.notebook).Show()
        finally:
            self.enable_all_optimize_buttons(True)

    def enable_all_optimize_buttons(self, enable=True):
        """Enable/disable the optimization buttons."""
        for button in self.all_optimize_buttons:
            button.Enable(enable)


class PsTestsTab(awx.Panel):
    def __init__(self, parent, oncv_dims):
        super(PsTestsTab, self).__init__(parent)
        self.notebook = parent

        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)
        add_opts = dict(proportion=0, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL | wx.EXPAND, border=5)

        layout = [
            LogDerField,
        ]

        self.wxctrls = {}

        for cls in layout:
            i = _FIELD_LIST.index(cls)
            f = empty_field(i, oncv_dims)
            wxctrl = f.make_wxctrl(self)
            self.wxctrls[cls] = wxctrl
            sbox_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, f.name + ":"), wx.VERTICAL)
            sbox_sizer.Add(wxctrl, **add_opts)
            main_sizer.Add(sbox_sizer, **add_opts)

        self.conf_txtctrl = wx.TextCtrl(self, -1, "#TEST CONFIGURATIONS\n  0", style=wx.TE_MULTILINE | wx.TE_LEFT)
        main_sizer.Add(self.conf_txtctrl, **add_opts)

        common_button = wx.Button(self, -1, label='Common Oxidation States')
        common_button.Bind(wx.EVT_BUTTON, self.addOxiStates)
        common_button.test_mode = "common_oxis"

        all_button = wx.Button(self, -1, label='All Oxidation States')
        all_button.Bind(wx.EVT_BUTTON, self.addOxiStates)
        all_button.test_mode = "all_oxis"

        hsz = wx.BoxSizer(wx.HORIZONTAL)
        hsz.Add(common_button, **add_opts)
        hsz.Add(all_button, **add_opts)
        main_sizer.Add(hsz, **add_opts)

        self.SetSizerAndFit(main_sizer)

    def conftests_str(self):
        return self.conf_txtctrl.GetValue()

    def addOxiStates(self, event):
        button = event.GetEventObject()
        test_mode = button.test_mode

        element = self.notebook.getElement()

        if test_mode == "common_oxis":
            oxi_states = element.common_oxidation_states
            title = "(common oxidation states)"
        elif test_mode == "all_oxis":
            oxi_states = element.oxidation_states
            title = "(all oxidation states)"
        else:
            raise ValueError("Wrong test_mode %s" % test_mode)

        self.conf_txtctrl.Clear()
        self.conf_txtctrl.WriteText("# TEST CONFIGURATIONS" + title + "\n")

        core_aeatom = self.notebook.ae_tab.get_core()
        val_aeatom = self.notebook.ae_tab.get_valence()
        print("core", core_aeatom)
        print("val", val_aeatom)

        """
        # TEST CONFIGURATIONS
        # ncnf
            3
        #
        #   nvcnf (repeated ncnf times)
        #   n, l, f  (nvcnf lines, repeated follwing nvcnf's ncnf times)
            2
            2    0    2.0
            2    1    3.0
        #
            2
            2    0    1.0
            2    1    4.0
        #
            2
            2    0    1.0
            2    1    3.0
        """
        test_confs = []

        for oxi in oxi_states:
            #self.conf_txtctrl.WriteText(str(oxi) + "\n")

            # Get the electronic configuration of atom with Z = Z + oxi
            if oxi == 0: continue

            # Here we found the valence configuration by comparing
            # the full configuration of oxiele and the one of the initial element.
            oxi_element = Element.from_Z(element.Z - oxi)
            oxi_estruct = oxi_element.full_electronic_structure
            oxi_estruct = [(t[0], char2l(t[1]), t[2]) for t in oxi_estruct]

            if oxi < 0:
                test_conf = [t for t in oxi_estruct if t not in core_aeatom]
            else:
                test_conf = [t for t in oxi_estruct if t not in core_aeatom]

            self.conf_txtctrl.WriteText(str(oxi) + "\n")
            self.conf_txtctrl.WriteText(str(oxi_element) + "\n")
            self.conf_txtctrl.WriteText(str(test_conf) + "\n")

        self.conf_txtctrl.WriteText("# ncnf\n" + str(len(test_confs)) + "\n")
        self.conf_txtctrl.WriteText("""\
#
#   nvcnf (repeated ncnf times)
#   n, l, f  (nvcnf lines, repeated follwing nvcnf's ncnf times)\n""")

        #for test in test_confs:
        #    self.conf_txtctrl.WriteText(len(test) + "\n")
        #    for row in test:
        #        self.conf_txtctrl.WriteText(str(row) + "\n")

        self.main_sizer.Layout()


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
    _COLUMNS = ["#", 'status', "max_ecut", "atan_logder_err", "max_psexc_abserr", "herm_err", "warnings"]

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
    def notebook(self):
        """Reference to the Wx window with the input file."""
        try:
            return self._notebook
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
            self._notebook = parent.notebook
            return self._notebook

    @staticmethod
    def make_entry(index, psgen):
        """Returns the entry associated to the generator psgen with the given index."""
        max_ecut = None
        max_atan_logder_l1err = None
        max_psexc_abserr = None
        herm_err = None

        if psgen.results is None and psgen.status == psgen.S_OK:
            psgen.check_status()

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
            "%s" % len(psgen.warnings),
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
        self.ID_POPUP_COMPUTE_HINTS = wx.NewId()
        self.ID_POPUP_COMPUTE_GBRV = wx.NewId()
        self.ID_POPUP_COMPUTE_PSPS = wx.NewId()
        self.ID_POPUP_SAVE_PSGEN = wx.NewId()

        menu = wx.Menu()

        # Make sub-menu with the list of supported quantities
        plot_submenu = wx.Menu()

        # TODO: this list could be taken from the class or from the plotter instance.
        all_keys = [
            "radial_wfs",
            "projectors",
            "densities",
            "potentials",
            "der_potentials",
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
        menu.Append(self.ID_POPUP_COMPUTE_HINTS, "Compute hints for ecut")
        #menu.Append(self.ID_POPUP_COMPUTE_GBRV, "Perform GBRV tests")
        menu.Append(self.ID_POPUP_COMPUTE_PSPS, "Get PSPS.nc file and plot data")
        menu.Append(self.ID_POPUP_SAVE_PSGEN, "Save PS generation")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_POPUP_STDIN, self.onShowStdin),
            (self.ID_POPUP_STDOUT, self.onShowStdout),
            (self.ID_POPUP_STDERR, self.onShowStderr),
            (self.ID_POPUP_CHANGE_INPUT, self.onChangeInput),
            (self.ID_POPUP_COMPUTE_HINTS, self.onComputeHints),
            #(self.ID_POPUP_COMPUTE_GBRV, self.onGBRV),
            (self.ID_POPUP_COMPUTE_PSPS, self.onPsps),
            (self.ID_POPUP_SAVE_PSGEN, self.onSavePsgen),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu

    def onChangeInput(self, event):
        """Change the template input file."""
        psgen = self.getSelectedPseudoGen()
        if psgen is None: return

        # Change the parameters in the notebook using those in psgen.stdin_path
        if not os.path.exists(psgen.stdin_path):
            return awx.showErrorMessage(self, "Input file %s does not exist" % psgen.stdin_path)

        main_frame = self.notebook.GetParent()
        new_notebook = OncvNotebook.from_file(main_frame, psgen.stdin_path)
        main_frame.BuildUI(notebook=new_notebook)

    def onComputeHints(self, event):
        psgen = self.getSelectedPseudoGen()
        if psgen is None: return
        print("workdir", psgen.workdir)

        # Change the parameters in the notebook using those in psgen.stdin_path
        #if not os.path.exists(psgen.stdin_path):
        #    return awx.showErrorMessage(self, "Input file %s does not exist" % psgen.stdin_path)
        from abipy import abilab
        from pseudo_dojo.dojo.dojo_workflows import PPConvergenceFactory
        factory = PPConvergenceFactory()

        print("pseudo", psgen.pseudo)
        workdir = os.path.join("HINTS_", psgen.pseudo.name)
        flow = abilab.AbinitFlow(workdir, pickle_protocol=0)

        work = factory.work_for_pseudo(psgen.pseudo, ecut_slice=slice(4, None, 1), nlaunch=2)
        flow.register_work(work)
        flow.allocate()
        flow.build_and_pickle_dump()
        scheduler = abilab.PyFlowScheduler.from_user_config()
        scheduler.add_flow(flow)
        scheduler.start()

    def onSavePsgen(self, event):
        # Update the notebook
        self.onChangeInput(event)

        psgen = self.getSelectedPseudoGen()
        if psgen is None: return

        main_frame = self.notebook.GetParent()
        input_file = main_frame.input_file

        if input_file is None:
            raise NotImplementedError()

        dirpath = os.path.dirname(input_file)
        basename = os.path.basename(input_file).replace(".in", "")

        ps_dest = os.path.join(dirpath, basename + ".psp8")
        upf_dest = os.path.join(dirpath, basename + ".upf")
        out_dest = os.path.join(dirpath, basename + ".out")
        djrepo_dest = os.path.join(dirpath, basename + ".djrepo")

        exists = []
        for f in [input_file, ps_dest, upf_dest, out_dest, djrepo_dest]:
            if os.path.exists(f): exists.append(os.path.basename(f))

        if exists:
            msg = "File(s):\n%s already exist.\nDo you want to owerwrite them?" % "\n".join(exists)
            answer = awx.askUser(self, msg)
            if not answer: return

        # Update the input file, then copy the pseudo file and the output file.
        with open(input_file, "wt") as fh:
            fh.write(self.notebook.makeInputString())

        print(psgen.pseudo.path)
        shutil.copy(psgen.pseudo.path, ps_dest)
        upf_src = psgen.pseudo.path.replace(".psp8", ".upf")
        if os.path.exists(upf_src):
            shutil.copy(upf_src, upf_dest)
        shutil.copy(psgen.stdout_path, out_dest)

        # Parse the output file
        onc_parser = OncvOutputParser(out_dest)
        onc_parser.scan()
        if not onc_parser.run_completed:
            raise RuntimeError("oncvpsp output is not complete. Exiting")

        # Build dojoreport
        pseudo = Pseudo.from_file(ps_dest)
        report = DojoReport.empty_from_pseudo(pseudo, onc_parser.hints, devel=False)
        report.json_write()

    def onGBRV(self, event):
        psgen = self.getSelectedPseudoGen()
        if psgen is None: return

        # Change the parameters in the notebook using those in psgen.stdin_path
        #if not os.path.exists(psgen.stdin_path):
        #    return awx.showErrorMessage(self, "Input file %s does not exist" % psgen.stdin_path)
        from abipy import abilab
        from pseudo_dojo.dojo.dojo_workflows import GbrvFactory
        factory = GbrvFactory()

        flow = abilab.AbinitFlow(workdir="GBRV", pickle_protocol=0)
        print("pseudo", psgen.pseudo)
        for struct_type in ["fcc", "bcc"]:
            work = factory.relax_and_eos_work(psgen.pseudo, struct_type)
            flow.register_work(work)

        flow.allocate()
        flow.build_and_pickle_dump()
        scheduler = abilab.PyFlowScheduler.from_user_config()
        scheduler.add_flow(flow)
        scheduler.start()

    def onPsps(self, event):
        psgen = self.getSelectedPseudoGen()
        if psgen is None: return

        with psgen.pseudo.open_pspsfile(ecut=30) as psps:
            print("Printing data from:", psps.filepath)
            psps.plot(ecut_ffnl=60)

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

    def launch_psgens(self):
        return self.psgens_wxlist.psgen_list_ctrl.launch_psges()

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
            (self.ID_SHOW_INPUTS, self.onShowInputs),
            (self.ID_SHOW_OUTPUTS, self.onShowOutputs),
            (self.ID_SHOW_ERRS, self.onShowErrors),
            (self.ID_CHECK_STATUS, self.onCheckStatus),
            (self.ID_MULTI_PLOTTER, self.onMultiPlotter),
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
        self.launch_psgens()

    def launch_psgens(self):
        max_nlaunch = int(self.max_nlaunch.GetValue())

        nlaunch = 0
        for psgen in self.psgens:
            nlaunch += psgen.start()
            if nlaunch == max_nlaunch:
                break

        self.statusbar.PushStatusText("Submitted %d tasks" % nlaunch)

    def onCheckStatus(self, event):
        """
        Callback triggered by the checkstatus button.
        Check the status of the `PseudoGenerators` and refresh the panel.
        """
        self.checkStatusAndRedraw()

    def checkStatusAndRedraw(self):
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
        main_sizer.Insert(0, new_psgen_wxlist, 1, wx.ALL | wx.EXPAND, 5)
        self.psgen_wxlist = new_psgen_wxlist

        self.panel.Layout()

        # Write number of jobs with given status.
        #message = ", ".join("%s: %s" % (k, v) for (k, v) in counter.items())
        #self.statusbar.PushStatusText(message)

    def OnIdle(self, event):
        """Function executed when the GUI is idle."""
        now = time.time()
        if (now - self.last_refresh) > self.REFRESH_INTERVAL:
            self.checkStatusAndRedraw()
            self.last_refresh = time.time()

    def onShowInputs(self, event):
        """Show all input files."""
        TextNotebookFrame(self, text_list=[psgen.get_stdin() for psgen in self.psgens],
                          page_names=["PSGEN # %d" % i for i in range(len(self.psgens))]).Show()

    def onShowOutputs(self, event):
        """Show all output files."""
        TextNotebookFrame(self, text_list=[psgen.get_stdout() for psgen in self.psgens],
                          page_names=["PSGEN # %d" % i for i in range(len(self.psgens))]).Show()

    def onShowErrors(self, event):
        """Show all error files."""
        TextNotebookFrame(self, text_list=[psgen.get_stderr() for psgen in self.psgens],
                          page_names=["PSGEN # %d" % i for i in range(len(self.psgens))]).Show()

    def onPlotColumns(self, event):
        self.psgens_wxlist.plot_columns()

    def onMultiPlotter(self, event):
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


def wxapp_oncvpsp(filepath=None):
    """Standalone application."""
    class OncvApp(awx.App):
        def OnInit(self):
            # Enforce WXAgg as matplotlib backend to avoid nasty SIGSEGV in the C++ layer
            # occurring when WX Guis produce plots with other backends.
            import matplotlib
            matplotlib.use('WXAgg')

            frame = WxOncvFrame(None, filepath)
            #frame = my_periodic_table(None)
            #frame = OncvParamsFrame(None, z=12)
            frame.Show(True)
            self.SetTopWindow(frame)
            return True

    app = OncvApp()
    return app

if __name__ == "__main__":
    import sys
    #onc_inp = OncvInput.from_file("08_O.dat")
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
    try:
        filepaths = sys.argv[1:]
    except IndexError:
        filepaths = None

    if filepaths is not None:
        if filepaths[0] == "table":
            for symbol in all_symbols():
                path = symbol + ".dat"
                if os.path.exists(path):
                    print("Will open file %s" % path)
                    wxapp_oncvpsp(path).MainLoop()
        else:
            for filepath in filepaths:
                #print(filepath)
                wxapp_oncvpsp(filepath).MainLoop()
    else:
        wxapp_oncvpsp().MainLoop()

    #app = awx.App()
    #frame = QcutOptimizationFrame(None, lmax=1)
    #frame.Show()
    #app.MainLoop()


from __future__ import print_function, division, unicode_literals, absolute_import

import wx
import collections
import fnmatch
import abipy.gui.awx as awx

from monty.string import list_strings
from abipy.electrons import ElectronBandsPlotter, ElectronDosPlotter, MdfPlotter, SigresPlotter
from abipy.gui.electronswx import ElectronDosDialog


class FileCheckBoxPanel(awx.Panel):
    """A panel with a list of filepaths and checkboxes."""
    def __init__(self, parent, filepaths, **kwargs):
        """
        Args:
            parent:
                Parent window.
            filepaths:
                String or List of strings with filepaths.
        """
        super(FileCheckBoxPanel, self).__init__(parent, -1, **kwargs)

        self.all_filepaths = list_strings(filepaths)

        self.BuildUi()

    def BuildUi(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        static_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Files"), wx.VERTICAL)

        self.check_boxes = collections.OrderedDict()

        self.wildcards = ["*"]

        for path in self.all_filepaths:
            if not self._SelectFilename(path): continue
            assert path not in self.check_boxes
            cbox = wx.CheckBox(self, -1, path, wx.DefaultPosition, wx.DefaultSize, 0)
            cbox.SetValue(True)
            static_sizer.Add(cbox, 0, wx.ALL | wx.EXPAND, 5)
            self.check_boxes[path] = cbox

        main_sizer.Add(static_sizer, 1, wx.EXPAND, 5)

        # Add buttons to (select|deselect) all checkboxes.
        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        all_button = wx.Button(self, -1, "Select all", wx.DefaultPosition, wx.DefaultSize, 0)
        all_button.Bind(wx.EVT_BUTTON, self.OnSelectAll)
        hsizer.Add(all_button, 0, wx.ALL, 5)

        deselect_button = wx.Button(self, -1, "Deselect all", wx.DefaultPosition, wx.DefaultSize, 0)
        deselect_button.Bind(wx.EVT_BUTTON, self.OnDeselectAll)
        hsizer.Add(deselect_button, 0, wx.ALL, 5)

        filter_label = wx.StaticText(self, -1, "Filter:", wx.DefaultPosition, wx.DefaultSize, 0)
        filter_label.Wrap(-1)
        hsizer.Add(filter_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        wildcard_choices = ["*", "*.nc"]
        self.filter_combobox = wx.ComboBox(self, wx.ID_ANY, "*", wx.DefaultPosition, wx.DefaultSize, wildcard_choices,
                                           0)
        self.filter_combobox.Bind(wx.EVT_COMBOBOX, self.OnFilterComboBox)
        self.filter_combobox.Bind(wx.EVT_TEXT_ENTER, self.OnFilterComboBox)
        hsizer.Add(self.filter_combobox, 0, wx.ALL, 5)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizerAndFit(main_sizer)

    def _SelectFilename(self, filename):
        for wcard in self.wildcards:
            if not fnmatch.fnmatch(filename, wcard):
                return False
        return True

    def OnFilterComboBox(self, event):
        self.wildcards = self.filter_combobox.GetValue().split("|")

        for (filepath, cbox) in self.check_boxes.items():
            if self._SelectFilename(filepath):
                cbox.SetValue(True)
            else:
                cbox.SetValue(False)

    def OnSelectAll(self, event):
        for cbox in self.check_boxes.values():
            cbox.SetValue(True)

    def OnDeselectAll(self, event):
        for cbox in self.check_boxes.values():
            cbox.SetValue(False)

    def GetSelectedFilepaths(self):
        """
        Return the list of filepaths selected by the user.
        Main entry point for client code.
        """
        selected = []

        for path, cbox in self.check_boxes.items():
            if cbox.GetValue():
                selected.append(path)

        return selected


class ComparisonFrame(awx.Frame):
    """
    This frame allows the user to select/deselect a list of files and to produce plots
    for all the files selected. Useful for convergence studies.
    """

    def __init__(self, parent, dirpaths=None, filepaths=None, wildcard=None, **kwargs):
        """
        Args:
            parent:
                parent window
            dirpaths
                List of directory names
            filepaths
                List of filepaths.
            wildcard
                Regular expression for selecting files.
        """
        super(ComparisonFrame, self).__init__(parent, -1, **kwargs)

        # TODO
        #self.dirpaths = dirpaths if dirpaths is not None else []
        #self.filepaths = filepaths if filepaths is not None else []
        #
        #self.dirpaths = map(os.path.abspath, self.dirpaths)
        #self.filepaths = map(os.path.abspath, self.filepaths)
        #
        #self.wildcard = wildcard if wildcard is not None else ""

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        st1 = wx.StaticText(self, -1, "Quantity:", wx.DefaultPosition, wx.DefaultSize, 0)
        st1.Wrap(-1)
        hsizer.Add(st1, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        plotter_choices = ["ebands", "edos", "mdf", "sigres"]
        self.plotter_cbox = wx.ComboBox(self, -1, "ebands", wx.DefaultPosition, wx.DefaultSize, plotter_choices, 0)
        hsizer.Add(self.plotter_cbox, 0, wx.ALL, 5)

        compare_button = wx.Button(self, -1, "Compare", wx.DefaultPosition, wx.DefaultSize, 0)
        compare_button.Bind(wx.EVT_BUTTON, self.OnCompareButton)
        hsizer.Add(compare_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)

        main_sizer.Add(hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.panel = FileCheckBoxPanel(self, filepaths)
        main_sizer.Add(self.panel, 1, wx.EXPAND, 5)

        self.SetSizerAndFit(main_sizer)

    def OnCompareButton(self, event):
        selected_files = self.panel.GetSelectedFilepaths()

        choice = self.plotter_cbox.GetValue()

        try:
            if choice == "ebands":
                plotter = ElectronBandsPlotter()

                for filepath in selected_files:
                    plotter.add_ebands_from_file(filepath)

                plotter.plot()

            elif choice == "edos":
                # Open dialog to ask the user the DOS parameters.
                dos_dialog = ElectronDosDialog(None)

                if dos_dialog.ShowModal() == wx.ID_OK:
                    p = dos_dialog.GetParams()

                    plotter = ElectronDosPlotter()

                    for filepath in selected_files:
                        plotter.add_edos_from_file(filepath, **p)

                    plotter.plot()

                dos_dialog.Destroy()

            elif choice == "mdf":
                plotter = MdfPlotter()

                for filepath in selected_files:
                    plotter.add_mdf_from_file(filepath, mdf_type="exc")

                plotter.plot()

            elif choice == "sigres":
                plotter = SigresPlotter()

                plotter.add_files(selected_files)

                plotter.plot_qpgaps()
                #plotter.plot_qpenes()

            else:
                awx.showErrorMessage(self, message="No function registered for choice %s" % choice)

        except Exception:
            awx.showErrorMessage(self)


def wxapp_comparison(dirpaths=None, filepaths=None, wildcard=None):
    app = awx.App()
    frame = ComparisonFrame(None, dirpaths=dirpaths, filepaths=filepaths, wildcard=wildcard)
    frame.Show()
    return app

#!/usr/bin/env python

import wx
import awx
import collections
import fnmatch

from abipy.electrons import EBandsPlotter

class FileCheckBoxPanel ( wx.Panel ):

    def __init__(self, parent, filepaths, **kwargs):
        wx.Panel.__init__ (self, parent, id=-1, **kwargs)

        if isinstance(filepaths, str): filepaths = [filepaths]
        self.all_filepaths = filepaths

        self.BuildUi()

    def BuildUi(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        static_sizer = wx.StaticBoxSizer( wx.StaticBox(self, -1, "Files"), wx.VERTICAL )

        self.check_boxes = collections.OrderedDict()

        self.wildcards = ["*"]

        for path in self.all_filepaths:
            if not self._SelectFilename(path):
                continue
            assert path not in self.check_boxes
            cbox = wx.CheckBox(self, -1, path, wx.DefaultPosition, wx.DefaultSize, 0)
            cbox.SetValue(True)
            static_sizer.Add(cbox, 0, wx.ALL|wx.EXPAND, 5)
            self.check_boxes[path] = cbox

        main_sizer.Add(static_sizer, 1, wx.EXPAND, 5)

        # Add buttons to (select|deselect) all checkboxes.
        hsizer = wx.BoxSizer( wx.HORIZONTAL )

        self.all_button = wx.Button( self, -1, "Select all", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.Bind(wx.EVT_BUTTON, self.OnSelectAll, self.all_button)
        hsizer.Add( self.all_button, 0, wx.ALL, 5 )

        self.deselect_button = wx.Button( self, -1, "Deselect all", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.Bind(wx.EVT_BUTTON, self.OnDeselectAll, self.deselect_button)
        hsizer.Add( self.deselect_button, 0, wx.ALL, 5 )

        filter_label = wx.StaticText( self, -1, "Filter:", wx.DefaultPosition, wx.DefaultSize, 0 )
        filter_label.Wrap( -1 )
        hsizer.Add(filter_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        wildcard_choices = ["*", "*.nc"]
        self.filter_combobox = wx.ComboBox( self, wx.ID_ANY, "*", wx.DefaultPosition, wx.DefaultSize, wildcard_choices, 0)
        self.filter_combobox.Bind(wx.EVT_COMBOBOX, self.OnFilterComboBox)
        self.filter_combobox.Bind(wx.EVT_TEXT_ENTER, self.OnFilterComboBox)
        hsizer.Add( self.filter_combobox, 0, wx.ALL, 5 )

        main_sizer.Add( hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.SetSizer( main_sizer )
        self.Fit()

    def _SelectFilename(self, filename):
        for wcard in self.wildcards:
            if not fnmatch.fnmatch(filename, wcard):
                return False
        return True

    def OnFilterComboBox(self, event):
        self.wildcards = self.filter_combobox.GetValue().split("|")
        for filepath, cbox in self.check_boxes.items():
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


class ComparisonFrame(wx.Frame):
    def __init__(self, parent, dirpaths=None, filepaths=None, wildcard=None, **kwargs):
        super(ComparisonFrame, self).__init__(parent, -1, **kwargs)

        #self.dirpaths = dirpaths if dirpaths is not None else []
        #self.filepaths = filepaths if filepaths is not None else []
        #
        #self.dirpaths = map(os.path.abspath, self.dirpaths)
        #self.filepaths = map(os.path.abspath, self.filepaths)
        #
        #self.wildcard = wildcard if wildcard is not None else ""

        main_sizer = wx.BoxSizer( wx.VERTICAL )

        hsizer = wx.BoxSizer( wx.HORIZONTAL )

        st1 = wx.StaticText(self, -1, "Quantity:", wx.DefaultPosition, wx.DefaultSize, 0 )
        st1.Wrap( -1 )
        hsizer.Add(st1, 0, wx.ALIGN_CENTER_VERTICAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        plotter_choices = ["bands",]
        self.plotter_cbox = wx.ComboBox( self, -1, "bands", wx.DefaultPosition, wx.DefaultSize, plotter_choices, 0 )
        hsizer.Add( self.plotter_cbox, 0, wx.ALL, 5 )

        compare_button = wx.Button( self, -1, "Compare", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.Bind(wx.EVT_BUTTON, self.OnCompare, compare_button)
        hsizer.Add(compare_button, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        main_sizer.Add( hsizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.panel = FileCheckBoxPanel(self, filepaths)
        main_sizer.Add(self.panel, 1, wx.EXPAND, 5)

        self.SetSizer( main_sizer )
        self.Layout()

    def OnCompare(self, event):
        selected = self.panel.GetSelectedFilepaths()
        awx.PRINT("selected",selected)

        choice = self.plotter_cbox.GetValue()

        try:
            if choice == "bands":
                plotter = EBandsPlotter()
                for file in selected:
                    plotter.add_bands_from_file(file)
                plotter.plot()

            else:
                awx.showErrorMessage(self, message="No fallback registered for choice %s" % choice)

        except Exception as exc:
            awx.showErrorMessage(self, message=str(exc))


def main():
    import sys
    import os
    app = wx.App()

    filepaths = ["ciao.nc", "hello.out", "a"]
    dir = sys.argv[1]
    filepaths = [os.path.join(dir, f) for f in os.listdir(sys.argv[1])]
    #frame = wx.Frame(None, -1)
    #panel = FileCheckBoxPanel(frame, filepaths)
    frame = ComparisonFrame(None, filepaths=filepaths)
    app.SetTopWindow(frame)
    frame.Show()
    app.MainLoop()

if __name__ == "__main__":
    main()

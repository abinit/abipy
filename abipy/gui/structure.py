from __future__ import print_function, division, unicode_literals, absolute_import

import wx
import abipy.gui.awx as awx

from monty.string import is_string
from abipy.core.structure import Structure
from abipy.gui.editor import SimpleTextViewer


class StructureConverterFrame(wx.Frame):
    """
    This frame allows the user to convert the structure to different formats (CIF, POSCAR, ...).
    """
    def __init__(self, parent, obj, **kwargs):
        """
        Args:
            parent:
                Parent window.
            obj:
                Structure object, filename or object with a structure attribute.
        """
        super(StructureConverterFrame, self).__init__(parent, -1, **kwargs)
        self._get_structure(obj)

        panel = wx.Panel(self, id=-1)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(panel, -1, "Convert to:")
        label.Wrap(-1)

        # list of supported formats.
        formats = ["cif", "POSCAR", "cssr", "json"]
        self.format_choice = wx.Choice(panel, -1, choices=formats)
        self.format_choice.SetSelection(0)

        show_button = wx.Button(panel, -1, "Show", wx.DefaultPosition, wx.DefaultSize, 0)
        show_button.Bind(wx.EVT_BUTTON, self.OnShow)

        save_button = wx.Button(panel, -1, "Save", wx.DefaultPosition, wx.DefaultSize, 0)
        save_button.Bind(wx.EVT_BUTTON, self.OnSave)

        hsizer.Add(label, 0, wx.ALL, 5)
        hsizer.Add(self.format_choice, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        hsizer.Add(show_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        hsizer.Add(save_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        main_sizer.Add(hsizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        str_structure = wx.TextCtrl(panel, -1, str(self.structure), style=wx.TE_MULTILINE|wx.TE_LEFT|wx.TE_READONLY)
        main_sizer.Add(str_structure, 1, wx.ALL | wx.EXPAND, 5)

        panel.SetSizerAndFit(main_sizer)

    def _get_structure(self, obj):
        """Extract the structure from the input object."""
        return Structure.as_structure(obj)

    @property
    def format(self):
        """The format specified by the user."""
        return self.format_choice.GetStringSelection()

    def _convert(self):
        """Returns string with the structure converted in the user-specified format."""
        return self.structure.convert(fmt=self.format)

    def OnShow(self, event):
        s = self._convert()
        SimpleTextViewer(self, text=s, title=self.structure.formula).Show()

    def OnSave(self, event):
        save_dialog = wx.FileDialog(self, "Save %s file" % self.format, "", "",
                                    wildcard="%{format} files (*.{format})|*.{format}".format(format=self.format),
                                    style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if save_dialog.ShowModal() == wx.ID_CANCEL: return

        with open(save_dialog.GetPath(), "w") as fh:
            fh.write(self._convert())


def wxapp_structure_converter(obj):
    """
    Standalong WX application for structure conversion.

    Args:
        obj:
            Structure object, filename or object with a structure attribute.
    """
    app = awx.App()
    frame = StructureConverterFrame(None, obj)
    app.SetTopWindow(frame)
    frame.Show()

    return app

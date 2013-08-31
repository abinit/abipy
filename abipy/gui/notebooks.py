from __future__ import print_function, division

import wx

import abipy.gui.awx as awx

__all__ = [
    "ReadOnlyTextNotebookFrame",
]


class ReadOnlyTextNotebookFrame(awx.Frame):
    """
    This frame receives a list of strings and displays them in a notebook (read-only mode)
    """
    def __init__(self, parent, text_list, page_names, **kwargs):
        """
        Args:
            parent:
                Parent Widget.
            text_list:
                List of strings. Each string is displayed in its own page in the notebook.
            page_names:
                List of strings giving the name of the tab for each string in text_list.
        """
        if "title" not in kwargs:
            kwargs["title"] = self.__class__.__name__

        super(ReadOnlyTextNotebookFrame, self).__init__(parent, **kwargs)

        # Here we create a panel and a notebook on the panel
        import wx.lib.agw.flatnotebook as fnb
        nb_panel = awx.Panel(self)
        nb = fnb.FlatNotebook(nb_panel)

        # Add the pages to the notebook with the name to show on the tab
        if not isinstance(text_list, (list, tuple)):
            text_list = [text_list]

        if not isinstance(page_names, (list, tuple)):
            page_names = [page_names]

        assert len(page_names) == len(text_list)

        for page_name, text in zip(page_names, text_list):
            page = wx.TextCtrl(nb_panel, -1, text, style=wx.TE_MULTILINE|wx.TE_LEFT|wx.TE_READONLY)
            nb.AddPage(page, text=page_name)

        # Finally, put the notebook in a sizer for the nb_panel to manage the layout
        sizer = wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND)

        nb_panel.SetSizerAndFit(sizer)


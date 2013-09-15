from __future__ import print_function, division

import os
import wx
import abipy.gui.awx as awx

from abipy.tools.text import WildCard
from wx.py.editor import EditorFrame, EditorNotebookFrame, EditorNotebook, Editor

__all__ = [
    "TextNotebookFrame",
]


class SimpleTextViewer(awx.Frame):
    """Very simple frame that displays text (string )in read-only mode."""
    def __init__(self, parent, filename, **kwargs):
        super(SimpleTextViewer, self).__init__(parent, **kwargs)
        wx.TextCtrl(self, -1, text, style=wx.TE_MULTILINE|wx.TE_LEFT|wx.TE_READONLY)


#class AbinitEditorFrame(awx.Frame):
class AbinitEditorFrame(EditorFrame):
    def __init__(self, parent, filename, **kwargs):
        if "size" not in kwargs:
            kwargs["size"] = awx.FRAME_SIZE

        super(AbinitEditorFrame, self).__init__(parent, filename=filename, **kwargs)

        #self.editor = Editor(self)
        #with open(filename , "r") as fh:
        #    self.editor.setText(fh.read())
        #EditorFrame.__init__(self, parent=parent, filename=filename, **kwargs)

class TextNotebookFrame(awx.Frame):
#class TextNotebookFrame(EditorNotebookFrame):
    """
    This frame receives a list of strings and displays them in a notebook (read-only mode)
    """
    def __init__(self, parent, text_list, page_names, num_dirs=2, **kwargs):
        """
        Args:
            parent:
                Parent Widget.
            text_list:
                List of strings. Each string is displayed in its own page in the notebook.
            page_names:
                List of strings giving the name of the tab for each string in text_list.
            num_dirs:
                Maximum number of directories that will be shown in the tab.
        """
        if "title" not in kwargs:
            kwargs["title"] = self.__class__.__name__

        super(TextNotebookFrame, self).__init__(parent, **kwargs)

        # Add the pages to the notebook with the name to show on the tab
        if not isinstance(text_list, (list, tuple)):
            text_list = [text_list]

        if not isinstance(page_names, (list, tuple)):
            page_names = [page_names]

        assert len(page_names) == len(text_list)

        # Here we create a panel and a notebook on the panel
        import wx.lib.agw.flatnotebook as fnb
        nb_panel = awx.Panel(self)

        style = fnb.FNB_X_ON_TAB | fnb.FNB_NAV_BUTTONS_WHEN_NEEDED 
        nb = fnb.FlatNotebook(nb_panel, style=style)

        for page_name, text in zip(page_names, text_list):
            page = wx.TextCtrl(nb, -1, text, style=wx.TE_MULTILINE|wx.TE_LEFT|wx.TE_READONLY)

            if num_dirs > 0:
                tokens = page_name.split(os.path.sep)
                page_name = os.path.join(*tokens[-num_dirs:])

            nb.AddPage(page, text=page_name)

        # Finally, put the notebook in a sizer for the nb_panel to manage the layout
        sizer = wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND)

        nb_panel.SetSizerAndFit(sizer)

    @classmethod
    def from_files_and_dir(cls, parent, filenames=None, dirpath=None, walk=True, wildcard=""):
        """
        Static constructure that reads the content of the files/directory specified in input.

        Args:
            filenames:
                List of files to show in the botebook. Defaults to an empty list.
            dirpath:
                Directory to scan for additional files.
            walk:
                Used only if dirpath is not None.
                If True, we scan all the files contained within dirpath and 
                we add them to the list if their name match the regular expression 
                given in wildcard.
            wildcard:
                String with regular expressions separated by |.
                Only the files matching one of the regular expressions will be showed.
                example: wildcard="*.nc|*.txt" shows only the files whose extension is in ["nc", "txt"].
        """
        wildcard = WildCard(wildcard)

        if filenames is None:
            filenames = []

        filenames = wildcard.filter(filenames)

        if dirpath is not None:
            if not walk:
                filenames += wildcard.filter(os.listdir(dirpath))

            else:
                for root, dirnames, fnames in os.walk(dirpath):
                    for fname in fnames:
                        if wildcard.match(fname):
                            filenames.append(os.path.join(root, fname))

        #frame = EditorNotebookFrame(parent)
        #frame.notebook.DeletePage(0)
        #for fname in filenames:
        #    frame.bufferCreate(filename=fname)
        #return frame

        # Open the files and read the content in a string
        text_list = []
        for fname in filenames:
            with open(fname, "r") as fh:
                text_list.append(fh.read())

        # Instanciate the frame.
        return cls(parent, text_list, page_names=filenames)


def wxapp_showfiles(filenames=None, dirpath=None, walk=True, wildcard=None):
    """
    Standalone applications that reads the content of the files specified 
    in input and show them in a noteboox.

    Args:
        filenames:
            List of files to show in the botebook. Defaults to an empty list.
        dirpath:
            Directory to scan for additional files.
        walk:
            Used only if dirpath is not None.
            If True, we scan all the files contained within dirpath and 
            we add them to the list if their name match the regular expression 
            given in wildcard.
        wildcard:
            String with regular expressions separated by |.
            Only the files matching one of the regular expressions will be showed.
            example: wildcard="*.nc|*.txt" shows only the files whose extension is in ["nc", "txt"].

    Returns:
        wxpython application.
    """
    app = wx.App()
    frame = TextNotebookFrame.from_files_and_dir(None, filenames=filenames, dirpath=dirpath, walk=walk, wildcard=wildcard)
    frame.Show()
    return app


if __name__ == "__main__":
    app = wx.App()
    frame = EditorNotebookFrame()
    #frame = wx.Frame(None, -1)
    #notebook = EditorNotebook(frame)
    for filename in ["editor.py", "__init__.py"]:
        frame.bufferCreate(filename=filename)

    #Editor(frame)
    frame.Show()
    #frame = EditorFrame()
    #frame.bufferCreate(file.path)
    #frame = EditorNotebookFrame(filename=file.path)
    app.MainLoop()

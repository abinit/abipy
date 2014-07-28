"""
Utilities and helper functions.
"""
from __future__ import division, print_function

import wx

__all__ = [

]

def getApplicationConfigDirectory():
    """Returns the directory for application configuration information."""
    app = getApplication()
    return app.getConfigDirectory()


def getApplicationConfiguration():
    """Returns the application configuration object."""
    app = getApplication()
    return app.getConfiguration()


def beautifySize(size):
    """Returns a string representing the size in bytes in a more human readable way."""
    if size / 1073741824:
        return "%(size)0.2f GBytes" % {'size': (float(size) / 1073741824.0)}
    elif size / 1048576:
        return "%(size)0.2f MBytes" % {'size': (float(size) / 1048576.0)}
    elif size / 1024:
        return "%(size)0.2f KBytes" % {'size': (float(size) / 1024.0)}
    else:
        return "%(size)s Bytes" % {'size': size}


def addHorizontalSpaceTool(toolbar, width):
    """Adds horizontal space in a portable manner."""
    space = wx.StaticText(toolbar, -1, "")
    space.SetSize((width, -1))
    toolbar.AddControl(space)


def addLineSeparator(toolbar, height):
    """Adds a line separator to a toolbar."""
    addHorizontalSpaceTool(toolbar, 3)
    line = wx.StaticLine(toolbar, -1, style=wx.LI_VERTICAL)
    line.SetSize((-1, height))
    toolbar.AddControl(line)
    addHorizontalSpaceTool(toolbar, 3)


def getColumnText(list, index, col):
    """Gets the text from the specified column entry in a list."""
    item = list.GetItem(index, col)
    return item.GetText()


def getSelected(list):
    """Gets the selected items from a list object."""
    selected = []
    item = -1
    while 1:
        item = list.GetNextItem(item, wx.LIST_NEXT_ALL, wx.LIST_STATE_SELECTED)
        if item == -1:
            break
        selected.append(item)
    return selected


def setListColumnAlignment(list, col, align):
    """Sets the column alignment for a column in a list."""
    item_info = list.GetColumn(col)
    item_info.SetAlign(align)
    list.SetColumn(col, item_info)


def colorToTuple(color):
    """Coverts a color object to a tuple RGB tuple representing the color."""
    return (color.Red(), color.Green(), color.Blue())


def isStrOrUnicode(value):
    """Returns true if the value is either a string or a unicode type."""
    return isinstance(value, str) or isinstance(value, unicode)


def is_string(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        dummy = s + " "
        return True
    except TypeError:
        return False


def straceback(color=None):
    """
    Returns a string with the traceback.

    Use ANSII color formatting for output in terminal if color is not None.
    """
    import traceback
    s = traceback.format_exc()

    if color is not None:
        try:
            from termcolor import colored
            return colored(s, color)
        except ImportError:
            return s
    else:
        return s

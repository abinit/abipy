from __future__ import print_function, division

import traceback
import wx

__all__ = [
    "showErrorMessage",
    "askUser",
]


def _straceback():
    """Returns a string with the traceback."""
    return traceback.format_exc()


def showErrorMessage(parent, message=None):
    """
    Open a `MessageDialog` with an error message.
    If message is None, the python traceback is used.
    """
    if message is None:
        message = _straceback()

    dlg = wx.MessageDialog(parent, message=message,
                           caption='Error Message',
                           style=wx.OK | wx.ICON_INFORMATION | wx.ICON_ERROR | wx.STAY_ON_TOP
    )
    dlg.ShowModal()
    dlg.Destroy()


def askUser(parent, message):
    ask = wx.MessageDialog(parent, message)
    answer = ask.ShowModal() == wx.ID_OK
    ask.Destroy()
    return answer

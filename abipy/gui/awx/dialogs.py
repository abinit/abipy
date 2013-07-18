from __future__ import print_function, division

import traceback
import wx

__all__ = [
    "showErrorMessage",
    "askUser",
    "showLicense",
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


def showLicense(parent=None, codename=None):
    codename = "Abipy" if codename is None else codename

    license_text = """%(codename)s is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

%(codename)s is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have
received a copy of the GNU General Public License along with the code;
if not, write to the Free Software Foundation, Inc., 59 Temple Place,
Suite 330, Boston, MA  02111-1307  USA""" % {"codename": codename}

    dialog = License(parent, license_text)
    dialog.ShowModal()
    dialog.Destroy()

class License(wx.Dialog):
    def __init__(self, parent, license_text,  **kwargs):
        wx.Dialog.__init__ (self, parent, id=-1, title="License")

        vsizer = wx.BoxSizer( wx.VERTICAL )

        text = wx.TextCtrl( self, -1, license_text, style=wx.TE_MULTILINE|wx.TE_READONLY )
        vsizer.Add(text, 0, wx.ALL|wx.EXPAND, 5 )

        vsizer.Add(wx.Button( self, wx.ID_OK), 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.SetSizerAndFit(vsizer)

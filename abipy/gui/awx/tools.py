from __future__ import print_function, division

import os
import wx

__all__ = [
    "path_img",
    "makeAboutBox",
]

#class Logger(object):
#    def WriteText(self, text):
#        if text[-1:] == '\n':
#            text = text[:-1]
#        wx.LogMessage(text)
#    write = WriteText

def path_img(filename):
    """Returns the absolute path of an image."""
    dirname = os.path.dirname(__file__)
    return os.path.join(dirname, "images", filename)


def makeAboutBox(codename, version, description, developers, website=None, icon_path=None):

    licence = """%(codename)s is free software; you can redistribute 
it and/or modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation; either version 2 of the License, 
or (at your option) any later version.

%(codename)s is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details. You should have 
received a copy of the GNU General Public License along with File Hunter; 
if not, write to the Free Software Foundation, Inc., 59 Temple Place, 
Suite 330, Boston, MA  02111-1307  USA""" % {"codename": codename}

    info = wx.AboutDialogInfo()

    if icon_path is not None:
        info.SetIcon(wx.Icon(icon_path, wx.BITMAP_TYPE_PNG))

    info.SetName(codename)
    info.SetVersion(version)
    info.SetDescription(description)

    if website is not None:
        info.SetWebSite(website)

    info.SetLicence(licence)

    if isinstance(developers, str): 
        developers = [developers,]

    for dev in developers:
        info.AddDeveloper(dev)

    wx.AboutBox(info)

"""Core objects and helper functions."""
from __future__ import print_function, division

import os
import sys
import wx

from monty.string import list_strings
#import abipy.tools.decorators as dec


#import logging
#logger = logging.getLogger(__name__)


__all__ = [
    "path_img",
    "makeAboutBox",
    "verbose",
    "Panel",
    "Frame",
    "get_width_height",
]

_DEBUG = True
_DEBUG = False


#if _DEBUG:
#    verbose = dec.verbose
#
#else:
def verbose(func):
    return func


#class Error(Exception):
#    """Base class for exceptions raised by awx library"""


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
received a copy of the GNU General Public License along with the code; 
if not, write to the Free Software Foundation, Inc., 59 Temple Place, 
Suite 330, Boston, MA  02111-1307  USA""" % {"codename": codename}

    # Make a template for the description 
    #desc = "\n".join(["\nwxPython Cookbook Chapter 5\n",
    #                  "Platform Info: (%s,%s)",
    #                  "License: Public Domain"])

    ## Get the platform information 
    #py_version = [sys.platform, ", python ", sys.version.split()[0]] 
    #platform = list(wx.PlatformInfo[1:])
    #platform[0] += (" " + wx.VERSION_STRING) 
    #wx_info = ", ".join(platform)

    #info.SetDescription(desc % (py_version, wx_info))

    info = wx.AboutDialogInfo()

    if icon_path is not None:
        info.SetIcon(wx.Icon(icon_path, wx.BITMAP_TYPE_PNG))

    info.SetName(codename)
    info.SetVersion(version)
    info.SetDescription(description)
    info.SetCopyright('(C) Abipy group')

    if website is not None:
        info.SetWebSite(website)

    info.SetLicence(licence)

    for dev in list_strings(developers):
        info.AddDeveloper(dev)

    wx.AboutBox(info)


def get_width_height(window, string, pads=(10, 10)):
    """
    Returns the width and the height (in pixels) of a string used in window.
    Returns are padded with pads

    See http://docs.wxwidgets.org/2.8/wx_wxdc.html#wxdcgettextextent
    """
    f = window.GetFont()
    dc = wx.WindowDC(window)
    dc.SetFont(f)
    width, height = dc.GetTextExtent(string)
    return width + pads[0], height + pads[1]


class MyWindow(object):
    """
    Mixin class providing helper functions and commonly used callbacks.

    Attributes:
        HELP_MSG:
            string with a short help (will be displayed in MessageDialog in onHelp callback)
    """
    HELP_MSG = "No help available"

    def getParentWithType(self, cls):
        """
        Returns the first parent window of type cls.

        Raises:
            RuntimeError if we have reached the head of the linked list.
        """
        parent = self.GetParent() 

        while True: 
            if parent is None:
                raise RuntimeError("Cannot find parent with class %s, reached None parent!" % cls)

            if isinstance(parent, cls): 
                return parent
            else:
                parent = parent.GetParent()

    def onHelp(self, event):
        """Short help."""
        dialog = wx.MessageDialog(self, self.HELP_MSG, " Quick Reference", wx.OK | wx.ICON_INFORMATION)
        dialog.ShowModal()
        dialog.Destroy()


class Panel(wx.Panel, MyWindow):
    def __init__(self, parent, *args, **kwargs):
        super(Panel, self).__init__(parent, *args, **kwargs)


class Frame(wx.Frame, MyWindow):
    """Base class for frames."""

    def __init__(self, parent, *args, **kwargs):
        if "title" not in kwargs:
            kwargs["title"] = self.__class__.__name__
        
        super(Frame, self).__init__(parent, *args, **kwargs)

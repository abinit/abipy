from __future__ import print_function, division

import os
import sys
import warnings
import wx

import abipy.tools.decorators as dec

from abipy.tools.text import list_strings

#import logging
#logger = logging.getLogger(__name__)


__all__ = [
    "path_img",
    "makeAboutBox",
    "Error",
    "verbose",
    "Panel",
    "Frame",
    "FRAME_SIZE",
    "get_width_height",
]

_DEBUG = True
_DEBUG = False

class Error(Exception):
    """Base class for exceptions"""

if _DEBUG:
    verbose = dec.verbose

else:
    def verbose(func):
        return func


FRAME_SIZE = (800, 600)

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

    info = wx.AboutDialogInfo()

    if icon_path is not None:
        info.SetIcon(wx.Icon(icon_path, wx.BITMAP_TYPE_PNG))

    info.SetName(codename)
    info.SetVersion(version)
    info.SetDescription(description)

    if website is not None:
        info.SetWebSite(website)

    info.SetLicence(licence)

    for dev in list_strings(developers):
        info.AddDeveloper(dev)

    wx.AboutBox(info)


class MyWindow(object):
    """Mixin class providing helper functions."""

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


class Panel(wx.Panel, MyWindow):
    def __init__(self, parent, *args, **kwargs):
        super(Panel, self).__init__(parent, *args, **kwargs)


class Frame(wx.Frame, MyWindow):
    def __init__(self, parent, *args, **kwargs):
        if "size" not in kwargs:
            kwargs["size"] = FRAME_SIZE

        if "title" not in kwargs:
            kwargs["title"] = self.__class__.__name__
        
        super(Frame, self).__init__(parent, *args, **kwargs)


def get_width_height(window, string, pads=(10,10)):
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




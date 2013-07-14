"""This module ..."""
from __future__ import print_function, division

import wx
import collections

import wx.lib.dialogs as wxdg 
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx
from abipy import ncfile_subclass_from_filename, abiopen

from abipy.iotools.files import NcDumper 

from abipy.iotools.files import AbinitNcFile
from abipy import WFK_File, SIGRES_File
#from abipy.phonons import PHBST_File


__all__ = [
    "popupmenu_for_filename",
]

def popupmenu_for_filename(filename, parent=None):
    """
    Factory function that returns the appropriate popup menu 

    Args:
        parent:
            The parent wx window. Used to connect the children created
            by the popupmenu to the calling window.
    """
    menu = PopupMenu.from_filename(filename)

    menu.set_target(filename)
    if parent is not None:
        menu.set_parent(parent)
    return menu

#--------------------------------------------------------------------------------------------------
# Callbacks

def showNcdumpMessage(parent, filepath):
    """Open a dialog with the output of ncdump."""
    caption = "ncdump output for file %s" % filepath
    msg = NcDumper().dump(filepath)
    style = wx.DEFAULT_FRAME_STYLE
    wxdg.ScrolledMessageDialog(parent, msg, caption=caption, size=(600, 600), style=style).Show()

#def showStructure(parent, filepath):
#    ncfile = abiopen(filepath)
#    structure = ncfile.get_structure() 
#    visu = structure.visualize("xcrysden")
#    visu()

#--------------------------------------------------------------------------------------------------

class PopupMenu(wx.Menu):
    """
    Base class for popup menus. `A PopupMenu` has a list of callback functions
    indexed by the menu title and a list of `AbinitNcFile` associated to it.
    The signature of the callback function is func(filename, parent) where 
    filename is the name of the file selected in the Widget and parent is the wx
    Window that will become the parent of the new frame created by the callback.

    How to subclass PopupMenu:

        1. Define a new class that inherits from PopupMenu. 

        2. Define the callbacks in the class variable MENU_TITLES.
           Use OrderedDict to have a fixed ordering of the labels.

        3. Define the class variable HANDLED_NCFILES with the list
           of `AbinitNcFile` subclasses associated to the popupmenu.

        4. Done (most of the work is indeed done in the base class and in 
           the factory function popupmenu_for_filename.
    """
    MENU_TITLES = collections.OrderedDict([
        ("ncdump", showNcdumpMessage),
        #("structure", showStructure),
        #"Properties":
    ])

    HANDLED_NCFILES = []

    def __init__(self):
        super(PopupMenu, self).__init__()
        self._make_menu()

    @staticmethod
    def from_filename(filename):
        """
        Static factory function that instanciates the appropriate subclass of `PopupMenu`
        """
        # Find the AbinitNcFile subclass associated to files.
        ncfile_class = ncfile_subclass_from_filename(filename)
        print("ncfile_class:",ncfile_class)

        # Check whether a subclass handles this file..
        # Fallback to a simple PopupMenu if no match.
        for popsubc in PopupMenu.__subclasses__():
            if popsubc.handle_ncfile_class(ncfile_class):
                return popsubc()
        else:
            return PopupMenu()

    @classmethod
    def handle_ncfile_class(cls, ncfile_class):
        """True if the popupmenu is associated to ncfile_class."""
        return ncfile_class in cls.HANDLED_NCFILES

    def _make_menu(self):
        """Build the menu taking into account the options of the superclasses."""
        base_classes = list(self.__class__.__bases__) + [self.__class__]

        assert not hasattr(self, "menu_title_by_id")
        assert not hasattr(self, "menu_titles")
        self.menu_title_by_id, self.menu_titles = {}, {}

        for cls in base_classes:
            print(cls)
            try:
                menus = cls.MENU_TITLES
            except AttributeError as exc:
                print("exc ",exc," for cls", cls)
                continue

            self.menu_titles.update(menus)
            for title in menus:
                self.menu_title_by_id[wx.NewId()] = title
                                                                  
        print(self.menu_title_by_id)
                                                                  
        for (id, title) in self.menu_title_by_id.items():
            self.Append(id, title)
            # registers menu handlers with EVT_MENU, on the menu.
            wx.EVT_MENU(self, id, self.OnMenuSelection)

    def set_parent(self, parent):
        """Set the parent window."""
        self._parent = parent

    @property
    def parent(self):
        """Returns the parent window."""
        try:
            return self._parent
        except AttributeError:
            print("Warning Popup menu has not parent")
            return None

    def set_target(self, target):
        """Set the target of the callback."""
        self._target = target

    @property
    def target(self):
        """The target of the callback."""
        try:
            return self._target
        except AttributeError:
            return None

    def _get_callback(self, title):
        return self.menu_titles[title]

    def OnMenuSelection(self, event):
        title = self.menu_title_by_id[event.GetId()]
        callback = self._get_callback(title)
        print("Calling callback %s on target %s" % (callback, self.target))

        try:
            callback(parent=self.parent, filepath=self.target)
        except:
            awx.showErrorMessage(parent=self.parent)


class EbandsPopupMenu(PopupMenu):
    """Popup menu for ncfiles that contain the electron band structure."""
    MENU_TITLES = collections.OrderedDict([
        ("eDos",  ewx.showElectronDosFrame),
        ("ePlot", ewx.showElectronBandsPlot),
    ])

    HANDLED_NCFILES = [WFK_File] 


def showQPData(parent, filepath):
    sigres = abiopen(filepath)
    qps_spin = sigres.get_allqps()
    assert len(qps_spin) == 1
    qps_spin[0].plot_qps_vs_e0()


class SigResPopupMenu(PopupMenu):
    """Popup menu for SIGRES files."""
    MENU_TITLES = collections.OrderedDict([
        ("qpDataPlot", showQPData),
    ])

    HANDLED_NCFILES = [SIGRES_File] 


#def showEXCMDF(parent, filepath):
#    mdf_file = abiopen(filepath)
#    mdf_file.read
#
#class MDFPopupMenu(PopupMenu):
#    """Popup menu for MDF files."""
#    MENU_TITLES = collections.OrderedDict([
#        ("excmdfPlot", showEXCMDF),
#    ])
#
#    HANDLED_NCFILES = [MDF_File] 

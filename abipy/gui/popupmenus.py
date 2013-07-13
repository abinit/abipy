from __future__ import print_function, division

import sys
import os
import wx

import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy.abifiles import abiopen
from abipy.electrons import ElectronBands
from abipy.waves import WFK_File
from wx.lib.agw.floatspin import FloatSpin
from wxmplot import PlotApp, PlotFrame

__all__ = [
    "popupmenu_from_ext",
]

def popupmenu_from_ext(ncfile):
    menu = BasePopupMenu()
    menu.add_target(ncfile)
    return menu

class BasePopupMenu(wx.Menu):
    MENU_TITLES = [ 
        "Properties",
        #"Ncdump",
        "DOS",
    ]

    def __init__(self, *args, **kwargs):
        super(BasePopupMenu, self).__init__()

        menu_title_by_id = {}
        for title in self.MENU_TITLES:
            menu_title_by_id[wx.NewId()] = title

        if not hasattr(self, "menu_title_by_id"):
            self.menu_title_by_id = {}

        self.menu_title_by_id.update(menu_title_by_id)

        for (id, title) in self.menu_title_by_id.items():
            self.Append(id, title)
            # registers menu handlers with EVT_MENU, on the menu.
            wx.EVT_MENU(self, id, self.OnMenuSelection)

    def add_target(self, target):
        self._target = target

    @property
    def target(self):
        try:
            return self._target
        except AttributeError:
            return None

    def OnMenuSelection(self, event):
        operation = self.menu_title_by_id[event.GetId()]
        print("Perform operation %s on target %s" % (operation, self.target))
        try:
            ewx.showElectronDosFrame(self.target)
        except:
            awx.showErrorMessage(self)


"""Collections of Popup menus associated to filenames."""
import wx
import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from collections import OrderedDict
from abipy.abilab import abifile_subclass_from_filename, abiopen
from abipy.core.mixins import NcDumper, get_filestat
from abipy.abio.outputs import AbinitLogFile, AbinitOutputFile
from abipy.iotools.visualizer import Visualizer
from abipy.waves import WfkFile
from abipy.electrons import SigresFile, GsrFile
from abipy.electrons.bse import MdfFile
from abipy.gui.events import AbinitEventsFrame
from abipy.gui.timer import AbinitTimerFrame
from abipy.gui.editor import SimpleTextViewer, MyEditorFrame
from abipy.gui.wxncview import NcViewerFrame


__all__ = [
    "popupmenu_for_filename",
]


def popupmenu_for_filename(parent, filename):
    """
    Factory function that returns the appropriate popup menu.

    Args:
        parent:
            The parent wx window. Used to connect the children created
            by the popupmenu to the calling window.
    """
    menu = PopupMenu.from_filename(filename)

    if menu is not None:
        menu.set_target(filename)
        menu.set_parent(parent)

    return menu

#--------------------------------------------------------------------------------------------------
# Callbacks


def showNcdumpMessage(parent, filepath):
    """Open a dialog with the output of ncdump."""
    title = "ncdump output for file %s" % filepath
    text = NcDumper().dump(filepath)
    # TODO: Get a decent wxpython editor somewhere
    #SimpleTextViewer(parent, text, title=title).Show()
    MyEditorFrame.from_text(parent, text, title=title).Show()

def openWxncview(parent, filepath):
    """Open a dialog with the output of ncdump."""
    title = "wxncview:  %s" % filepath
    text = NcDumper().dump(filepath)
    NcViewerFrame(parent, filepath, title=title).Show()


def showFileStat(parent, filepath):
    """Open a dialog reporting file stats."""
    caption = "Info on file %s" % filepath
    stat_dict = get_filestat(filepath)
    msg = str(stat_dict)
    style = wx.DEFAULT_FRAME_STYLE
    wxdg.ScrolledMessageDialog(parent, msg, caption=caption, size=(600, 600), style=style).Show()


def showAbinitEventsFrame(parent, filepath):
    """Open a dialog reporting file stats."""
    AbinitEventsFrame(parent, filepath).Show()


def showAbinitTimerFrame(parent, filepath):
    """Open a dialog reporting file stats."""
    try:
        frame = AbinitTimerFrame(parent, filepath)
        frame.Show()
    except awx.Error as exc:
        awx.showErrorMessage(parent, str(exc))


def showStructure(parent, filepath):
    ncfile = abiopen(filepath)
    visu_classes = Visualizer.get_available(ext="xsf")
    if not visu_classes:
        print("Not visualizer found for extension xsf")
        return
    vname = visu_classes[0].name

    visu = ncfile.structure.visualize(vname)

    thread = awx.WorkerThread(parent, target=visu)
    thread.start()



class PopupMenu(wx.Menu):
    """
    Base class for popup menus. `A PopupMenu` has a list of callback functions
    indexed by the menu title. The signature of the callback function is func(parent, filename) where
    filename is the name of the file selected in the Widget and parent is the wx
    Window that will become the parent of the new frame created by the callback.
    """
    MENU_TITLES = OrderedDict([
    ])

    HANDLED_FILES = []

    def __init__(self):
        super(PopupMenu, self).__init__()
        self._make_menu()

    @staticmethod
    def from_filename(filename):
        """
        Static factory function that instanciates the appropriate subclass of `NcFilePopupMenu`
        Returns None if the extesion of filename is not supported.
        """
        # Find the AbinitNcFile subclass associated to files.
        try:
            file_class = abifile_subclass_from_filename(filename)
        except KeyError:
            if filename.endswith(".nc"): NcFilePopupMenu()
            return None

        # Check whether a subclass handles this file..
        # Fallback to a simple PopupMenu if no match.
        def allsubclasses(cls):
            """Returns the set of subclasses of cls."""
            children = [cls]
            for sc in cls.__subclasses__():
                if sc.__subclasses__():
                    for k in sc.__subclasses__():
                        children.extend(allsubclasses(k))
                else:
                    children.append(sc)
            return set(children)

        for cls in allsubclasses(PopupMenu):
            if cls.handle_file_class(file_class):
                return cls()
        else:
            if filename.endswith(".nc"): NcFilePopupMenu()
            return PopupMenu()

    @classmethod
    def handle_file_class(cls, file_class):
        """True if the popupmenu is associated to file_class."""
        return file_class in cls.HANDLED_FILES

    def _make_menu(self):
        """Build the menu taking into account the options of the superclasses."""
        base_classes = list(self.__class__.__bases__) + [self.__class__]
        base_classes.reverse()

        assert not hasattr(self, "menu_title_by_id")
        assert not hasattr(self, "menu_titles")
        self.menu_title_by_id, self.menu_titles = OrderedDict(), OrderedDict()

        for cls in base_classes:
            try:
                menus = cls.MENU_TITLES
            except AttributeError as exc:
                awx.WARNING("exc ",exc," for cls", cls)
                continue

            self.menu_titles.update(menus)

            for title in menus:
                self.menu_title_by_id[wx.NewId()] = title

            # Add sentinel for Menu separator.
            self.menu_title_by_id["separator_" + str(len(self.menu_titles))] = None

        for (id, title) in self.menu_title_by_id.items():
            if title is None:
                sepid = int(id.split("_")[-1])
                if sepid != len(self.menu_titles):
                    self.AppendSeparator()
            else:
                # Register menu handlers with EVT_MENU, on the menu.
                self.Append(id, title)
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
            awx.WARNING("Popup menu doesn't have parent")
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
        #print("Calling callback %s on target %s" % (callback, self.target))
        try:
            callback(parent=self.parent, filepath=self.target)
        except:
            awx.showErrorMessage(parent=self.parent)


class AbinitTextFilePopupMenu(PopupMenu):
    """
    """
    MENU_TITLES = OrderedDict([
        ("events",     showAbinitEventsFrame),
        ("properties", showFileStat),
        ("timer",      showAbinitTimerFrame),
    ])

    HANDLED_FILES = [AbinitLogFile, AbinitOutputFile]


class NcFilePopupMenu(PopupMenu):
    """
    Base class for popup menus. `A PopupMenu` has a list of callback functions
    indexed by the menu title and a list of `AbinitNcFile` associated to it.
    The signature of the callback function is func(filename, parent) where
    filename is the name of the file selected in the Widget and parent is the wx
    Window that will become the parent of the new frame created by the callback.

    How to subclass PopupMenu:

        1. Define a new class that inherits from NcFilePopupMenu.

        2. Define the callbacks in the class variable MENU_TITLES.
           Use OrderedDict to have a fixed ordering of the labels.

        3. Define the class variable HANDLED_FILES with the list of
           `AbinitNcFile` subclasses associated to the popupmenu.

        4. Done (most of the work is indeed done in the base class and in
           the factory function popupmenu_for_filename.
    """
    MENU_TITLES = OrderedDict([
        ("structure",  showStructure),
        ("ncdump",     showNcdumpMessage),
        ("wxncview",   openWxncview),
        ("properties", showFileStat),
    ])

    HANDLED_FILES = []


class EbandsPopupMenu(NcFilePopupMenu):
    """Popup menu for Netcdf files that contain the electron band structure."""
    MENU_TITLES = OrderedDict([
        ("ePlot", ewx.showElectronBandsPlot),
        ("eDos",  ewx.showElectronDosFrame),
        ("eJdos", ewx.showElectronJdosFrame),
    ])

    HANDLED_FILES = [WfkFile, GsrFile]


def showQPData(parent, filepath):
    sigres = abiopen(filepath)
    qps_spin = sigres.qplist_spin
    assert len(qps_spin) == 1
    qps_spin[0].plot_qps_vs_e0()


class SigResPopupMenu(NcFilePopupMenu):
    """Popup menu for SIGRES files."""
    MENU_TITLES = OrderedDict([
        ("qpDataPlot", showQPData),
    ])

    HANDLED_FILES = [SigresFile]


def showEXCMDF(parent, filepath):
    mdf_file = MdfFile(filepath)
    mdf_file.plot_mdfs()


class MDFPopupMenu(NcFilePopupMenu):
    """Popup menu for MDF files."""
    MENU_TITLES = OrderedDict([
        ("mdfPlot", showEXCMDF),
    ])

    HANDLED_FILES = [MdfFile]


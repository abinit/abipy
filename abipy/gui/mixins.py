from __future__ import print_function, division

import abc
import os
import wx
import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from pymatgen.util.io_utils import which 
from abipy.iotools.visualizer import Visualizer
from abipy.iotools.files import NcDumper
from abipy.gui.structure import StructureConverterFrame
from abipy.gui.converter import ConverterFrame


class Has_Structure(object):
    __metaclass__ = abc.ABCMeta

    # Structure Menu ID's
    ID_STRUCT_CONVERT = wx.NewId()
    ID_STRUCT_VISUALIZE = wx.NewId()
    ID_STRUCT_SHOWBZ = wx.NewId()

    @abc.abstractproperty
    def structure(self):
        """Structure object."""

    def CreateStructureMenu(self):
        """Creates the structure menu."""
        menu = wx.Menu()
        menu.Append(self.ID_STRUCT_CONVERT, "Convert", "Convert structure data to cif, POSCAR ...")
        self.Bind(wx.EVT_MENU, self.OnStructureConvert, id=self.ID_STRUCT_CONVERT)

        menu.Append(self.ID_STRUCT_SHOWBZ, "Show BZ", "Visualize the first Brillouin zone with matplotlib.")
        self.Bind(wx.EVT_MENU, self.OnStructureShowBz, id=self.ID_STRUCT_SHOWBZ)

        # Make sub-menu with the list of supported visualizers.
        visu_menu = wx.Menu()
        self._id2visuname = {}

        available_visus = [visu.name for visu in Visualizer.get_available()]

        for visu_name in available_visus:
            _id =  wx.NewId()
            visu_menu.Append(_id, visu_name)
            self._id2visuname[_id] = visu_name
            self.Bind(wx.EVT_MENU, self.OnStructureVisualize, id=_id)
                                                                                                                     
        menu.AppendMenu(-1, 'Visualize', visu_menu)
                                                                                                                     
        return menu

    def OnStructureConvert(self, event):
        """Processes a connect initiation event. is initiated."""
        StructureConverterFrame(self, self.structure).Show()

    def OnStructureVisualize(self, event):
        """"Call visualizer to visualize the crystalline structure."""
        visualizer = self._id2visuname[event.GetId()]
        print("eventID", event.GetId(), "map", self._id2visuname, "visualizer", visualizer)

        try:
            visu = self.structure.visualize(visualizer)
                                                                                            
            thread = awx.WorkerThread(self, target=visu)
            thread.start()
                                                                                            
        except:
            awx.showErrorMessage(self)

    def OnStructureShowBz(self, event):
        """"Visualize the Brillouin zone with matplotlib."""
        self.structure.show_bz()


class Has_Ebands(object):
    __metaclass__ = abc.ABCMeta

    # Ebands Menu ID's
    ID_EBANDS_PLOT = wx.NewId()
    ID_EBANDS_DOS = wx.NewId()
    ID_EBANDS_JDOS = wx.NewId()

    @abc.abstractproperty
    def ebands(self):
        """`Electron Bands object."""

    def CreateEbandsMenu(self):
        """Creates the ebands menu."""
        menu = wx.Menu()
        menu.Append(self.ID_EBANDS_PLOT, "Plot ebands", "Plot electron bands with matplotlib")
        self.Bind(wx.EVT_MENU, self.OnEbandsPlot, id=self.ID_EBANDS_PLOT)
        menu.Append(self.ID_EBANDS_DOS, "DOS", "Compute the electron DOS")
        self.Bind(wx.EVT_MENU, self.OnEbandsDos, id=self.ID_EBANDS_DOS)
        menu.Append(self.ID_EBANDS_JDOS, "JDOS", "Compute the electron Joint DOS")
        self.Bind(wx.EVT_MENU, self.OnEbandsJdos, id=self.ID_EBANDS_JDOS)
                                                                                                    
        return menu

    def OnEbandsPlot(self, event):
        """Plot band energies with matplotlib."""
        self.ebands.plot()
                                                                  
    def OnEbandsDos(self, event):
        """Open Frame for the computation of the DOS."""
        ewx.ElectronDosFrame(self, bands=self.ebands).Show()
                                                                  
    def OnEbandsJdos(self, event):
        """Open Frame for the computation of the JDOS."""
        ewx.ElectronJdosFrame(self, bands=self.ebands).Show()


#class Has_Kpoints(object):

class Has_Tools(object):

    # Tools Menu ID's
    ID_TOOLS_UNIT_CONVERTER = wx.NewId()

    def CreateToolsMenu(self):
        """Creates the ebands menu."""
        menu = wx.Menu()
        menu.Append(self.ID_TOOLS_UNIT_CONVERTER, "Unit converter", "Unit Converter")
        self.Bind(wx.EVT_MENU, self.OnTools_UnitConverter, id=self.ID_TOOLS_UNIT_CONVERTER)

        return menu

    def OnTools_UnitConverter(self, event):
        ConverterFrame(self).Show()


class Has_Netcdf(object):
    __metaclass__ = abc.ABCMeta

    # Netcdf Menu ID's
    ID_NETCDF_NCDUMP = wx.NewId()
    ID_NETCDF_NCVIEW = wx.NewId()
    ID_NETCDF_WXNCVIEW = wx.NewId()

    @abc.abstractproperty
    def nc_filepath(self):
        """String with the absolute path of the netcdf file."""

    def CreateNetcdfMenu(self):
        """Creates the ebands menu."""
        menu = wx.Menu()
        menu.Append(self.ID_NETCDF_NCDUMP, "ncdump", "Show the output of ncdump")
        self.Bind(wx.EVT_MENU, self.OnNetcdf_NcDump, id=self.ID_NETCDF_NCDUMP)

        menu.Append(self.ID_NETCDF_NCVIEW, "ncview", "Call ncview")
        self.Bind(wx.EVT_MENU, self.OnNetcdf_NcView, id=self.ID_NETCDF_NCVIEW)

        menu.Append(self.ID_NETCDF_WXNCVIEW, "wxncview", "Call wxncview")
        self.Bind(wx.EVT_MENU, self.OnNetcdf_WxNcView, id=self.ID_NETCDF_WXNCVIEW)
                                                                                            
        return menu

    def OnNetcdf_NcDump(self, event):
        """Call ncdump and show results in a dialog."""
        s = NcDumper().dump(self.nc_filepath)
        caption = "ncdump output for %s" % self.nc_filepath
        wxdg.ScrolledMessageDialog(self, s, caption=caption, style=wx.MAXIMIZE_BOX).Show()

    def OnNetcdf_NcView(self, event):
        """Call ncview in an subprocess."""
        if which("ncview") is None:
            awx.showErrorMessage(self, "Cannot find ncview in $PATH")
            return 

        def target():
            os.system("ncview %s" % self.nc_filepath)

        try:
            thread = awx.WorkerThread(self, target=target)
            thread.start()
        except:
            awx.showErrorMessage(self)

    def OnNetcdf_WxNcView(self, event):
        """Open wxncview frame."""
        from abipy.gui.wxncview import NcViewFrame
        NcViewFrame(None, filepaths=self.nc_filepath).Show()


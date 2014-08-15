from __future__ import print_function, division

import abc
import os
import wx
import wx.lib.dialogs as wxdg
import abipy.gui.awx as awx
import abipy.gui.electronswx as ewx

from abipy.tools import which 
from abipy.iotools.visualizer import Visualizer
from abipy.iotools.files import NcDumper
from abipy.electrons.ebands import ElectronBandsPlotter, ElectronDosPlotter
from abipy.gui.structure import StructureConverterFrame
from abipy.gui.converter import ConverterFrame
from abipy.gui.wxncview import NcViewerFrame


class Has_Structure(object):
    """
    Mixin class that provides a menu and callbacks 
    for analyzing the crystalline structure.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def structure(self):
        """Structure object."""

    def CreateStructureMenu(self):
        """Creates the structure menu."""
        # Structure Menu ID's
        self.ID_STRUCT_CONVERT = wx.NewId()
        self.ID_STRUCT_VISUALIZE = wx.NewId()
        self.ID_STRUCT_SHOWBZ = wx.NewId()

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
            _id = wx.NewId()
            visu_menu.Append(_id, visu_name)
            self._id2visuname[_id] = visu_name
            self.Bind(wx.EVT_MENU, self.OnStructureVisualize, id=_id)
                                                                                                                     
        menu.AppendMenu(-1, 'Visualize', visu_menu)
        return menu

    def OnStructureConvert(self, event):
        """Open new frame that allows the user to convert the structure."""
        StructureConverterFrame(self, self.structure).Show()

    def OnStructureVisualize(self, event):
        """"Call the visualizer to visualize the crystalline structure."""
        visu_name = self._id2visuname[event.GetId()]

        try:
            visu = self.structure.visualize(visu_name)
                                                                                            
            thread = awx.WorkerThread(self, target=visu)
            thread.start()

        except Exception:
            awx.showErrorMessage(self)

    def OnStructureShowBz(self, event):
        """"Visualize the Brillouin zone with matplotlib."""
        self.structure.show_bz()


class Has_Ebands(object):
    """
    Mixin class that provides a menu and callbacks for analyzing electron bands.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def ebands(self):
        """`ElectronBands` object."""

    def CreateEbandsMenu(self):
        """Creates the ebands menu."""
        # Ebands Menu ID's
        self.ID_EBANDS_GETINFO = wx.NewId()
        self.ID_EBANDS_PLOT = wx.NewId()
        self.ID_EBANDS_DOS = wx.NewId()
        self.ID_EBANDS_JDOS = wx.NewId()
        #self.ID_EBANDS_FSURF = wx.NewId()
        #self.ID_EBANDS_SCISSORS = wx.NewId()

        menu = wx.Menu()

        menu.Append(self.ID_EBANDS_GETINFO, "Get Info", "Show info on the band structure")
        self.Bind(wx.EVT_MENU, self.OnEbandsGetInfo, id=self.ID_EBANDS_GETINFO)

        menu.Append(self.ID_EBANDS_PLOT, "Plot ebands", "Plot electron bands with matplotlib")
        self.Bind(wx.EVT_MENU, self.OnEbandsPlot, id=self.ID_EBANDS_PLOT)

        menu.Append(self.ID_EBANDS_DOS, "DOS", "Compute the electron DOS")
        self.Bind(wx.EVT_MENU, self.OnEbandsDos, id=self.ID_EBANDS_DOS)

        menu.Append(self.ID_EBANDS_JDOS, "JDOS", "Compute the electron Joint DOS")
        self.Bind(wx.EVT_MENU, self.OnEbandsJdos, id=self.ID_EBANDS_JDOS)

        # TODO
        #menu.Append(self.ID_EBANDS_FSURF, "Fermi surface", "Visualize the Fermi surface with Xcrysden")
        #self.Bind(wx.EVT_MENU, self.OnFermiSurface, id=self.ID_EBANDS_FSURF)

        # TODO
        #menu.Append(self.ID_EBANDS_SCISSORS, "Apply scissors", "Apply a scissors operator")
        #self.Bind(wx.EVT_MENU, self.OnApplyScissors, id=self.ID_EBANDS_SCISSORS)

        return menu

    def OnEbandsGetInfo(self, event):
        """Shows info on the bandstructure."""
        s = self.ebands.info
        caption = "Ebands info" 
        wxdg.ScrolledMessageDialog(self, s, caption=caption, style=wx.MAXIMIZE_BOX).Show()

    def OnEbandsPlot(self, event):
        """Plot band energies with matplotlib."""
        self.ebands.plot()
                                                                  
    def OnEbandsDos(self, event):
        """Open Frame for the computation of the DOS."""
        ewx.ElectronDosFrame(self, bands=self.ebands).Show()
                                                                  
    def OnEbandsJdos(self, event):
        """Open Frame for the computation of the JDOS."""
        ewx.ElectronJdosFrame(self, bands=self.ebands).Show()

    def OnFermiSurface(self, event):
        """Visualize the Fermi surface with Xcrysden."""
        try:
            visu = self.ebands.export_bxsf(".bxsf")
                                                                                            
            thread = awx.WorkerThread(self, target=visu)
            thread.start()
                                                                                            
        except Exception:
            awx.showErrorMessage(self)

    #def OnApplyScissors(self, event):
    #    """
    #    Read the scissors operator from a pickle file, apply it to the electron bands and save the results.
    #    """
    #    Get the scissors operator from file
    #    Apply the scissors.
    #    new_ebands = self.ebands.apply_scissors(self, scissors)
    #    new_ebands.plot()
    #    new_ebands.pickle_dump()


class Has_MultipleEbands(Has_Ebands):
    """
    Mixin class that provides a menu and callbacks 
    for analyzing and comparing multiple electron bands.
    """
    __metaclass__ = abc.ABCMeta

    def CreateEbandsMenu(self):
        """Creates the ebands menu."""
        menu = super(Has_MultipleEbands, self).CreateEbandsMenu()
        menu.AppendSeparator()

        # Multiple Ebands Menu ID's
        self.ID_MULTI_EBANDS_PLOT = wx.NewId()
        self.ID_MULTI_EBANDS_DOS = wx.NewId()
        #self.ID_MULTI_EBANDS_JDOS = wx.NewId()
        self.ID_MULTI_EBANDS_BANDSWITHDOS = wx.NewId()
                                                                                               
        menu.Append(self.ID_MULTI_EBANDS_PLOT, "Compare ebands", "Plot multiple electron bands")
        self.Bind(wx.EVT_MENU, self.OnCompareEbands, id=self.ID_MULTI_EBANDS_PLOT)
        menu.Append(self.ID_MULTI_EBANDS_DOS, "Compare DOSes", "Compare multiple electron DOSes")
        self.Bind(wx.EVT_MENU, self.OnCompareEdos, id=self.ID_MULTI_EBANDS_DOS)
        #menu.Append(self.ID_MULTI_EBANDS_JDOS, "Compare JDOSes", "Compare multiple electron JDOSes")
        #self.Bind(wx.EVT_MENU, self.OnCompareJdos, id=self.ID_MULTI_EBANDS_JDOS)

        menu.Append(self.ID_MULTI_EBANDS_BANDSWITHDOS, "Plot bands and DOS", "Plot electron bands and DOS on the same figure")
        self.Bind(wx.EVT_MENU, self.onPlotEbandsWithDos, id=self.ID_MULTI_EBANDS_BANDSWITHDOS)

        return menu

    @abc.abstractproperty
    def ebands_list(self):
        """Return a list of `ElectronBands`."""

    @abc.abstractproperty
    def ebands_filepaths(self):
        """
        Return a list with the absolute paths of the files 
        from which the `ElectronBands` have been read.
        """

    def OnCompareEbands(self, event):
        """Plot multiple electron bands"""
        plotter = ElectronBandsPlotter()
        for path, ebands in zip(self.ebands_filepaths, self.ebands_list):
            label = os.path.relpath(path)
            plotter.add_ebands(label, ebands)

        try:
            print(plotter.bands_statdiff())
        except:
            pass
        plotter.plot()

    def OnCompareEdos(self, event):
        """Plot multiple electron DOSes"""
        # Open dialog to get DOS parameters.
        dialog = ewx.ElectronDosDialog(self)
        if dialog.ShowModal() == wx.ID_CANCEL: return 
        dos_params = dialog.GetParams()

        plotter = ElectronDosPlotter()
        for path, ebands in zip(self.ebands_filepaths, self.ebands_list):
            try:
                edos = ebands.get_edos(**dos_params)
                label = os.path.relpath(path)
                plotter.add_edos(label, edos)
            except:
                awx.showErrorMessage(self)

        plotter.plot()

    #def OnCompareJdos(self, event):
        #"""Plot multiple electron JDOSes"""
        # Open dialog to get DOS parameters.
        #dialog = ElectronJdosDialog(self, nsppol, mband)
        #jdos_params = dialog.GetParams()

        #plotter = ElectronBandsPlotter()
        #for ebands in self.ebands_list:
        #    jos = ebands.get_jdos(**jdos_params)
        #    plotter.add_edos(label, edos)
        #                                      
        #plotter.plot()

    def onPlotEbandsWithDos(self, event):
        """Plot electron bands with DOS. Requires the specification of two files."""
        # Open dialog to get files and DOS parameters.
        dialog = ewx.EbandsDosDialog(self, self.ebands_filepaths)
        if dialog.ShowModal() == wx.ID_CANCEL: return 

        try:
            dos_params = dialog.getEdosParams()
            ipath, idos = dialog.getBandsDosIndex()

            ebands_path = self.ebands_list[ipath]
            ebands_mesh = self.ebands_list[idos]

            edos = ebands_mesh.get_edos(**dos_params)
            ebands_path.plot_with_edos(edos)
        except:
            awx.showErrorMessage(self)


#class Has_GsResults(object):
#    """
#    Mixin class for GUIs with ground-state results (etotal, forces, stresses...)
#    """
#    __metaclass__ = abc.ABCMeta

#    def CreateToolsMenu(self):
#        """Create the tools menu."""
#        # Tools Menu ID's
#        self.ID_GSRESULTS_EOSFIT = wx.NewId()
#                                                                                            
#        menu = wx.Menu()
#        menu.Append(self.ID_GSRESULTS_EOSFIT, "Fit E(V)", "Equation of State")
#        self.Bind(wx.EVT_MENU, self.onEosFit, id=self.ID_GSRESULTS_EOSFIT)
#                                                                                            
#        return menu
#
#    def onEosFit(self, event):
#        EosFrame(self, volumes, energies, vol_unit="ang^3", ene_unit="eV").Show()

class Has_Tools(object):
    """
    Mixin class that provides a menu with external tools.
    """
    def CreateToolsMenu(self):
        """Create the tools menu."""
        # Tools Menu ID's
        self.ID_TOOLS_UNIT_CONVERTER = wx.NewId()
        self.ID_TOOLS_PERIODIC_TABLE = wx.NewId()

        menu = wx.Menu()
        menu.Append(self.ID_TOOLS_PERIODIC_TABLE, "Periodic table", "Periodic Table")
        self.Bind(wx.EVT_MENU, self.OnTools_PeriodicTable, id=self.ID_TOOLS_PERIODIC_TABLE)

        menu.Append(self.ID_TOOLS_UNIT_CONVERTER, "Unit converter", "Unit Converter")
        self.Bind(wx.EVT_MENU, self.OnTools_UnitConverter, id=self.ID_TOOLS_UNIT_CONVERTER)

        return menu

    def OnTools_PeriodicTable(self, event):
        """Open new frame with the periodic table."""
        from awx.elements_gui import WxPeriodicTable
        WxPeriodicTable(self).Show()

    def OnTools_UnitConverter(self, event):
        """Open new frame with the unit converter."""
        ConverterFrame(self).Show()


class Has_NetcdfFiles(object):
    __metaclass__ = abc.ABCMeta
    """
    Mixin class that provides a menu and callbacks 
    for analyzing and comparing netcdf files.
    """
    @abc.abstractproperty
    def nc_filepaths(self):
        """List of absolute paths of the netcdf file."""

    def CreateNetcdfMenu(self):
        """Creates the ebands menu."""
        # Netcdf Menu ID's
        self.ID_NETCDF_WXNCVIEW = wx.NewId()
        self.ID_NETCDF_NCDUMP = wx.NewId()
        self.ID_NETCDF_NCVIEW = wx.NewId()

        menu = wx.Menu()

        menu.Append(self.ID_NETCDF_WXNCVIEW, "wxncview", "Call wxncview")
        self.Bind(wx.EVT_MENU, self.OnNetcdf_WxNcView, id=self.ID_NETCDF_WXNCVIEW)

        menu.Append(self.ID_NETCDF_NCDUMP, "ncdump", "Show the output of ncdump")
        self.Bind(wx.EVT_MENU, self.OnNetcdf_NcDump, id=self.ID_NETCDF_NCDUMP)

        menu.Append(self.ID_NETCDF_NCVIEW, "ncview", "Call ncview")
        self.Bind(wx.EVT_MENU, self.OnNetcdf_NcView, id=self.ID_NETCDF_NCVIEW)

        return menu

    def OnNetcdf_NcDump(self, event):
        """Call ncdump and show results in a dialog."""
        for path in self.nc_filepaths:
            s = NcDumper().dump(path)
            caption = "ncdump output for %s" % path
            wxdg.ScrolledMessageDialog(self, s, caption=caption, style=wx.MAXIMIZE_BOX).Show()

    def OnNetcdf_NcView(self, event):
        """Call ncview in an subprocess."""
        if which("ncview") is None:
            return awx.showErrorMessage(self, "Cannot find ncview in $PATH")

        for path in self.nc_filepaths:
            def target():
                os.system("ncview %s" % path)

            thread = awx.WorkerThread(self, target=target)
            thread.start()

    def OnNetcdf_WxNcView(self, event):
        """Open wxncview frame."""
        NcViewerFrame(self, filepaths=self.nc_filepaths).Show()


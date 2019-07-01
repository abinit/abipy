# coding: utf-8
"""TRANSPORT.nc file."""

import numpy as np
import pymatgen.core.units as units
import abipy.core.abinit_units as abu

from monty.functools import lazy_property
from monty.string import marquee
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons.ebands import ElectronsReader
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt

__all__ = [
    "TransportFile",
]

class TransportFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands):
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(TransportFile, self).__init__(filepath)
        self.reader = TransportReader(filepath)

        ebands = self.reader.read_ebands()
        self.fermi = ebands.fermie*abu.eV_Ha
        self.volume = ebands.structure.volume*abu.Ang_Bohr**3

        self.tmesh = self.reader.tmesh

    @property
    def ntemp(self):
        return len(self.tmesh)

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    def get_transport_result(self,itemp=None,el_temp=300):
        """
        Get one instance of TransportResult according to a temperature

        Args:
            itemp: the index of the temperature from which to create the TransportResult class
                   if None, read the transport value without tau
            el_temp: temperature to use in the Fermi integrations
        """
        reader = self.reader
        wmesh, dos, idos = reader.read_dos()

        if itemp is None:
            wmesh, vvdos = reader.read_vvdos()
            tau_temp = None
        else:
            wmesh, vvdos = reader.read_vvdos_tau(itemp)
            tau_temp = reader.tmesh[itemp]
            el_temp = tau_temp

        #todo spin
        #self.nsppol
        dos = dos[0]
        vvdos = vvdos[:,:,0]

        tr = TransportResult(wmesh,wmesh,dos,vvdos,self.fermi,el_temp,self.volume,self.nelect,tau_temp=tau_temp,nsppol=reader.nsppol)
        tr._L0,tr._L1,tr._L2 = reader.read_onsager(itemp)
        tr._sigma,tr._seebeck,tr._kappa,_ = reader.read_transport(itemp)
        return tr

    def get_all_transport_results(self):
        """
        Return multiple instances of TransportResults from the data in the TRANSPORT.nc file
        """
        results = [self.get_transport_result(itemp=itemp) for itemp in list(range(self.ntemp))]
        return TransportResultRobot(results)

    def plot_dos(self,ax=None):
        """
        Plot the density of states
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        wmesh, dos, idos = self.reader.read_dos()
        ax.plot(wmesh,dos)
        ax.set_xlabel('Fermi level (eV)')
        ax.set_ylabel('DOS')
        return fig

    def plot_vvdos(self,component='xx',ax=None,colormap='jet',**kwargs):
        """
        Plot velocity * lifetime density of states
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        for itemp in range(self.ntemp):
            temp = self.tmesh[itemp]
            wmesh, vvdos = self.reader.read_vvdos_tau(itemp,component=component)
            ax.plot(wmesh,vvdos,c=cmap(itemp/self.ntemp),label='T = %dK'%temp)
        ax.set_xlabel('Fermi level (eV)')
        ax.set_ylabel('VVDOS')
        ax.set_yscale('log')
        plt.legend()
        return fig

    def plot_mobility(self,ax=None,component='xx',colormap='jet',**kwargs):
        """
        Read the Mobility from the netcdf file and plot it
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        for itemp in range(self.ntemp):
            temp = self.tmesh[itemp]
            wmesh, mu = self.reader.read_mobility(0,itemp,component,0)
            ax.plot(wmesh,mu,c=cmap(itemp/self.ntemp),label='T = %dK'%temp)
        ax.set_xlabel('Fermi level (eV)')
        ax.set_ylabel('mobility $\mu(\epsilon_F)$ [cm$^2$/Vs]')
        ax.set_yscale('log')
        plt.legend()
        return fig

    def get_mobility_mu(self,eh,itemp,component='xx',ef=None,spin=0):
        """
        Get the value of the mobility at a chemical potential ef

          Arguments:
            itemp : Index of the temperature.
            ef : Value of the doping in eV. The default None uses the chemical potential at the temperature
        """
        from scipy import interpolate
        if ef is None: ef = self.reader.read_value('transport_mu_e')[itemp]
        wmesh, mobility = self.reader.read_mobility(eh,itemp,component,spin)
        f = interpolate.interp1d(wmesh,mobility)
        return f(ef)

    def __str__(self):
        """String representation"""
        return self.to_string()

    def to_string(self,verbose=0):
        """String representation"""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="KS Electron Bands"))
        app("")

        # Transport section.
        app(marquee("Transport calculation", mark="="))
        app("Number of temperatures: %d"%self.ntemp)
        app("Mobility:")
        app("Temperature [K]     Electrons [cm^2/Vs]     Holes [cm^2/Vs]")
        for itemp in range(self.ntemp):
            temp = self.tmesh[itemp]
            mobility_mu_e = self.get_mobility_mu(0,itemp)
            mobility_mu_h = self.get_mobility_mu(1,itemp)
            app("%14.1lf %18.6lf %18.6lf"%(temp,mobility_mu_e,mobility_mu_h))
        return "\n".join(lines)

    def yield_figs(self):
        """
        Return figures plotting the transport data
        """
        yield self.plot_dos()
        yield self.plot_vvdos()
        yield self.plot_mobility()

    def close(self):
        """Close the file."""
        self.reader.close()

class TransportReader(ElectronsReader):
    """
    This class reads the results stored in the TRANSPORT.nc file
    It provides helper function to access the most important quantities.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        super(TransportReader, self).__init__(filepath)
        ktmesh = self.read_value("kTmesh")
        self.tmesh = ktmesh / abu.kb_HaK
        self.nsppol = self.read_dimvalue('nsppol')

    def read_vvdos(self,component='xx',spin=1):
        """
        Read the group velocity density of states
        The vvdos_vals array has 3 dimensions (9,nsppolplus1,nw)
          1. 3x3 components of the tensor
          2. the spin polarization + 1 for the sum
          3. the number of frequencies
        """
        i,j = abu.s2itup(component)
        wmesh = self.read_variable("vvdos_mesh")[:]*abu.Ha_eV
        vals = self.read_variable("vvdos_vals")
        vvdos = vals[i,j,spin,:]
        return wmesh, vvdos

    def read_vvdos_tau(self,itemp,component='xx',spin=1):
        """
        Read the group velocity density of states times lifetime for different temperatures
        The vvdos_tau array has 4 dimensions (ntemp,9,nsppolplus1,nw)
          1. the number of temperatures
          2. 3x3 components of the tensor
          3. the spin polarization + 1 for the sum
          4. the number of frequencies
        """
        i,j = abu.s2itup(component)
        wmesh = self.read_variable("vvdos_mesh")[:]*abu.Ha_eV
        vals = self.read_variable("vvdos_tau")
        vvdos_tau = vals[itemp,i,j,spin,:]/(2*abu.Ha_s)
        return wmesh, vvdos_tau

    def read_dos(self,spin=0):
        """
        Read the density of states
        """
        vals = self.read_variable("edos_mesh")
        wmesh = vals[:]
        vals = self.read_variable("edos_dos")
        dos = vals[spin,:]
        vals = self.read_variable("edos_idos")
        idos = vals[spin,:]
        return wmesh, dos, idos

    def read_onsager(self,itemp):
        """
        Read the onsager coeficients computed in the transport driver in Abinit
        """
        L0 = np.moveaxis(self.read_variable("L0")[itemp,:],[0,1,2,3],[3,2,0,1])
        L1 = np.moveaxis(self.read_variable("L1")[itemp,:],[0,1,2,3],[3,2,0,1])
        L2 = np.moveaxis(self.read_variable("L2")[itemp,:],[0,1,2,3],[3,2,0,1])
        return L0,L1,L2

    def read_transport(self,itemp):
        sigma   = np.moveaxis(self.read_variable("sigma")[itemp,:],   [0,1,2,3],[3,2,0,1])
        kappa   = np.moveaxis(self.read_variable("kappa")[itemp,:],   [0,1,2,3],[3,2,0,1])
        seebeck = np.moveaxis(self.read_variable("seebeck")[itemp,:], [0,1,2,3],[3,2,0,1])
        pi      = np.moveaxis(self.read_variable("pi")[itemp,:],      [0,1,2,3],[3,2,0,1])
        return sigma, kappa, seebeck, pi

    def read_mobility(self,eh,itemp,component,spin):
        """
        Read mobility from the TRANSPORT.nc file
        The mobility is computed separately for electrons and holes.
        """
        i,j = abu.s2itup(component)
        wvals = self.read_variable("vvdos_mesh")
        mobility = self.read_variable("mobility")[eh,itemp,i,j,spin,:]
        return wvals,mobility

    def read_evk_diagonal(self):
        """
        Read the group velocities i.e the diagonal matrix elements.
        Return (nsppol, nkpt) |numpy-array| of real numbers.
        """
        vels = self.read_variable("vred_diagonal")
        # Cartesian? Ha --> eV?
        return vels * (units.Ha_to_eV / units.bohr_to_ang)

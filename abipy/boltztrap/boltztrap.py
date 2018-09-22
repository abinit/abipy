# coding: utf-8
"""
This module containes a Bolztrap2 class to interpolate and analyse the results
It also provides interfaces with Abipy objects allowing to
initialize the Boltztrap2 calculation from Abinit files
"""
import pickle
import numpy as np
from abipy.tools.plotting import add_fig_kwargs
from abipy.tools import duck
import abipy.core.abinit_units as abu
import time

def timeit(method):
    """
    timeit decorator adapted from:
    https://medium.com/pythonhive/python-decorator-to-measure-the-execution-time-of-methods-fa04cb6bb36d
    sets the timing of the routine as an attribute of the class
    """
    def timed(self, *args, **kw):
        ts = time.time()
        result = method(self, *args, **kw)
        te = time.time()

        setattr(self,"time_"+method.__name__, (te - ts) * 1000)
        return result
    return timed

def xyz_comp(component):
    """
    Take a string in the form xx, xy or xz and 
    return indexes (0,0) (0,1) and (0,2) respectively  
    """
    index = {'x':0,'y':1,'z':2} 
    i = index[component[0]]
    j = index[component[0]]
    return i,j

class AbipyBoltztrap():
    """
    Wrapper to Boltztrap2 interpolator
    This class contains the same quantities as the Loader classes from dft.py in Boltztrap2
    Addifionally it has methods to call the Boltztrap2 interpolator.
    It creates an instance of Bolztrap2Results to save the data
    Enter with quantities in the IBZ and interpolate to a fine BZ mesh
    """
    def __init__(self,fermi,atoms,nelect,kpoints,eig,volume,linewidths=None,tmesh=None,mumesh=None,
                 lpratio=1,nworkers=1):
        self.fermi = fermi
        self.atoms = atoms
        self.nelect = nelect
        self.kpoints = kpoints
        self.eig = eig
        self.volume = volume
        self.linewidths = linewidths
        self.tmesh = tmesh
        self.mumesh = mumesh
        self.mommat = None
        self.nworkers = nworkers
        self.lpratio = lpratio

    @property
    def equivalences(self):
        if not hasattr(self,'_equivalences'):
            self.compute_equivalences()
        return self._equivalences

    @property
    def coefficients(self):
        if not hasattr(self,'_coefficients'):
            self.compute_coefficients()
        return self._coefficients

    @property
    def rmesh(self):
        if not hasattr(self,'_rmesh'):
            self.get_interpolation_mesh()
        return self._rmesh

    @property
    def nequivalences(self):
        return len(self.equivalences)   

    @property
    def ncoefficients(self):
        return len(self.coefficients)
    
    @property
    def ntemps(self):
        return len(self.linewidths)

    @classmethod
    def from_sigeph(cls,sigeph,itemp_list=None,bstart=None,bstop=None,lpratio=1):
        """Initialize interpolation of the bands and lifetimes from a sigeph object"""

        #units conversion
        eV_Ry = 2 * abu.eV_Ha
        eV_s = abu.eV_to_THz*1e12 * 2*np.pi

        #get the lifetimes as an array
        qpes = sigeph.get_qp_array(mode='ks+lifetimes')

        #get other dimensions
        if bstart is None: bstart = sigeph.reader.max_bstart
        if bstop is None:  bstop  = sigeph.reader.min_bstop
        fermi  = sigeph.ebands.fermie*eV_Ry
        atoms  = sigeph.ebands.structure.to_ase_atoms()
        volume = sigeph.ebands.structure.volume
        nelect = sigeph.ebands.nelect
        kpoints = [k.frac_coords for k in sigeph.sigma_kpoints]

        #TODO handle spin
        eig = qpes[0,:,bstart:bstop,0].real.T*eV_Ry

        itemp_list = list(range(sigeph.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)
        linewidths = []
        tmesh = []
        mumesh = []
        for itemp in itemp_list:
            tmesh.append(sigeph.tmesh[itemp])
            mumesh.append(sigeph.mu_e[itemp])
            #TODO handle spin
            linewidth = qpes[0, :, bstart:bstop, itemp].imag.T*eV_Ry
            linewidths.append(linewidth)

        return cls(fermi, atoms, nelect, kpoints, eig, volume, linewidths, 
                   tmesh, mumesh, lpratio=lpratio)

    def get_lattvec(self):
        """this method is required by Bolztrap"""
        return self.lattvec

    @property
    def lattvec(self):
        if not hasattr(self,"_lattvec"):
            self._lattvec = self.atoms.get_cell().T / abu.Bohr_Ang
        return self._lattvec

    def get_interpolation_mesh(self):
        """From the array of equivalences determine the mesh that was used"""
        max1, max2, max3 = 0,0,0
        for equiv in self.equivalences:
            max1 = max(np.max(equiv[:,0]),max1)
            max2 = max(np.max(equiv[:,1]),max2)
            max3 = max(np.max(equiv[:,2]),max3)
        self._rmesh = (2*max1+1,2*max2+1,2*max3+1)
        return self._rmesh

    @timeit
    def compute_equivalences(self):
        """Compute equivalent k-points"""
        from BoltzTraP2 import sphere
        self._equivalences = sphere.get_equivalences(self.atoms, self.lpratio)

    @timeit
    def compute_coefficients(self):
        """Call fitde3D routine from Boltztrap2"""
        from BoltzTraP2 import fite
        #we will set ebands to compute teh coefficients
        self.ebands = self.eig
        self._coefficients = fite.fitde3D(self, self.equivalences, nworkers=self.nworkers)

        if self.linewidths:
            self._linewidth_coefficients = []
            for itemp in range(self.ntemps):
                self.ebands = self.linewidths[itemp]
                coeffs = fite.fitde3D(self, self.equivalences, nworkers=self.nworkers)
                self._linewidth_coefficients.append(coeffs)

        #at the end we always unset ebands
        delattr(self,"ebands")

    @timeit
    def run(self,npts=500,dos_method='gaussian:0.1 eV',erange=None,verbose=True):
        """
        Interpolate the eingenvalues This part is quite memory intensive
        """
        boltztrap_results = []; app = boltztrap_results.append
        eV_s = abu.eV_to_THz*1e12 * 2*np.pi
        from BoltzTraP2 import fite
        import BoltzTraP2.bandlib as BL

        #TODO change this!
        if erange is None: erange = (np.min(self.eig),np.max(self.eig))

        #interpolate the electronic structure
        results = fite.getBTPbands(self.equivalences, self.coefficients, 
                                   self.lattvec, nworkers=self.nworkers)
        eig_fine, vvband, cband = results

        #calculate DOS and VDOS without lifetimes
        wmesh,dos,vvdos,_ = BL.BTPDOS(eig_fine, vvband, erange=erange, npts=npts, mode=dos_method) 
        app(BoltztrapResult(self,wmesh,dos,vvdos,self.fermi,self.mumesh,self.tmesh,self.volume))
    
        #if we have linewidths
        if self.linewidths:
            for itemp in range(self.ntemps):
                #calculate the lifetimes on the fine grid
                results = fite.getBTPbands(self.equivalences, self._linewidth_coefficients[itemp], 
                                           self.lattvec, nworkers=self.nworkers)
                linewidth_fine, vvband, cband = results
                tau_fine = 1.0/np.abs(2*linewidth_fine*eV_s) 
                
                #calculate vvdos with the lifetimes
                wmesh, dos_tau, vvdos_tau, _ = BL.BTPDOS(eig_fine, vvband, erange=erange, npts=npts,
                                                         scattering_model=tau_fine, mode=dos_method)
                #store results
                app(BoltztrapResult(self,wmesh,dos_tau,vvdos_tau,self.fermi,self.mumesh,self.tmesh,
                                    self.volume,tau_temp=self.tmesh[itemp]))
 
        return BoltztrapResultRobot(boltztrap_results)

    def __str__(self):
        lines = []; app = lines.append
        app("nequiv: {}".format(self.nequivalences))
        app("rmesh:  {}".format(self.rmesh))
        return "\n".join(lines)

class BoltztrapResult():
    """
    Container for BoltztraP2 results
    Provides a object oriented interface to BoltztraP2 for plotting, 
    storing and analysing the results
    """
    def __init__(self,abipyboltztrap,wmesh,dos,vvdos,fermi,mumesh,tmesh,volume,tau_temp=None):
        self.abipyboltztrap = abipyboltztrap

        self.fermi  = fermi
        self.mumesh = mumesh
        self.volume = volume
        self.wmesh  = np.array(wmesh)
        self.tmesh  = np.array(tmesh)

        self.tau_temp = tau_temp

        self.dos = dos
        self.vvdos = vvdos

    @property
    def ntemp(self):
        return len(self.tmesh)

    @property
    def L0(self):
        if not hasattr(self,'_L0'):
            self.compute_fermiintegrals()
        return self._L0
 
    @property
    def L1(self):
        if not hasattr(self,'_L1'):
            self.compute_fermiintegrals()
        return self._L1

    @property
    def L2(self):
        if not hasattr(self,'_L2'):
            self.compute_fermiintegrals()
        return self._L2
 
    @property
    def sigma(self):
        if not hasattr(self,'_sigma'):
            self.compute_onsager_coefficients()
        return self._sigma
 
    @property
    def seebeck(self):
        if not hasattr(self,'_seebeck'):
            self.compute_onsager_coefficients()
        return self._seebeck

    @property
    def powerfactor(self):
        return self.sigma * self.seebeck**2
     
    @property
    def kappa(self):
        if not hasattr(self,'_kappa'):
            self.compute_onsager_coefficients()
        return self._kappa
   
    def compute_fermiintegrals(self):
        """Compute and store the results of the Fermi integrals""" 
        import BoltzTraP2.bandlib as BL
        results = BL.fermiintegrals(self.wmesh, self.dos, self.vvdos, mur=self.wmesh, Tr=self.tmesh)
        _, self._L0, self._L1, self._L2, self._Lm11 = results

    def compute_onsager_coefficients(self):
        """Compute Onsager coefficients"""
        import BoltzTraP2.bandlib as BL
        L0,L1,L2 = self.L0,self.L1,self.L2
        results = BL.calc_Onsager_coefficients(L0,L1,L2,self.wmesh,self.tmesh,self.volume)
        self._sigma, self._seebeck, self._kappa, self._hall = results

    @staticmethod
    def from_pickle(filename):
        """load results from file"""
        with open(filename,'rb') as f:
            instance = pickle.load(f)
        return instance
 
    def pickle(self,filename):
        """Write a file with the results from the calculation"""
        with open(filename,'wb') as f:
            pickle.dump(self,f)

    def istensor(self,what):
        """check if a certain quantity is a tensor"""
        if not hasattr(self,what): return None
        return len(getattr(self,what).shape) > 2

    def get_component(self,what,component,itemp):
        i,j = xyz_comp(component)
        return getattr(self,what)[itemp,:,i,j]
   
    def plot_vvdos_ax(self,ax1,components):
        """ Plot components of vvdos"""
        wmesh = (self.wmesh-self.fermi) * abu.Ha_eV
        for component in components:
            i,j = xyz_comp(component)
            label = "vvdos %s"%component
            if self.tau_temp: label += " $\\tau_T$ = %dK"%self.tau_temp
            ax1.plot(wmesh,self.vvdos[i,j,:],label=label)
        ax1.set_xlabel('Energy (eV)')

    def plot_ax(self,ax1,what,components=['xx'],itemp_list=None,colormap='viridis',**kwargs):
        """Plot the DOS for all the dopings as a function of temperature for an axes"""
        from matplotlib import pyplot as plt

        if what == 'vvdos':
            self.plot_vvdos_ax(ax1,components)
            return

        itemp_list = list(range(self.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)

        cmap = plt.get_cmap(colormap)
        color = None
        wmesh = (self.wmesh-self.fermi) * abu.Ha_eV

        if self.istensor(what):
            for itemp in itemp_list:
                for component in components:
                    y = self.get_component(what,component,itemp)
                    if len(itemp_list) > 1: color=cmap(itemp/len(itemp_list))
                    label = "%s %s"%(what,component)
                    if self.tau_temp: label += " $\\tau_T$ = %dK"%self.tau_temp
                    ax1.plot(wmesh,y,label=label,c=color,**kwargs)
        else:
            ax1.plot(wmesh,getattr(self,what),label=what)
        ax1.set_xlabel('Energy (eV)')

    @add_fig_kwargs
    def plot(self,what,colormap='viridis',directions=['xx']):
        """Plot the DOS for all the dopings as a function of temperature"""
        from matplotlib import pyplot as plt

        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        if what == 'dos': self.plot_dos_ax(ax1)
        if what == 'powerfactor': self.plot_powerfactor_ax(ax1)
        fig.legend()
        return fig

    def __str__(self):
        lines = []; app = lines.append
        return "".join(lines)


class BoltztrapResultRobot():
    """
    Robot to analyse multiple Boltztrap calculations
    """
    def __init__(self,results):
        if not all([isinstance(r,BoltztrapResult) for r in results]):
            raise ValueError('Must provide BolztrapResult instances')
        self.results = results

    @staticmethod
    def from_pickle(filename):
        """load results from file"""
        with open(filename,'rb') as f:
            instance = pickle.load(f)
        return instance
 
    def pickle(self,filename):
        """Write a file with the results from the calculation"""
        with open(filename,'wb') as f:
            pickle.dump(self,f)
 
    def plot_ax(self,ax,what,components=['xx'],itemp_list=None,with_tau=None):
        """
        Plot a quantity in an axis for all the results
        """
        for result in self.results:
            if with_tau is True and result.tau_temp == None: continue
            result.plot_ax(ax,what,components,itemp_list)
 
    @add_fig_kwargs
    def plot(self,what,itemp_list=None,components=['xx'],with_tau=None):
        """
        Plot all the boltztrap results
        """    
        from matplotlib import pyplot as plt

        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        self.plot_ax(ax1,what,components=components,itemp_list=itemp_list,with_tau=with_tau)
        fig.legend()
        return fig 

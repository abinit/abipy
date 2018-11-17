# coding: utf-8
"""
This module containes a Bolztrap2 class to interpolate and analyse the results
It also provides interfaces with Abipy objects allowing to
initialize the Boltztrap2 calculation from Abinit files

Warning:

    Work in progress
"""
import pickle
import numpy as np
import abipy.core.abinit_units as abu
from monty.string import marquee
from monty.termcolor import cprint
from monty.dev import deprecated
from abipy.tools.plotting import add_fig_kwargs
from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands
from abipy.core.kpoints import Kpath
from abipy.core.structure import Structure
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt #, set_axlims, set_visible, set_ax_xylabels
from abipy.tools.decorators import timeit

class AbipyBoltztrap():
    """
    Wrapper to Boltztrap2 interpolator
    This class contains the same quantities as the Loader classes from dft.py in Boltztrap2
    Additionally it has methods to call the Boltztrap2 interpolator.
    It creates multiple instances of BolztrapResult storing the results of the interpolation
    Enter with quantities in the IBZ and interpolate to a fine BZ mesh
    """
    def __init__(self,fermi,structure,nelect,kpoints,eig,volume,linewidths=None,tmesh=None,
                 mommat=None,magmom=None,lpratio=5):
        #data needed by boltztrap
        self.fermi = fermi
        self.atoms = structure.to_ase_atoms()
        self.nelect = nelect
        self.kpoints = np.array(kpoints)
        self.volume = volume
        self.mommat = mommat
        self.magmom = magmom

        #additional parameters
        self.eig = eig
        self.structure = structure
        self.linewidths = linewidths
        self.tmesh = tmesh
        self.lpratio = lpratio

    @property
    def nkpoints(self):
        return len(self.kpoints)

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
    def linewidth_coefficients(self):
        if not hasattr(self,'_linewidth_coefficients'):
            self.compute_coefficients()
        return self._linewidth_coefficients

    @property
    def linewidth_coefficients(self):
        if not hasattr(self,'_linewidth_coefficients'):
            self.compute_coefficients()
        return self._linewidth_coefficients

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

    def pickle(self,filename):
        with open(filename,'wb') as f:
            pickle.dump(self,f)

    @classmethod
    def from_pickle(cls,filename):
        with open(filename,'rb') as f:
            cls = pickle.load(f)
        return cls

    @classmethod
    def from_ebands(cls):
        """Initialize from an ebands object"""
        raise NotImplementedError('TODO')

    @classmethod
    def from_evk(cls):
        """Intialize from a EVK file"""
        raise NotImplementedError('TODO')

    @classmethod
    def from_dftdata(cls,dftdata,tmesh,lpratio=5):
        """
        Initialize an instance of this class from a DFTData instance from Boltztrap

        Args:
            dftdata: DFTData
            tmesh: a list of temperatures to use in the fermi integrations
            lpratio: ratio to multiply by the number of k-points in the IBZ and give the
                     number of real space points inside a sphere
        """
        structure = Structure.from_ase_atoms(dftdata.atoms)
        return cls(dftdata.fermi,structure,dftdata.nelect,dftdata.kpoints,dftdata.ebands,
                   dftdata.get_volume(),linewidths=None,tmesh=tmesh,
                   mommat=dftdata.mommat,magmom=None,lpratio=lpratio)

    @classmethod
    def from_sigeph(cls, sigeph, itemp_list=None, bstart=None, bstop=None, lpratio=5):
        """
        Initialize interpolation of the bands and lifetimes from a SigEphFile object

        Args:
            sigeph: |SigEphFile| instance
            itemp_list: list of the temperature indexes to consider
            bstart, bstop: only consider bands between bstart and bstop
            lpratio: ratio to multiply by the number of k-points in the IBZ and give the
                     number of real space points inside a sphere
        """
        #get the lifetimes as an array
        qpes = sigeph.get_qp_array(mode='ks+lifetimes')

        #get other dimensions
        if bstart is None: bstart = sigeph.reader.max_bstart
        if bstop is None:  bstop  = sigeph.reader.min_bstop
        fermi  = sigeph.ebands.fermie*abu.eV_Ha
        structure = sigeph.ebands.structure
        volume = sigeph.ebands.structure.volume*abu.Ang_Bohr**3
        nelect = sigeph.ebands.nelect
        kpoints = [k.frac_coords for k in sigeph.sigma_kpoints]

        if sigeph.nsppol == 2:
            raise NotImplementedError("nsppol 2 not implemented")

        #TODO handle spin
        eig = qpes[0,:,bstart:bstop,0].real.T*abu.eV_Ha

        itemp_list = list(range(sigeph.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)
        linewidths = []
        tmesh = []
        for itemp in itemp_list:
            tmesh.append(sigeph.tmesh[itemp])
            fermi = sigeph.mu_e[itemp]*abu.eV_Ha
            #TODO handle spin
            linewidth = qpes[0, :, bstart:bstop, itemp].imag.T*abu.eV_Ha
            linewidths.append(linewidth)

        return cls(fermi, structure, nelect, kpoints, eig, volume, linewidths=linewidths,
                   tmesh=tmesh, lpratio=lpratio)

    def get_lattvec(self):
        """this method is required by Bolztrap"""
        return self.lattvec

    @property
    def nbands(self):
        nbands, rpoints = self.coefficients.shape
        return nbands

    @property
    def lattvec(self):
        if not hasattr(self,"_lattvec"):
            self._lattvec = self.atoms.get_cell().T / abu.Bohr_Ang
        return self._lattvec

    def get_ebands(self,kpath=None,line_density=20,vertices_names=None,linewidth_itemp=False):
        """
        Compute the band-structure using the computed coefficients

        Args:
            kpath: |Kpath| instance where to interpolate the eigenvalues and linewidths
            line_density: Number of points used to sample the smallest segment of the path
            vertices_names:  List of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path.
            linewith_itemp: list of indexes refering to the temperatures where the linewidth will be interpolated
        """
        from BoltzTraP2 import fite

        if kpath is None:
            if vertices_names is None:
                vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]

            kpath = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)

        #call boltztrap to interpolate
        coeffs = self.coefficients
        eigens_kpath, vvband = fite.getBands(kpath.frac_coords, self.equivalences, self.lattvec, coeffs)

        linewidths_kpath = None
        if linewidth_itemp is not False:
            coeffs = self.linewidth_coefficients[linewidth_itemp]
            linewidths_kpath, vvband = fite.getBands(kpath.frac_coords, self.equivalences, self.lattvec, coeffs)
            linewidths_kpath = linewidths_kpath.T[np.newaxis,:,:]*abu.Ha_eV

        #convert units and shape
        eigens_kpath   = eigens_kpath.T[np.newaxis,:,:]*abu.Ha_eV
        occfacts_kpath = np.zeros_like(eigens_kpath)
        nspinor1 = 1
        nspden1 = 1

        #return a ebands object
        return ElectronBands(self.structure, kpath, eigens_kpath, self.fermi*abu.Ha_eV, occfacts_kpath,
                             self.nelect, nspinor1, nspden1, linewidths=linewidths_kpath)

    @deprecated(message="get_bands is deprecated, use get_ebands")
    def get_bands(self, **kwargs):
        return self.get_ebands(**kwargs)

    def get_interpolation_mesh(self):
        """From the array of equivalences determine the mesh that was used"""
        max1, max2, max3 = 0,0,0
        for equiv in self.equivalences:
            max1 = max(np.max(equiv[:,0]),max1)
            max2 = max(np.max(equiv[:,1]),max2)
            max3 = max(np.max(equiv[:,2]),max3)
        self._rmesh = (2*max1+1,2*max2+1,2*max3+1)
        return self._rmesh

    def dump_rsphere(self,filename):
        """ Write a file with the real space points"""
        with open(filename, 'wt') as f:
            for iband in range(self.nbands):
                for ie,equivalence in enumerate(self.equivalences):
                    coeff = self.coefficients[iband,ie]
                    for ip,point in enumerate(equivalence):
                        f.write("%5d %5d %5d "%tuple(point)+"%lf\n"%((abs(coeff))**(1./3)))
                f.write("\n\n")

    @timeit
    def compute_equivalences(self):
        """Compute equivalent k-points"""
        from BoltzTraP2 import sphere
        try:
            self._equivalences = sphere.get_equivalences(self.atoms, self.magmom, self.lpratio*self.nkpoints)
        except TypeError:
            self._equivalences = sphere.get_equivalences(self.atoms, self.lpratio*self.nkpoints)

    @timeit
    def compute_coefficients(self):
        """Call fitde3D routine from Boltztrap2"""
        from BoltzTraP2 import fite
        #we will set ebands to compute teh coefficients
        self.ebands = self.eig
        self._coefficients = fite.fitde3D(self, self.equivalences)

        if self.linewidths:
            self._linewidth_coefficients = []
            for itemp in range(self.ntemps):
                self.ebands = self.linewidths[itemp]
                coeffs = fite.fitde3D(self, self.equivalences)
                self._linewidth_coefficients.append(coeffs)

        #at the end we always unset ebands
        delattr(self,"ebands")

    @timeit
    def run(self,npts=500,dos_method='gaussian:0.05 eV',erange=None,margin=0.1,nworkers=1,verbose=0):
        """
        Interpolate the eingenvalues to compute dos and vvdos
        This part is quite memory intensive

        Args:
            npts: number of frequency points
            dos_method: when using a patched version of Boltztrap
        """
        boltztrap_results = []; app = boltztrap_results.append
        import inspect
        from BoltzTraP2 import fite
        import BoltzTraP2.bandlib as BL

        def BTPDOS(eband,vvband,cband=None,erange=None,npts=None,scattering_model="uniform_tau",mode=dos_method):
            """
            This is a small wrapper for Boltztrap2 to use the official version or a modified
            verison using gaussian or lorentzian smearing
            """
            try:
                return BL.BTPDOS(eband, vvband, erange=erange, npts=npts, scattering_model=scattering_model, mode=dos_method)
            except TypeError:
                return BL.BTPDOS(eband, vvband, erange=erange, npts=npts, scattering_model=scattering_model)


        #TODO change this!
        if erange is None: erange = (np.min(self.eig),np.max(self.eig))
        else: erange = np.array(erange)/abu.Ha_eV+self.fermi

        #interpolate the electronic structure
        if verbose: print('interpolating bands')
        results = fite.getBTPbands(self.equivalences, self.coefficients,
                                   self.lattvec, nworkers=nworkers)
        eig_fine, vvband, cband = results

        #calculate DOS and VDOS without lifetimes
        if verbose: print('calculating dos and vvdos without lifetimes')
        wmesh,dos,vvdos,_ = BTPDOS(eig_fine, vvband, erange=erange, npts=npts, mode=dos_method)
        app(BoltztrapResult(self,wmesh,dos,vvdos,self.fermi,self.tmesh,self.volume,margin=margin))

        #if we have linewidths
        if self.linewidths:
            for itemp in range(self.ntemps):
                if verbose: print('itemp %d\ninterpolating bands')
                #calculate the lifetimes on the fine grid
                results = fite.getBTPbands(self.equivalences, self._linewidth_coefficients[itemp],
                                           self.lattvec, nworkers=nworkers)
                linewidth_fine, vvband, cband = results
                tau_fine = 1.0/np.abs(2*linewidth_fine*abu.eV_s)

                #calculate vvdos with the lifetimes
                if verbose: print('calculating dos and vvdos with lifetimes')
                wmesh, dos_tau, vvdos_tau, _ = BTPDOS(eig_fine, vvband, erange=erange, npts=npts,
                                                      scattering_model=tau_fine, mode=dos_method)
                #store results
                app(BoltztrapResult(self,wmesh,dos_tau,vvdos_tau,self.fermi,self.tmesh,
                                    self.volume,tau_temp=self.tmesh[itemp],margin=margin))

        return BoltztrapResultRobot(boltztrap_results)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=2):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__,mark="="))
        app("equivalent points: {}".format(self.nequivalences))
        app("real space mesh:   {}".format(self.rmesh))
        app("lpratio:           {}".format(self.lpratio))
        return "\n".join(lines)


class BoltztrapResult():
    """
    Container for BoltztraP2 results
    Provides a object oriented interface to BoltztraP2
    for plotting, storing and analysing the results
    """
    _attrs = ['_L0','_L1','_L2','_sigma','_seebeck','_kappa']

    def __init__(self,abipyboltztrap,wmesh,dos,vvdos,fermi,tmesh,volume,tau_temp=None,margin=0.1):
        self.abipyboltztrap = abipyboltztrap

        self.fermi  = fermi
        self.volume = volume
        self.wmesh  = np.array(wmesh)
        idx_margin = int(margin*len(wmesh))
        self.mumesh = self.wmesh[idx_margin:-(idx_margin+1)]
        self.tmesh  = np.array(tmesh)

        #Temperature fix
        if any(self.tmesh < 1):
            cprint("Boltztrap does not handle 0K well.\n"
                   "I avoid potential problems by setting all T<1K to T=1K",color="yellow")
            self.tmesh[self.tmesh < 1] = 1

        self.tau_temp = tau_temp

        self.dos = dos
        self.vvdos = vvdos

    @property
    def has_tau(self):
        return self.tau_temp is not None

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

    def set_tmesh(self,tmesh):
        """ Set the temperature mesh"""
        self.tmesh = tmesh

    def set_tmesh(self,tmesh):
        """ Set the temperature mesh"""
        self.tmesh = tmesh

    def del_attrs(self):
        """ Remove all the atributes so they are recomputed """
        for attr in self._attrs:
            delattr(attr)

    def set_mumesh(self,emin,emax):
        """
        Set the range in which to plot the change of the doping

        Args:
            emin: minimun energy in eV
            emax: maximum energy in eV
        """
        start_idx = np.abs(self.wmesh - emin*abu.eV_Ha - self.fermi).argmin()
        stop_idx  = np.abs(self.wmesh - emax*abu.eV_Ha - self.fermi).argmin()
        self.mumesh = self.wmesh[start_idx:stop_idx]

    def compute_fermiintegrals(self):
        """Compute and store the results of the Fermi integrals"""
        import BoltzTraP2.bandlib as BL
        results = BL.fermiintegrals(self.wmesh, self.dos, self.vvdos, mur=self.mumesh, Tr=self.tmesh)
        _, self._L0, self._L1, self._L2, self._Lm11 = results

    def compute_onsager_coefficients(self):
        """Compute Onsager coefficients"""
        import BoltzTraP2.bandlib as BL
        L0,L1,L2 = self.L0,self.L1,self.L2
        results = BL.calc_Onsager_coefficients(L0,L1,L2,mur=self.mumesh,Tr=self.tmesh,vuc=self.volume)
        self._sigma, self._seebeck, self._kappa, self._hall = results

    @staticmethod
    def from_pickle(filename):
        """Load BoltztrapResult from a pickle file"""
        with open(filename,'rb') as f:
            instance = pickle.load(f)
        return instance

    def pickle(self,filename):
        """Write a file with the results from the calculation"""
        with open(filename,'wb') as f:
            pickle.dump(self,f)

    def istensor(self,what):
        """Check if a certain quantity is a tensor"""
        if not hasattr(self,what): return None
        return len(getattr(self,what).shape) > 2

    def get_component(self,what,component,itemp):
        i,j = abu.s2itup(component)
        return getattr(self,what)[itemp,:,i,j]

    def plot_dos_ax(self,ax,fontsize=8,**kwargs):
        """
        Plot the density of states on axis ax.

        Args:
            ax: |matplotlib-Axes|.
            kwargs: Passed to ax.plot
        """
        wmesh = (self.wmesh-self.fermi) * abu.Ha_eV
        ax.plot(wmesh,self.dos,label=self.get_letter('dos'),**kwargs)
        ax.set_xlabel('Energy (eV)',fontsize=fontsize)

    def plot_vvdos_ax(self,ax,components=('xx',),fontsize=8,**kwargs):
        """
        Plot components of vvdos on the axis ax.

        Args:
            ax: |matplotlib-Axes|.
            components: Choose the components of the tensor to plot ['xx','xy','xz','yy',(...)]
            kwargs: Passed to ax.plot
        """
        wmesh = (self.wmesh-self.fermi) * abu.Ha_eV

        for component in components:
            i,j = abu.s2itup(component)
            label = "%s $_{%s}$" % (self.get_letter('vvdos'),component)
            if self.tau_temp: label += r" $\tau_T$ = %dK" % self.tau_temp
            ax.plot(wmesh,self.vvdos[i,j,:],label=label,**kwargs)
        ax.set_xlabel('Energy (eV)',fontsize=fontsize)

    def plot_ax(self, ax, what, components=('xx',), itemp_list=None, fontsize=8, **kwargs):
        """
        Plot a quantity for all the dopings as a function of temperature on the axis ax.

        Args:
            ax: |matplotlib-Axes|.
            what: choose the quantity to plot can be: ['sigma','kappa','powerfactor']
            components: Choose the components of the tensor to plot ['xx','xy','xz','yy',(...)]
            itemp_list: list of indexes of the tempratures to plot
            colormap: Colormap used to plot the results
            kwargs: Passed to ax.plot
        """
        from matplotlib import pyplot as plt
        colormap = kwargs.pop('colormap','plasma')
        cmap = plt.get_cmap(colormap)
        color = None

        itemp_list = list(range(self.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)
        maxitemp = max(itemp_list)
        minitemp = min(itemp_list)
        if maxitemp > self.ntemp or minitemp < 0:
            raise ValueError('Invalid itemp_list, should be between 0 and %d. Got %d.'%(self.ntemp,maxitemp))

        mumesh = (self.mumesh-self.fermi) * abu.Ha_eV

        if self.istensor(what):
            kwargs.pop('c',None)
            for itemp in itemp_list:
                for component in components:
                    y = self.get_component(what,component,itemp)
                    if len(itemp_list) > 1: color=cmap(itemp/len(itemp_list))
                    label = "%s $_{%s}$ $b_T$ = %dK" % (self.get_letter(what),component,self.tmesh[itemp])
                    if self.has_tau: label += r" $\tau_T$ = %dK" % self.tau_temp
                    ax.plot(mumesh,y,label=label,c=color,**kwargs)
        else:
            ax.plot(mumesh,getattr(self,what), label=what, **kwargs)

        ax.set_ylabel(self.get_ylabel(what), fontsize=fontsize)
        ax.set_xlabel('Energy (eV)', fontsize=fontsize)

    def get_ylabel(self,what):
        """
        Get a label with units for the quntities stores in this object.
        """
        if self.has_tau: tau = ''
        else: tau = 's^{-1}'
        if what == 'sigma':       return r'$\sigma$ [$Sm^{-1}%s$]'%tau
        if what == 'seebeck':     return r'$S$ [$VSm^{-1}%s$]'%tau
        if what == 'kappa':       return r'$\kappa_e$ [$VJSm^{-1}%s$]'%tau
        if what == 'powerfactor': return r'$S^2\sigma$ [$VJSm^{-1}%s$]'%tau
        return ''

    def get_letter(self,what):
        letters = {'sigma':      r'$\sigma$',
                   'seebeck':    r'$S$',
                   'kappa':      r'$\kappa_e$',
                   'powerfactor':r'$S^2\sigma$',
                   'vvdos':      r'$v\otimes v$',
                   'dos':        r'$n(\epsilon)$'}
        return letters[what]

    @add_fig_kwargs
    def plot(self, what, colormap='plasma', directions=('xx'), ax=None, fontsize=8, **kwargs):
        """
        Plot the qantity for all the temperatures as a function of the doping
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        self.plot_ax(ax, what, colormap=colormap, directions=directions, **kwargs)
        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def to_string(self, title=None, mark="=", verbose=0):
        """
        String representation of the class
        """
        lines = []; app = lines.append
        if title is None: app(marquee(self.__class__.__name__,mark=mark))
        app("fermi:    %8.5lf eV"%(self.fermi*abu.Ha_eV))
        app("mumesh:   %8.5lf <-> %8.5lf eV"%(self.mumesh[0]*abu.Ha_eV,self.mumesh[-1]*abu.Ha_eV))
        app("tmesh:    %s K"%self.tmesh)
        app("has_tau:  %s"%self.has_tau)
        if self.tau_temp: app("tau_temp: %.1lf K"%self.tau_temp)
        return "\n".join(lines)

    def __str__(self):
        return self.to_string()

class BoltztrapResultRobot():
    """
    Robot to analyse multiple Boltztrap calculations
    Behaves as a list of BoltztrapResult
    Provides methods to plot multiple results on a single figure
    """
    def __init__(self,results,erange=None):
        if not all([isinstance(r,BoltztrapResult) for r in results]):
            raise ValueError('Must provide BolztrapResult instances.')

        #consistency check in the temperature meshes
        res0 = results[0]
        if np.any([res0.tmesh != res.tmesh for res in results]):
            cprint("Comparing BoltztrapResults with different temperature meshes.", color="yellow")

        #consistency check in chemical potential meshes
        if np.any([res0.wmesh != res.wmesh for res in results]):
            cprint("Comparing BoltztrapResults with different energy meshes.", color="yellow")

        #store the results
        self.results = results
        self.erange = erange

        if not all([np.allclose(results[0].mumesh,result.mumesh) for result in results[1:]]):
            raise ValueError('The doping meshes of the results differ, cannot continue')
        self.mumesh = results[0].mumesh

        if not all([np.allclose(results[0].tmesh,result.tmesh) for result in results[1:]]):
            raise ValueError('The temperature meshes of the results differ, cannot continue')
        self.tmesh = results[0].tmesh

    def __getitem__(self,index):
        """Access the results stored in the class as a list"""
        return self.results[index]

    @property
    def ntemp(self):
        return len(self.tmesh)

    @property
    def tau_list(self):
        """Get all the results with tau included"""
        return [ res.tau_temp for res in self.results if res.tau_temp is not None ]

    @property
    def notau_results(self):
        """Get all the results without the tau included"""
        instance = self.__class__([ res for res in self.results if res.tau_temp is None ])
        if self.erange: instance.erange = self.erange
        return instance

    @property
    def tau_results(self):
        """Return all the results that have temperature dependence"""
        instance = self.__class__([ res for res in self.results if res.tau_temp ])
        if self.erange: instance.erange = self.erange
        return instance

    @property
    def nresults(self):
        return len(self.results)

    @staticmethod
    def from_pickle(filename):
        """
        Load results from file
        """
        with open(filename,'rb') as f:
            instance = pickle.load(f)
        return instance

    def pickle(self,filename):
        """
        Write a file with the results from the calculation
        """
        with open(filename,'wb') as f:
            pickle.dump(self,f)

    def plot_vvdos_ax(self,ax,legend=True,components=('xx',),itau_list=None,fontsize=8,erange=None,**kwargs):
        """
        Plot the vvdos for all the results in the robot
        """
        from matplotlib import pyplot as plt
        colormap = kwargs.pop('colormap','plasma')
        cmap = plt.get_cmap(colormap)

        #set erange
        erange = erange or self.erange
        if erange is not None: ax.set_xlim(erange)

        if itau_list:
            #filter results by temperature
            tau_list = [self.tmesh[itau] for itau in itau_list]
            filtered_results = sorted([res for res in self.results if res.tau_temp in tau_list],key=lambda x: x.tau_temp)

            for itemp,result in enumerate(filtered_results):
                color = kwargs.pop('c',cmap(itemp/len(filtered_results)))
                result.plot_vvdos_ax(ax,fontsize=fontsize,c=color,components=components,**kwargs)
            ax.set_ylabel(r'with $\tau$',fontsize=fontsize)
            if legend: ax.legend(loc="best", shadow=True, fontsize=fontsize)
        else:
            #results without temperature
            for result in self.notau_results:
                result.plot_vvdos_ax(ax,fontsize=fontsize,components=components,**kwargs)
            ax.set_ylabel(r'without $\tau$',fontsize=fontsize)
            if legend: ax.legend(loc="best", shadow=True, fontsize=fontsize)

    def plot_dos_ax(self, ax1, legend=True, fontsize=8, erange=None, **kwargs):
        """
        Plot the dos for all the results in the robot
        """
        #set erange
        erange = erange or self.erange
        if erange is not None: ax1.set_xlim(erange)

        for result in self.results:
            result.plot_dos_ax(ax1,fontsize=fontsize,**kwargs)
        if legend: ax1.legend(loc="best", shadow=True, fontsize=fontsize)

    def plot_ax(self,ax1,what,components=('xx',),itemp_list=None,itau_list=None,fontsize=8,erange=None,**kwargs):
        """
        Plot the same quantity for all the results on axis ax1

        Args:
            ax1: |matplotlib-Axes|.
            what: choose the quantity to plot can be: ['sigma','kappa','powerfactor']
            itemp_list: list of indexes of the tempratures to plot
            itau_list: list of indexes of the tempratures at which the lifetimes were computed
            components: Choose the components of the tensor to plot ['xx','xy','xz','yy',(...)]
            erange: choose energy range of the plot
            kwargs: Passed to ax.plot
       """
        from matplotlib import pyplot as plt
        colormap = kwargs.pop('colormap','plasma')
        cmap = plt.get_cmap(colormap)

        #set erange
        erange = erange or self.erange
        if erange is not None: ax1.set_xlim(erange)

        if itau_list:
            #filter results by temperature
            tau_list = self.tmesh if itau_list is None else [self.tmesh[itau] for itau in itau_list]
            filtered_results = [res for res in self.results if res.tau_temp in tau_list]

            #plot the results
            for itemp,result in enumerate(filtered_results):
                color = kwargs.pop('c',cmap(itemp/len(filtered_results)))
                result.plot_ax(ax1,what,components,itemp_list,fontsize=fontsize,c=color,**kwargs)
        else:
            #plot result without tau
            for result in self.notau_results:
                result.plot_ax(ax1,what,components,itemp_list,fontsize=fontsize,**kwargs)


    @add_fig_kwargs
    def plot_transport(self, itemp_list=None, itau_list=None, components=('xx',),
                       erange=None, ax_array=None, fontsize=8, legend=True, **kwargs):
        """
        Plot the different quantities relevant for transport for all the results in the robot
        """
        ax_array, fig, plt = get_axarray_fig_plt(ax_array,nrows=2,ncols=2)
        self.plot_ax(ax_array[0,0],'sigma',      itemp_list=itemp_list,itau_list=itau_list,fontsize=fontsize,**kwargs)
        self.plot_ax(ax_array[0,1],'seebeck',    itemp_list=itemp_list,itau_list=itau_list,fontsize=fontsize,**kwargs)
        self.plot_ax(ax_array[1,0],'kappa',      itemp_list=itemp_list,itau_list=itau_list,fontsize=fontsize,**kwargs)
        self.plot_ax(ax_array[1,1],'powerfactor',itemp_list=itemp_list,itau_list=itau_list,fontsize=fontsize,**kwargs)

        if legend:
            for ax in ax_array.flatten(): ax.legend(loc="best", shadow=True, fontsize=fontsize)

        #fig.tight_layout()
        return fig

    @add_fig_kwargs
    def plot(self,what,itemp_list=None,itau_list=None,components=('xx',),
             erange=None,fontsize=8,legend=True,**kwargs):
        """
        Plot all the boltztrap results in the Robot

        Args:
            what: choose the quantity to plot can be: ['sigma','kappa','powerfactor']
            itemp_list: list of indexes of the tempratures to plot
            itau_list: list of indexes of the tempratures at which the lifetimes were computed
            components: Choose the components of the tensor to plot ['xx','xy','xz','yy',(...)]
            erange: choose energy range of the plot
            kwargs: Passed to ax.plot
        """
        ax1, fig, plt = get_ax_fig_plt(ax=None)
        self.plot_ax(ax1,what,components=components,itemp_list=itemp_list,itau_list=itau_list,
                     fontsize=fontsize,erange=erange,**kwargs)
        if legend: ax1.legend(loc="best", shadow=True, fontsize=fontsize)
        return fig

    @add_fig_kwargs
    def plot_dos_vvdos(self,dos_color=None,erange=None,ax_array=None,components=('xx',),fontsize=8,legend=True,**kwargs):
        """
        Plot dos and vvdos on the same figure
        """
        ax_array, fig, plt = get_axarray_fig_plt(ax_array,nrows=3)
        self.plot_dos_ax(ax_array[0],erange=erange,legend=legend,fontsize=fontsize,**kwargs)
        self.plot_vvdos_ax(ax_array[1],components=components,erange=erange,fontsize=fontsize,legend=legend)
        self.plot_vvdos_ax(ax_array[2],itau_list=range(self.ntemp),components=components,erange=erange,
                           fontsize=fontsize,legend=legend)
        return fig

    @add_fig_kwargs
    def plot_dos(self,ax=None,erange=None,fontsize=8,legend=True,**kwargs):
        """
        Plot dos for the results in the Robot
        """
        ax1, fig, plt = get_ax_fig_plt(ax=ax)
        self.plot_dos_ax(ax1,erange=erange,legend=legend,fontsize=fontsize,**kwargs)
        return fig

    @add_fig_kwargs
    def plot_vvdos(self,ax_array=None,itau_list=None,components=('xx',),erange=None,fontsize=8,legend=True,**kwargs):
        """
        Plot vvdos for all the results in the Robot
        """
        ax_array, fig, plt = get_axarray_fig_plt(ax_array=ax_array,sharex=True,nrows=2)

        self.plot_vvdos_ax(ax_array[0],components=components,erange=erange,fontsize=fontsize,legend=legend)
        self.plot_vvdos_ax(ax_array[1],itau_list=range(self.ntemp),components=components,erange=erange,
                           fontsize=fontsize,legend=legend)
        return fig

    def set_erange(self,emin,emax):
        """ Get an energy range based on an energy margin above and bellow the fermi level"""
        self.erange = (emin,emax)

    def unset_erange(self):
        """ Unset the energy range"""
        self.erange = None

    def to_string(self, verbose=0):
        """
        Return a string representation of the data in this class
        """
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__,mark="="))
        app('nresults: %d'%self.nresults)
        for result in self.results:
            app(result.to_string(mark='-'))
        return "\n".join(lines)

    def set_mumesh(self,emin,emax):
        """
        Set the range in which to plot the change of the doping
        for all the results

        Args:
            emin: minimun energy in eV
            emax: maximum energy in eV
        """
        for result in self.results:
            result.set_mumesh(emin,emax)

    def set_tmesh(self,tmesh):
        """
        Set the temperature mesh of all the results

        Args:
            tmesh: array with temperatures at which to compute the Fermi integrals
        """
        for result in self.results:
            result.set_tmesh(tmesh)

    def __str__(self):
        return self.to_string()

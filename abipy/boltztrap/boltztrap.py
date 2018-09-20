# coding: utf-8
"""
This module containes a Bolztrap2 class to interpolate and analyse the results
It also provides interfaces with Abipy objects allowing to
initialize the Boltztrap2 calculation from Abinit files
"""
import numpy as np
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

class AbipyBoltztrap():
    """
    Wrapper to Boltztrap2 interpolator
    This class contains the same quantities as the Loader classes from dft.py in Boltztrap2
    Addifionally it has methods to call the Boltztrap2 interpolator.
    It creates an instance of Bolztrap2Results to save the data
    Enter with quantities in the IBZ and interpolate to a fine BZ mesh
    """
    def __init__(self,fermi,atoms,nelect,kpoints,eig,tau=None,bstart=None,bstop=None,lpratio=1,nworkers=1):
        self.fermi = fermi
        self.atoms = atoms
        self.nelect = nelect
        self.kpoints = kpoints
        self.ebands = eig
        self.tau = tau
        self.mommat = None
        self.bstart = bstart
        self.bstop  = bstop
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

    @classmethod
    def from_sigeph(cls,sigeph):
        """Initialize interpolation of the bands and lifetimes from a sigeph object"""

        #units conversion
        eV_Ry = 2 * abu.eV_Ha
        eV_s = abu.eV_to_THz*1e12 * 2*np.pi

        #get the lifetimes as an array
        qpes = sigeph.get_qp_array(mode='ks+lifetimes')

        #get other dimensions
        bstart = sigeph.reader.max_bstart
        bstop  = sigeph.reader.min_bstop
        fermie = sigeph.ebands.fermie*eV_Ry
        atoms  = sigeph.ebands.structure.to_ase_atoms()
        volume = sigeph.ebands.structure.volume
        nelect = sigeph.ebands.nelect
        kpoints = [k.frac_coords for k in sigeph.sigma_kpoints]

        #TODO handle spin
        eig = qpes[0,:,bstart:bstop,0].T*eV_Ry

        return cls(fermie, atoms, nelect, kpoints, eig)

    def get_lattvec(self):
        import abipy.core.abinit_units as abu
        try:
            self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T / abu.Bohr_Ang
        return self.lattvec

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
        self._coefficients = fite.fitde3D(self, self.equivalences, nworkers=self.nworkers)

    def __str__(self):
        lines = []; app = lines.append
        app("nequiv: {}".format(self.nequivalences))
        app("rmesh: {}".format(self.rmesh))
        return "\n".join(lines)

class Boltztrap2Results():
    """
    Container for BoltztraP2 results
    Provides a Object oriented interface to BoltztraP2 for plotting, storing and analysing the results
    """
    def __init__(self,structure,wmesh,dos,vvdos,mumesh,tmesh):
        self.structure = structure
        self.wmesh = wmesh
        self.mumesh = mumesh
        self.tmesh = tmesh
        self.dos = dos
        self.vvdos = vvdos

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
    def kappa(self):
        if not hasattr(self,'_kappa'):
            self.compute_onsager_coefficients()
        return self._kappa
   
    def compute_fermiintegrals(self):
        """Compute and store the results of the Fermi integrals""" 
        import BoltzTraP2.bandlib as BL
        results = BL.fermiintegrals(self.wmesh, self.dos, self.vvdos, self.mumesh, self.tmesh)
        _, self._L0, self._L1, self._L2, self._Lm11 = results

    def compute_onsager_coefficients(self):
        """Compute Onsager coefficients"""
        import BoltzTraP2.bandlib as BL
        L0,L1,L2 = self.L0,self.L1,self.L2
        results = BL.calc_Onsager_coefficients(L0,L1,L2,self.mumesh,self.tmesh,self.volume)
        self._sigma, self._seebeck, self._kappa, self._hall = results

    @classmethod
    def from_file(self):
        """load results from file"""
        return cls()
 
    def write_file(self):
        """Write a file with the results from the calculation"""
        return
 
    def plot(self):
        """Plot for all the dopings as a function of temperature"""
        return

    def __str__(self):
        lines = []; app = lines.append
        return "".join(lines)

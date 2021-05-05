mport numpy as np
import json
from abipy.core.structure import Structure
from abipy.abilab import abiopen
import pandas as pd


class DeltaSCF():
    """
    Object to post-process the results from a LumiWork.
    """

    @classmethod
    def from_four_points_file(cls,filepaths):
        """ Create the object from a list of netcdf files in the order (Ag,Agstar,Aestar,Ae)"""
        energies=[]
        structures=[]
        for path in filepaths:
            with abiopen(path) as gsr_file:
                energies.append(gsr_file.energy)
                structures.append(gsr_file.structure)

        return cls(structuregs=structures[0],
                   structureex=structures[2],
                   AgEnergy=energies[0],
                   AgstarEnergy=energies[1],
                   AestarEnergy=energies[2],
                   AeEnergy=energies[3],)

    @classmethod
    def from_relax_file(cls,filepaths):
        """ Create the object from the two relaxation files (relax_gs, relax_ex).
            Give only acccess to structural relaxation induced by the transition
            and to E_zpl
        """
        energies=[]
        structures=[]
        for path in filepaths:
            with abiopen(path) as gsr_file:
                energies.append(gsr_file.energy)
                structures.append(gsr_file.structure)

        return cls(structuregs=structures[0],
                   structureex=structures[1],
                   AgEnergy=energies[0],
                   AgstarEnergy=None,
                   AestarEnergy=energies[1],
                   AeEnergy=None,)


    def __init__(self,structuregs,structureex,AgEnergy,AgstarEnergy,AestarEnergy,AeEnergy):
        """
        :param structuregs: relaxed ground state structure
        :param structureex: relaxed excited state structure
        :param AgEnergy
        :param AgstarEnergy
        :param AestarEnergy
        :param AeEnergy
        """

        self.structuregs=structuregs
        self.structureex=structureex
        self.AgEnergy=AgEnergy
        self.AgstarEnergy=AgstarEnergy
        self.AestarEnergy=AestarEnergy
        self.AeEnergy=AeEnergy

    def structure_gs(self):
        return self.structuregs

    def structure_ex(self):
        return self.structureex

    def natom(self):
        """Number of atoms in structure."""
        return len(self.structuregs)

    def symbols(self):
        symbols=[]
        elements=self.structuregs.species_by_znucl
        for element in elements:
            symbols.append(element.name)
        return symbols

    def diff_pos(self):
        """"
        Difference between gs and ex structures in Angström, (n_atoms,3) shape
        """
        return (self.structureex.cart_coords - self.structuregs.cart_coords)

    def diff_pos_sq_per_atom(self,index):
        return sum((self.diff_pos()[index])**2)


    def diff_dq_sq_per_atom(self,index):
        return self.amu_list()[index]*sum((self.diff_pos()[index])**2)

    def diff_pos_sq_per_element(self,symbol):
        s2i = self.structuregs.get_symbol2indices()
        indices=s2i[symbol]
        total=0
        for i in indices:
            total += self.diff_pos_sq_per_atom(index=i)
        return total

    def diff_dq_sq_per_element(self,symbol):
        s2i = self.structuregs.get_symbol2indices()
        indices=s2i[symbol]
        total=0
        for i in indices:
            total += self.diff_dq_sq_per_atom(index=i)
        return total

    def delta_r(self):
        d_r_squared=np.sum(self.diff_pos()**2)
        return(np.sqrt(d_r_squared))

    def amu_list(self):
        amu_list=[]
        for atom in self.structuregs.species:
            amu_list.append(atom.atomic_mass)
        return(amu_list)

    def delta_q(self,unit='atomic'):
        """
        :return: total Delta_Q in unit of (amu^1/2.Angstrom)
        """
        sq_Q_matrix = np.zeros((self.natom(), 3))
        for a in np.arange(self.natom()):
            for i in np.arange(3):
                sq_Q_matrix[a, i] = self.amu_list()[a] * self.diff_pos()[a, i] ** 2
        delta_Q = np.sqrt(np.sum(sq_Q_matrix))

        if unit == 'SI':
            return (delta_Q* 1e-10 * np.sqrt(1.66053892173E-27))
        else:
            return(delta_Q)

    def effective_mass(self):
        M = self.delta_q() ** 2 / self.delta_r() ** 2
        return (M)

    def E_zpl(self):
        return (self.AestarEnergy - self.AgEnergy)

    def E_em(self):
        return(self.AestarEnergy-self.AeEnergy)

    def E_abs(self):
        return(self.AgstarEnergy-self.AgEnergy)

    def E_FC_ex(self,unit='eV'):
        e_fc=self.AgstarEnergy - self.AestarEnergy
        if unit == 'SI':
            return 1.602176565e-19*e_fc
        else:
            return e_fc

    def E_FC_gs(self,unit='eV'):
        e_fc=self.AeEnergy - self.AgEnergy
        if unit == 'SI':
            return 1.602176565e-19*e_fc
        else:
            return e_fc

    def Stoke_shift(self):
        return(self.E_FC_ex()+self.E_FC_gs())

    def eff_freq_gs(self):
        """
        :return: effective frequency of the ground state in eV
        """
        hbar_eV = 6.582119570e-16  # in eV*s
        omega_g=np.sqrt(self.E_FC_gs(unit='SI')/(self.delta_q(unit='SI'))**2)
        return(hbar_eV*omega_g)

    def eff_freq_ex(self):
        """
        :return: effective frequency of the ground state in eV
        """
        hbar_eV = 6.582119570e-16  # in eV*s
        omega_e=np.sqrt(self.E_FC_ex(unit='SI')/(self.delta_q(unit='SI'))**2)
        return(hbar_eV*omega_e)

    def S_em(self):
        S = (self.E_FC_gs()) / (self.eff_freq_gs())
        return (S)

    def S_abs(self):
        S = (self.E_FC_ex()) / (self.eff_freq_ex())
        return (S)

    def FWHM_1D(self,T=0):
        from mpmath import coth
        w_0 = np.sqrt(8*np.log(2))*(self.S_em()/np.sqrt(self.S_abs()))*self.eff_freq_gs()
        if T==0:
            return w_0
        else :
            k_b = 8.6173303e-5# in eV/K
            w_T = w_0 * np.sqrt(coth(self.eff_freq_ex() / (2 * k_b * T)))
            return w_T


    def get_dataframe(self):
        """
        :return: pd dataframe with the main results of a lumi work : transition energies, delta Q, Huang Rhys factor,...

        """
        rows=[]
        index=[]

        units=dict([
            (r'$E_{em}$','eV'),
            (r'E$_{abs}$','eV'),
            (r'E$_{zpl}$','eV'),
            (r'E$_{FC,g}$','eV'),
            (r'E$_{FC,e}$','eV'),
            (r'$\Delta S$','eV'),
            (r'$\Delta R$','Ang'),
            (r'$\Delta Q $','(amu^{1/2}.Ang)'),
            (r'$\hbar\Omega_g$','eV'),
            (r'$\hbar\Omega_e$','eV'),
            (r'$S_{em}$','/'),
            (r'$S_{abs}$','/'),
        ])

        d=dict([
            (r'$E_{em}$',self.E_em()),
            (r'$E_{abs}$' ,self.E_abs()),
            (r'$E_{zpl}$',self.E_zpl()),
            (r'$E_{FC,g}$',self.E_FC_gs()),
            (r'$E_{FC,e}$',self.E_FC_ex()),
            (r'$\Delta S$',self.Stoke_shift()),
            (r'$\Delta R $',self.delta_r()),
            (r'$\Delta Q$',self.delta_q()),
            (r'$\hbar\Omega_g$',self.eff_freq_gs()),
            (r'$\hbar\Omega_e$',self.eff_freq_ex()),
            (r'$S_{em}$',self.S_em()),
            (r'$S_{abs}$',self.S_abs()),
        ])

        rows.append(d)
        index.append('test')

        df=pd.DataFrame(rows,index=index)

        return df

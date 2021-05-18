import numpy as np
import json
from abipy.core.structure import Structure
from abipy.abilab import abiopen
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class DeltaSCF():
    """
    Object to post-process the results from a LumiWork.
    """

    @classmethod
    def from_json_file(cls,json_path):
        """ Create the object from a json file containing the path to netcdf files)"""

        with open(json_path) as f:
            data = json.load(f)

        if 'meta' in data:
            meta = data['meta']
        else:
            meta=None

        gs_relax_path=data["gs_relax_filepath"]
        ex_relax_path=data["ex_relax_filepath"]

        with abiopen(gs_relax_path) as gsr_file:
            structure_gs=gsr_file.structure
        with abiopen(ex_relax_path) as gsr_file:
            structure_ex=gsr_file.structure

        include_four_points='Ag_gsr_filepath' in data # True if the json file contains the four points paths

        if include_four_points:
            Ag_path=data['Ag_gsr_filepath']
            Agstar_path=data['Agstar_gsr_filepath']
            Aestar_path=data['Aestar_gsr_filepath']
            Ae_path = data['Ae_gsr_filepath']

            with abiopen(Ag_path) as gsr_file:
                Ag_energy = gsr_file.energy
            with abiopen(Agstar_path) as gsr_file:
                Agstar_energy = gsr_file.energy
            with abiopen(Aestar_path) as gsr_file:
                Aestar_energy = gsr_file.energy
            with abiopen(Ae_path) as gsr_file:
                Ae_energy = gsr_file.energy

        else:
            with abiopen(gs_relax_path) as gsr_file:
                Ag_energy = gsr_file.energy
            with abiopen(ex_relax_path) as gsr_file:
                Agstar_energy = gsr_file.energy

            Ae_energy=None
            Agstar_energy=None


        return cls(structuregs=structure_gs,
                   structureex=structure_ex,
                   AgEnergy=Ag_energy,
                   AgstarEnergy=Agstar_energy,
                   AestarEnergy=Aestar_energy,
                   AeEnergy=Ae_energy,
                   meta=meta)


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


    def __init__(self,structuregs,structureex,AgEnergy,AgstarEnergy,AestarEnergy,AeEnergy,meta=None):
        """
        :param structuregs: relaxed ground state structure
        :param structureex: relaxed excited state structure
        :param AgEnergy
        :param AgstarEnergy
        :param AestarEnergy
        :param AeEnergy
        :param meta : dict. of meta data of the lumiwork (can be the supercell size, ecut,...)
        """

        self.structuregs=structuregs
        self.structureex=structureex
        self.AgEnergy=AgEnergy
        self.AgstarEnergy=AgstarEnergy
        self.AestarEnergy=AestarEnergy
        self.AeEnergy=AeEnergy
        self.meta=meta


    def get_meta(self):
        return self.meta

    def structure_gs(self):
        return self.structuregs

    def structure_ex(self):
        return self.structureex

    def natom(self):
        """Number of atoms in structure."""
        return len(self.structuregs)

    def diff_pos(self):
        """"
        Difference between gs and ex structures in Angström, (n_atoms,3) shape
        """
        return (self.structureex.cart_coords - self.structuregs.cart_coords)

    def defect_index(self,defect_symbol):
        index=self.structuregs.get_symbol2indices()[defect_symbol][0]
        return index

    def get_dict_per_atom(self,index,defect_symbol):
        stru=self.structuregs
        def_index=self.defect_index(defect_symbol=defect_symbol)
        d={}
        d["symbol"]=stru.species[index].name
        d["mass"]=self.amu_list()[index]
        d[r"$\Delta R^2$"]=sum(self.diff_pos()[index]**2)
        d[r"$\Delta Q^2$"]=self.amu_list()[index]*sum(self.diff_pos()[index]**2)
        d["dist. from defect"]=stru[index].distance(other=stru[def_index])

        return d


    def get_dataframe_atoms(self,defect_symbol):
        list_of_dict=[]
        for index,atom in enumerate(self.structuregs):
            d=self.get_dict_per_atom(index,defect_symbol)
            list_of_dict.append(d)

        return pd.DataFrame(list_of_dict)

    def get_dict_per_element(self,element):
        stru=self.structuregs
        indices=stru.indices_from_symbol(element.name)
        dr_el=[]
        for index in indices:
            dr_el.append(sum(self.diff_pos()[index]**2))
        d={}
        d["symbol"]=element.name
        d["mass"]=element.atomic_mass
        d[r"$\Delta R^2$"]=sum(dr_el)
        d[r"$\Delta Q^2$"]=element.atomic_mass*sum(dr_el)

        return d


    def get_dataframe_element(self):
        list_of_dict=[]
        for index,element in enumerate(self.structuregs.types_of_species):
            d=self.get_dict_per_element(element)
            list_of_dict.append(d)

        return pd.DataFrame(list_of_dict)


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
        omega_g=np.sqrt(2*self.E_FC_gs(unit='SI')/(self.delta_q(unit='SI'))**2)
        return(hbar_eV*omega_g)

    def eff_freq_ex(self):
        """
        :return: effective frequency of the ground state in eV
        """
        hbar_eV = 6.582119570e-16  # in eV*s
        omega_e=np.sqrt(2*self.E_FC_ex(unit='SI')/(self.delta_q(unit='SI'))**2)
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

    def get_dict_results(self):
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
        return d


    def get_dataframe(self,label=None):
        """
        :return: pd dataframe with the main results of a lumi work : transition energies, delta Q, Huang Rhys factor,...

        """
        rows=[]
        index=[]
        d=self.get_dict_results()

 #       units=dict([
 #           (r'$E_{em}$','eV'),
 #           (r'E$_{abs}$','eV'),
 #           (r'E$_{zpl}$','eV'),
 #          (r'E$_{FC,g}$','eV'),
 #           (r'E$_{FC,e}$','eV'),
 #           (r'$\Delta S$','eV'),
 #           (r'$\Delta R$','Ang'),
 #           (r'$\Delta Q $','(amu^{1/2}.Ang)'),
 #           (r'$\hbar\Omega_g$','eV'),
 #           (r'$\hbar\Omega_e$','eV'),
 #           (r'$S_{em}$','/'),
 #           (r'$S_{abs}$','/'),
 #       ])

        rows.append(d)
        index.append(label)

        df=pd.DataFrame(rows,index=index)

        return df

    def displacements_visu(self):
        pos_gs=self.structuregs.cart_coords
        pos_ex=self.structureex.cart_coords
        # make the grid
        x = pos_ex[:, 0]
        y = pos_ex[:, 1]
        z = pos_ex[:, 2]

        # Make the direction data for the arrows
        u =(pos_ex[:, 0]-pos_gs[:, 0])
        v =(pos_ex[:, 1]-pos_gs[:, 1])
        w =(pos_ex[:, 2]-pos_gs[:, 2])

        M = self.amu_list()*(u**2+v**2+w**2)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        a_g=10 ### vector size increase coef
        ax.quiver(x, y, z,u*a_g,v*a_g,w*a_g, color='k',linewidths=1)
        sc = ax.scatter(x, y, z, c=M, marker='o', s=60, cmap="jet")
        #ax.scatter(self.pos_gs[0, 0], self.pos_gs[0, 1], self.pos_gs[0, 2], marker='o', s=200, color='fuchsia')

        #set plot to scale
       # a=0.37
       # ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1*a,1*a,(90/15)*a,1]))
        clb=plt.colorbar(sc)
        clb.set_label(r'$\Delta Q^2$ per atom')


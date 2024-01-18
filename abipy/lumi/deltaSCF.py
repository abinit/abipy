from __future__ import annotations

import numpy as np
import json
from abipy.core.structure import Structure
from abipy.abilab import abiopen
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from mpmath import coth
from scipy.integrate import simps
from scipy.special import factorial
from abipy.tools.plotting import get_ax_fig_plt, add_fig_kwargs,get_axarray_fig_plt
import abipy.core.abinit_units as abu
import os,shutil 

class DeltaSCF():
    """
    Object to post-process the results from a LumiWork, following a one-effective phonon mode model (1D-CCM).
    For equations, notations and formalism, please refer to :
     https://doi.org/10.1103/PhysRevB.96.125132
     https://doi.org/10.1002/adom.202100649
    """

    @classmethod
    def from_json_file(cls,json_path):
        """ Create the object from a json file containing the path to netcdf files, produced at the end of a LumiWork"""

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
                forces_ex=gsr_file.cart_forces
            with abiopen(Aestar_path) as gsr_file:
                Aestar_energy = gsr_file.energy
            with abiopen(Ae_path) as gsr_file:
                Ae_energy = gsr_file.energy
                forces_gs=gsr_file.cart_forces


        else:
            with abiopen(gs_relax_path) as gsr_file:
                Ag_energy = gsr_file.energy
            with abiopen(ex_relax_path) as gsr_file:
                Agstar_energy = gsr_file.energy

            Ae_energy=None
            Agstar_energy=None


        return cls(structuregs=structure_gs,
                   structureex=structure_ex,
                   forces_gs=forces_gs,
                   forces_ex=forces_ex,
                   ag_energy=Ag_energy,
                   ag_star_energy=Agstar_energy,
                   ae_star_energy=Aestar_energy,
                   ae_energy=Ae_energy,
                   meta=meta)

    @classmethod
    def from_four_points_file(cls,filepaths):
        """ 
            Create the object from a list of netcdf files in the order (Ag,Agstar,Aestar,Ae).
            Ag: Ground state at relaxed ground state atomic positions. 
            Agstar: Excited state at relaxed ground state atomic positions.
            Aestar: Excited state at relaxed excited state atomic positions.
            Ae: Ground state at relaxed excited state atomic positions.
    
        Args:
            filepaths: list of netcdf files in the order [Ag,Agstar,Aestar,Ae]
        
        Returns:
            A DeltaSCF object

        """
        energies=[]
        structures=[]
        forces=[]
        for path in filepaths:
            with abiopen(path) as gsr_file:
                energies.append(gsr_file.energy)
                structures.append(gsr_file.structure)
                forces.append(gsr_file.cart_forces)
        return cls(structuregs=structures[0],
                   structureex=structures[2],
                   forces_gs=forces[3],
                   forces_ex=forces[1],
                   ag_energy=energies[0],
                   ag_star_energy=energies[1],
                   ae_star_energy=energies[2],
                   ae_energy=energies[3],)

    @classmethod
    def from_relax_file(cls,filepaths):
        """ Create the object from the two relaxation files (relax_gs, relax_ex).
            Give only acccess to structural relaxation induced by the transition
            and to E_zpl
        """
        energies=[]
        structures=[]
        forces=[]
        for path in filepaths:
            with abiopen(path) as gsr_file:
                energies.append(gsr_file.energy)
                structures.append(gsr_file.structure)
                forces.append(gsr_file.cart_forces)

        return cls(structuregs=structures[0],
                   structureex=structures[1],
                   forces_gs=None,
                   forces_ex=None,
                   ag_energy=energies[0],
                   ag_star_energy=None,
                   ae_star_energy=energies[1],
                   ae_energy=None,)


    def __init__(self,structuregs,structureex,forces_gs,forces_ex,
                 ag_energy,ag_star_energy,ae_star_energy,ae_energy,meta=None):
        """
        :param structuregs: relaxed ground state structure
        :param structureex: relaxed excited state structure
        :param forces_gs: forces in the gs
        :param forces_ex: forces in the ex
        :param ag_energy
        :param ag_star_energy
        :param ae_star_energy
        :param ae_energy
        :param meta : dict. of meta data of the lumiwork (can be the supercell size, ecut,...)
        """

        self.structuregs=structuregs
        self.structureex=structureex
        self.forces_gs=forces_gs
        self.forces_ex=forces_ex
        self.ag_energy=ag_energy
        self.ag_star_energy=ag_star_energy
        self.ae_star_energy=ae_star_energy
        self.ae_energy=ae_energy
        self.meta=meta


    def structure_gs(self):
        """ Ground state relaxed structure """

        return self.structuregs

    def structure_ex(self):
        """ Excited state relaxed structure """

        return self.structureex

    def natom(self):
        """Number of atoms in the structure."""
        return len(self.structuregs)

    def diff_pos(self):
        """
        Difference between gs and ex structures in Angström, (n_atoms,3) shape
        """
        return (self.structureex.cart_coords - self.structuregs.cart_coords)
    
    def diff_pos_mass_weighted(self):
        """
        Difference between gs and ex structures in Angström, weighted by the squared atomic masses (n_atoms,3) shape
        """
        return np.einsum('i,ij->ij',np.array(self.amu_list()),self.diff_pos())

    def defect_index(self,defect_symbol):
        """
        Defect index in the structure from its symbol, ex defect_index("Eu").
        """
        index=self.structuregs.get_symbol2indices()[defect_symbol][0]
        return index

    def get_dict_per_atom(self,index,defect_symbol):
        """"
        Dict. with relevant properties per atom.
        """
        stru=self.structuregs
        def_index=self.defect_index(defect_symbol=defect_symbol)
        d={}
        d["symbol"]=stru.species[index].name
        d["mass "]=self.amu_list()[index]
        d[r"$\Delta R$"]=(sum(self.diff_pos()[index]**2))**(0.5)
        d[r"$\Delta Q^2$"]=self.amu_list()[index]*sum(self.diff_pos()[index]**2)
        d[r"$\Delta F$"]=(sum(self.forces_gs[index]**2))**(0.5)
        d["dist. from defect"]=stru[index].distance(other=stru[def_index])

        return d

    def get_dataframe_atoms(self,defect_symbol):
        """
        Panda dataframe with relevant properties per atom.
        Units : [ (mass,amu), (deltaR,Angstrom), (DeltaQ^2,amu.Angstrom^2), (DeltaF,eV/Angstrom) ] 
        """
        list_of_dict=[]
        for index,atom in enumerate(self.structuregs):
            d=self.get_dict_per_atom(index,defect_symbol)
            list_of_dict.append(d)

        return pd.DataFrame(list_of_dict)

    def get_dict_per_specie(self,specie):
        stru=self.structuregs
        indices=stru.indices_from_symbol(specie.name)
        dr_sp=[]
        for index in indices:
            dr_sp.append(sum(self.diff_pos()[index]**2))
        d={}
        d["symbol"]=specie.name
        d["mass"]=specie.atomic_mass
        d[r"$\Delta R^2$"]=sum(dr_sp)
        d[r"$\Delta Q^2$"]=specie.atomic_mass*sum(dr_sp)

        return d

    def get_dataframe_species(self):
        """
        Panda dataframe with relevant properties per species.
        """
        list_of_dict=[]
        for index,specie in enumerate(self.structuregs.types_of_species):
            d=self.get_dict_per_specie(specie)
            list_of_dict.append(d)

        return pd.DataFrame(list_of_dict)


    def delta_r(self):
        """
        Total Delta_R (Angstrom)
        """
        d_r_squared=np.sum(self.diff_pos()**2)
        return(np.sqrt(d_r_squared))

    def amu_list(self):
        """
        List of masses (amu)
        """
        amu_list=[]
        for atom in self.structuregs.species:
            amu_list.append(atom.atomic_mass)
        return(amu_list)

    def delta_q(self,unit='atomic'):
        """
        Total Delta_Q

        Args:
            unit: amu^1/2.Angstrom if unit = 'atomic', kg^1/2.m if 'SI'

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
        """
        Effective mass M
        """
        M = self.delta_q() ** 2 / self.delta_r() ** 2
        return (M)

    def E_zpl(self):
        """
        Zero-phonon line energy (eV)
        """
        return (self.ae_star_energy - self.ag_energy)

    def E_em(self):
        """
        Emission energy(eV)
        """
        return(self.ae_star_energy-self.ae_energy)

    def E_abs(self):
        """
        Absorption energy(eV)
        """
        return(self.ag_star_energy-self.ag_energy)

    def E_FC_ex(self,unit='eV'):
        """
        Franck condon energy in excited state (eV)
        = Relaxation energy between Ag* and Ae* states
        """
        e_fc=self.ag_star_energy - self.ae_star_energy
        if unit == 'SI':
            return 1.602176565e-19*e_fc
        else:
            return e_fc

    def E_FC_gs(self,unit='eV'):
        """
        Franck condon energy in ground state (eV)
        = Relaxation energy between Ae and Ag states
        """
        e_fc=self.ae_energy - self.ag_energy
        if unit == 'SI':
            return 1.602176565e-19*e_fc
        else:
            return e_fc

    def Stoke_shift(self):
        """
        Stokes shift (eV)
        """
        return(self.E_FC_ex()+self.E_FC_gs())

    def eff_freq_gs(self):
        """
        Phonon effective frequency of the ground state (eV)
        """
        omega_g=np.sqrt(2*self.E_FC_gs(unit='SI')/(self.delta_q(unit='SI'))**2)
        return(abu.hbar_eVs*omega_g)

    def eff_freq_ex(self):
        """
        Phonon effective frequency of the excited state (eV)
        """
        omega_e=np.sqrt(2*self.E_FC_ex(unit='SI')/(self.delta_q(unit='SI'))**2)
        return(abu.hbar_eVs*omega_e)

    def S_em(self):
        """
        Total Huang-Rhys factor for emission following the 1D-CCM
        """
        S = (self.E_FC_gs()) / (self.eff_freq_gs())
        return (S)

    def S_abs(self):
        """
        Total Huang-Rhys factor for absorption following the 1D-CCM
        """
        S = (self.E_FC_ex()) / (self.eff_freq_ex())
        return (S)

    def FWHM_1D(self,T=0):
        """
        Full width at half-maximum following a semi-classical approx in the 1D-CCM
        (eq.20-21 of https://doi.org/10.1103/PhysRevB.96.125132)

        Args:
            T: Temperature
        """
        w_0 = np.sqrt(8*np.log(2))*(self.S_em()/np.sqrt(self.S_abs()))*self.eff_freq_gs()
        if T==0:
            return w_0
        else :
            k_b = abu.kb_eVK
            w_T = w_0 * np.sqrt(coth(self.eff_freq_ex() / (2 * k_b * T)))
            return w_T

    def FC_factor_approx(self,n):
        """
        FC factor between initial vib. state m=0 to final vib. state n
        approx of same eff frequency in gs and ex and T = 0K.
        See eq. (9) of https://doi.org/10.1002/adom.202100649
        """
        FC = np.exp(self.S_em()) * self.S_em() ** n / math.factorial(n)
        return FC

    def lineshape_1D_zero_temp(self,energy_range=[0.5,5],max_m=25,phonon_width=0.01,with_omega_cube=True,normalized='Area'):
        """
        Compute the emission lineshape following the effective phonon 1D-CCM at T=0K.
        See eq. (9) of  https://doi.org/10.1002/adom.202100649.
        
        Args:
            energy_range:  Energy range at which the intensities are computed, ex : [0.5,5]
            max_m: Maximal vibrational state m considered
            phonon_width: fwhm of each phonon peak, in eV
            with_omega_cube: Considered or not the omega^3 dependence of the intensity
            normalized: Normalisation procedure. 'Area' if Area under the curve = 1. 'Sum' if maximum of the curve = 1.
                       
        Returns:
            E_x = Energies at which the intensities are computed
            I   = Intensities
        """
        n_x = 10000  #
        E_x = np.linspace(energy_range[0], energy_range[1], n_x)

        list_n = np.arange(0, max_m)
        A = np.zeros(n_x)
        sigma = phonon_width / (2.35482)
        for n in list_n:
            #gaussian_1D = np.zeros(n_x)
            fc_factor=self.FC_factor_approx(n)
            #f=np.exp(self.S_em())*self.S_em()**m[i]/math.factorial(m[i])
            arg_exp=-((self.E_zpl() - self.eff_freq_gs() * n - E_x) ** 2 / (2 * (sigma) ** 2))
            gaussian_1D = fc_factor * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(arg_exp)
            A += gaussian_1D

        if with_omega_cube==True:
            A=A*E_x**3

        if normalized=="Area":
            C = 1 / (simps(A, E_x))
        if normalized=="Sum":
            C=1/(max(A))

        return E_x, C*A

    @add_fig_kwargs
    def plot_lineshape_1D_zero_temp(self,energy_range=[0.5,5],max_m=25,phonon_width=0.01,with_omega_cube="True",
                                    normalized='Area', ax=None, **kwargs):
        """
        Plot the the emission lineshape following the effective phonon 1D-CCM at T=0K.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            energy_range:  Energy range at which the intensities are computed, ex : [0.5,5]
            max_m: Maximal vibrational state m considered
            phonon_width: fwhm of each phonon peak, in eV
            with_omega_cube: Considered or not the omega^3 dependence of the intensity
            normlized: Normalisation procedure. 'Area' if Area under the curve = 1. 'Sum' if maximum of the curve = 1.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        x,y=self.lineshape_1D_zero_temp(energy_range,max_m,phonon_width,with_omega_cube,
                                        normalized)
        ax.plot(x,y, **kwargs)
        ax.set_xlabel(r'Energy (eV)')
        ax.set_ylabel(r'Intensity ')
        return fig


    def get_dict_results(self):
        d=dict([
            (r'E_em',self.E_em()),
            (r'E_abs' ,self.E_abs()),
            (r'E_zpl',self.E_zpl()),
            (r'E_fc_gs',self.E_FC_gs()),
            (r'E_fc_ex',self.E_FC_ex()),
            (r'Delta_S',self.Stoke_shift()),
            (r'Delta_R ',self.delta_r()),
            (r'Delta_Q',self.delta_q()),
            (r'Eff_freq_gs',self.eff_freq_gs()),
            (r'Eff_freq_ex',self.eff_freq_ex()),
            (r'S_em',self.S_em()),
            (r'S_abs',self.S_abs()),
        ])
        return d

    def get_dataframe(self,label=None):
        """
        Panda dataframe with the main results of a LumiWork : transition energies, delta Q, Huang Rhys factor,...
        Units used are Angstrom, eV, amu.
        DeltaSCF object should be instantiated with the four points files, not with relax files only.
        """
        rows=[]
        index=[]
        d=self.get_dict_results()
        rows.append(d)
        index.append(label)

        df=pd.DataFrame(rows,index=index)

        return df
    
    def draw_displacements_vesta(self,in_path, mass_weighted = False,
                 scale_vector=20,width_vector=0.3,color_vector=[255,0,0],centered=True,
                 factor_keep_vectors=0.1,
                 out_path="VESTA_FILES",out_filename="gs_ex_relaxation"):
        """
        Draw the ground state to excited state atomic relaxation on a vesta structure. 

        Args:
            in_path : path where the initial .vesta structure in stored, should correspond to the ground state relaxed structure. 
            mass_weighted : If True, weight the displacements by the atomic masses. Draw the \Delta Q in that case.  
            scale_vector : scaling factor of the vector modulus
            width_vector : vector width
            color_vector : color in rgb format
            centered : center the vesta structure around [0,0,0]
            factor_keep_vectors : draw only the eigenvectors with magnitude > factor_keep_vectors * max(magnitude)
            out_path : path where .vesta files with vector are stored
        """

        vesta = open(in_path,'r').read()
        natoms = len(self.structure_gs())


        if os.path.isdir(out_path): 
            shutil.rmtree(out_path)
            os.mkdir(out_path)
        else : 
            os.mkdir(out_path)

        path=out_path

        towrite = vesta.split('VECTR')[0]
        towrite += 'VECTR\n'

        magnitudes=[]
        displacements=self.diff_pos()
        if mass_weighted == True:
            displacements=self.diff_pos_mass_weighted()

        for iatom in range(natoms) :
            magnitudes.append(np.sqrt(displacements[iatom][0]**2 + displacements[iatom][1]**2 + displacements[iatom][2]**2))

        for iatom in range(natoms) :
            if magnitudes[iatom] > factor_keep_vectors * max(np.real(magnitudes)):
                towrite += '%5d' %(iatom + 1)
                towrite += '%10.5f' %(displacements[iatom][0] * (scale_vector))
                towrite += '%10.5f' %(displacements[iatom][1] * (scale_vector))
                towrite += '%10.5f' %(displacements[iatom][2] * (scale_vector))
                towrite += '\n'
                towrite += '%5d' %(iatom + 1)  +  ' 0 0 0 0\n  0 0 0 0 0\n'
    
        towrite += '0 0 0 0 0\n' 
        towrite += 'VECTT\n'

        for atom in range(natoms) :
            towrite += '%5d' %(atom + 1)
            towrite += f'  {width_vector} {color_vector[0]}   {color_vector[1]}   {color_vector[2]} 0\n'

        towrite += '0 0 0 0 0\n' 
        towrite += 'SPLAN'
        towrite += vesta.split('SPLAN')[1]
        towrite += 'VECTS 1.00000'


        filename = path + '/'+out_filename
        filename += '.vesta'

        open(filename, 'w').write(towrite)

        if centered==True:

            with open(filename, 'r') as file:
                file_contents = file.read()
                search_word="BOUND\n       0        1         0        1         0        1\n  0   0   0   0  0"
                replace_word="BOUND\n       -0.5        0.5         -0.5        0.5         -0.5        0.5\n  0   0   0   0  0"

                updated_contents = file_contents.replace(search_word, replace_word)

            with open(filename, 'w') as file:
                file.write(updated_contents)
                
        print(f"Vesta files created and stored in : \n {os.getcwd()}/{out_path}")

        return 
     


    @add_fig_kwargs
    def displacements_visu(self,a_g=10,**kwargs):
        """
        Make a 3d visualisation of the displacements induced by the electronic transition =
        Difference between ground state and excited state atomic positions.
        The colors of the atoms are based on Delta_Q_^2 per atom.

        Args:
            a_g : coefficient that multiplies the displacement magnitudes

        Returns: |matplotlib-Figure|
        """
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

        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.quiver(x, y, z,u*a_g,v*a_g,w*a_g, color='k',linewidths=1,**kwargs)
        sc = ax.scatter(x, y, z, c=M, marker='o', s=60, cmap="jet",**kwargs)

        clb=plt.colorbar(sc)
        clb.set_label(r'$\Delta Q^2$ per atom')
        return fig


    @add_fig_kwargs
    def plot_delta_R_distance(self, defect_symbol,colors=["k","r","g","b","c","m"],ax=None, **kwargs):
        """
        Plot \DeltaR vs distance from defect for each atom, colored by species.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            defect_symbol:  defect_symbol, defect location will be the reference
            colors: list of colors for the species

        Returns: |matplotlib-Figure|
        """

        symbols=self.structuregs.symbol_set
        dfs = []
        xs = []
        ys = []
        df=self.get_dataframe_atoms(defect_symbol=defect_symbol)

        for i, symbol in enumerate(symbols):
            dfs.append(df.loc[df['symbol'] == symbol])
            xs.append(dfs[i]["dist. from defect"])
            ys.append(dfs[i]["$\\Delta R$"])

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for i, symbol in enumerate(symbols):
            ax.stem(xs[i], ys[i], label=symbol, linefmt=colors[i], markerfmt="o" + colors[i],**kwargs)
            ax.set_xlabel(r'Distance from defect ($\AA$)')
            ax.set_ylabel(r'$\Delta R $ ($\AA$)')
            ax.legend()

        return fig

    @add_fig_kwargs
    def plot_delta_F_distance(self, defect_symbol,colors=["k","r","g","b","c","m"],ax=None, **kwargs):
        """
        Plot \DeltaF vs distance from defect for each atom, colored by species.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            defect_symbol:  defect_symbol, defect location will be the reference
            colors: list of colors for the species

        Returns: |matplotlib-Figure|
        """

        symbols=self.structuregs.symbol_set
        dfs = []
        xs = []
        ys = []
        df=self.get_dataframe_atoms(defect_symbol=defect_symbol)

        for i, symbol in enumerate(symbols):
            dfs.append(df.loc[df['symbol'] == symbol])
            xs.append(dfs[i]["dist. from defect"])
            ys.append(dfs[i]["$\\Delta F$"])

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for i, symbol in enumerate(symbols):
            ax.stem(xs[i], ys[i], label=symbol, linefmt=colors[i], markerfmt="o" + colors[i],**kwargs)
            ax.set_xlabel(r'Distance from defect ($\AA$)')
            ax.set_ylabel(r'$\Delta F$ ($eV/\AA$)')
            ax.legend()

        return fig

    @add_fig_kwargs
    def plot_four_BandStructures(self,nscf_files,ax_mat=None,ylims=[-5,5],**kwargs):
        """
        plot the 4 band structures.
        nscf_files is the list of Ag, Agstar, Aestar, Ae nscf gsr file paths.
        """
        ebands = []
        for file in nscf_files:
            with abiopen(file) as f:
                ebands.append(f.ebands)

        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=1, ncols=4,
                                               sharex=True, sharey=True, squeeze=False)

        titles = [r'$A_g$', r'$A_g^*$', r'$A_e^*$', r'$A_e$']
        e0 = ebands[0].fermie

        for i,eband in enumerate(ebands):
            eband.plot_ax(ax=ax_mat[0,i],spin=0, e0=e0,color="k",**kwargs)
            eband.plot_ax(ax=ax_mat[0,i],spin=1, e0=e0,color="r",**kwargs)
            eband.decorate_ax(ax=ax_mat[0,i],title=titles[i])

        ax_mat[0,0].set_ylim(ylims)
        ax_mat[0,1].set_ylabel("")
        ax_mat[0,2].set_ylabel("")
        ax_mat[0,3].set_ylabel("")

        return fig

    @add_fig_kwargs
    def draw_displaced_parabolas(self,ax=None,scale_eff_freq=4,font_size=8, **kwargs):
        """
        Draw the four points diagram with relevant transition energies.

        Args:
        ax: |matplotlib-Axes| or None if a new figure should be created.
        scale_eff_freq:  scaling factor to adjust the parabolas curvatures.
        font_size: font size for the annotations

        Returns: |matplotlib-Figure|
        """
        ax,fig,plt=get_ax_fig_plt(ax=ax)

        delta_Q=self.delta_q()
        E_zpl=self.E_zpl()
        omega_gs_sq=scale_eff_freq*2*self.E_FC_gs()/self.delta_q()**2
        omega_ex_sq=scale_eff_freq*2*self.E_FC_ex()/self.delta_q()**2

        new_FC_gs=omega_gs_sq*delta_Q**2*0.5
        new_FC_ex=omega_ex_sq*delta_Q**2*0.5

        Qs=np.linspace(-delta_Q*0.2,delta_Q*1.5,1000)

        E_gs=0.5*omega_gs_sq.real*(Qs)**2+0 # ref at (0,0)
        E_ex=0.5*omega_ex_sq.real*(Qs-delta_Q)**2+ self.E_zpl()# min at (delta_Q,ae_energy)


        #  parabolas
        ax.plot(Qs,E_gs,'k',zorder=1)
        ax.plot(Qs,E_ex,'k',zorder=1)

        #  points
        xs=np.array([0,0,delta_Q,delta_Q])
        ys=np.array([0,E_zpl+new_FC_ex,E_zpl,new_FC_gs])

        ax.scatter(xs,ys,s=50,color='k',zorder=2)

        # arrows

        ax.annotate("", xy=(0, E_zpl+0.95*new_FC_ex), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->",color="b",lw=1))
        ax.annotate(r' $E_{abs}$='+format(self.E_abs(),".2f")+' eV  ', xy=(0,(E_zpl+new_FC_ex)/2),ha='left',fontsize=font_size)

        ax.annotate("", xy=(delta_Q, new_FC_gs*1.05), xytext=(delta_Q, E_zpl),
            arrowprops=dict(arrowstyle="->",color="r",lw=1))
        ax.annotate(r' $E_{em}$='+format(self.E_em(),".2f")+' eV  ', xy=(delta_Q,E_zpl-(E_zpl-new_FC_gs)/2),ha='left',fontsize=font_size)

        ax.annotate("", xy=(delta_Q, E_zpl), xytext=(delta_Q, E_zpl+new_FC_ex*1.5),
            arrowprops=dict(arrowstyle="-",color="k",lw=0.3,ls='--'))
        ax.annotate(r' $E_{FC,e}$='+format(self.E_FC_ex(),".2f")+' eV  ', xy=(delta_Q,E_zpl+new_FC_ex/2),ha='left',fontsize=font_size)

        ax.annotate("", xy=(delta_Q, new_FC_gs), xytext=(delta_Q, -new_FC_gs*0.5),
            arrowprops=dict(arrowstyle="-",color="k",lw=0.3,ls='--'))
        ax.annotate(r' $E_{FC,g}$='+format(self.E_FC_gs(),".2f")+' eV  ', xy=(delta_Q,new_FC_gs/2),ha='left',fontsize=font_size)

        ax.annotate("", xy=(0, 0), xytext=(delta_Q*1.1, 0),
            arrowprops=dict(arrowstyle="-",color="k",lw=0.3,ls='--'))
        ax.annotate("", xy=(0, E_zpl+new_FC_ex), xytext=(delta_Q*1.1, E_zpl+new_FC_ex),
            arrowprops=dict(arrowstyle="-",color="k",lw=0.3,ls='--'))

        ax.annotate("", xy=(0, -new_FC_gs*0.2), xytext=(delta_Q, -new_FC_gs*0.2),
            arrowprops=dict(arrowstyle="<->",color="k",lw=0.6))
        ax.annotate(r'$\Delta Q$ ='+format(self.delta_q(),".2f"), xy=(delta_Q/2, -new_FC_gs*0.4),ha='center',fontsize=font_size)

        ax.set_ylim(-new_FC_gs*1.5,E_zpl+2*new_FC_ex)
        ax.set_xlim(-0.5*delta_Q,2*delta_Q)

        ax.annotate("", xy=(-0.4*delta_Q, -new_FC_gs), xytext=(-0.4*delta_Q, E_zpl+2*new_FC_ex),
            arrowprops=dict(arrowstyle="<-",color="k",lw=1.5))
        ax.text(x=-0.45*delta_Q,y=(E_zpl+new_FC_ex)/2, s='Energy (eV)',fontsize=10,rotation=90,ha='center')

        ax.annotate("", xy=(-0.5*delta_Q, -new_FC_gs*0.6), xytext=(+1.4*delta_Q, -new_FC_gs*0.6),
            arrowprops=dict(arrowstyle="<-",color="k",lw=1.5))
        ax.text(x=0.7*delta_Q,y=-new_FC_gs, s='Configuration coordinate Q',fontsize=10,ha='center')


        ax.axis('off')

        return fig


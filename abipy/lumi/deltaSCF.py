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
from abipy.tools.plotting import get_ax_fig_plt, add_fig_kwargs


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
                   AgEnergy=Ag_energy,
                   AgstarEnergy=Agstar_energy,
                   AestarEnergy=Aestar_energy,
                   AeEnergy=Ae_energy,
                   meta=meta)


    @classmethod
    def from_four_points_file(cls,filepaths):
        """ Create the object from a list of netcdf files in the order (Ag,Agstar,Aestar,Ae)
            Ag     -->  Ground state at relaxed ground state atomic positions
            Agstar -->  Excited state at relaxed ground state atomic positions
            Aestar -->  Excited state at relaxed excited state atomic positions
            Ae     -->  Ground state at relaxed excited state atomic positions
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
                   AgEnergy=energies[0],
                   AgstarEnergy=None,
                   AestarEnergy=energies[1],
                   AeEnergy=None,)


    def __init__(self,structuregs,structureex,forces_gs,forces_ex,
                 AgEnergy,AgstarEnergy,AestarEnergy,AeEnergy,meta=None):
        """
        :param structuregs: relaxed ground state structure
        :param structureex: relaxed excited state structure
        :param forces_gs: forces in the gs
        :param forces_ex: forces in the ex
        :param AgEnergy
        :param AgstarEnergy
        :param AestarEnergy
        :param AeEnergy
        :param meta : dict. of meta data of the lumiwork (can be the supercell size, ecut,...)
        """

        self.structuregs=structuregs
        self.structureex=structureex
        self.forces_gs=forces_gs
        self.forces_ex=forces_ex
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
        """"
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
        units used are 
        mass     - amu
        DeltaR   - Angstrom
        DeltaQ^2 - amu.Angstrom^2
        DeltaF   - eV/Angstrom
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
        units used are 
        mass     - amu
        DeltaR   - Angstrom
        DeltaQ^2 - amu.Angstrom^2
        DeltaF   - eV/Angstrom
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
        return (self.AestarEnergy - self.AgEnergy)

    def E_em(self):
        """
        Emission energy(eV)
        """
        return(self.AestarEnergy-self.AeEnergy)

    def E_abs(self):
        """
        Absorption energy(eV)
        """
        return(self.AgstarEnergy-self.AgEnergy)

    def E_FC_ex(self,unit='eV'):
        """
        Franck condon energy in excited state (eV)
        = Relaxation energy between Ag* and Ae* states
        """
        e_fc=self.AgstarEnergy - self.AestarEnergy
        if unit == 'SI':
            return 1.602176565e-19*e_fc
        else:
            return e_fc

    def E_FC_gs(self,unit='eV'):
        """
        Franck condon energy in ground state (eV)
        = Relaxation energy between Ae and Ag states
        """
        e_fc=self.AeEnergy - self.AgEnergy
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
        hbar_eV = 6.582119570e-16  # in eV*s
        omega_g=np.sqrt(2*self.E_FC_gs(unit='SI')/(self.delta_q(unit='SI'))**2)
        return(hbar_eV*omega_g)

    def eff_freq_ex(self):
        """
        Phonon effective frequency of the excited state (eV)
        """
        hbar_eV = 6.582119570e-16  # in eV*s
        omega_e=np.sqrt(2*self.E_FC_ex(unit='SI')/(self.delta_q(unit='SI'))**2)
        return(hbar_eV*omega_e)

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
            k_b = 8.6173303e-5# in eV/K
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

    def lineshape_1D_zero_temp(self,energy_range=[0.5,5],max_m=25,phonon_width=0.01,with_omega_cube="True",normalized='Area'):
        """
        Compute the emission lineshape following the effective phonon 1D-CCM at T=0K.
        See eq. (9) of  https://doi.org/10.1002/adom.202100649 
        Args:
            energy_range:  Energy range at which the intensities are computed, ex : [0.5,5]
            max_m: Maximal vibrational state m considered
            phonon_width: fwhm of each phonon peak, in eV
            with_omega_cube: Considered or not the omega^3 dependence of the intensity
            normlized: Normalisation procedure. 'Area' if Area under the curve = 1
                       'Sum' if maximum of the curve = 1.
        Returns:
            E_x = Energies at which the intensities are computed, np.array
            I   = Intensities, np.array
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

        if with_omega_cube=="True":
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
        normlized: Normalisation procedure. 'Area' if Area under the curve = 1
                   'Sum' if maximum of the curve = 1.
        
        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        x,y=self.lineshape_1D_zero_temp(energy_range,max_m,phonon_width,with_omega_cube,
                                        normalized)
        ax.plot(x,y, **kwargs)
        return fig

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
    
    @add_fig_kwargs
    def displacements_visu(self,a_g=10,**kwargs):
        """
        Make a 3d visualisation of the displacements induced by the electronic transition =
        Difference between ground state and excited state atomic positions.
        The colors of the atoms are based on Delta_Q_^2 per atom.
        Args:
            a_g = coefficient that multiplies the displacement magnitudes
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

        ax.quiver(x, y, z,u*a_g,v*a_g,w*a_g, color='k',linewidths=1)
        sc = ax.scatter(x, y, z, c=M, marker='o', s=60, cmap="jet")

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
    def plot_four_BandStructures(self,nscf_files):
        """"
        plot the 4 band structures
        nscf_files is the list of Ag, Agstar, Aestar, Ae nscf file paths.
        """
        ebands = []
        for file in nscf_files:
            with abiopen(file) as f:
                ebands.append(f.ebands)

        fig, axs = plt.subplots(1, 4) #figsize=(9, 5))
        titles = [r'$A_g$', r'$A_g^*$', r'$A_e^*$', r'$A_e$']
        e0 = ebands[0].fermie
        for i,eband in enumerate(ebands):
            eband.plot(ax=axs[i], e0=e0, ylims=(-5, 5))
            eband.decorate_ax(ax=axs[i], title=titles[i])

        for ax in axs.flat:
            ax.label_outer()

        plt.subplots_adjust(hspace=0.5, wspace=0.1)
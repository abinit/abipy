
from phonopy import Phonopy
import numpy as np

from pymatgen.io.phonopy import get_pmg_structure,get_phonopy_structure

from abipy.core.abinit_units import eV_to_THz
from abipy.dfpt.converters import phonopy_to_abinit
import os,shutil
from abipy.embedding.utils_ifc import accoustic_sum,map_two_structures_coords,clean_structure
from abipy.core.structure import Structure


class Embedded_phonons(Phonopy):
    """
    Defines a Phonopy object implementing the Interatomic Force Constants embedding method defined in:
    https://pubs.acs.org/doi/10.1021/acs.chemmater.3c00537    
    https://iopscience.iop.org/article/10.1088/1367-2630/16/7/073026/meta
    https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.084603

    """
    
    # Inherits from Phonopy class
    def __init__(self,stru_pristine,stru_defect,stru_emb,ifc_emb,nac_params) :
        """
        :param stru_pristine: Supercell pristine structure
        :param stru_defect: Supercell defect structure
        :param stru_emb: Supercell embedded structure
        :param ifc_emb: Interatomic force constant associated to the supercell embedded structure
        :param nac_params : Non-analytical parameters associated to the supercell embedded structure, in phonopy format. 
        """
        super().__init__(unitcell=get_phonopy_structure(stru_emb),
                       supercell_matrix=[1,1,1],
                       primitive_matrix=np.eye(3),)
        self.force_constants=ifc_emb
        self.nac_params=nac_params
        self.stru_pristine=stru_pristine
        self.stru_defect=stru_defect
        self.stru_emb=stru_emb



    @classmethod
    def from_phonopy_instances(cls,
                               phonopy_pristine,
                               phonopy_defect,
                               structure_defect_wo_relax,# Should corresponds to order of phonopy_defect structure!
                               main_defect_coords_in_pristine,
                               main_defect_coords_in_defect, 
                               substitutions_list:list=None,  #index in pristine,specie ex: [0,"Eu"]
                               vacancies_list:list=None,       # index in pristine ex: [13,14]
                               interstitial_list:list=None,  # species, cart_coord ex: [['Eu',[0,0,0]],['Ce','[0,0,3]']]
                               tol_mapping=0.01,
                               cut_off_mode='auto',rc_1=None,rc_2=None,
                               factor_ifc=1,
                               verbose=0,
                               asr=True) -> Phonopy :
        """
        :param phonopy_pristine: Phonopy object of the pristine structure
        :param phonopy_defect  : Phonopy object of the defect structure
        :param structure_defect_wo_relax: Supercell structure associated to the defect structure, but without relaxation. Needed for an easier mapping
        :param main_defect_coords_in_pristine: Coordinates of th
        :param main_defect_coords_in_defect : Non-analytical parameters associated to the supercell embedded structure, in phonopy format. 
        :param substitutions_list:
        :param vacancies_list
        :param interstitial_list
        :param tol_mapping
        :param cut_off_mode
        :param rc_1
        :param rc_2
        :param factor_ifc
        :param asr
        """
        
        
        ########################
        # Structures manipulations and mapping
        ########################

        stru_defect=get_pmg_structure(phonopy_defect.supercell)
        stru_pristine=get_pmg_structure(phonopy_pristine.supercell)

        stru_emb=stru_pristine.copy()


        if vacancies_list is not None:
            stru_emb.remove_sites(indices=vacancies_list)

        if substitutions_list is not None:
            for sub in substitutions_list: 
                stru_emb.replace(sub[0],sub[1])

        if interstitial_list is not None:
            for inter in interstitial_list: 
                stru_emb.append(species=inter[0],coords=inter[1],coords_are_cartesian=True)

        stru_emb=clean_structure(stru_emb,defect_coord=main_defect_coords_in_pristine)
        structure_defect_wo_relax=clean_structure(structure_defect_wo_relax,defect_coord=main_defect_coords_in_defect)
        


        mapping=map_two_structures_coords(structure_defect_wo_relax,stru_emb,tol=tol_mapping)

        ########################
        # Init of the IFCs
        ########################

        ifc_pristine=phonopy_pristine.force_constants
        ifc_defect=phonopy_defect.force_constants

        # in case of IFC with vacancy, remove ifcs of the vac site
        if vacancies_list is not None:
            ifc_pristine=np.delete(ifc_pristine,vacancies_list,0)

            ifc_pristine=np.delete(ifc_pristine,vacancies_list,1)

        # interstitial
        if interstitial_list is not None:
            for inter in interstitial_list:
                ifc_pristine=np.append(ifc_pristine,np.zeros(shape=[1,len(ifc_pristine),3,3]),axis=0)

                ifc_pristine=np.append(ifc_pristine,np.zeros(shape=[len(ifc_pristine),1,3,3]),axis=1)

    

        ########################
        # Print infos
        ########################

        print(f"Number of atoms in the pristine supercell      : {len(stru_pristine)}")
        print(f"Number of atoms in the defective supercell     : {len(stru_defect)}")

        print(f"Defect infos")
        if substitutions_list is not None:
            print(f"    Substitutions :")
            for sub in substitutions_list:
                print(f"       {sub[0]}, {stru_pristine[sub[0]].coords}, {stru_pristine[sub[0]].species} replaced by {sub[1]}")

        if vacancies_list is not None:
            print(f"    Vacancies :" )
            for vac in vacancies_list:
                print(f"       {vac}, {stru_pristine[vac].coords}, {stru_pristine[vac].species} removed")

        if interstitial_list is not None:
            print(f"    Interstitials :" )
            for inter in interstitial_list:
                print(f"       {inter[0]}, {inter[1]} added")


        print(f"Mapping after structure manipulation           : {len(mapping)}/{len(stru_defect)}")

        ########################
        # change the IFCs
        ########################

        ifc_emb=np.zeros(shape=np.shape(ifc_pristine))

        if cut_off_mode=='auto':
            rc_1=100000 # very large value to include all the ifcs, no sparse matrix. 
            rc_2=0.99*min(np.array(structure_defect_wo_relax.lattice.abc)/2) # largest sphere inscribed in defect supercell
        
        if cut_off_mode=='manual':
            rc_1=rc_1
            rc_2=rc_2

        # Set to zero if 2 atoms distance is bigger than rc_1

        for i,atom1 in enumerate(stru_emb):
            for j,atom2 in enumerate(stru_emb):
                dist_ij=np.sqrt(sum((atom1.coords-atom2.coords)**2))

                if dist_ij > rc_1:
                    ifc_emb[i][j]=np.zeros(shape=(3,3))

                else:
                    ifc_emb[i][j]=ifc_pristine[i][j]
                    
        # Set to doped phonons if 2 atoms are separated from defect by distance < R_c2
        print(f"\n Set IFC to explicit defect phonons calculations if both atoms are separated from defect by a distance < R_c2 = {round(rc_2,3)}")

        for i,atom1 in enumerate(stru_emb):
            for j,atom2 in enumerate(stru_emb):
                # structure centered around defect!!!, main defect is at [0,0,0]
                dist_1_from_defect=np.sqrt(sum((atom1.coords)**2))
                dist_2_from_defect=np.sqrt(sum((atom2.coords)**2))

                if dist_1_from_defect < rc_2 and dist_2_from_defect < rc_2:

                    if verbose > 0:
                        print(f"\n \n Atomic pair : {i,atom1} - {j,atom2}  \n")
                        print(f"atom1 : Dist. from. defect = {dist_1_from_defect}")
                        print(f"atom2 : Dist. from. defect = {dist_2_from_defect}")

                        print(f"Replacing pristine cell IFC = \n {ifc_pristine[i][j]}" )
                        print(f"by defect cell IFC = \n {ifc_defect[mapping.index(i)][mapping.index(j)]}" )
                        print(f"Diff IFC = \n {ifc_pristine[i][j]-ifc_defect[mapping.index(i)][mapping.index(j)]}" )


                    ifc_emb[i][j]=factor_ifc*ifc_defect[mapping.index(i)][mapping.index(j)]
                    
                # enforce ASR
        if asr==True:
            print(f"\n Enforce ASR")
            for alpha in [0,1,2]:
                for beta in [0,1,2]:
                    for i,atom1 in enumerate(stru_emb):
                        ifc_emb[i][i][alpha,beta] = -(accoustic_sum(ifc_emb,i)[alpha,beta]-ifc_emb[i][i][alpha,beta])
        
        ########################
        # change the nac params
        ########################
        
        if phonopy_pristine.nac_params is not None:

            nac_params_emb=phonopy_pristine.nac_params.copy()

            if vacancies_list is not None:
                nac_params_emb["born"]=np.delete(nac_params_emb["born"],vacancies_list,0)
        else:
            nac_params_emb=None

        ########################

        print(f"\n Embedding procedure done")

        return cls(stru_pristine,stru_defect,stru_emb,ifc_emb,nac_params_emb)
    

    def get_gamma_freq_with_vec_abipy_fmt(self):
    # compute with phonopy the Gamma phonons freq and vectors and returns them in abipy format
        ph_freq_phonopy,ph_vec_phonopy=self.get_frequencies_with_eigenvectors(q=[0,0,0])

        ph_freq = ph_freq_phonopy / (eV_to_THz)  # put it in eV
        ph_vec = ph_vec_phonopy.transpose().reshape(3*len(self.supercell),len(self.supercell),3)  
        # such that ph_vec[iband][iatom] gives vector of iband,iatom
        
        return ph_freq, ph_vec
    
    
    def to_ddb(self,embedded_ddb_path='out_DDB',workdir=None):
        # Call the converter to go from phonopy to a DDB
        ddb_sc = phonopy_to_abinit(unit_cell=get_pmg_structure(self.supercell), supercell_matrix=[1,1,1], qpt_list=[[0,0,0]],
                                   out_ddb_path=embedded_ddb_path, force_constants=self.force_constants,
                                   born=self.nac_params, primitive_matrix=np.eye(3), symprec=1e-5,
                                   tolsym=None,workdir=workdir,nsym=1)

        return ddb_sc





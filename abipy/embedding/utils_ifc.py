from phonopy import Phonopy
import numpy as np

from pymatgen.io.phonopy import get_pmg_structure,get_phonopy_structure

from abipy.core.abinit_units import eV_to_THz
from abipy.dfpt.converters import phonopy_to_abinit
import os,shutil
## Usefull functions to be used in embedding_ifc ###

def stru_0_1_to_minus_05_05(structure):
    """
    translate the structure so that the frac_coords go from [0,1] to [-0.5,0,5]
    """

    new_stru = structure.copy()
    for site in new_stru:
        if site.frac_coords[0] >= 0.5:
            site.frac_coords[0] = site.frac_coords[0] - 1

        if site.frac_coords[1] >= 0.5:
            site.frac_coords[1] = site.frac_coords[1] - 1

        if site.frac_coords[2] >= 0.5:
            site.frac_coords[2] = site.frac_coords[2] - 1
    return new_stru

def frac_coords_05_to_minus05(stru):
    """
    Move atoms that are close to frac_coord = 0.5 to -0.5
    Needed for atoms that are at -0.5 in the perfect supercell but move to for instance 0.498 after relaxation
    """
    new_stru = stru.copy()
    for i, atom in enumerate(new_stru):
        if atom.frac_coords[0] > 0.495:
            atom.frac_coords[0] = -0.5
        if atom.frac_coords[1] > 0.495:
            atom.frac_coords[1] = -0.5
        if atom.frac_coords[2] > 0.495:
            atom.frac_coords[2] = -0.5

    return new_stru

def center_wrt_defect(structure,defect_coord):
    """
    Center the structure around defect_coord. Defect is now at [0,0,0]
    """

    new_stru = structure.copy()
#    site = structure[index_in_structure]
    new_stru.translate_sites(indices=np.arange(0, len(structure)),
                            vector=-defect_coord,frac_coords=False,
                            to_unit_cell=False)
    return new_stru


def clean_structure(structure,defect_coord):
    """
    Apply successively :
    center_wrt_defect(), stru_0_1_to_minus_05_05(), frac_coords_05_to_minus05()
    Usefull to match structures after clean_structure()
    """
    stru=structure.copy()
    stru=center_wrt_defect(stru,defect_coord)
    stru=stru_0_1_to_minus_05_05(stru)
    stru=frac_coords_05_to_minus05(stru)
    return stru


def map_two_structures_coords(stru_1, stru_2,tol=0.1):
    """
    Returns a mapping between two structures based on 
    the coordinates, within a given tolerance in Angstrom.
    stru_1 should be a subset of the superset structure stru_2
    """

    cart_1 = stru_1.cart_coords
    cart_2 = stru_2.cart_coords
    mapping = []
    for i, site_1 in enumerate(stru_1):  # subset structure
        for j, site_2 in enumerate(stru_2):  # superset structure
            if max(abs(cart_1[i] - cart_2[j])) < tol:
                mapping.append(j)
    
    if len(mapping)!= len(stru_1):
        print(f"Mapping incomplete : only {len(mapping)}/{len(stru_1)} atoms mapped.")
    return mapping

def accoustic_sum(Hessian_matrix, atom_label):
    Sum_hessian=np.zeros((3,3))
    for m in range(len(Hessian_matrix[0])):
            Sum_hessian += Hessian_matrix[m][atom_label]
        
    return Sum_hessian


def inverse_participation_ratio(eigenvectors):
    """
    Compute equation (10) of https://pubs.acs.org/doi/10.1021/acs.chemmater.3c00537
    for a given q-point, eigenvectors shape=(nbands,natoms,3)
    """
    # for a given q-point,  eigenvectors shape=(nbands,natoms,3)
    ipr=[]
    for iband in range(len(eigenvectors)):    
        sum_atoms=0
        for iatom in range(len(eigenvectors[iband])):
            sum_atoms += np.dot(eigenvectors[iband,iatom],eigenvectors[iband,iatom])**2
        ipr.append(1/sum_atoms)
    return np.array(ipr).real

def localization_ratio(eigenvectors):
    """
    See equation (10) of https://pubs.acs.org/doi/10.1021/acs.chemmater.3c00537   
    for a given q-point, eigenvectors shape=(nbands,natoms,3)
    """
    ipr=inverse_participation_ratio(eigenvectors)
    return len(eigenvectors[0])/ipr


def vesta_phonon(eigenvectors,in_path,ibands=None,
                 scale_vector=20,width_vector=0.3,color_vector=[255,0,0],centered=True,
                 factor_keep_vectors=0.1,
                 out_path="VESTA_FILES"):
    """
    Draw the phonons eigenvectors on a vesta file. 
    Inspired from https://github.com/AdityaRoy-1996/Phonopy_VESTA/tree/master

    Args:
        eigenvectors: phonons eigenvectors for given q-point with shape (nbands,natom,3)
        in_path : path where the initial vesta structure in stored, the structure should be consistent 
        with the eigenvectors provided.
        ibands: list of indices of the bands to include, if not provided, all the bands
        scale_vector : scaling factor of the vector modulus
        width_vector : vector width
        color_vector : color in rgb format
        centered : center the vesta structure around [0,0,0]
        factor_keep_vectors : draw only the eigenvectors with magnitude > factor_keep_vectors * max(magnitude)
        out_path : path where .vesta files with vector are stored
    """

    vesta = open(in_path,'r').read()
    nbands = len(eigenvectors)
    natoms = len(eigenvectors[0])

    if ibands is None :
        ibands = range(nbands)

    if os.path.isdir(out_path): 
        shutil.rmtree(out_path)
        os.mkdir(out_path)
    else : 
        os.mkdir(out_path)

    path=out_path

    for iband in ibands : 
        towrite = vesta.split('VECTR')[0]
        towrite += 'VECTR\n'

        magnitudes=[]
        for iatom in range(natoms) :
            magnitudes.append(np.sqrt(eigenvectors[iband][iatom][0]**2+eigenvectors[iband][iatom][1]**2+eigenvectors[iband][iatom][2]**2))
        for iatom in range(natoms) :
            if magnitudes[iatom] > factor_keep_vectors * max(np.real(magnitudes)):
                towrite += '%5d' %(iatom + 1)
                towrite += '%10.5f' %(eigenvectors[iband][iatom][0] * int(scale_vector))
                towrite += '%10.5f' %(eigenvectors[iband][iatom][1] * int(scale_vector))
                towrite += '%10.5f' %(eigenvectors[iband][iatom][2] * int(scale_vector))
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


        filename = path + f'/{iband:05}_'
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

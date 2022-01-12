from phonopy import Phonopy
from phonopy.harmonic import force_constants

import numpy as np

from pymatgen.io.phonopy import get_pmg_structure,get_phonopy_structure
from pymatgen.analysis import structure_matcher
import abipy.core.kpoints as kpoints


def DDB_to_phonopy_Gamma(DDB):
    # DDB_file of a unit cell with a q-grid corresponding to the supercell size
    # (first phonon supercell size= displacement supercell)

    # return a phonopy object with supercell size = unitcell size (Gamma mapping)

    ngqpt = DDB.guessed_ngqpt

    phonopy_unit = DDB.anaget_phonopy_ifc(dipdip=1, asr=2, set_masses=True)
    #    nac_unit=phonopy_unit.get_nac_params()# extract non anal infos

    # Phonopy structures
    #    stru_phonopy_unit=get_pmg_structure(phonopy_unit.unitcell)
    #    stru_phonopy_supercell=get_pmg_structure(phonopy_unit.supercell)

    correct_stru = DDB.structure.copy()
    correct_stru.make_supercell(ngqpt)

    phonopy_stru = get_pmg_structure(phonopy_unit.get_supercell())

    mapping = structure_matcher.StructureMatcher(primitive_cell=False, ).get_mapping(superset=phonopy_stru,
                                                                                     subset=correct_stru)

    # usefull maps
    #    s2u_map=phonopy_unit.supercell.get_supercell_to_unitcell_map()
    #    u2u_map=phonopy_unit.supercell.get_unitcell_to_unitcell_map()

    # map the born charges to the supercell
    #    born_unit=nac_unit['born']
    #    born_supercell=np.zeros(shape=(len(stru_phonopy_supercell),3,3))

    #    for i,site in enumerate(stru_phonopy_supercell):
    #        index_in_unit=u2u_map[s2u_map[i]] # find the corresponding index in unit cell
    #        born_supercell[i]=born_unit[index_in_unit]  # fill the new Born

    #    nac_supercell=nac_unit
    #    nac_supercell['born']=born_supercell
    # extract full force constants size = (n_atom_supercell,n_atom_supercell,3,3)
    full_fc = force_constants.compact_fc_to_full_fc(phonon=phonopy_unit,
                                                    compact_fc=phonopy_unit.force_constants)

    ordered_full_fc = np.zeros(shape=np.shape(full_fc))

    for i in range(len(phonopy_stru)):
        for j in range(len(phonopy_stru)):
            ordered_full_fc[i, j] = full_fc[mapping[i], mapping[j]]

            # create the phonopy supercell
    phonopy_supercell = Phonopy(unitcell=get_phonopy_structure(correct_stru),
                                # the new unit cell is the 'old' supercell
                                supercell_matrix=[1, 1, 1],  # sup_size= unit_size
                                primitive_matrix=np.eye(3), )

    phonopy_supercell.set_force_constants(ordered_full_fc)
    #    phonopy_supercell.set_nac_params(nac_supercell)

    return phonopy_supercell


def get_phonopy_phonons_supercell(phonopy_supercell_instance):
    # return the phonons freq and vectors from the phonopy supercell instance
    ph_modes = phonopy_supercell_instance.get_frequencies_with_eigenvectors(q=[0, 0, 0])
    ph_freq_phonopy, ph_vec_phonopy = ph_modes

    ph_freq = ph_freq_phonopy * 0.00413566553853599  # put it in eV
    ph_vec = ph_vec_phonopy.transpose()  # such that ph_vec[i] gives the i-th mode (as in abipy)

    return ph_freq, ph_vec

def get_interpolated_phonopy_instance(supercell_size,prim_DDB):
    inter_DDB=prim_DDB.anaget_interpolated_ddb(kpoints.kmesh_from_mpdivs(mpdivs=supercell_size,shifts=[0,0,0]))
    phono_instance=DDB_to_phonopy_Gamma(inter_DDB)
    return phono_instance

def get_supercell_structure_from_phonopy_instance(phonopy_instance):
    phonon_supercell=get_pmg_structure(phonopy_instance.get_supercell())
    return phonon_supercell

def get_phonons_mapped_and_spcell_structure(supercell_size,prim_DDB):
    phonopy_instance=get_interpolated_phonopy_instance(supercell_size,prim_DDB)
    ph_freq,ph_vec=get_phonopy_phonons_supercell(phonopy_instance)
    structure=get_supercell_structure_from_phonopy_instance(phonopy_instance)
    return ph_freq,ph_vec,structure

def get_matching_dSCF_phonon_spcell(dSCF_specell,supercell_size_dSCF,phonon_spcell,supercell_size_phonon):

    prim_matrix=phonon_spcell.lattice.matrix/supercell_size_phonon
    inv_prim_matrix=np.linalg.inv(prim_matrix)

    super_fracs_phonon=np.dot(phonon_spcell.cart_coords,inv_prim_matrix)
    super_fracs_dSCF=np.dot(dSCF_specell.cart_coords,inv_prim_matrix)

    super_fracs_phonon_center = super_fracs_phonon.copy()
    super_fracs_dSCF_center = super_fracs_dSCF.copy()

    # centering both structures

    for i in range(len(super_fracs_dSCF_center)):
        if super_fracs_dSCF[i][0] > supercell_size_dSCF[0] / 2:
            super_fracs_dSCF_center[i][0] = super_fracs_dSCF[i][0] - supercell_size_dSCF[0]

        if super_fracs_dSCF[i][1] > supercell_size_dSCF[1] / 2:
            super_fracs_dSCF_center[i][1] = super_fracs_dSCF[i][1] - supercell_size_dSCF[1]

        if super_fracs_dSCF[i][2] > supercell_size_dSCF[2] / 2:
            super_fracs_dSCF_center[i][2] = super_fracs_dSCF[i][2] - supercell_size_dSCF[2]

    for i in range(len(super_fracs_phonon_center)):
        if super_fracs_phonon[i][0] > supercell_size_phonon[0] / 2:
            super_fracs_phonon_center[i][0] = super_fracs_phonon[i][0] - supercell_size_phonon[0]

        if super_fracs_phonon[i][1] > supercell_size_phonon[1] / 2:
            super_fracs_phonon_center[i][1] = super_fracs_phonon[i][1] - supercell_size_phonon[1]

        if super_fracs_phonon[i][2] > supercell_size_phonon[2] / 2:
            super_fracs_phonon_center[i][2] = super_fracs_phonon[i][2] - supercell_size_phonon[2]

    # perform the matching
    mapping = []
    for i, site_1 in enumerate(dSCF_specell):  # subset structure
        for j, site_2 in enumerate(phonon_spcell):  # superset structure
            if max(abs(super_fracs_dSCF_center[i] - super_fracs_phonon_center[j])) < 0.01:
                mapping.append(j)
    print(mapping)
    print(len(mapping))
    return mapping


def get_forces_on_phonon_supercell(dSCF_supercell,dSCF_size,phonon_supercell,phonon_size,forces_dSCF):
    forces_in_supercell = np.zeros(shape=(len(phonon_supercell), 3))
    mapping=get_matching_dSCF_phonon_spcell(dSCF_supercell,dSCF_size,phonon_supercell,phonon_size)
    for i in range(len(mapping)):
        forces_in_supercell[mapping[i]]=forces_dSCF[i]

    return forces_in_supercell







from phonopy import Phonopy
from phonopy.harmonic import force_constants

import numpy as np

from pymatgen.io.phonopy import get_pmg_structure
from pymatgen.analysis import structure_matcher

def DDB_to_phonopy_Gamma(DDB):
    # DDB_file of a unit cell with a q-grid corresponding to the supercell size
    # (first phonon supercell size= displacement supercell)

    # return a phonopy object with supercell size = unitcell size (Gamma mapping)

    ngqpt = DDB.guessed_ngqpt

    phonopy_unit = DDB.anaget_phonopy_ifc(dipdip=1)
    nac_unit = phonopy_unit.get_nac_params()  # extract non anal infos

    # Phonopy structures
    stru_phonopy_unit = get_pmg_structure(phonopy_unit.unitcell)
    stru_phonopy_supercell = get_pmg_structure(phonopy_unit.supercell)

    # usefull maps
    s2u_map = phonopy_unit.supercell.get_supercell_to_unitcell_map()
    u2u_map = phonopy_unit.supercell.get_unitcell_to_unitcell_map()

    # map the born charges to the supercell
    born_unit = nac_unit['born']
    born_supercell = np.zeros(shape=(len(stru_phonopy_supercell), 3, 3))

    for i, site in enumerate(stru_phonopy_supercell):
        index_in_unit = u2u_map[s2u_map[i]]  # find the corresponding index in unit cell
        born_supercell[i] = born_unit[index_in_unit]  # fill the new Born

    nac_supercell = nac_unit
    nac_supercell['born'] = born_supercell

    # extract full force constants size = (n_atom_supercell,n_atom_supercell,3,3)
    full_fc = force_constants.compact_fc_to_full_fc(phonon=phonopy_unit,
                                                    compact_fc=phonopy_unit.force_constants)

    # create the phonopy supercell
    phonopy_supercell = Phonopy(unitcell=phonopy_unit.supercell,  # the new unit cell is the 'old' supercell
                                supercell_matrix=[1, 1, 1],  # sup_size= unit_size
                                primitive_matrix=np.eye(3), )

    phonopy_supercell.set_force_constants(full_fc)
    phonopy_supercell.set_nac_params(nac_supercell)

    return phonopy_supercell


def get_phonopy_phonons_supercell(phonopy_supercell_instance):
    # return the phonons freq and vectors from the phonopy supercell instance
    ph_modes = phonopy_supercell_instance.get_frequencies_with_eigenvectors(q=[0, 0, 0])
    ph_freq_phonopy, ph_vec_phonopy = ph_modes

    ph_freq = ph_freq_phonopy * 0.00413566553853599  # put it in eV
    ph_vec = ph_vec_phonopy.transpose()  # such that ph_vec[i] gives the i-th mode (as in abipy)

    return ph_freq, ph_vec


def get_phonons_abipy_order(phonopy_supercell_instance):
    # the atoms ordering is different between phonopy supercell and abipy supercell
    # put the phonons eigenvectors in the abipy ordering

    structure_phonopy = get_pmg_structure(phonopy_supercell_instance.unitcell)  # phonopy structure
    structure_abipy = structure_phonopy.get_sorted_structure()  # to be checked twice
    # reasoning is that structure.get_sorted_structure() gives the same order than make_supercell() of abipy
    mapping = structure_matcher.StructureMatcher(primitive_cell=False).get_mapping(superset=structure_phonopy,
                                                                                   subset=structure_abipy)

    ph_freq, ph_vec_phonopy = get_phonopy_phonons_supercell(phonopy_supercell_instance)

    n_atoms = len(structure_phonopy)
    n_modes = 3 * n_atoms

    reshaped_vec = ph_vec_phonopy.reshape((n_modes, n_atoms, 3))  # reshaped for easier manipulation

    final_vectors = np.zeros(shape=(np.shape(reshaped_vec)))

    for j in range(n_modes):
        for i in range(n_atoms):
            final_vectors[j, i] = reshaped_vec[j, mapping[i]]

    final_vectors = final_vectors.reshape(n_modes, 3 * n_atoms)  # reshape to the initial format

    return ph_freq, final_vectors

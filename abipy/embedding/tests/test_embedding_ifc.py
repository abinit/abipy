"""Tests for embedding.embedding_ifc module"""
import abipy.data as abidata
from abipy.core.testing import AbipyTest
from abipy.abilab import abiopen

from abipy.dfpt.ddb import DdbFile
import tempfile
import os

from abipy.dfpt.converters import ddb_ucell_to_phonopy_supercell
from abipy.core.kpoints import kmesh_from_mpdivs

from abipy.embedding.utils_ifc import localization_ratio
from abipy.embedding.embedding_ifc import Embedded_phonons
from pymatgen.io.phonopy import get_pmg_structure,get_phonopy_structure


try:
    import phonopy 
except ImportError:
    Phonopy = None


class Embedding_ifcTest(AbipyTest):
    
    def test_embedding_vacancy(self):

        self.skip_if_not_phonopy()
        tmp_dir = tempfile.mkdtemp()

        # open a 4x4x4 q-mesh ddb of CaO
        ddb_pristine = DdbFile(abidata.ref_file("refs/embedding_ifc/CaO_444_DDB"))

        ## interpolate on a 3x3x3 q-mesh

        qpts=kmesh_from_mpdivs(mpdivs=[3,3,3],shifts=[0,0,0],order="unit_cell")
        ddb_pristine_333=ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts)



        ph_pristine=ddb_ucell_to_phonopy_supercell(unit_ddb=ddb_pristine_333,nac=True)


        ddb_defect = DdbFile(abidata.ref_file("refs/embedding_ifc/CaO_16at_vacancy_DDB"))
        ph_defect=ddb_defect.anaget_phonopy_ifc()

        ########

        # index of the vac. = 8 (in defect structure), this is found manually
        idefect_defect_stru=8

        # We need first to create the defect structure without relax
        structure_defect_wo_relax=ddb_pristine.structure.copy()
        structure_defect_wo_relax.make_supercell(2)

        main_defect_coords_in_defect=structure_defect_wo_relax.cart_coords[idefect_defect_stru]

        structure_defect_wo_relax.remove_sites(indices=[idefect_defect_stru])
        structure_defect_wo_relax.sort()

        # 
        idefect_pristine_stru=27
        main_defect_coords_in_pristine=get_pmg_structure(ph_pristine.supercell).cart_coords[idefect_pristine_stru]


        ph_emb=Embedded_phonons.from_phonopy_instances(
                               phonopy_pristine=ph_pristine,
                               phonopy_defect=ph_defect,
                               structure_defect_wo_relax=structure_defect_wo_relax,
                               main_defect_coords_in_pristine=main_defect_coords_in_pristine,
                               main_defect_coords_in_defect=main_defect_coords_in_defect, 
                               substitutions_list=None, 
                               vacancies_list=[27],       
                               interstitial_list=None,  
                               factor_ifc=3, #  increasing ifc, this will induce a fictious local mode above the bulk frequencies
                               cut_off_mode='auto',
                               verbose=0,
                               asr=True,)
        
        #test the conversion to ddb
        ph_emb.to_ddb(tmp_dir+"out_DDB")
        ddb_emb = DdbFile(tmp_dir+"out_DDB")

        modes=ddb_emb.anaget_phmodes_at_qpoint([0,0,0])
        freqs_anaddb=modes.phfreqs[0]
        freqs_phonopy,vecs=ph_emb.get_gamma_freq_with_vec_abipy_fmt()

        self.assert_almost_equal(freqs_anaddb[-1],7.31410150e-02,decimal=5)
        self.assert_almost_equal(freqs_phonopy[-1],7.31411352e-02,decimal=5)
        self.assert_almost_equal(max(abs(freqs_anaddb-freqs_phonopy)),0,decimal=5)

        ratio=localization_ratio(vecs)

        self.assert_almost_equal(ratio[-1],7.41692860,decimal=5)



    def test_embedding_substitution(self):

        self.skip_if_not_phonopy()
        tmp_dir = tempfile.mkdtemp()

        ddb_pristine = DdbFile(abidata.ref_file("refs/embedding_ifc/SrCl2_DDB"))

        ## interpolate on a 6x6x6 q-mesh

        qpts=kmesh_from_mpdivs(mpdivs=[6,6,6],shifts=[0,0,0],order="unit_cell")
        ddb_pristine_666=ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts)



        ph_pristine=ddb_ucell_to_phonopy_supercell(unit_ddb=ddb_pristine_666,nac=False)

        # test with phonopy loading 
        ph_defect = phonopy.load(supercell_filename=abidata.ref_file("refs/embedding_ifc/SrCl2_Eu_POSCAR"),
                         force_sets_filename=abidata.ref_file("refs/embedding_ifc/SrCl2_Eu_FORCE_SETS"))


        ########
        # We need first to create the defect structure without relax
        structure_defect_wo_relax=ddb_pristine.structure.copy()
        structure_defect_wo_relax.make_supercell(3)
        structure_defect_wo_relax.replace(0,'Eu')
        structure_defect_wo_relax.sort()

        # index of the sub. = 26 (in defect structure), this is found manually
        idefect_defect_stru=26
        main_defect_coords_in_defect=structure_defect_wo_relax.cart_coords[idefect_defect_stru]

        # index of the sub. = 0 (in pristine structure), this is found manually
        idefect_pristine_stru=0
        main_defect_coords_in_pristine=get_pmg_structure(ph_pristine.supercell).cart_coords[idefect_pristine_stru]


        ph_emb=Embedded_phonons.from_phonopy_instances(
                               phonopy_pristine=ph_pristine,
                               phonopy_defect=ph_defect,
                               structure_defect_wo_relax=structure_defect_wo_relax,
                               main_defect_coords_in_pristine=main_defect_coords_in_pristine,
                               main_defect_coords_in_defect=main_defect_coords_in_defect, 
                               substitutions_list=[[idefect_pristine_stru,"Eu"]], 
                               vacancies_list=None,       
                               interstitial_list=None,  
                               cut_off_mode='auto',
                               verbose=0,
                               asr=True,)


        freqs,vecs=ph_emb.get_gamma_freq_with_vec_abipy_fmt()

        self.assert_almost_equal(freqs[-1],0.026440349386134484,decimal=5)

        ratio=localization_ratio(vecs)

        self.assert_almost_equal(ratio[446],21.65624596304,decimal=5)


    def test_embedding_interstitial(self):

        self.skip_if_not_phonopy()
        tmp_dir = tempfile.mkdtemp()

        ddb_pristine = DdbFile(abidata.ref_file("refs/embedding_ifc/C_conv_DDB"))


        qpts=kmesh_from_mpdivs(mpdivs=[3,3,3],shifts=[0,0,0],order="unit_cell")
        ddb_unit_333=ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts,)

        ph_pristine=ddb_ucell_to_phonopy_supercell(unit_ddb=ddb_unit_333,
                                                   nac=False,)


        # test with phonopy loading 
        ddb_defect = DdbFile(abidata.ref_file("refs/embedding_ifc/C_interstitial_DDB"))
        ph_defect=ddb_defect.anaget_phonopy_ifc()

        ########

        ### model of the split-interstitial as one vacancy + 2 inter.

        stru=abiopen(abidata.ref_file("refs/embedding_ifc/C_sc.cif"))
        stru.remove_sites(indices=[14])
        stru.append(species="C",coords=[0.600, 0.5000, 0.2500])
        stru.append(species="C",coords=[0.400, 0.5000, 0.2500])
        structure_defect_wo_relax=stru
        
        # main defect is 
        main_defect_coords_in_pristine=get_pmg_structure(ph_pristine.supercell)[31].coords
        main_defect_coords_in_defect=main_defect_coords_in_pristine

        ph_emb=Embedded_phonons.from_phonopy_instances(phonopy_pristine=ph_pristine,
                                               phonopy_defect=ph_defect,
                                               structure_defect_wo_relax=structure_defect_wo_relax,
                                               main_defect_coords_in_pristine=main_defect_coords_in_pristine,
                                               main_defect_coords_in_defect=main_defect_coords_in_defect, 
                                               vacancies_list=[31],
                                               interstitial_list=[['C',[4.2885, 3.5737, 1.7869]],
                                                                  ['C',[2.8590, 3.5737, 1.7869]]
                                                                 ],
                                               verbose=1,
                                               cut_off_mode='auto',
                                               factor_ifc=1,
                                               asr=True)



        freqs,vecs=ph_emb.get_gamma_freq_with_vec_abipy_fmt()

        self.assert_almost_equal(freqs[-1],0.212826358519,decimal=5)

        ratio=localization_ratio(vecs)

        self.assert_almost_equal(ratio[650],78.76445591,decimal=5)



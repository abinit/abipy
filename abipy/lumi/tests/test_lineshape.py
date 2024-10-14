"""Tests for lumi.lineshape module"""
import abipy.data as abidata
import phonopy

from pymatgen.io.phonopy import get_pmg_structure
from abipy.core.testing import AbipyTest
from abipy.lumi.deltaSCF import DeltaSCF
from abipy.embedding.embedding_ifc import Embedded_phonons
from abipy.dfpt.ddb import DdbFile
from abipy.lumi.lineshape import Lineshape
from abipy.dfpt.converters import ddb_ucell_to_phonopy_supercell
from abipy.core.kpoints import kmesh_from_mpdivs

class DeltaSCFTest(AbipyTest):

    def test_deltaSCF(self):
        """Testing DeltaSCF"""

        Delta_333=DeltaSCF.from_four_points_file([abidata.ref_file("A_g_out_GSR.nc"),
                                                  abidata.ref_file("A_g_starout_GSR.nc"),
                                                  abidata.ref_file("A_e_starout_GSR.nc"),
                                                  abidata.ref_file("A_e_out_GSR.nc"),])

        # open initial DDB file
        ddb_pristine = DdbFile(abidata.ref_file("refs/embedding_ifc/SrCl2_DDB"))

        ## interpolate on a different q-meshes
        qgrids=[
                [3,3,3], # test same size dSCF scell than phonon scell
                [4,4,4], # test same size dSCF scell smaller phonon scell
            ]
        qpts_s=[kmesh_from_mpdivs(mpdivs=qgrid, shifts=[0,0,0], order="unit_cell") for qgrid in qgrids]

        ddb_pristine_s=[ddb_pristine.anaget_interpolated_ddb(qpt_list=qpts) for qpts in qpts_s]

        ## folding procedure
        ph_pristine=[ddb_ucell_to_phonopy_supercell(unit_ddb=ddb_pristine,nac=False) for ddb_pristine in ddb_pristine_s]

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
        main_defect_coords_in_pristine=get_pmg_structure(ph_pristine[0].supercell).cart_coords[idefect_pristine_stru]


        ph_emb_s=[Embedded_phonons.from_phonopy_instances(
                            phonopy_pristine=ph,
                            phonopy_defect=ph_defect,
                            structure_defect_wo_relax=structure_defect_wo_relax,
                            main_defect_coords_in_pristine=main_defect_coords_in_pristine,
                            main_defect_coords_in_defect=main_defect_coords_in_defect,
                            substitutions_list=[[idefect_pristine_stru,"Eu"]],
                            vacancies_list=None,
                            interstitial_list=None,
                            cut_off_mode='auto',
                            rc_2=0,
                            rc_1=100000,
                            verbose=0,
                            asr=True,) for ph in ph_pristine ]

        idefect_dSCF=0
        coords_defect_dSCF=Delta_333.structuregs.cart_coords[idefect_dSCF]

        lineshapes=[Lineshape.from_phonopy_phonons(E_zpl=Delta_333.E_zpl(),
                                                phonopy_ph=ph,
                                                dSCF_structure=Delta_333.structure_gs(),
                                                use_forces=True,
                                                dSCF_displacements=Delta_333.diff_pos(),
                                                dSCF_forces=Delta_333.forces_gs,
                                                coords_defect_dSCF=coords_defect_dSCF) for ph in ph_emb_s]


        self.assert_almost_equal(lineshapes[0].S_tot(),4.697242297808644, decimal=5)
        self.assert_almost_equal(lineshapes[1].S_tot(),5.030138903884942, decimal=5)

        if self.has_matplotlib():
            assert lineshapes[0].plot_emission_spectrum(show=False)


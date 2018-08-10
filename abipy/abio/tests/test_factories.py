from __future__ import unicode_literals, division, print_function

import abipy.data as abidata
import abipy.abilab as abilab

from abipy.flowtk import Flow, RelaxWork, G0W0Work
from abipy.core.testing import AbipyTest
from abipy.abio.inputs import AbinitInput
from abipy.abio.factories import *
from abipy.abio.factories import BandsFromGsFactory, IoncellRelaxFromGsFactory, HybridOneShotFromGsFactory
from abipy.abio.factories import ScfForPhononsFactory, PhononsFromGsFactory
from abipy.abio.factories import PiezoElasticFactory, PiezoElasticFromGsFactory
import json

write_inputs_to_json = False


class ShiftModeTest(AbipyTest):

    def test_shiftmode(self):
        """Testing shiftmode"""
        from abipy.abio.factories import ShiftMode
        gamma = ShiftMode.GammaCentered
        assert ShiftMode.from_object("G") == gamma
        assert ShiftMode.from_object(gamma) == gamma
        with self.assertRaises(TypeError):
            ShiftMode.from_object({})


class FactoryTest(AbipyTest):

    def setUp(self):
        # Si ebands
        self.si_structure = abidata.structure_from_cif("si.cif")
        self.si_pseudo = abidata.pseudos("14si.pspnc")

        self.gan_structure = abidata.structure_from_cif("gan.cif")
        self.n_pseudo = abidata.pseudos("7n.pspnc")
        self.ga_pseudo = abidata.pseudos("31ga.pspnc")

    def test_gs_input(self):
        """Testing gs_input factory."""
        inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        self.abivalidate_input(inp)

        if False:
            print("Updating json reference files.")
            #with open('gs_input.json', mode='w') as fp:
            #    json.dump(inp.as_dict(), fp, indent=2)

        self.assertIn('scf', inp.runlevel)
        self.assertIn('ground_state', inp.runlevel)
        self.assert_input_equality('gs_input.json', inp)

    def test_ebands_input(self):
        """Testing ebands_input factory."""
        multi = ebands_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)

        scf_inp, nscf_inp = multi.split_datasets()

        if False:
            print("Updating json reference files.")
            #with open('scf_input.json', mode='w') as fp:
            #    json.dump(scf_inp.as_dict(), fp, indent=2)
            #with open('nscf_input.json', mode='w') as fp:
            #    json.dump(nscf_inp.as_dict(), fp, indent=2)

        self.assertIn('bands', nscf_inp.runlevel)
        self.assertIn('nscf', nscf_inp.runlevel)
        self.assertIn('ground_state', nscf_inp.runlevel)
        self.abivalidate_multi(multi)
        self.assert_input_equality('scf_input.json', scf_inp)
        self.assert_input_equality('nscf_input.json', nscf_inp)

        # Test dos_kppa and other options.
        multi_dos = ebands_input(self.si_structure, self.si_pseudo, nscf_nband=10, kppa=10, ecut=2,
                                 spin_mode="unpolarized", smearing=None, charge=2.0, dos_kppa=50)
        assert len(multi_dos) == 3
        assert all(i["charge"] == 2 for i in multi_dos)
        self.assert_equal(multi_dos.get("nsppol"), [1, 1, 1])
        self.assert_equal(multi_dos.get("iscf"), [None, -2, -2])
        self.abivalidate_multi(multi_dos)

        multi_dos = ebands_input(self.si_structure, self.si_pseudo, nscf_nband=10, kppa=10, ecut=2,
                                 spin_mode="unpolarized", smearing=None, charge=2.0, dos_kppa=[50, 100])
        assert len(multi_dos) == 4
        self.assert_equal(multi_dos.get("iscf"), [None, -2, -2, -2])



    def test_ion_ioncell_relax_input(self):
        """Testing ion_ioncell_relax_input factory."""
        multi = ion_ioncell_relax_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)
        # scf_kppa, scf_nband #accuracy="normal", spin_mode="polarized",
        # smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

        ion_inp, ioncell_inp = multi.split_datasets()
        self.assertIn('ion_relax', ion_inp.runlevel)
        self.assertIn('relax', ion_inp.runlevel)
        self.assertIn('ground_state', ion_inp.runlevel)
        flow = Flow.temporary_flow()
        flow.register_work(RelaxWork(ion_inp, ioncell_inp))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_ion_ioncell_relax_and_ebands_input(self):
        """Testing ion_ioncell_relax_ands_ebands_input factory."""
        multi = ion_ioncell_relax_and_ebands_input(self.si_structure, self.si_pseudo, ecut=4)
        # kppa=None, nband=None,
        # pawecutdg=None, accuracy="normal", spin_mode="polarized",
        # smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

        self.abivalidate_multi(multi)
        # ion_inp, ioncell_inp = multi.split_datasets()

    def test_g0w0_with_ppmodel_inputs(self):
        """Testing g0w0_with_ppmodel_input factory."""
        scf_kppa, scf_nband, nscf_nband = 10, 10, 10
        ecuteps, ecutsigx = 2, 2

        multi = g0w0_with_ppmodel_inputs(self.si_structure, self.si_pseudo,
                                         scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                         shifts=(0.5, 0.5, 0.5), ecut=2)

        scf_input, nscf_input, scr_input, sigma_input = multi.split_datasets()

        self.abivalidate_multi(multi)

        self.assertIn('many_body', scr_input.runlevel)
        self.assertIn('screening', scr_input.runlevel)
        self.abivalidate_input(sigma_input)
        self.assertIn('many_body', sigma_input.runlevel)
        self.assertIn('sigma', sigma_input.runlevel)
        self.assertNotIn('hybrid', sigma_input.runlevel)

        if write_inputs_to_json:
            with open('g0w0_with_ppmodel_scf_input.json', mode='w') as fp:
                json.dump(scf_input.as_dict(), fp, indent=2)
            with open('g0w0_with_ppmodel_nscf_input.json', mode='w') as fp:
                json.dump(nscf_input.as_dict(), fp, indent=2)
            with open('g0w0_with_ppmodel_scr_input.json', mode='w') as fp:
                json.dump(scr_input.as_dict(), fp, indent=2)
            with open('g0w0_with_ppmodel_sigma_input.json', mode='w') as fp:
                json.dump(sigma_input.as_dict(), fp, indent=2)

        self.assert_input_equality('g0w0_with_ppmodel_scf_input.json', scf_input)
        self.assert_input_equality('g0w0_with_ppmodel_nscf_input.json', nscf_input)
        self.assert_input_equality('g0w0_with_ppmodel_scr_input.json', scr_input)
        self.assert_input_equality('g0w0_with_ppmodel_sigma_input.json', sigma_input)

        flow = Flow.temporary_flow()
        flow.register_work(G0W0Work(scf_input, nscf_input, scr_input, sigma_input))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0


        # The default value of `shifts` changed in v0.3 from (0.5, 0.5, 0.5) to (0.0, 0.0, 0.0)
        multi = g0w0_with_ppmodel_inputs(self.si_structure, self.si_pseudo,
                                        scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                        ecut=2)
        for inp in multi:
            self.assert_equal(inp["shiftk"].flatten(), (0, 0, 0))

    def test_convergence_inputs_single(self):
        """Testing g0w0_convergence_input factory single calculation."""
        scf_kppa, scf_nband, nscf_nband = 10, 10, [10]
        ecuteps, ecutsigx = [2], 2

        inputs = g0w0_convergence_inputs(self.si_structure, self.si_pseudo, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                         extra_abivars={'ecut_s': [2]}, scf_nband=scf_nband, ecut=2)
        # accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
        # ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
        # sigma_nband=None, gw_qprange=1):

        # one scf, one nscf and one screening / sigma multi
        assert len(inputs) == 4
        self.assertIsInstance(inputs[0][0], AbinitInput)
        self.assertIsInstance(inputs[1][0], AbinitInput)
        self.assertIsInstance(inputs[2][0], AbinitInput)
        self.assertIsInstance(inputs[3][0], AbinitInput)

        if write_inputs_to_json:
            for t in ['00', '10', '20', '30']:
                input_dict = inputs[int(t[0])][int(t[1])].as_dict()
                with open('convergence_inputs_single_factory_' + t + '.json', mode='w') as fp:
                    json.dump(input_dict, fp, indent=2)

        for t in ['00', '10', '20', '30']:
            ref_file = 'convergence_inputs_single_factory_' + t + '.json'
            self.assert_input_equality(ref_file, inputs[int(t[0])][int(t[1])])

        for inp in [item for sublist in inputs for item in sublist]:
            self.abivalidate_input(inp)

        self.assertEqual(inputs[3][0]['gwpara'], 2)
        self.assertEqual(inputs[3][0]['gwmem'], '10')
        self.assertEqual(inputs[2][0]['optdriver'], 3)
        self.assertEqual(inputs[3][0]['optdriver'], 4)
        self.assertNotIn('hybrid', inputs[3][0].runlevel)
        self.assertIn('many_body', inputs[3][0].runlevel)

    def test_convergence_inputs_conv(self):
        """Testing g0w0_convergence_input factory convergence calculation."""
        scf_kppa, scf_nband, nscf_nband = 10, 10, [10, 12, 14]
        ecuteps, ecutsigx = [2, 3, 4], 2

        inputs = g0w0_convergence_inputs(self.si_structure, self.si_pseudo, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                         extra_abivars={'ecut_s': [6, 4, 2]}, scf_nband=scf_nband, ecut=2, nksmall=20)

        inputs_flat = [item for sublist in inputs for item in sublist]

        for i, inp in enumerate(inputs_flat):
            ref_file = "gw_convergence_full_" + str(i) + ".json"
            if write_inputs_to_json:
                with open(ref_file, mode='w') as fp:
                    json.dump(inp.as_dict(), fp, indent=2)
            # self.assert_input_equallity(ref_file, inp)

        for inp in [item for sublist in inputs for item in sublist]:
            self.abivalidate_input(inp)

        # the rest is redundant now..
        self.assertEqual(len(inputs_flat), 24)
        nbands = [inp['nband'] for inp in inputs_flat]
        print(nbands)
        ecuteps = [inp.get('ecuteps', None) for inp in inputs_flat]
        print(ecuteps)
        ecuts = [inp.get('ecut', None) for inp in inputs_flat]
        print(ecuts)

        self.assertEqual(nbands, [10, 10, 10, 14, 14, 14, 10, 12, 14, 10, 12, 14, 10, 12, 14, 10, 12, 14, 10, 12, 14,
                                  10, 12, 14])
        self.assertEqual(ecuteps, [None, None, None, None, None, None, 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4,
                                   4, 4])
        self.assertEqual(ecuts, [6, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

        self.assertEqual(inputs_flat[-1]['ecuteps'], 4)
        self.assertEqual(inputs_flat[-1]['nband'], 14)

    def test_bse_with_mdf(self):
        """Testing bse_with_mdf input factory."""
        scf_kppa, scf_nband, nscf_nband, dos_kppa = 10, 10, 10, 4
        # ecuteps, ecutsigx = 3, 2
        nscf_ngkpt, nscf_shiftk = [2, 2, 2], [[0, 0, 0]]

        multi = bse_with_mdf_inputs(self.si_structure, self.si_pseudo, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                                    ecuteps=2, bs_loband=1, bs_nband=2, mbpt_sciss="0.1 eV", mdf_epsinf=12, ecut=2)
        # exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized",
        # smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
        scf_input, nscf_input, bse_input = multi.split_datasets()
        self.assertIn('many_body', bse_input.runlevel)
        self.assertIn('bse', bse_input.runlevel)
        self.assertNotIn('hybrid', bse_input.runlevel)

        self.abivalidate_multi(multi)

    def test_scf_phonons_inputs(self):
        """Testing scf_phonons_inputs."""
        scf_kppa, scf_nband, nscf_nband, dos_kppa = 10, 10, 10, 4
        ecut = 4
        multi = scf_phonons_inputs(self.si_structure, self.si_pseudo, scf_kppa, ecut=ecut)
        self.abivalidate_multi(multi)
        self.assertIn('dfpt', multi[-1].runlevel)
        self.assertIn('ph_q_pert', multi[-1].runlevel)

    def test_phonons_from_gsinput(self):
        """Testing phonons_from_gsinput"""
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        multi = phonons_from_gsinput(gs_inp, ph_ngqpt=[4, 4, 4], with_ddk=True, with_dde=True,
                                     with_bec=False, ph_tol=None, ddk_tol=None, dde_tol=None)
        self.abivalidate_multi(multi)

        from abipy.abio.input_tags import DDK, DDE, PH_Q_PERT
        inp_ddk = multi.filter_by_tags(DDK)[0]
        inp_dde = multi.filter_by_tags(DDE)[0]
        inp_ph_q_pert_1 = multi.filter_by_tags(PH_Q_PERT)[0]
        inp_ph_q_pert_2 = multi.filter_by_tags(PH_Q_PERT)[10]

        if False:
            print("Updating json reference files.")
            #with open('phonons_from_gsinput_ddk.json', mode='w') as fp:
            #    json.dump(inp_ddk.as_dict(), fp, indent=2)
            #with open('phonons_from_gsinput_dde.json', mode='w') as fp:
            #    json.dump(inp_dde.as_dict(), fp, indent=2)
            #with open('phonons_from_gsinput_ph_q_pert_1.json', mode='w') as fp:
            #    json.dump(inp_ph_q_pert_1.as_dict(), fp, indent=2)
            #with open('phonons_from_gsinput_ph_q_pert_2.json', mode='w') as fp:
            #    json.dump(inp_ph_q_pert_2.as_dict(), fp, indent=2)

        self.assert_input_equality('phonons_from_gsinput_ddk.json', inp_ddk)
        self.assert_input_equality('phonons_from_gsinput_dde.json', inp_dde)
        self.assert_input_equality('phonons_from_gsinput_ph_q_pert_1.json', inp_ph_q_pert_1)
        self.assert_input_equality('phonons_from_gsinput_ph_q_pert_2.json', inp_ph_q_pert_2)

        factory_obj = PhononsFromGsFactory(ph_ngqpt=[4, 4, 4], with_ddk=True, with_dde=True, with_bec=False,
                                           ph_tol=None, ddk_tol=None, dde_tol=None)
        self.assertMSONable(factory_obj)
        multi_obj = factory_obj.build_input(gs_inp)

        inp_ddk_obj = multi_obj.filter_by_tags(DDK)[0]
        inp_dde_obj = multi_obj.filter_by_tags(DDE)[0]
        inp_ph_q_pert_1_obj = multi_obj.filter_by_tags(PH_Q_PERT)[0]
        inp_ph_q_pert_2_obj = multi_obj.filter_by_tags(PH_Q_PERT)[10]

        self.assert_input_equality('phonons_from_gsinput_ddk.json', inp_ddk_obj)
        self.assert_input_equality('phonons_from_gsinput_dde.json', inp_dde_obj)
        self.assert_input_equality('phonons_from_gsinput_ph_q_pert_1.json', inp_ph_q_pert_1_obj)
        self.assert_input_equality('phonons_from_gsinput_ph_q_pert_2.json', inp_ph_q_pert_2_obj)


    def test_elastic_inputs_from_gsinput(self):
        """Testing elastic_inputs_from_gsinput."""
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        multi = piezo_elastic_inputs_from_gsinput(gs_inp, ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False)
        self.abivalidate_multi(multi)
        factory_obj = PiezoElasticFromGsFactory(ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False)
        self.assertMSONable(factory_obj)
        multi_obj = factory_obj.build_input(gs_inp)

    def test_scf_piezo_elastic_inputs(self):
        """Testing scf_piezo_elastic_inputs."""
        kppa = 800
        multi = scf_piezo_elastic_inputs(self.si_structure, self.si_pseudo, kppa, ecut=3, pawecutdg=None,
                                         scf_nband=None, accuracy="normal", spin_mode="polarized",
                                         smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                                         ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False)
        self.abivalidate_multi(multi)

        factory_obj = PiezoElasticFactory(self.si_structure, self.si_pseudo, kppa, ecut=3, pawecutdg=None,
                                          scf_nband=None, accuracy="normal", spin_mode="polarized",
                                          smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                                          ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False)
        self.assertMSONable(factory_obj)
        multi_obj = factory_obj.build_input()

    def test_scf_input(self):
        """Testing scf_input"""
        from abipy.abio.factories import scf_input
        inp = scf_input(self.si_structure, self.si_pseudo, kppa=None, ecut=None, pawecutdg=None, nband=None,
                        accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0,
                        scf_algorithm=None, shift_mode="Monkhorst-Pack")

        with self.assertRaises(AssertionError):
            self.abivalidate_input(inp)
        inp["ecut"] = 2
        self.abivalidate_input(inp)

    def test_ebands_dos_from_gsinput(self):
        """Testing ebands_from_gsinput and dos_from_gsinput"""
        from abipy.abio.factories import ebands_from_gsinput, dos_from_gsinput
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        ebands_inp = ebands_from_gsinput(gs_inp, nband=None, ndivsm=15, accuracy="normal")
        self.abivalidate_input(ebands_inp)

        dos_kppa = 3000
        edos_inp = dos_from_gsinput(gs_inp, dos_kppa, nband=None, accuracy="normal", pdos=False)
        self.abivalidate_input(edos_inp)

        factory_obj = BandsFromGsFactory(nband=None, ndivsm=15, accuracy="normal")
        self.assertMSONable(factory_obj)
        ebands_input_obj = factory_obj.build_input(gs_inp)

    def test_ioncell_relax_from_gsinput(self):
        """Testing ioncell_relax_from_gsinput"""
        from abipy.abio.factories import ioncell_relax_from_gsinput
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=100, ecut=2, spin_mode="polarized")
        icrelax_input = ioncell_relax_from_gsinput(gs_inp)
        self.abivalidate_input(icrelax_input)

        factory_obj = IoncellRelaxFromGsFactory()
        self.assertMSONable(factory_obj)
        icrelax_input_obj = factory_obj.build_input(gs_inp)


    def test_hybrid_oneshot_input(self):
        """Testing hybrid_oneshot_input."""
        from abipy.abio.factories import hybrid_oneshot_input
        ecut = 2
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=100, ecut=ecut, spin_mode="polarized")
        hyb_inp = hybrid_oneshot_input(gs_inp, functional="hse06", ecutsigx=None, gw_qprange=1)
        assert "ecutsigx" in hyb_inp and hyb_inp["ecutsigx"] == ecut * 2
        self.abivalidate_input(hyb_inp)
        self.assertIn('hybrid', hyb_inp.runlevel)
        self.assertNotIn('many_body', hyb_inp.runlevel)

        factory_obj = HybridOneShotFromGsFactory(functional="hse06", ecutsigx=None, gw_qprange=1)
        self.assertMSONable(factory_obj)
        hyb_inp_obj = factory_obj.build_input(gs_inp)

    def test_scf_for_phonons(self):
        """Testing scf_for_phonons."""
        from abipy.abio.factories import scf_for_phonons
        scf_inp = scf_for_phonons(self.si_structure, self.si_pseudo, kppa=1000, ecut=3)
        self.abivalidate_input(scf_inp)

        if False:
            print("Updating json reference files.")
            #with open('scf_for_phonons.json', mode='w') as fp:
            #    json.dump(scf_inp.as_dict(), fp, indent=2)

        self.assert_input_equality('scf_for_phonons.json', scf_inp)

        factory_obj = ScfForPhononsFactory(self.si_structure, self.si_pseudo, kppa=1000, ecut=3)
        self.assertMSONable(factory_obj)
        scf_inp_obj = factory_obj.build_input(scf_inp)

        self.assert_input_equality('scf_for_phonons.json', scf_inp_obj)

    def test_dte_from_gsinput(self):
        """Testing for dte_from_gsinput"""
        #self.skip_if_not_abinit(version='8.3.2')
        from abipy.abio.factories import dte_from_gsinput
        pseudos = [self.ga_pseudo.pseudo_with_symbol('Ga'), self.n_pseudo.pseudo_with_symbol('N')]
        gs_inp = gs_input(self.gan_structure, pseudos, kppa=None, ecut=2,
                          spin_mode="unpolarized", smearing=None)
        # dte calculations only work with selected values of ixc
        gs_inp['ixc'] = 7
        multi = dte_from_gsinput(gs_inp, use_phonons=False, skip_dte_permutations=False)
        self.abivalidate_multi(multi)
        multi = dte_from_gsinput(gs_inp, use_phonons=True, skip_dte_permutations=True)
        self.abivalidate_multi(multi)

        inp1 = multi[0]
        inp2 = multi[5]

        if False:
            print("Updating json reference files.")
            #with open('dte_from_gsinput_1.json', mode='w') as fp:
            #    json.dump(inp1.as_dict(), fp, indent=2)
            #with open('dte_from_gsinput_2.json', mode='w') as fp:
            #    json.dump(inp2.as_dict(), fp, indent=2)

        self.assert_input_equality('dte_from_gsinput_1.json', inp1)
        self.assert_input_equality('dte_from_gsinput_2.json', inp2)

    def test_dfpt_from_gsinput(self):
        """Testing phonons_from_gsinput"""
        pseudos = [self.ga_pseudo.pseudo_with_symbol('Ga'), self.n_pseudo.pseudo_with_symbol('N')]
        gs_inp = gs_input(self.gan_structure, pseudos, kppa=None, ecut=2,
                          spin_mode="unpolarized", smearing=None)
        # dte calculations only work with selected values of ixc
        gs_inp['ixc'] = 7
        multi = dfpt_from_gsinput(gs_inp, ph_ngqpt=[4, 4, 4], do_ddk=True, do_dde=True, do_strain=True,
                                  do_dte=True, ph_tol=None, ddk_tol=None, dde_tol=None)
        self.abivalidate_multi(multi)

        from abipy.abio.input_tags import DDK, DDE, PH_Q_PERT, STRAIN, DTE
        inp_ddk = multi.filter_by_tags(DDK)[0]
        inp_dde = multi.filter_by_tags(DDE)[0]
        inp_ph_q_pert_1 = multi.filter_by_tags(PH_Q_PERT)[0]
        inp_ph_q_pert_2 = multi.filter_by_tags(PH_Q_PERT)[10]
        inp_strain = multi.filter_by_tags(STRAIN)[0]
        inp_dte = multi.filter_by_tags(DTE)[0]

        if False:
            print("Updating json reference files.")
            #with open('dfpt_from_gsinput_ddk.json', mode='w') as fp:
            #    json.dump(inp_ddk.as_dict(), fp, indent=2)
            #with open('dfpt_from_gsinput_dde.json', mode='w') as fp:
            #    json.dump(inp_dde.as_dict(), fp, indent=2)
            #with open('dfpt_from_gsinput_ph_q_pert_1.json', mode='w') as fp:
            #    json.dump(inp_ph_q_pert_1.as_dict(), fp, indent=2)
            #with open('dfpt_from_gsinput_ph_q_pert_2.json', mode='w') as fp:
            #    json.dump(inp_ph_q_pert_2.as_dict(), fp, indent=2)
            #with open('dfpt_from_gsinput_strain.json', mode='w') as fp:
            #    json.dump(inp_strain.as_dict(), fp, indent=2)
            #with open('dfpt_from_gsinput_dte.json', mode='w') as fp:
            #    json.dump(inp_dte.as_dict(), fp, indent=2)

        self.assert_input_equality('dfpt_from_gsinput_ddk.json', inp_ddk)
        self.assert_input_equality('dfpt_from_gsinput_dde.json', inp_dde)
        self.assert_input_equality('dfpt_from_gsinput_ph_q_pert_1.json', inp_ph_q_pert_1)
        self.assert_input_equality('dfpt_from_gsinput_ph_q_pert_2.json', inp_ph_q_pert_2)
        self.assert_input_equality('dfpt_from_gsinput_strain.json', inp_strain)
        self.assert_input_equality('dfpt_from_gsinput_dte.json', inp_dte)


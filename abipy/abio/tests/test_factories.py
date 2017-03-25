from __future__ import unicode_literals, division, print_function

import abipy.data as abidata
import abipy.abilab as abilab

from abipy.flowtk import Flow, BandStructureWork, RelaxWork, G0W0Work, BseMdfWork
from abipy.core.testing import AbipyTest
from abipy.abio.inputs import AbinitInput
from abipy.abio.factories import *
import json


class ShiftModeTest(AbipyTest):
    def test_shiftmode(self):
        """Testing shifmode"""
        from abipy.abio.factories import ShiftMode
        gamma = ShiftMode.GammaCentered
        assert ShiftMode.from_object("G") == gamma
        assert ShiftMode.from_object(gamma) == gamma


class FactoryTest(AbipyTest):

    def setUp(self):
        # Si ebands
        self.si_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        self.si_pseudo = abidata.pseudos("14si.pspnc")

    @staticmethod
    def validate_multi(multi):
        """Test validity of MultiDataset or a list of input objects."""
        if hasattr(multi, "split_datasets"):
            dtlist = multi.split_datasets()
        else:
            dtlist = multi

        rcode = 0
        for dtset in dtlist:
            v = dtset.abivalidate()
            if v.retcode != 0: print("Validation error in %s" % str(v))
            rcode += v.retcode
        assert rcode == 0

    def test_gs_input(self):
        """Testing gs_input factory."""
        inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        inp.abivalidate()

    def test_ebands_input(self):
        """Testing ebands_input factory."""
        multi = ebands_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)

        scf_inp, nscf_inp = multi.split_datasets()
        self.validate_multi(multi)

    def test_ion_ioncell_relax_input(self):
        """Testing ion_ioncell_relax_input factory."""
        multi = ion_ioncell_relax_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)
        # scf_kppa, scf_nband #accuracy="normal", spin_mode="polarized",
        # smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

        ion_inp, ioncell_inp = multi.split_datasets()

        flow = Flow.temporary_flow()
        flow.register_work(RelaxWork(ion_inp, ioncell_inp))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_ion_ioncell_relax_and_ebands_input(self):
        """Testing ion_ioncell_relax_ands_ebands_input factory."""
        multi = ion_ioncell_relax_and_ebands_input(self.si_structure, self.si_pseudo, ecut=4)
        # kppa=None, nband=None,
        # pawecutdg=None, accuracy="normal", spin_mode="polarized",
        # smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

        self.validate_multi(multi)
        # ion_inp, ioncell_inp = multi.split_datasets()

    def test_g0w0_with_ppmodel_inputs(self):
        """Testing g0w0_with_ppmodel_input factory."""
        scf_kppa, scf_nband, nscf_nband = 10, 10, 10
        ecuteps, ecutsigx = 2, 2

        multi = g0w0_with_ppmodel_inputs(self.si_structure, self.si_pseudo,
                                         scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                         ecut=2)

        scf_input, nscf_input, scr_input, sigma_input = multi.split_datasets()

        flow = Flow.temporary_flow()
        flow.register_work(G0W0Work(scf_input, nscf_input, scr_input, sigma_input))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

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
        self.assertEqual(len(inputs), 4)
        self.assertIsInstance(inputs[0][0], AbinitInput)
        self.assertIsInstance(inputs[1][0], AbinitInput)
        self.assertIsInstance(inputs[2][0], AbinitInput)
        self.assertIsInstance(inputs[3][0], AbinitInput)

        if False:
            for t in ['00', '10', '20', '30']:
                input_dict = inputs[int(t[0])][int(t[1])].as_dict()
                with open('convergence_inputs_single_factory_' + t + '.json', mode='w') as fp:
                    json.dump(input_dict, fp, indent=2)

        for t in ['00', '10', '20', '30']:
            ref_file = 'convergence_inputs_single_factory_' + t + '.json'
            self.assert_input_equallity(ref_file, inputs[int(t[0])][int(t[1])])

        for inp in [item for sublist in inputs for item in sublist]:
            val = inp.abivalidate()
            if val.retcode != 0:
                print(inp)
                print(val.log_file.read())
                self.assertEqual(val.retcode, 0)

        self.assertEqual(inputs[3][0]['gwpara'], 2)
        self.assertEqual(inputs[3][0]['gwmem'], '10')
        self.assertEqual(inputs[2][0]['optdriver'], 3)
        self.assertEqual(inputs[3][0]['optdriver'], 4)

    def test_convergence_inputs_conv(self):
        """Testing g0w0_convergence_input factory convergence calculation."""
        scf_kppa, scf_nband, nscf_nband = 10, 10, [10, 12, 14]
        ecuteps, ecutsigx = [2, 3, 4], 2

        inputs = g0w0_convergence_inputs(self.si_structure, self.si_pseudo, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                         extra_abivars={'ecut_s': [6, 4, 2]}, scf_nband=scf_nband, ecut=2, nksmall=20)

        inputs_flat = [item for sublist in inputs for item in sublist]

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

        for inp in [item for sublist in inputs for item in sublist]:
            val = inp.abivalidate()
            if val.retcode != 0:
                print(inp)
                print(val.log_file.read())
                self.assertEqual(val.retcode, 0)

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
        self.validate_multi(multi)

    def test_scf_phonons_inputs(self):
        """Testing scf_phonons_inputs."""
        scf_kppa, scf_nband, nscf_nband, dos_kppa = 10, 10, 10, 4
        ecut = 4
        multi = scf_phonons_inputs(self.si_structure, self.si_pseudo, scf_kppa, ecut=ecut)
        # , pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
        self.validate_multi(multi)

    def test_phonons_from_gsinput(self):
        """Testing phonons_from_gsinput"""
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        multi = phonons_from_gsinput(gs_inp, ph_ngqpt=[4, 4, 4], with_ddk=True, with_dde=True,
                                     with_bec=False, ph_tol=None, ddk_tol=None, dde_tol=None)
        self.validate_multi(multi)

    def test_elastic_inputs_from_gsinput(self):
        """Testing elastic_inputs_from_gsinput."""
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        multi = piezo_elastic_inputs_from_gsinput(gs_inp, ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False)
        self.validate_multi(multi)

    def test_scf_piezo_elastic_inputs(self):
        """Testing scf_piezo_elastic_inputs."""
        kppa = 800
        multi = scf_piezo_elastic_inputs(self.si_structure, self.si_pseudo, kppa, ecut=3, pawecutdg=None,
                                         scf_nband=None, accuracy="normal", spin_mode="polarized",
                                         smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                                         ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False)
        self.validate_multi(multi)

    def test_scf_input(self):
        """Testing scf_input"""
        from abipy.abio.factories import scf_input
        inp = scf_input(self.si_structure, self.si_pseudo, kppa=None, ecut=None, pawecutdg=None, nband=None,
                        accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0,
                        scf_algorithm=None, shift_mode="Monkhorst-Pack")
        inp.abivalidate()

    def test_ebands_dos_from_gsinput(self):
        """Testing ebands_from_gsinput and dos_from_gsinput"""
        from abipy.abio.factories import ebands_from_gsinput, dos_from_gsinput
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=None, ecut=2, spin_mode="unpolarized")
        ebands_inp = ebands_from_gsinput(gs_inp, nband=None, ndivsm=15, accuracy="normal")
        ebands_inp.abivalidate()

        dos_kppa = 3000
        edos_inp = dos_from_gsinput(gs_inp, dos_kppa, nband=None, accuracy="normal", pdos=False)
        edos_inp.abivalidate()

    def test_ioncell_relax_from_gsinput(self):
        """Testing ioncell_relax_from_gsinput"""
        from abipy.abio.factories import ioncell_relax_from_gsinput
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=100, ecut=2, spin_mode="polarized")
        icrelax_input = ioncell_relax_from_gsinput(gs_inp)
        icrelax_input.abivalidate()

    def test_hybrid_oneshot_input(self):
        """Testing hybrid_oneshot_input."""
        from abipy.abio.factories import hybrid_oneshot_input
        ecut = 2
        gs_inp = gs_input(self.si_structure, self.si_pseudo, kppa=100, ecut=ecut, spin_mode="polarized")
        hyb_inp = hybrid_oneshot_input(gs_inp, functional="hse06", ecutsigx=None, gw_qprange=1)
        assert "ecutsigx" in hyb_inp and hyb_inp["ecutsigx"] == ecut * 2
        hyb_inp.abivalidate()

    def test_scf_for_phonons(self):
        """Testing scf_for_phonons."""
        from abipy.abio.factories import scf_for_phonons
        scf_inp = scf_for_phonons(self.si_structure, self.si_pseudo, kppa=1000, ecut=3)
        scf_inp.abivalidate()

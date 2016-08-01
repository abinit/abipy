from __future__ import unicode_literals, division, print_function

import abipy.data as abidata
import abipy.abilab as abilab
from abipy.core.testing import AbipyTest
from abipy.abio.inputs import AbinitInput
from abipy.abio.factories import *


class FactoryTest(AbipyTest):

    def setUp(self):
        # Si ebands
        self.si_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        self.si_pseudo = abidata.pseudos("14si.pspnc")

    def validate_multi(self, multi):
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

        flow = abilab.Flow.temporary_flow()
        flow.register_scf_task(inp)
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_ebands_input(self):
        """Testing ebands_input factory."""
        multi = ebands_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)

        scf_inp, nscf_inp = multi.split_datasets()

        flow = abilab.Flow.temporary_flow()
        flow.register_work(abilab.BandStructureWork(scf_inp, nscf_inp))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_ion_ioncell_relax_input(self):
        """Testing ion_ioncell_relax_input factory."""
        multi = ion_ioncell_relax_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)
                            #scf_kppa, scf_nband #accuracy="normal", spin_mode="polarized",
                            #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

        ion_inp, ioncell_inp = multi.split_datasets()

        flow = abilab.Flow.temporary_flow()
        flow.register_work(abilab.RelaxWork(ion_inp, ioncell_inp))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    #def test_ion_ioncell_relax_and_ebands_input(self):
    #    """Testing ion_ioncell_relax_ands_ebands_input factory."""
    #    multi = ion_ioncell_relax_and_ebands_input(structure, pseudos,
    #                                   kppa=None, nband=None,
    #                                   ecut=None, pawecutdg=None, accuracy="normal", spin_mode="polarized",
    #                                   smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

    #    ion_inp, ioncell_inp = multi.split_datasets()

    #    flow = abilab.Flow.temporary_flow()
    #    flow.register_work(abilab.RelaxWork(ion_inp, ioncell_inp))
    #    assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_g0w0_with_ppmodel_inputs(self):
        """Testing g0w0_with_ppmodel_input factory."""
        scf_kppa, scf_nband, nscf_nband = 10, 10, 10
        ecuteps, ecutsigx = 2, 2

        multi = g0w0_with_ppmodel_inputs(self.si_structure, self.si_pseudo,
                                         scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                         ecut=2)

        scf_input, nscf_input, scr_input, sigma_input = multi.split_datasets()

        flow = abilab.Flow.temporary_flow()
        flow.register_work(abilab.G0W0Work(scf_input, nscf_input, scr_input, sigma_input))
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

        self.assertEqual(inputs[0][0].variable_checksum(), "1f51104b0dac945bd669d7f363692baf2ced4695")
        self.assertEqual(inputs[1][0].variable_checksum(), "2397edaa6748216e14877140ec70f1d3774b5646")
        self.assertEqual(inputs[2][0].variable_checksum(), "b12bb64fb2e7aca84d13d6c0467f79715cf7ed0e")
        self.assertEqual(inputs[3][0].variable_checksum(), "7b2e23a0b622595de7b3a5d5bcb0f464a4152103")

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

    def test_convergence_inputs_conve(self):
        """Testing g0w0_convergence_input factory convergence calculation."""
        scf_kppa, scf_nband, nscf_nband = 10, 10, [10, 12, 14]
        ecuteps, ecutsigx = [2, 3, 4], 2

        inputs = g0w0_convergence_inputs(self.si_structure, self.si_pseudo, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                                         extra_abivars={'ecut_s': [2, 4, 6]}, scf_nband=scf_nband, ecut=2, nksmall=20)

        inputs_flat = [item for sublist in inputs for item in sublist]

        self.assertEqual(len(inputs_flat), 24)

        for inp in [item for sublist in inputs for item in sublist]:
            val = inp.abivalidate()
            if val.retcode != 0:
                print(inp)
                print(val.log_file.read())
                self.assertEqual(val.retcode, 0)

    def test_bse_with_mdf(self):
        """Testing bse_with_mdf input factory."""
        scf_kppa, scf_nband, nscf_nband, dos_kppa = 10, 10, 10, 4
        ecuteps, ecutsigx = 3, 2
        nscf_ngkpt, nscf_shiftk = [2,2,2], [[0,0,0]]

        multi = bse_with_mdf_inputs(self.si_structure, self.si_pseudo, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                                    ecuteps=2, bs_loband=1, bs_nband=2, mbpt_sciss="0.1 eV", mdf_epsinf=12, ecut=2)
                                    #exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized",
                                    #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

        scf_input, nscf_input, bse_input = multi.split_datasets()

        flow = abilab.Flow.temporary_flow()
        flow.register_work(abilab.BseMdfWork(scf_input, nscf_input, bse_input))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_scf_phonons_inputs(self):
        """Testing scf_phonons_inputs."""
        scf_kppa, scf_nband, nscf_nband, dos_kppa = 10, 10, 10, 4
        ecut = 4
        inps = scf_phonons_inputs(self.si_structure, self.si_pseudo, scf_kppa,
                                  ecut=ecut) #, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
        self.validate_multi(inps)

    #def test_phonons_from_gsinput(self):
    #    """Testing phonons_from_gsinput"""
    #    phonons_from_gsinput(gs_inp, ph_ngqpt=None, with_ddk=True, with_dde=True,
    #                        with_bec=False, ph_tol=None, ddk_tol=None, dde_tol=None)

    #def test_elastic_inputs_from_gsinput(self)
        #piezo_elastic_inputs_from_gsinput(gs_inp, ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False):

    #def test_scf_piezo_elastic_inputs(self):
    #    scf_piezo_elastic_inputs(structure, pseudos, kppa, ecut=None, pawecutdg=None, scf_nband=None,
    #                             accuracy="normal", spin_mode="polarized",
    #                             smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
    #                             ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False):

    #def test_scf_input(self):
    #    scf_input(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
    #              spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
    #              shift_mode="Monkhorst-Pack"):


    #def test_ebands_from_gsinput(self):
    #    ebands_from_gsinput(gsinput, nband=None, ndivsm=15, accuracy="normal"):

    #def test_ioncell_relax_from_gsinput(self):
    #    ioncell_relax_from_gsinput(gsinput, accuracy="normal"):

    #def test_hybrid_oneshot_input(self):
    #    def hybrid_oneshot_input(gsinput, functional="hse06", ecutsigx=None, gw_qprange=1):

    #def test_scf_for_phonons(self):
    #    scf_for_phonons(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",

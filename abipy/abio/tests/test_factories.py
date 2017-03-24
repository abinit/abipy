from __future__ import unicode_literals, division, print_function

import abipy.data as abidata
import abipy.abilab as abilab

from abipy.flowapi import Flow, BandStructureWork, RelaxWork, G0W0Work, BseMdfWork
from abipy.core.testing import AbipyTest
from abipy.abio.inputs import AbinitInput
from abipy.abio.factories import *
import numpy
import json


def input_equality_check(ref_file, input2, rtol=1e-05, atol=1e-08, equal_nan=False):
    """
    function to compare two inputs
    ref_file takes the path to reference input in json: json.dump(input.as_dict(), fp, indent=2)
    input2 takes an AbinintInput object
    tol relative tolerance for floats
    we check if all vars are uniquely present in both inputs and if the values are equal (integers, strings)
    or almost equal (floats)
    """

    def check_int(i, j):
        return i != j

    def check_float(x, y):
        return not numpy.isclose(x, y, rtol=rtol, atol=atol, equal_nan=equal_nan)

    def check_str(s, t):
        return s != t

    def check_var(v, w):
        _error = False
        if isinstance(v, int):
            _error = check_int(v, w)
        elif isinstance(v, float):
            _error = check_float(v, w)
        elif isinstance(v, (str, unicode)):
            _error = check_str(v, w)
        return _error

    def flatten_var(o, tree_types=(list, tuple, numpy.ndarray)):
        flat_var = []
        if isinstance(o, tree_types):
            for value in o:
                for sub_value in flatten_var(value, tree_types):
                    flat_var.append(sub_value)
        else:
            flat_var.append(o)
        return flat_var

    with open(ref_file) as fp:
        input_ref = AbinitInput.from_dict(json.load(fp))

    errors = []
    diff_in_ref = [var for var in input_ref.vars if var not in input2.vars]
    diff_in_actual = [var for var in input2.vars if var not in input_ref.vars]
    if len(diff_in_ref) > 0 or len(diff_in_actual) > 0:
        error_description = 'not the same input parameters:\n' \
                            '     %s were found in ref but not in actual\n' \
                            '     %s were found in actual but not in ref\n' % \
                            (diff_in_ref, diff_in_actual)
        errors.append(error_description)

    for var, val_r in input_ref.vars.items():
        try:
            val_t = input2.vars[var]
        except KeyError:
            errors.append('variable %s from the reference is not in the actual input\n' % str(var))
            continue
        val_list_t = flatten_var(val_t)
        val_list_r = flatten_var(val_r)
        error = False
        print(var)
        print(val_list_r, type(val_list_r[0]))
        print(val_list_t, type(val_list_t[0]))
        for k, var_item in enumerate(val_list_r):
            try:
                error = error or check_var(val_list_t[k], val_list_r[k])
            except IndexError:
                print(val_list_t, type(val_list_t[0]))
                print(val_list_r, type(val_list_r[0]))
                raise RuntimeError('two value lists were not flattened in the same way, try to add the collection'
                                   'type to the tree_types tuple in flatten_var')

        if error:
            error_description = 'var %s differs: %s (reference) != %s (actual)' % \
                                (var, val_r, val_t)
            errors.append(error_description)

    if input2.structure != input_ref.structure:
        errors.append('Structures are not the same.\n')
        print(input2.structure, input_ref.structure)

    if len(errors) > 0:
        msg = 'Two inputs were found to be not equal:\n'
        for err in errors:
            msg += '   ' + err + '\n'
        raise AssertionError(msg)


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

        flow = Flow.temporary_flow()
        flow.register_scf_task(inp)
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_ebands_input(self):
        """Testing ebands_input factory."""
        multi = ebands_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)

        scf_inp, nscf_inp = multi.split_datasets()

        flow = Flow.temporary_flow()
        flow.register_work(BandStructureWork(scf_inp, nscf_inp))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    def test_ion_ioncell_relax_input(self):
        """Testing ion_ioncell_relax_input factory."""
        multi = ion_ioncell_relax_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)
                            #scf_kppa, scf_nband #accuracy="normal", spin_mode="polarized",
                            #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

        ion_inp, ioncell_inp = multi.split_datasets()

        flow = Flow.temporary_flow()
        flow.register_work(RelaxWork(ion_inp, ioncell_inp))
        assert flow.build_and_pickle_dump(abivalidate=True) == 0

    #def test_ion_ioncell_relax_and_ebands_input(self):
    #    """Testing ion_ioncell_relax_ands_ebands_input factory."""
    #    multi = ion_ioncell_relax_and_ebands_input(structure, pseudos,
    #                                   kppa=None, nband=None,
    #                                   ecut=None, pawecutdg=None, accuracy="normal", spin_mode="polarized",
    #                                   smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

    #    ion_inp, ioncell_inp = multi.split_datasets()

    #    flow = Flow.temporary_flow()
    #    flow.register_work(RelaxWork(ion_inp, ioncell_inp))
    #    assert flow.build_and_pickle_dump(abivalidate=True) == 0

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
            with open('convergence_inputs_single_factory_00.json', mode='w') as fp:
                json.dump(inputs[0][0].as_dict(), fp, indent=2)
            with open('convergence_inputs_single_factory_10.json', mode='w') as fp:
                json.dump(inputs[1][0].as_dict(), fp, indent=2)
            with open('convergence_inputs_single_factory_20.json', mode='w') as fp:
                json.dump(inputs[2][0].as_dict(), fp, indent=2)
            with open('convergence_inputs_single_factory_30.json', mode='w') as fp:
                json.dump(inputs[3][0].as_dict(), fp, indent=2)

        for t in ['00', '10', '20', '30']:
            ref_file = 'convergence_inputs_single_factory_' + t + '.json'
            input_equality_check(ref_file, inputs[int(t[0])][int(t[1])])

        self.assertTrue(False)

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
        ecuteps, ecutsigx = 3, 2
        nscf_ngkpt, nscf_shiftk = [2,2,2], [[0,0,0]]

        multi = bse_with_mdf_inputs(self.si_structure, self.si_pseudo, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                                    ecuteps=2, bs_loband=1, bs_nband=2, mbpt_sciss="0.1 eV", mdf_epsinf=12, ecut=2)
                                    #exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized",
                                    #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

        scf_input, nscf_input, bse_input = multi.split_datasets()

        flow = Flow.temporary_flow()
        flow.register_work(BseMdfWork(scf_input, nscf_input, bse_input))
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

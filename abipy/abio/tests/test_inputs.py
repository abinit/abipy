"""Tests for input  module"""
from __future__ import print_function, division, unicode_literals

import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import *
from abipy.abio.inputs import *
from abipy.abio.input_tags import *

import abipy.abio.decorators as ideco
from abipy.abio.factories import *


class TestAbinitInput(AbipyTest):
    """Unit tests for AbinitInput."""

    def test_api(self):
        """Test AbinitInput API."""
        # Build simple input with structure and pseudos
        unit_cell = {
            "acell": 3*[10.217],       
            'rprim': [[.0, .5, .5],
                      [.5, .0, .5],
                      [.5, .5, .0]],
            'ntypat': 1,
            'znucl': [14,],
            'natom': 2,
            'typat': [1, 1],
            'xred': [[.0, .0, .0],
                     [.25,.25,.25]]
        }

        inp = AbinitInput(structure=unit_cell, pseudos=abidata.pseudos("14si.pspnc"))

        print(repr(inp))
        assert len(inp) == 0 and not inp 
        assert inp.get("foo", "bar") == "bar" and inp.pop("foo", "bar") == "bar"
        assert inp.comment is None
        inp.set_comment("This is a comment")
        assert inp.comment == "This is a comment"
        assert inp.isnc and not inp.ispaw
        assert not inp.decorators
        assert len(inp.structure) == 2 and inp.num_valence_electrons == 8

        # foo is not a valid Abinit variable
        with self.assertRaises(inp.Error): inp["foo"] = 1
        inp["ecut" ] = 1
        assert inp.get("ecut") == 1 and len(inp) == 1 and "ecut" in inp.keys() and "foo" not in inp

        assert inp.mnemonics == False
        inp.set_mnemonics(True)
        assert inp.mnemonics == True

        # Test to_string
        inp.to_string(sortmode="a", with_structure=True, with_pseudos=True)
        inp.to_string(sortmode="section", with_structure=True, with_pseudos=True)

        inp.set_vars(ecut=5, toldfe=1e-6)

        # Cannot change structure variables directly.
        with self.assertRaises(inp.Error):
            inp.set_vars(unit_cell)

        # Test deepcopy and remove_vars.
        inp["bdgw"] = [1, 2]
        inp_copy = inp.deepcopy()
        inp_copy["bdgw"][1] = 3
        assert inp["bdgw"] == [1, 2]
        assert inp.remove_vars("bdgw") and "bdgw" not in inp

        removed = inp.pop_tolerances()
        assert len(removed) == 1 and removed["toldfe"] == 1e-6

        # Test set_structure 
        new_structure = inp.structure.copy() 
        new_structure.perturb(distance=0.1)
        inp.set_structure(new_structure)
        assert inp.structure == new_structure

        # Compatible with Pickle and MSONable?
        self.serialize_with_pickle(inp, test_eq=False)
        self.assertMSONable(inp)

        # Test tags
        assert isinstance(inp.tags, set)
        assert len(inp.tags) == 0
        inp.add_tags([GROUND_STATE, SCF])
        assert len(inp.tags) == 2
        inp.remove_tags([GROUND_STATE])
        assert len(inp.tags) == 1

    def test_input_errors(self):
        """Testing typical AbinitInput Error"""
        si_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

        # Ambiguous list of pseudos.
        with self.assertRaises(AbinitInput.Error):
            AbinitInput(si_structure, pseudos=abidata.pseudos("14si.pspnc", "Si.oncvpsp"))

        # Pseudos do not match structure.
        with self.assertRaises(AbinitInput.Error):
            AbinitInput(si_structure, pseudos=abidata.pseudos("13al.981214.fhi"))

        # Negative triple product.
        with self.assertRaises(AbinitInput.Error):
            s = abidata.structure_from_ucell("Al-negative-volume")
            AbinitInput(s, pseudos=abidata.pseudos("13al.981214.fhi"))

    def test_helper_functions(self):
        """Testing AbinitInput helper functions."""
        inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))

        inp.set_kmesh(ngkpt=(1, 2, 3), shiftk=(1, 2, 3, 4, 5, 6))
        assert inp["kptopt"] == 1 and inp["nshiftk"] == 2

        inp.set_autokmesh(nksmall=2)
        assert inp["kptopt"] == 1 and np.all(inp["ngkpt"] == [2, 2, 2]) and inp["nshiftk"] == 4

        inp.set_kpath(ndivsm=3, kptbounds=None)
        assert inp["iscf"] == -2 and len(inp["kptbounds"]) == 12

        inp.set_kptgw(kptgw=(1, 2, 3, 4, 5, 6), bdgw=(1, 2))
        assert inp["nkptgw"] == 2 and np.all(inp["bdgw"].ravel() == np.array(len(inp["kptgw"]) * [1,2]).ravel())

        linps = inp.linspace("ecut", 2, 6, num=3, endpoint=True)
        assert len(linps) == 3 and (linps[0]["ecut"] == 2 and (linps[-1]["ecut"] == 6))

        ranginps = inp.arange("ecut", start=3, stop=5, step=1)
        assert len(ranginps) == 2 and (ranginps[0]["ecut"] == 3 and (ranginps[-1]["ecut"] == 4))

        prod_inps = inp.product("ngkpt", "tsmear", [[2,2,2], [4,4,4]], [0.1, 0.2, 0.3])
        #prod_inps = inp.product([("ngkpt", [[2,2,2], [4,4,4]]), ("tsmear", [0.1, 0.2, 0.3])])
        assert len(prod_inps) == 6 
        assert prod_inps[0]["ngkpt"] == [2,2,2] and prod_inps[0]["tsmear"] == 0.1
        assert prod_inps[-1]["ngkpt"] ==  [4,4,4] and prod_inps[-1]["tsmear"] == 0.3

    def test_abinit_calls(self):
        """Testing AbinitInput methods invoking Abinit."""
        inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))
        inp.set_kmesh(ngkpt=(2, 2, 2), shiftk=(0, 0, 0))

        # The code below invokes Abinit.
        if self.has_abinit():

            # Test validate with wrong input
            inp.set_vars(ecut=-1)
            v = inp.abivalidate()
            assert v.retcode != 0 and v.log_file.read()

            # Test validate with correct input
            inp.set_vars(ecut=2, toldfe=1e-6)
            v = inp.abivalidate()
            assert v.retcode == 0

            # Test abiget_ibz
            ibz = inp.abiget_ibz()
            assert np.all(ibz.points == [[ 0. ,  0. ,  0. ], [ 0.5,  0. ,  0. ], [ 0.5,  0.5,  0. ]])
            assert np.all(ibz.weights == [0.125,  0.5,  0.375])

            # Test abiget_irred_phperts
            # [{'idir': 1, 'ipert': 1, 'qpt': [0.0, 0.0, 0.0]}]
            irred_perts = inp.abiget_irred_phperts(qpt=(0, 0, 0))
            assert len(irred_perts) == 1
            pert = irred_perts[0]
            assert pert.idir == 1 and (pert.idir, pert.ipert) == (1, 1) and all(c == 0 for c in pert.qpt)

            # Test abiget_autoparal_pconfs
            inp["paral_kgb"] = 0
            pconfs = inp.abiget_autoparal_pconfs(max_ncpus=5)
            inp["paral_kgb"] = 1
            pconfs = inp.abiget_autoparal_pconfs(max_ncpus=5)

    def test_dict_methods(self):
        """ Testing AbinitInput dict methods """
        inp = ebands_input(abidata.cif_file("si.cif"), abidata.pseudos("14si.pspnc"), kppa=10, ecut=2)[0]
        inp = ideco.SpinDecorator("spinor")(inp)
        inp_dict = inp.as_dict()
        #self.assertIsInstance(inp_dict['abi_kwargs'], collections.OrderedDict)
        assert "abi_args" in inp_dict and len(inp_dict["abi_args"]) == len(inp)
        assert all(k in inp for k, _ in inp_dict["abi_args"])
        self.assertMSONable(inp)

    def test_dfpt_methods(self):
        """Testing DFPT methods."""
        gs_inp = AbinitInput(structure=abidata.structure_from_ucell("AlAs"),
                             pseudos=abidata.pseudos("13al.981214.fhi", "33as.pspnc"))

        gs_inp.set_vars(
            nband=4,             
            ecut=2,         
            ngkpt=[4, 4, 4],
            nshiftk=4,
            shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
                    0.0, 0.5, 0.0,
                    0.5, 0.0, 0.0,
                    0.5, 0.5, 0.5],
            #shiftk=[0, 0, 0],
            paral_kgb=1,
            nstep=25,
            tolvrs=1.0e-10,
        )

        ################
        # Phonon methods
        ################
        with self.assertRaises(gs_inp.Error):
            ddk_inputs = gs_inp.make_ddk_inputs(tolerance={"tolfoo": 1e10})

        phg_inputs = gs_inp.make_ph_inputs_qpoint(qpt=(0, 0, 0), tolerance=None)
        print("phonon inputs at Gamma\n", phg_inputs)
        assert len(phg_inputs) == 2
        assert np.all(phg_inputs[0]["rfatpol"] == [1, 1])
        assert np.all(phg_inputs[1]["rfatpol"] == [2, 2])
        assert all(np.all(inp["rfdir"] == [1, 0, 0] for inp in phg_inputs))
        assert all(np.all(inp["kptopt"] == 2 for inp in phg_inputs))

        if self.has_abinit():
            # Validate
            vs = phg_inputs.abivalidate()
            assert all(v.retcode == 0 for v in vs)

        #############
        # DDK methods
        #############
        with self.assertRaises(gs_inp.Error):
            ddk_inputs = gs_inp.make_ddk_inputs(tolerance={"tolvrs": 1e10})

        ddk_inputs = gs_inp.make_ddk_inputs(tolerance=None)
        print("DDK inputs\n", ddk_inputs)
        assert len(ddk_inputs) == 3
        assert all(inp["iscf"] == -3 for inp in ddk_inputs)
        assert all(inp["rfelfd"] == 2 for inp in ddk_inputs)

        if self.has_abinit():
            # Validate
            vs = ddk_inputs.abivalidate()
            assert all(v.retcode == 0 for v in vs)
            #assert 0

    def TestInputCheckSum(self):
        """Testing the hash method of AbinitInput"""
        inp = ebands_input(abidata.cif_file("si.cif"), abidata.pseudos("14si.pspnc"), kppa=10, ecut=2)[0]
        inp_cs = inp.variable_checksum()
        ecut = inp.pop('ecut')
        inp.set_vars({'ecut': ecut})
        self.assertEqual(inp_cs, inp.variable_checksum())


class TestMultiDataset(AbipyTest):
    """Unit tests for MultiDataset."""
    def test_api(self):
        """Test MultiDataset API."""
        structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        multi = MultiDataset(structure=structure, pseudos=abidata.pseudos("14si.pspnc"))

        assert len(multi) == 1 and multi.ndtset == 1 
        for i, inp in enumerate(multi):
            assert inp.keys() == multi[i].keys()

        multi.addnew_from(0)
        assert multi.ndtset == 2 and multi[0] is not multi[1]
        assert multi[0].structure ==  multi[1].structure
        assert multi[0].structure is not multi[1].structure

        multi.set_vars(ecut=2)
        assert all(inp["ecut"] == 2 for inp in multi)
        multi[1].set_vars(ecut=1)
        assert multi[0]["ecut"] == 2 and multi[1]["ecut"] == 1

        pert_structure = structure.copy()
        pert_structure.perturb(distance=0.1)
        assert structure != pert_structure

        assert multi.set_structure(structure) == multi.ndtset * [structure]
        assert all(s == structure for s in multi.structure)
        assert multi.has_same_structures
        multi[1].set_structure(pert_structure)
        assert multi[0].structure != multi[1].structure and multi[1].structure == pert_structure
        assert not multi.has_same_structures

        split = multi.split_datasets()
        assert len(split) == 2 and all(split[i] == multi[i] for i in range(multi.ndtset))

        print(multi)
        #print(dir(multi))
        #assert 0

        new_multi = MultiDataset.from_inputs([inp for inp in multi])
        assert new_multi.ndtset == multi.ndtset
        assert new_multi.structure == multi.structure

        for old_inp, new_inp in zip(multi, new_multi):
            assert old_inp is not new_inp
            self.assertDictEqual(old_inp.as_dict(), new_inp.as_dict())

        ref_input = multi[0]
        new_multi = MultiDataset.replicate_input(input=ref_input, ndtset=4)

        assert new_multi.ndtset == 4
        for inp in new_multi:
            assert ref_input is not inp
            self.assertDictEqual(ref_input.as_dict(), inp.as_dict())

        # Compatible with Pickle and MSONable?
        #self.serialize_with_pickle(multi, test_eq=False)
        #self.assertMSONable(multi)

        # Test tags
        new_multi.add_tags([GROUND_STATE, RELAX], [0,2])
        assert len(new_multi[0].tags) == 2
        sub_multi = new_multi.filter_by_tags(GROUND_STATE)
        assert len(sub_multi) == 2
        new_multi.remove_tags([RELAX])
        assert len(new_multi[0].tags) == 1



class AnaddbInputTest(AbipyTest):
    """Tests for AnaddbInput."""

    def setUp(self):
        self.structure = abidata.structure_from_ucell("Si")

    def test_phbands_and_dos(self):
        """Test phbands_and_dos constructor."""
        inp = AnaddbInput(self.structure, comment="hello anaddb", anaddb_kwargs={"brav": 1})
        self.assertTrue("brav" in inp)
        self.assertEqual(inp["brav"], 1)
        self.assertEqual(inp.get("brav"), 1)

        # Unknown variable.
        with self.assertRaises(AnaddbInput.Error):
            AnaddbInput(self.structure, anaddb_kwargs={"foo": 1})

        ndivsm = 1
        nqsmall = 3
        ngqpt = (4, 4, 4)

        inp2 = AnaddbInput.phbands_and_dos(self.structure, ngqpt, ndivsm, nqsmall, asr=0, dos_method="tetra")
        self.assertEqual(inp2['ifcflag'], 1)

        s2 = inp2.to_string(sortmode="a")
        print(s2)

        inp3 = AnaddbInput.phbands_and_dos(self.structure, ngqpt, ndivsm, nqsmall,
                                           qptbounds=[0,0,0,1,1,1], dos_method="gaussian:0.001 eV")
        self.assertEqual(inp3['ifcflag'], 1)
        self.assertEqual(inp3['prtdos'], 1)
        s3 = inp3.to_string(sortmode="a")
        print(s3)

        # Compatible with deepcopy and Pickle?
        for i in (inp, inp2, inp3):
            self.serialize_with_pickle(i, test_eq=False)
            i.deepcopy()

    def test_thermo(self):
        """Test the thermodynamics constructor"""
        anaddb_input = AnaddbInput.thermo(self.structure, ngqpt=(40, 40, 40), nqsmall=20)
        self.assertTrue(str(anaddb_input))
        for var in ('thmtol', 'ntemper', 'temperinc', 'thmtol'):
            self.assertTrue(anaddb_input[var] >= 0)
        for flag in ('ifcflag', 'thmflag'):
            self.assertEqual(anaddb_input[flag], 1)

        self.serialize_with_pickle(anaddb_input, test_eq=False)
        anaddb_input.deepcopy()

    def test_modes(self):
        """Test the modes constructor"""
        anaddb_input = AnaddbInput.modes(self.structure)
        self.assertTrue(str(anaddb_input))
        for flag in ('ifcflag', 'dieflag'):
            self.assertEqual(anaddb_input[flag], 1)

        self.serialize_with_pickle(anaddb_input, test_eq=False)
        anaddb_input.deepcopy()


if __name__ == "__main__":
    import unittest
    unittest.main()

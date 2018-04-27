"""Tests for input  module"""
from __future__ import print_function, division, unicode_literals

import os
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.abio.inputs import *
from abipy.abio.input_tags import *

import abipy.abio.decorators as ideco
from abipy.abio.factories import *


class TestAbinitInput(AbipyTest):
    """Unit tests for AbinitInput."""

    def test_api(self):
        """Testing AbinitInput API."""
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

        repr(inp), str(inp)
        assert len(inp) == 0 and not inp
        assert inp.get("foo", "bar") == "bar" and inp.pop("foo", "bar") == "bar"
        assert inp.comment is None
        inp.set_comment("This is a comment")
        assert inp.comment == "This is a comment"
        assert inp.isnc and not inp.ispaw
        assert not inp.decorators
        assert len(inp.structure) == 2 and inp.num_valence_electrons == 8

        # foo is not a valid Abinit variable
        with self.assertRaises(inp.Error):
            inp["foo"] = 1

        # unless we deactive spell_check
        assert inp.spell_check
        inp.set_spell_check(False)
        inp["foo"] = 1
        assert inp["foo"] == 1
        inp.pop("foo")
        assert "foo" not in inp
        inp.set_spell_check(True)

        inp["ecut" ] = 1
        assert inp.get("ecut") == 1 and len(inp) == 1 and "ecut" in inp.keys() and "foo" not in inp

        # Default is kptopt 1
        assert inp.uses_ktimereversal

        assert inp.mnemonics == False
        inp.set_mnemonics(True)
        assert inp.mnemonics == True

        # Test to_string
        assert inp.to_string(sortmode="a", with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode="section", with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode=None, with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode=None, with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode=None, with_structure=True, with_pseudos=False, mode="html")
        assert inp.to_string(sortmode="a", with_structure=False, with_pseudos=False, mode="html")
        assert inp._repr_html_()

        inp.set_vars(ecut=5, toldfe=1e-6)
        assert inp["ecut"] == 5
        inp.set_vars_ifnotin(ecut=-10)
        assert inp["ecut"] == 5
        assert inp.scf_tolvar == ("toldfe", inp["toldfe"])

        inp.write(filepath=self.get_tmpname(text=True))

        # Cannot change structure variables directly.
        with self.assertRaises(inp.Error):
            inp.set_vars(unit_cell)

        with self.assertRaises(TypeError):
            inp.add_abiobjects({})

        with self.assertRaises(KeyError):
            inp.remove_vars("foo", strict=True)
        assert not inp.remove_vars("foo", strict=False)

        # Test deepcopy and remove_vars.
        inp["bdgw"] = [1, 2]
        inp_copy = inp.deepcopy()
        inp_copy["bdgw"][1] = 3
        assert inp["bdgw"] == [1, 2]
        assert inp.remove_vars("bdgw") and "bdgw" not in inp

        removed = inp.pop_tolerances()
        assert len(removed) == 1 and removed["toldfe"] == 1e-6
        assert inp.scf_tolvar == (None, None)

        # Test set_spin_mode
        old_vars = inp.set_spin_mode("polarized")
        assert "nsppol" in inp and inp["nspden"] == 2 and inp["nspinor"] == 1
        inp.set_vars(old_vars)

        # Test set_structure
        new_structure = inp.structure.copy()
        new_structure.perturb(distance=0.1)
        inp.set_structure(new_structure)
        assert inp.structure == new_structure

        # Compatible with Pickle and MSONable?
        self.serialize_with_pickle(inp, test_eq=False)
        self.assertMSONable(inp)

        # Test generate method.
        ecut_list = [10, 20]
        for i, ginp in enumerate(inp.generate(ecut=ecut_list)):
            assert ginp["ecut"] == ecut_list[i]

        inp_list = list(inp.generate(ecut=[10, 20], nsppol=[1, 2]))
        assert len(inp_list) == 4

        # Test tags
        assert isinstance(inp.tags, set)
        assert len(inp.tags) == 0
        inp.add_tags([GROUND_STATE, SCF])
        assert len(inp.tags) == 2
        inp.remove_tags([GROUND_STATE])
        assert len(inp.tags) == 1

        new_inp = AbinitInput(structure=inp.structure, pseudos=inp.pseudos, abi_kwargs=inp.vars)
        new_inp.pop_par_vars(new_inp)
        new_inp["npband"] = 2
        new_inp["npfft"] = 3
        popped = new_inp.pop_par_vars(new_inp)
        assert popped["npband"] == 2 and "npband" not in new_inp
        assert popped["npfft"] == 3 and "npfft" not in new_inp

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
        pseudo = abidata.pseudo("14si.pspnc")
        pseudo_dir = os.path.dirname(pseudo.filepath)
        inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=pseudo.basename, pseudo_dir=pseudo_dir)

        nval_atoms = inp.valence_electrons_per_atom
        assert len(nval_atoms) == 2
        assert nval_atoms == [4, 4]

        inp.set_kmesh(ngkpt=(1, 2, 3), shiftk=(1, 2, 3, 4, 5, 6))
        assert inp["kptopt"] == 1 and inp["nshiftk"] == 2
        assert inp.uses_ktimereversal

        inp.set_gamma_sampling()
        assert inp["kptopt"] == 1 and inp["nshiftk"] == 1
        assert np.all(inp["shiftk"] == 0)

        inp.set_autokmesh(nksmall=2)
        assert inp["kptopt"] == 1 and np.all(inp["ngkpt"] == [2, 2, 2]) and inp["nshiftk"] == 4

        inp.set_kpath(ndivsm=3, kptbounds=None)
        assert inp["ndivsm"] == 3 and inp["iscf"] == -2 and len(inp["kptbounds"]) == 12

        inp.set_qpath(ndivsm=3, qptbounds=None)
        assert len(inp["ph_qpath"]) == 12 and inp["ph_nqpath"] == 12 and inp["ph_ndivsm"] == 3

        inp.set_phdos_qmesh(nqsmall=16, method="tetra")
        assert inp["ph_intmeth"] == 2 and np.all(inp["ph_ngqpt"] == 16) and np.all(inp["ph_qshift"] == 0)

        inp.set_kptgw(kptgw=(1, 2, 3, 4, 5, 6), bdgw=(1, 2))
        assert inp["nkptgw"] == 2 and np.all(inp["bdgw"].ravel() == np.array(len(inp["kptgw"]) * [1,2]).ravel())

        linps = inp.linspace("ecut", 2, 6, num=3, endpoint=True)
        assert len(linps) == 3 and (linps[0]["ecut"] == 2 and (linps[-1]["ecut"] == 6))

        ranginps = inp.arange("ecut", start=3, stop=5, step=1)
        assert len(ranginps) == 2 and (ranginps[0]["ecut"] == 3 and (ranginps[-1]["ecut"] == 4))

        with self.assertRaises(inp.Error):
            inp.product("ngkpt", "tsmear", [[2, 2, 2], [4, 4, 4]])

        prod_inps = inp.product("ngkpt", "tsmear", [[2, 2, 2], [4, 4, 4]], [0.1, 0.2, 0.3])
        assert len(prod_inps) == 6
        assert prod_inps[0]["ngkpt"] == [2, 2, 2] and prod_inps[0]["tsmear"] == 0.1
        assert prod_inps[-1]["ngkpt"] ==  [4, 4, 4] and prod_inps[-1]["tsmear"] == 0.3

        inp["kptopt"] = 4
        assert not inp.uses_ktimereversal

    # TODO
    def test_new_with_structure(self):
        """Testing new_with_structure."""

        si2_inp = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"),
                abi_kwargs={"ecut": 4, "toldfe": 1e-10, "nband": 6, "ngkpt": [12, 12, 12]})

        # Build new input with same parameters and compressed unit cell.
        new_vol = 0.9 * si2_inp.structure.volume
        new_lattice = si2_inp.structure.lattice.scale(new_vol)
        new_structure = abilab.Structure(new_lattice, si2_inp.structure.species, si2_inp.structure.frac_coords)

        new_inp = si2_inp.new_with_structure(new_structure, scdims=None)
        assert new_inp.structure.volume == new_vol
        assert new_inp["ecut"] == 4

        #al2_structure = abilab.Structure.from_file(abidata.cif_file("al.cif"))

        # Let's build an input file for a (2, 3, 4) supercell
        #super_structure = si2_inp.structure.make_supercell
        scdims = np.array((2, 3, 4))
        # Note: This will return a pymatgen structure, not an Abipy structure.
        super_structure = si2_inp.structure * scdims
        sc_inp = si2_inp.new_with_structure(super_structure, scdims=scdims)
        assert isinstance(sc_inp.structure, abilab.Structure)
        assert len(sc_inp.structure) == 2 * scdims.prod()
        self.assert_almost_equal(sc_inp.structure.volume, si2_inp.structure.volume * scdims.prod())
        assert sc_inp["chkprim"] == 0
        assert sc_inp["nband"] == si2_inp["nband"] * scdims.prod()
        self.assert_equal(sc_inp["ngkpt"], [6, 4, 3])
        self.abivalidate_input(sc_inp)
        #sc_inp.write(filepath=self.tmpfileindir("run.abi"))

        # TODO
        # Test new_with_structure with nsppol == 2
        # In this case, we have the `spinat` array that depends on natom
        #si2_inp["nsppol"] = 2
        #si2_inp["spinat"] = [[0, 0, 1], [0, 0 -1]]
        #scdims = np.array((1, 1, 2))
        #super_structure = si2_inp.structure * scdims
        #with self.assertRaises(ValueError):
        #    si2_inp.new_with_structure(new_structure, scdims=scdims)

        #new_inp = si2_inp.new_with_structure(super_structure, scdims=scdims)
        #self.abivalidate_input(new_inp)

    def test_abinit_calls(self):
        """Testing AbinitInput methods invoking Abinit."""
        inp_si = AbinitInput(structure=abidata.cif_file("si.cif"), pseudos=abidata.pseudos("14si.pspnc"))
        inp_si.set_kmesh(ngkpt=(2, 2, 2), shiftk=(0, 0, 0))

        inp_gan = AbinitInput(structure=abidata.cif_file("gan.cif"),
                              pseudos=abidata.pseudos("31ga.pspnc", "7n.pspnc"))
        inp_gan.set_kmesh(ngkpt=(2, 2, 2), shiftk=(0, 0, 0))

        # The code below invokes Abinit (the test must fail because of wrong input)
        inp_si.set_vars(ecut=-1)
        self.abivalidate_input(inp_si, must_fail=True)

        # Test validate with correct input
        inp_si.set_vars(ecut=2, toldfe=1e-6)
        self.abivalidate_input(inp_si)

        # TODO: Here spglib and abinit do not agree.
        # Test abiget_spacegroup
        #structure_with_abispg = inp_gan.abiget_spacegroup()
        #assert structure_with_abispg.abispg is not None
        #assert structure_with_abispg.abispg.spgid == 227

        # Test abiget_spacegroup for Si
        structure_with_abispg = inp_si.abiget_spacegroup()
        assert structure_with_abispg.abi_spacegroup is not None
        assert structure_with_abispg.abi_spacegroup.spgid == 227

        # Test abiget_ibz
        ibz = inp_si.abiget_ibz()
        assert np.all(ibz.points == [[ 0.,  0.,  0.], [0.5,  0.,  0.], [0.5, 0.5, 0.]])
        assert np.all(ibz.weights == [0.125,  0.5,  0.375])

        # This to test what happes with wrong inputs and Abinit errors.
        wrong = inp_si.deepcopy()
        removed = wrong.pop_vars("ecut")
        assert "ecut" not in wrong
        assert "ecut" in removed
        assert removed["ecut"] == 2
        with self.assertRaises(wrong.Error):
            #wrong["ecut"] = -12.0
            wrong.abiget_ibz(ngkpt=[-1, -1, -1])

        # Test abiget_irred_phperts
        # [{'idir': 1, 'ipert': 1, 'qpt': [0.0, 0.0, 0.0]}]
        irred_perts = inp_si.abiget_irred_phperts(qpt=(0, 0, 0))
        assert len(irred_perts) == 1
        pert = irred_perts[0]
        assert pert.idir == 1 and (pert.idir, pert.ipert) == (1, 1) and all(c == 0 for c in pert.qpt)

        irred_perts = inp_gan.abiget_irred_phperts(qpt=(0.5, 0, 0))
        print(irred_perts)
        assert len(irred_perts) == 6
        irred_perts_values = [{'idir': 1, 'ipert': 1, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 2, 'ipert': 1, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 3, 'ipert': 1, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 1, 'ipert': 3, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 2, 'ipert': 3, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 3, 'ipert': 3, 'qpt': [0.5, 0.0, 0.0]}]
        for a, b in zip(irred_perts, irred_perts_values):
            self.assertDictEqual(a, b)

        # Test abiget_autoparal_pconfs
        inp_si["paral_kgb"] = 0
        pconfs = inp_si.abiget_autoparal_pconfs(max_ncpus=5)
        inp_si["paral_kgb"] = 1
        pconfs = inp_si.abiget_autoparal_pconfs(max_ncpus=5)

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

        # qpt is not in gs_inp and not passed to method.
        with self.assertRaises(ValueError):
            gs_inp.abiget_irred_phperts()

        ################
        # Phonon methods
        ################
        with self.assertRaises(gs_inp.Error):
            try:
                ddk_inputs = gs_inp.make_ddk_inputs(tolerance={"tolfoo": 1e10})
            except Exception as exc:
                print(exc)
                raise

        phg_inputs = gs_inp.make_ph_inputs_qpoint(qpt=(0, 0, 0), tolerance=None)
        #print("phonon inputs at Gamma\n", phg_inputs)
        assert len(phg_inputs) == 2
        assert np.all(phg_inputs[0]["rfatpol"] == [1, 1])
        assert np.all(phg_inputs[1]["rfatpol"] == [2, 2])
        assert all(np.all(inp["rfdir"] == [1, 0, 0] for inp in phg_inputs))
        assert all(np.all(inp["kptopt"] == 2 for inp in phg_inputs))

        # Validate with Abinit
        self.abivalidate_multi(phg_inputs)

        ##############
        # BEC methods
        ##############
        with self.assertRaises(gs_inp.Error):
            gs_inp.make_bec_inputs(tolerance={"tolfoo": 1e10})

        bec_inputs = gs_inp.make_bec_inputs(tolerance=None)
        print("BEC inputs at Gamma\n", bec_inputs)
        assert len(bec_inputs) == 2
        assert np.all(bec_inputs[0]["rfatpol"] == [1, 1])
        assert np.all(bec_inputs[1]["rfatpol"] == [2, 2])
        assert all(np.all(inp["rfelfd"] == 3 for inp in bec_inputs))
        assert all(np.all(inp["tolvrs"] == 1.0e-10 for inp in bec_inputs))
        assert all(np.all(inp["kptopt"] == 2 for inp in bec_inputs))

        # Validate with Abinit
        self.abivalidate_multi(bec_inputs)

        #############
        # DDK methods
        #############
        with self.assertRaises(gs_inp.Error):
            gs_inp.make_ddk_inputs(tolerance={"tolvrs": 1e10})

        ddk_inputs = gs_inp.make_ddk_inputs(tolerance=None)
        print("DDK inputs\n", ddk_inputs)
        assert len(ddk_inputs) == 3
        assert all(inp["iscf"] == -3 for inp in ddk_inputs)
        assert all(inp["rfelfd"] == 2 for inp in ddk_inputs)

        # Validate with Abinit
        self.abivalidate_multi(ddk_inputs)

        #############
        # DDE methods
        #############
        with self.assertRaises(gs_inp.Error):
            gs_inp.make_dde_inputs(tolerance={"tolfoo": 1e10})

        dde_inputs = gs_inp.make_dde_inputs(tolerance=None, use_symmetries=False)
        assert len(dde_inputs) == 3
        assert all(inp["prepanl"] == 1 for inp in dde_inputs)
        for i, dde_rfdirs in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
            assert np.all(dde_inputs[i]["rfdir"] == dde_rfdirs)
        self.abivalidate_multi(dde_inputs)

        dde_inputs = gs_inp.make_dde_inputs(tolerance=None, use_symmetries=True)
        print("DDE inputs\n", dde_inputs)
        assert len(dde_inputs) == 1
        assert all(np.all(inp["tolvrs"] == 1.0e-22 for inp in dde_inputs))
        assert all(inp["rfelfd"] == 3 for inp in dde_inputs)
        assert all(inp["kptopt"] == 2 for inp in dde_inputs)

        # Validate with Abinit
        self.abivalidate_multi(dde_inputs)

        #################
        # Strain methods
        #################
        strain_inputs = gs_inp.make_strain_perts_inputs(tolerance=None)
        print("STRAIN inputs\n", strain_inputs)
        assert all(inp["tolvrs"] == 1e-12 for inp in strain_inputs)

        # Validate with Abinit
        self.abivalidate_multi(strain_inputs)

        #####################
        # Non-linear methods
        ####################
        if self.has_abinit(version='8.3.2'):
            dte_inputs = gs_inp.make_dte_inputs(phonon_pert=True, skip_permutations=True)
            print("dte inputs\n", dte_inputs)
            assert len(dte_inputs) == 8
            assert np.all(dte_inputs[0]["d3e_pert2_dir"] == [1, 0, 0])
            assert np.all(dte_inputs[3]["d3e_pert1_atpol"] == [2, 2])
            assert all(np.all(inp["optdriver"] == 5 for inp in dte_inputs))

            # Validate with Abinit
            self.abivalidate_multi(dte_inputs)

    def TestInputCheckSum(self):
        """Testing the hash method of AbinitInput"""
        inp = ebands_input(abidata.cif_file("si.cif"), abidata.pseudos("14si.pspnc"), kppa=10, ecut=2)[0]
        inp_cs = inp.variable_checksum()
        ecut = inp.pop('ecut')
        inp.set_vars({'ecut': ecut})
        assert inp_cs == inp.variable_checksum()


class TestMultiDataset(AbipyTest):
    """Unit tests for MultiDataset."""

    def test_api(self):
        """Testing MultiDataset API."""
        structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        pseudo = abidata.pseudo("14si.pspnc")
        pseudo_dir = os.path.dirname(pseudo.filepath)
        multi = MultiDataset(structure=structure, pseudos=pseudo)
        with self.assertRaises(ValueError):
            MultiDataset(structure=structure, pseudos=pseudo, ndtset=-1)

        multi = MultiDataset(structure=structure, pseudos=pseudo.basename, pseudo_dir=pseudo_dir)

        assert len(multi) == 1 and multi.ndtset == 1
        assert multi.isnc
        for i, inp in enumerate(multi):
            assert list(inp.keys()) == list(multi[i].keys())

        multi.addnew_from(0)
        assert multi.ndtset == 2 and multi[0] is not multi[1]
        assert multi[0].structure ==  multi[1].structure
        assert multi[0].structure is not multi[1].structure

        multi.set_vars(ecut=2)
        assert all(inp["ecut"] == 2 for inp in multi)
        self.assert_equal(multi.get("ecut"), [2, 2])

        df = multi.get_vars_dataframe("ecut", "foobar")
        assert "ecut" in df
        self.assert_equal(df["ecut"].values, [2, 2])

        multi[1].set_vars(ecut=1)
        assert multi[0]["ecut"] == 2 and multi[1]["ecut"] == 1
        self.assert_equal(multi.get("ecut"), [2, 1])

        self.assert_equal(multi.get("foo", "default"), ["default", "default"])

        multi[1].set_vars(paral_kgb=1)
        assert "paral_kgb" not in multi[0]
        self.assert_equal(multi.get("paral_kgb"), [None, 1])

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
        repr(multi); str(multi)
        assert multi.to_string(mode="text")
        assert multi._repr_html_()

        inp.write(filepath=self.tmpfileindir("run.abi"))
        multi.write(filepath=self.tmpfileindir("run.abi"))

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
        self.serialize_with_pickle(multi, test_eq=False)
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
        """Testing phbands_and_dos constructor."""
        inp = AnaddbInput(self.structure, comment="hello anaddb", anaddb_kwargs={"brav": 1})
        repr(inp); str(inp)
        assert inp.to_string(sortmode="a")
        assert inp._repr_html_()
        assert "brav" in inp
        assert inp["brav"] == 1
        assert inp.get("brav") == 1

        self.serialize_with_pickle(inp, test_eq=False)
        #self.assertMSONable(inp, test_if_subclass=False)

        # Unknown variable.
        with self.assertRaises(AnaddbInput.Error):
            AnaddbInput(self.structure, anaddb_kwargs={"foo": 1})

        # unless we deactive spell_check
        assert inp.spell_check
        inp.set_spell_check(False)
        inp["foo"] = 1
        assert inp["foo"] == 1

        inp.pop("foo")
        assert "foo" not in inp
        inp.set_spell_check(True)

        ndivsm = 1
        nqsmall = 3
        ngqpt = (4, 4, 4)

        inp2 = AnaddbInput.phbands_and_dos(self.structure, ngqpt, ndivsm, nqsmall, asr=0, dos_method="tetra")
        assert inp2['ifcflag'] == 1
        print(inp2.to_string(sortmode="a"))
        self.abivalidate_input(inp2)

        inp3 = AnaddbInput.phbands_and_dos(self.structure, ngqpt, ndivsm, nqsmall,
                                           qptbounds=[0,0,0,1,1,1], dos_method="gaussian:0.001 eV")
        assert inp3['ifcflag'] == 1
        assert inp3['prtdos'] == 1
        print(inp3.to_string(sortmode="a"))
        self.abivalidate_input(inp3)

        # Compatible with deepcopy and Pickle?
        for i in (inp, inp2, inp3):
            self.serialize_with_pickle(i, test_eq=False)
            i.deepcopy()

        # Test AnaddbInput with lo_to_splitting.
        inp_loto = AnaddbInput.phbands_and_dos(self.structure, ngqpt, ndivsm, nqsmall,
                                               lo_to_splitting=True)
        #print(inp_loto)
        #print(inp_loto["qpath"], inp_loto["qph2l"])
        assert "qpath" in inp_loto
        assert inp_loto["nph2l"] == 3
        self.assert_almost_equal(inp_loto["qph2l"],
            [[ 0.        , 0.184959  , 0.       , 0.],
             [ 0.13871925, 0.13871925, 0.       , 0.],
             [ 0.0924795 , 0.0924795 , 0.0924795, 0.]])
        self.abivalidate_input(inp_loto)

    def test_thermo(self):
        """Testing the thermodynamics constructor"""
        anaddb_input = AnaddbInput.thermo(self.structure, ngqpt=(40, 40, 40), nqsmall=20)
        assert str(anaddb_input)
        for var in ('thmtol', 'ntemper', 'temperinc', 'thmtol'):
            assert anaddb_input[var] >= 0
        for flag in ('ifcflag', 'thmflag'):
            assert anaddb_input[flag] == 1

        self.serialize_with_pickle(anaddb_input, test_eq=False)
        anaddb_input.deepcopy()

    def test_modes(self):
        """Testing modes constructor"""
        anaddb_input = AnaddbInput.modes(self.structure)
        assert str(anaddb_input)
        for flag in ('ifcflag', 'dieflag'):
            assert anaddb_input[flag] == 1
        self.abivalidate_input(anaddb_input)

        self.serialize_with_pickle(anaddb_input, test_eq=False)
        anaddb_input.deepcopy()

    def test_ifc(self):
        """Testing ifc constructor"""
        anaddb_input = AnaddbInput.ifc(self.structure, ngqpt=[4,4,4])
        assert str(anaddb_input)
        for flag in ('ifcflag', 'dipdip'):
            assert anaddb_input[flag] == 1

        self.serialize_with_pickle(anaddb_input, test_eq=False)
        anaddb_input.deepcopy()
        self.abivalidate_input(anaddb_input)

    def test_piezo_elastic(self):
        """Testing piezo_elastic constructor."""
        anaddb_input = AnaddbInput.piezo_elastic(self.structure, stress_correction=True)
        assert anaddb_input["elaflag"] == 5 and anaddb_input["piezoflag"] == 3 and anaddb_input["asr"] == 0
        self.abivalidate_input(anaddb_input)

        anaddb_input = AnaddbInput.piezo_elastic(self.structure, stress_correction=False)
        assert anaddb_input["elaflag"] == 3 and anaddb_input["piezoflag"] == 3 and anaddb_input["chneut"] == 1
        self.abivalidate_input(anaddb_input)


class TestCut3DInput(AbipyTest):
    """Unit tests for AbinitInput."""

    def setUp(self):
        self.structure = abidata.structure_from_ucell("Si")

    def test_dict_methods(self):
        cut3d_input = Cut3DInput.den_to_cube('/path/to/den', 'outfile_name')
        repr(cut3d_input); str(cut3d_input)
        cut3d_input.write(self.get_tmpname(text=True))

        self.serialize_with_pickle(cut3d_input, test_eq=False)
        self.assertMSONable(cut3d_input, test_if_subclass=False)

    def test_generation_methods(self):
        cut3d_input = Cut3DInput.den_to_cube('/path/to/den', 'outfile_name')
        cut3d_input = Cut3DInput.den_to_xsf('/path/to/den', 'outfile_name', shift=[2, 2, 2])
        cut3d_input = Cut3DInput.den_to_3d_indexed('/path/to/den', 'outfile_name')
        cut3d_input = Cut3DInput.den_to_3d_formatted('/path/to/den', 'outfile_name')
        cut3d_input = Cut3DInput.den_to_tecplot('/path/to/den', 'outfile_name')
        cut3d_input = Cut3DInput.den_to_molekel('/path/to/den', 'outfile_name')

        # TODO
        #cut3d_input = Cut3DInput.hirshfeld(density_filepath, all_el_dens_paths)
        #cut3d_input = Cut3DInput.hirshfeld_from_fhi_path(density_filepath, structure, fhi_all_el_path)


class OpticInputTest(AbipyTest):
    """Tests for OpticInput."""

    def test_optic_input_api(self):
        """Testing OpticInput API."""
        optic_input = OpticInput()
        repr(optic_input); str(optic_input)
        #assert optic_input._repr_html_()

        for var in OpticInput._VARIABLES:
            repr(var); str(var)
            assert str(var.help) and var.group
            assert str(var.html_link(label="foo"))

        with self.assertRaises(optic_input.Error):
            optic_input["foo"] = 23

        with self.assertRaises(optic_input.Error):
            optic_input.get_default("foo")

        optic_input.set_vars(num_lin_comp=0)
        assert optic_input["num_lin_comp"] == 0
        assert optic_input["broadening"] == 0.01

    def test_real_optic_input(self):
        """Testing real optic input."""
        # Optic does not support MPI with ncpus > 1.
        optic_input = OpticInput(
            broadening=0.002,          # Value of the smearing factor, in Hartree
            domega=0.0003,             # Frequency mesh.
            maxomega=0.3,
            scissor=0.000,             # Scissor shift if needed, in Hartree
            tolerance=0.002,           # Tolerance on closeness of singularities (in Hartree)
            num_lin_comp=1,            # Number of components of linear optic tensor to be computed
            lin_comp=11,               # Linear coefficients to be computed (x=1, y=2, z=3)
            num_nonlin_comp=2,         # Number of components of nonlinear optic tensor to be computed
            nonlin_comp=(123, 222),    # Non-linear coefficients to be computed
        )

        repr(optic_input); str(optic_input)
        assert optic_input.to_string(verbose=2)
        # TODO
        #assert optic_input._repr_html_()
        assert optic_input.vars

        # Compatible with Pickle and MSONable?
        self.serialize_with_pickle(optic_input, test_eq=True)

        # TODO: But change function that build namelist to ignore @module ...
        #self.assertMSONable(optic_input)

        self.abivalidate_input(optic_input)

        # Test helper functions
        si_structure = abidata.structure_from_ucell("Si")
        comps = optic_input.only_independent_chi_components(si_structure)
        assert comps == ["xx"]
        assert optic_input["num_lin_comp"] == 1
        assert optic_input["lin_comp"] == [11]
        self.abivalidate_input(optic_input)
        #print(optic_input)

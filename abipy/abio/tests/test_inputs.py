"""Tests for inputs module"""
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
        assert inp.comment == "# This is a comment"
        assert inp.isnc and not inp.ispaw
        assert not inp.decorators
        assert len(inp.structure) == 2 and inp.num_valence_electrons == 8

        # foo is not a valid Abinit variable
        with self.assertRaises(inp.Error):
            inp["foo"] = 1

        with self.assertRaises(inp.Error): inp._check_nsppol_nspinor(3, 1)
        with self.assertRaises(inp.Error): inp._check_nsppol_nspinor(2, 2)
        with self.assertRaises(inp.Error): inp._check_nsppol_nspinor(1, 4)
        with self.assertRaises(inp.Error): inp.set_cutoffs_for_accuracy("normal")

        # unless we deactivate spell_check
        assert inp.spell_check
        inp.set_spell_check(False)
        inp["foo"] = 1
        assert inp["foo"] == 1
        inp.pop("foo")
        assert "foo" not in inp
        inp.set_spell_check(True)

        inp["ecut"] = 1
        assert inp.get("ecut") == 1 and len(inp) == 1 and "ecut" in inp.keys() and "foo" not in inp

        # Default is kptopt 1
        assert inp.uses_ktimereversal

        assert not inp.mnemonics
        inp.set_mnemonics(True)
        assert inp.mnemonics

        # Test to_string
        assert inp.to_string(sortmode="a", with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode="section", with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode=None, with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode=None, with_structure=True, with_pseudos=True)
        assert inp.to_string(sortmode=None, with_structure=True, with_pseudos=False, mode="html")
        assert inp.to_string(sortmode="a", with_structure=False, with_pseudos=False, mode="html")
        assert inp._repr_html_()

        inp.set_vars(ecut=5, toldfe=1e-6, comment="hello")
        assert inp["ecut"] == 5
        assert inp.comment == "# hello"
        inp.set_vars_ifnotin(ecut=-10)
        assert inp["ecut"] == 5
        assert inp.scf_tolvar == ("toldfe", inp["toldfe"])

        inp.write(filepath=self.get_tmpname(text=True))
        assert inp.make_targz().endswith(".tar.gz")

        # Cannot change structure variables directly.
        with self.assertRaises(inp.Error):
            inp.set_vars(unit_cell)

        with self.assertRaises(TypeError):
            inp.add_abiobjects({})

        with self.assertRaises(KeyError):
            inp.remove_vars("foo", strict=True)
        assert not inp.remove_vars("foo", strict=False)

        with self.assertRaises(inp.Error):
            # Pseudos do not provide hints
            inp.set_cutoffs_for_accuracy("normal")

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

        new_inp.set_vars(shiftk=[0, 0, 0, 0.5, 0, 0, 0, 0, 0.5])
        assert new_inp["nshiftk"] == 3
        other_inp = new_inp.new_with_vars(ph_qpath=[0, 0, 0, 0.5, 0, 0])
        assert other_inp["ph_nqpath"] == 2

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
        ngkpt, shiftk = inp.get_ngkpt_shiftk()
        assert ngkpt.tolist() == [1, 2, 3]
        assert len(shiftk) == 2 and shiftk.ravel().tolist() == [1, 2, 3, 4, 5, 6]

        inp.pop("ngkpt")
        kptrlatt = [1, 0, 0, 0, 4, 0, 0, 0, 8]
        shiftk = (0.5, 0.0, 0.0)
        inp.set_vars(kptrlatt=kptrlatt, nshiftk=1, shiftk=shiftk)
        ngkpt, shiftk = inp.get_ngkpt_shiftk()
        assert ngkpt.tolist() == [1, 4, 8]
        assert len(shiftk) == 1 and shiftk.ravel().tolist() == [0.5, 0.0, 0.0]

        inp.set_vars(kptrlatt=[1, 2, 0, 0, 4, 0, 0, 0, 8], nshiftk=1, shiftk=shiftk)
        ngkpt, shiftk = inp.get_ngkpt_shiftk()
        assert ngkpt is None

        inp.set_gamma_sampling()
        assert inp["kptopt"] == 1 and inp["nshiftk"] == 1
        assert np.all(inp["shiftk"] == 0)

        inp.set_autokmesh(nksmall=2)
        assert inp["kptopt"] == 1 and np.all(inp["ngkpt"] == [2, 2, 2]) and inp["nshiftk"] == 4

        inp.set_kpath(ndivsm=3, kptbounds=None)
        assert inp["ndivsm"] == 3 and inp["iscf"] == -2 and len(inp["kptbounds"]) == 12

        inp.set_kpath(ndivsm=-20)
        assert inp["nkpt"] == 156 and inp["iscf"] == -2

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
        assert prod_inps[-1]["ngkpt"] == [4, 4, 4] and prod_inps[-1]["tsmear"] == 0.3

        inp["kptopt"] = 4
        assert not inp.uses_ktimereversal

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
        sc_inp["chksymtnons"] = 0
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
        inp_gan["ecut"] = 2

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

        scr_ibz = inp_si.abiget_scr_ibz()
        assert np.all(scr_ibz.points == [[ 0.,  0.,  0.], [0.5,  0.,  0.], [0.5, 0.5, 0.]])
        assert np.all(scr_ibz.weights == [0.125,  0.5,  0.375])

        # This to test what happens with wrong inputs and Abinit errors.
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
        #print(irred_perts)
        assert len(irred_perts) == 6
        irred_perts_values = [{'idir': 1, 'ipert': 1, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 2, 'ipert': 1, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 3, 'ipert': 1, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 1, 'ipert': 3, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 2, 'ipert': 3, 'qpt': [0.5, 0.0, 0.0]},
                              {'idir': 3, 'ipert': 3, 'qpt': [0.5, 0.0, 0.0]}]
        for a, b in zip(irred_perts, irred_perts_values):
            # The nkpt_rbz entry was added in Abinit v9.5 but it's not used by AbiPy.
            if "nkpt_rbz" in a: b["nkpt_rbz"] = a["nkpt_rbz"]
            self.assertDictEqual(a, b)

        # Test abiget_autoparal_pconfs
        inp_si["paral_kgb"] = 0
        pconfs = inp_si.abiget_autoparal_pconfs(max_ncpus=5)
        inp_si["paral_kgb"] = 1
        pconfs = inp_si.abiget_autoparal_pconfs(max_ncpus=5)

        # Test run_in_shell
        assert inp_si.run_in_shell()

    def test_dict_methods(self):
        """ Testing AbinitInput dict methods """
        inp = ebands_input(abidata.cif_file("si.cif"), abidata.pseudos("14si.pspnc"), kppa=10, ecut=2)[0]
        inp = ideco.SpinDecorator("spinor")(inp)
        inp_dict = inp.as_dict()
        #self.assertIsInstance(inp_dict['abi_kwargs'], collections.OrderedDict)
        assert "abi_args" in inp_dict and len(inp_dict["abi_args"]) == len(inp)
        assert all(k in inp for k, _ in inp_dict["abi_args"])
        self.assertMSONable(inp)

    def test_enforce_znucl_and_typat(self):
        """
        Test the order of typat and znucl in the Abinit input when enforce_typat, enforce_znucl are used.
        """
        # Ga  Ga1  1  0.33333333333333  0.666666666666667  0.500880  1.0
        # Ga  Ga2  1  0.66666666666667  0.333333333333333  0.000880  1.0
        # N  N3  1  0.333333333333333  0.666666666666667  0.124120  1.0
        # N  N4  1  0.666666666666667  0.333333333333333  0.624120  1.0
        gan2 = Structure.from_file(abidata.cif_file("gan2.cif"))
        pseudos = abidata.pseudos("31ga.pspnc", "7n.pspnc")

        def_znucl = [31, 7]
        def_typat = [1, 1, 2, 2]

        # Build AbinitInput with default sorting algorithm for znucl and typat
        def_inp = AbinitInput(gan2, pseudos)
        def_dict = def_inp.as_dict()
        assert def_dict["enforce_znucl"] is None
        assert def_dict["enforce_typat"] is None

        # Build AbinitInput with specific znucl and typat.
        enforce_znucl = [7 ,31]
        enforce_typat = [2, 2, 1, 1]
        enf_inp = AbinitInput(gan2, pseudos, enforce_znucl=enforce_znucl, enforce_typat=enforce_typat)

        enf_dict = enf_inp.as_dict()
        self.assert_equal(enf_dict["enforce_znucl"], enforce_znucl)
        self.assert_equal(enf_dict["enforce_typat"], enforce_typat)

        # Parser the input string and make sure znucl and typat are what we expect.
        from abipy.abio.abivars import AbinitInputFile
        def_string = def_inp.to_string()
        enf_string = enf_inp.to_string()
        #print("\ndef_string\n", def_string, "\nenf_string:\n", enf_string)
        assert def_string != enf_string
        def_inpfile = AbinitInputFile.from_string(def_string)
        enf_inpfile = AbinitInputFile.from_string(enf_string)
        assert def_inpfile.structure == enf_inpfile.structure

        self.assert_equal(def_inpfile.datasets[0]["znucl"], " ".join(str(e) for e in def_znucl))
        self.assert_equal(def_inpfile.datasets[0]["typat"], " ".join(str(t) for t in def_typat))
        self.assert_equal(enf_inpfile.datasets[0]["znucl"], " ".join(str(e) for e in enforce_znucl))
        self.assert_equal(enf_inpfile.datasets[0]["typat"], " ".join(str(t) for t in enforce_typat))

        other_enf_dict = AbinitInput.from_dict(enf_dict).as_dict()
        self.assert_equal(other_enf_dict["enforce_znucl"], enforce_znucl)
        self.assert_equal(other_enf_dict["enforce_typat"], enforce_typat)

        # Make sure we detect wrong calls.
        with self.assertRaises(ValueError):
            AbinitInput(gan2, pseudos, enforce_znucl=[1,], enforce_typat=enforce_typat)
        with self.assertRaises(ValueError):
            AbinitInput(gan2, pseudos, enforce_znucl=enforce_znucl, enforce_typat=[1, 2])

        # Make sure enforce flags are propagated in MultiDataset and other methods.
        inputs = enf_inp.product("ngkpt", "tsmear", [[2, 2, 2], [4, 4, 4]], [0.1, 0.2, 0.3])
        for inp in inputs:
            self.assert_equal(inp.enforce_znucl, enforce_znucl)
            self.assert_equal(inp.enforce_typat, enforce_typat)

        multi = MultiDataset.from_inputs(inputs)
        for znucl in multi.enforce_znucl:
            self.assert_equal(znucl, enforce_znucl)
        for typat in multi.enforce_typat:
            self.assert_equal(typat, enforce_typat)

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
            paral_kgb=1,
            nstep=25,
            tolvrs=1.0e-10,
        )

        assert gs_inp.is_input and not gs_inp.is_multidataset
        multi = gs_inp.replicate(2)
        assert multi.ndtset == 2 and multi.is_multidataset and not multi.is_input

        # Test make_nscf_kptopt0_input
        nscf_inp = gs_inp.make_nscf_kptopt0_input(kpts=[1, 2, 3, 4, 5, 6])
        assert "ngkpt" not in nscf_inp and "shiftk" not in nscf_inp
        assert nscf_inp["kptopt"] == 0 and nscf_inp["nkpt"] == 2 and nscf_inp["iscf"] == -2

        # Test make_ebands_input
        nscf_inp = gs_inp.make_ebands_input(ndivsm=3, tolwfr=1e-5)
        assert "ngkpt" not in nscf_inp and "shiftk" not in nscf_inp
        assert nscf_inp["iscf"] == -2 and nscf_inp["tolwfr"] == 1e-5
        assert nscf_inp["nband"] == gs_inp["nband"] + 10

        nscf_inp = gs_inp.make_ebands_input(ndivsm=-20, tolwfr=1e-5)
        assert "ngkpt" not in nscf_inp and "shiftk" not in nscf_inp
        assert nscf_inp["iscf"] == -2 and nscf_inp["tolwfr"] == 1e-5
        assert nscf_inp["nkpt"] == 149
        self.abivalidate_input(nscf_inp)

        # Test make_ebands_input with nband = "*91" (occopt 2)
        gs_occopt2 = gs_inp.new_with_vars(nband="*91", occopt=2)
        nscf_inp = gs_occopt2.make_ebands_input(nb_extra=10)
        assert nscf_inp["nband"] == "*101"
        self.abivalidate_input(nscf_inp)

        # Test make_edos_input
        ngkpt = [4, 4, 4]
        shiftk = [(0.5, 0.5, 0.5)]
        dos_input = gs_inp.make_edos_input(ngkpt=ngkpt, shiftk=shiftk, nscf_nband=9)
        self.assert_equal(dos_input["ngkpt"], ngkpt)
        self.assert_equal(dos_input["shiftk"], np.array(shiftk))
        assert dos_input["nshiftk"] == 1
        assert dos_input["iscf"] == -2 and dos_input["tolwfr"] == 1e-20
        assert dos_input["nband"] == 9

        # Test make_dfpt_effmass_input
        multi =  gs_inp.make_dfpt_effmass_inputs(kpts=[0, 0, 0, 0.5, 0, 0], effmass_bands_f90=[1, 4, 5, 5])
        assert len(multi) == 3
        assert all(inp["kptopt"] == 0 for inp in multi)
        assert all(inp["nkpt"] == 2 for inp in multi)
        assert multi.is_multidataset and not multi.is_input
        assert multi.make_targz().endswith(".tar.gz")

        inp0, inp1, inp2 = multi
        assert inp0["iscf"] == -2
        assert inp1["rfelfd"] == 2 and inp1["efmas"] == 1 and inp1["efmas_calc_dirs"] == 1 and inp1["efmas_n_dirs"] == 7
        assert inp2["eph_frohlichm"] == 1 and inp2["eph_task"] == 6 and inp2["asr"] == 2 and inp2["chneut"] == 1

        # Validate with Abinit
        self.abivalidate_multi(multi)

        ################
        # Phonon methods
        ################
        # qpt is not in gs_inp and not passed to method.
        with self.assertRaises(ValueError):
            gs_inp.abiget_irred_phperts()

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

        #######################
        # dynamical quadrupoles
        #######################
        quad_input = gs_inp.make_quad_input()
        assert quad_input["kptopt"] == 2
        assert quad_input["optdriver"] == 10
        assert quad_input["lw_qdrpl"] == 1
        assert quad_input["useylm"] == 1
        assert quad_input["nstep"] == 100

        # Validate with Abinit
        # This is expected to fail since we have NC pseudos with NLCC.
        #self.abivalidate_input(quad_input)

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

        oneshot_ddk_inputs = gs_inp.make_ddk_inputs(kptopt=3, only_vk=True)
        for inp in oneshot_ddk_inputs:
            assert inp["kptopt"] == 3
            assert inp["nstep"] == 1
            assert inp["nline"] == 1

        #############
        # DKDK methods
        #############
        dkdk_input = gs_inp.make_dkdk_input()
        assert dkdk_input["kptopt"] == 2
        assert dkdk_input["rf2_dkdk"] == 1
        assert dkdk_input["useylm"] == 1
        # Validate with Abinit
        self.abivalidate_input(dkdk_input)

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
        #print("DDE inputs\n", dde_inputs)
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

        ###############
        # EPH mobility
        ###############

        nscf_inp = gs_inp.new_with_vars(iscf=-2, tolwfr=1e-5, comment="this is just to have a NSCF input.")
        sigma_erange = sigma_kerange = [0, 1.0]
        sigma_ngkpt = [4, 4, 4]
        kerange_inputs = nscf_inp.make_wfk_kerange_inputs(sigma_kerange, sigma_ngkpt, einterp=(1, 5, 0, 0))
        inp_0, inp_1 = kerange_inputs
        assert inp_0["wfk_task"] == '"wfk_kpts_erange"'
        assert "iscf" not in inp_0
        assert inp_0["sigma_erange"] == sigma_erange
        assert inp_0["kptopt"] == 1
        assert inp_0["sigma_ngkpt"] == sigma_ngkpt

        assert inp_1["iscf"] == -2
        assert inp_1["kptopt"] == 0
        assert inp_1["ngkpt"] == sigma_ngkpt

        self.abivalidate_multi(kerange_inputs)

        ddb_ngqpt = [4, 4, 4]
        tmesh = [0, 300, 10]

        trans_inp = nscf_inp.make_eph_transport_input(ddb_ngqpt, sigma_erange, tmesh, kptopt=2, eph_ngqpt_fine=None,
                                                      mixprec=1, boxcutmin=1.1, ibte_prep=0, ibte_niter=200,
                                                      ibte_abs_tol=1e-3)
        assert trans_inp["optdriver"] == 7
        assert trans_inp["eph_task"] == -4
        assert trans_inp["kptopt"] == 2
        assert trans_inp["eph_ngqpt_fine"] == trans_inp["ngkpt"]

        self.abivalidate_input(trans_inp)

        #####################
        # Non-linear methods
        ####################
        #if self.has_abinit(version='8.3.2'):
        dte_inputs = gs_inp.make_dte_inputs(phonon_pert=True, skip_permutations=True, ixc=3)
        #print("dte inputs\n", dte_inputs)
        assert len(dte_inputs) == 8
        assert np.all(dte_inputs[0]["d3e_pert2_dir"] == [1, 0, 0])
        assert np.all(dte_inputs[3]["d3e_pert1_atpol"] == [2, 2])
        assert all(np.all(inp["optdriver"] == 5 for inp in dte_inputs))

        # Validate with Abinit
        self.abivalidate_multi(dte_inputs)

    def test_input_check_sum(self):
        """Testing the hash method of AbinitInput"""
        inp = ebands_input(abidata.cif_file("si.cif"), abidata.pseudos("14si.pspnc"), kppa=10, ecut=2)[0]
        inp_cs = inp.variable_checksum()
        ecut = inp.pop('ecut')
        inp.set_vars({'ecut': ecut})
        assert inp_cs == inp.variable_checksum()

    def test_explicit_occ(self):
        """Testing helper function to set occupancies when occopt == 2"""

        abi_kwargs = dict(ecut=15, pawecutdg=30, tolvrs=1e-12, chksymbreak=0)

        inp = AbinitInput(structure=abidata.cif_file("SrO_Eu_222.cif"),
                          pseudos=abidata.pseudos("Sr.xml", "o.paw", "Eu.xml"),
                          #pseudos=abidata.pseudos("Sr.xml", "O.xml", "Eu.xml"))
                          abi_kwargs=abi_kwargs)

        assert inp.ispaw and len(inp.structure) == 16
        #abivars = inp.to_abivars()
        symb2luj = {"Eu": {"lpawu": 3, "upawu": 7, "jpawu": 0.7}}
        usepawu = 1
        inp.set_usepawu(usepawu, symb2luj, units="eV")

        assert inp["usepawu"] == usepawu
        assert inp["lpawu"] == [-1, 3, -1]
        assert inp["upawu"] == [0, 7, 0, "eV"]
        assert inp["jpawu"] == [0, 0.7, 0, "eV"]

        symb2spinat = {"Eu": [0, 0, 7]}
        inp.set_spinat_from_symbols(symb2spinat, default=(0, 0, 0))
        assert inp["spinat"] == "21*0 0 0 7 24*0 "

        nsppol = 2
        ngkpt = [2, 2, 2]
        shiftk = [0.5, 0.5, 0.5]
        n_cond = round(20)
        n_val = inp.num_valence_electrons

        spin_up_gs = f"\n{int((n_val - 7) / 2)}*1 7*1 {n_cond}*0"
        spin_up_ex = f"\n{int((n_val - 7) / 2)}*1 6*1 0 1/3 1/3 1/3 {n_cond - 3}*0" # Fractional occ works better for SrO
        spin_dn = f"\n{int((n_val - 7) / 2)}*1 7*0 {n_cond}*0"

        #occ1k_spin = [spin_up_gs, spin_dn]
        occ1k_spin = [spin_up_ex, spin_dn]

        inp.set_kmesh_nband_and_occ(ngkpt, shiftk, nsppol, occ1k_spin,
                                    # nspinor=1, kptopt=1, occopt=2)
                                    )

        assert inp["kptopt"] == 1 and inp["occopt"] == 2

        # nband
        # shiftk    0.5    0.5    0.5
        # ngkpt 2 2 2
        assert inp["nsppol"] == 2 and inp["occopt"] == 2 # and inp["nspden"] == 2

        # occ
        assert inp["nband"] == "*91"
        assert inp["occ"] == """
64*1 6*1 0 1/3 1/3 1/3 17*0
64*1 6*1 0 1/3 1/3 1/3 17*0

64*1 7*0 20*0
64*1 7*0 20*0
"""

        # Call Abinit in dry run mode to validate input.
        self.abivalidate_input(inp)


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
        assert multi[0].structure == multi[1].structure
        assert multi[0].structure is not multi[1].structure

        multi.set_vars(ecut=2, comment="hello")
        assert all(inp["ecut"] == 2 for inp in multi)
        assert all(inp.comment == "# hello" for inp in multi)
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
        multi.write(filepath=self.tmpfileindir("run.abi"), split=True)
        multi.write(filepath=self.tmpfileindir("run.abi"), split=False)

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
        new_multi.add_tags([GROUND_STATE, RELAX], [0, 2])
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
        self.assertMSONable(inp)

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
            [[0.        , 0.184959  , 0.       , 0.],
             [0.13871925, 0.13871925, 0.       , 0.],
             [0.0924795 , 0.0924795 , 0.0924795, 0.]])
        self.abivalidate_input(inp_loto)

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
        assert anaddb_input["elaflag"] == 5 and anaddb_input["piezoflag"] == 3 and anaddb_input["asr"] == 2
        assert anaddb_input["instrflag"] == 1 and len(anaddb_input.comment) > 0
        self.abivalidate_input(anaddb_input)

        anaddb_input = AnaddbInput.piezo_elastic(self.structure, stress_correction=False, asr=0)
        assert anaddb_input["elaflag"] == 3 and anaddb_input["piezoflag"] == 3 and anaddb_input["chneut"] == 1
        assert anaddb_input["instrflag"] == 1 and anaddb_input["asr"] == 0
        self.abivalidate_input(anaddb_input)

    def test_dfpt(self):
        """Testing dfpt constructor."""
        anaddb_input = AnaddbInput.dfpt(self.structure, stress_correction=True, relaxed_ion=True, piezo=True, dde=True,
                                        strain=True, dte=False)
        assert anaddb_input["elaflag"] == 5
        assert anaddb_input["dieflag"] == 3
        assert anaddb_input["piezoflag"] == 7
        self.abivalidate_input(anaddb_input)

        anaddb_input = AnaddbInput.dfpt(self.structure, stress_correction=True, relaxed_ion=False, piezo=True, dde=True,
                                        strain=True, dte=False)
        assert anaddb_input["elaflag"] == 1
        assert anaddb_input["dieflag"] == 2
        assert anaddb_input["piezoflag"] == 1
        self.abivalidate_input(anaddb_input)

        anaddb_input = AnaddbInput.dfpt(self.structure, stress_correction=False, relaxed_ion=True, piezo=True, dde=False,
                                        strain=True, dte=False)
        assert anaddb_input["elaflag"] == 3
        assert anaddb_input["dieflag"] == 0
        assert anaddb_input["piezoflag"] == 3
        self.abivalidate_input(anaddb_input)

        ndivsm = 1
        nqsmall = 3
        ngqpt = (4, 4, 4)
        anaddb_input = AnaddbInput.dfpt(self.structure, ngqpt=ngqpt, ndivsm=ndivsm, nqsmall=nqsmall, asr=0, dos_method="tetra")
        assert anaddb_input['ifcflag'] == 1
        self.abivalidate_input(anaddb_input)

        anaddb_input = AnaddbInput.dfpt(self.structure, dte=True)
        assert anaddb_input['nlflag'] == 3
        assert anaddb_input['alphon'] == 0

        anaddb_input = AnaddbInput.dfpt(self.structure, raman=True)
        assert anaddb_input['nlflag'] == 1
        assert anaddb_input['ramansr'] == 1


class TestCut3DInput(AbipyTest):
    """Unit tests for Cut3dInput."""

    def setUp(self):
        self.structure = abidata.structure_from_ucell("Si")

    def test_dict_methods(self):
        cut3d_input = Cut3DInput.den_to_cube('/path/to/den', 'outfile_name')
        repr(cut3d_input); str(cut3d_input)
        cut3d_input.write(self.get_tmpname(text=True))

        self.serialize_with_pickle(cut3d_input, test_eq=False)
        self.assertMSONable(cut3d_input)

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

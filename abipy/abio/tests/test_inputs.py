"""Tests for input  module"""
from __future__ import print_function, division, unicode_literals

import numpy as np
import abipy.data as abidata

from abipy.core.testing import *
from abipy.abio.inputs import *


class TestAbinitInput(AbipyTest):
    """Unit tests for AbinitInput."""

    def test_api(self):
        """Test AbinitInput API."""
        # Build empty input
        inp = AbinitInput(pseudos=abidata.pseudos("14si.pspnc"))

        print(repr(inp))
        assert len(inp) == 0 and not inp 
        assert inp.get("foo", "bar") == "bar" and inp.pop("foo", "bar") == "bar"
        assert inp.comment is None
        inp.set_comment("This is a comment")
        assert inp.comment == "This is a comment"
        assert inp.isnc and not inp.ispaw
        assert not inp.decorators

        # foo is not a valid Abinit variable
        with self.assertRaises(inp.Error): inp["foo"] = 1
        inp["ecut" ] = 1
        assert inp.get("ecut") == 1 and len(inp) == 1 and "ecut" in inp.keys() and "foo" not in inp

        assert inp.mnemonics == False
        inp.set_mnemonics(True)
        assert inp.mnemonics == True

        inp.set_vars(ecut=5, toldfe=1e-6)
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
        inp.set_vars(unit_cell)

        # Now we have a structure
        assert len(inp.structure) == 2 and inp.num_valence_electrons == 8

        inp_copy = inp.deepcopy()
        inp_copy["typat"][1] = 2
        assert inp["typat"] == [1, 1]

        # Test set_structure 
        new_structure = inp.structure.copy() 
        new_structure.perturb(distance=0.1)
        inp.set_structure(new_structure)
        assert inp.structure == new_structure

        # Compatible with deepcopy, Pickle and PMGSONable?
        self.serialize_with_pickle(inp, test_eq=False)
        #self.assertPMGSONable(inp)

    def test_helper_functions(self):
        """Testing AbinitInput helper functions."""
        inp = AbinitInput(pseudos=abidata.pseudos("14si.pspnc"), structure=abidata.cif_file("si.cif"))

        inp.set_kmesh(ngkpt=(1, 2, 3), shiftk=(1, 2, 3, 4, 5, 6))
        assert inp["kptopt"] == 1 and inp["nshiftk"] == 2

        inp.set_autokmesh(nksmall=2)
        assert inp["kptopt"] == 1 and np.all(inp["ngkpt"] == [2, 2, 2]) and inp["nshiftk"] == 4

        inp.set_kpath(ndivsm=3, kptbounds=None)
        assert inp["iscf"] == -2 and len(inp["kptbounds"]) == 12

        inp.set_kptgw(kptgw=(1, 2, 3, 4, 5, 6), bdgw=(1, 2))
        assert inp["nkptgw"] == 2 and np.all(inp["bdgw"].ravel() == np.array(len(inp["kptgw"]) * [1,2]).ravel())

    def test_abinit_calls(self):
        """Testing AbinitInput methods invoking Abinit."""
        inp = AbinitInput(pseudos=abidata.pseudos("14si.pspnc"))
        inp.set_structure(abidata.cif_file("si.cif"))

        inp.set_kmesh(ngkpt=(2, 2, 2), shiftk=(0, 0, 0))

        # Test validate with wrong input
        inp.set_vars(ecut=-1)
        v = inp.validate()
        assert v.retcode != 0 and v.log_file.read()

        # Test validate with correct input
        inp.set_vars(ecut=2, toldfe=1e-6)
        v = inp.validate()
        assert v.retcode == 0
                                                    
        # Test get_ibz
        ibz = inp.get_ibz()
        assert np.all(ibz.points == [[ 0. ,  0. ,  0. ], [ 0.5,  0. ,  0. ], [ 0.5,  0.5,  0. ]])
        assert np.all(ibz.weights == [0.125,  0.5,  0.375])

        # [{'idir': 1, 'ipert': 1, 'qpt': [0.0, 0.0, 0.0]}]
        irred_perts = inp.get_irred_phperts(qpt=(0, 0, 0))
        assert len(irred_perts) == 1
        pert = irred_perts[0]
        assert pert.idir == 1 and (pert.idir, pert.ipert) == (1, 1) and all(c == 0 for c in pert.qpt)

        #assert 0


if __name__ == "__main__":
    import unittest
    unittest.main()

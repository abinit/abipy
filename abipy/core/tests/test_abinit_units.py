"""Tests for core.abinit_units"""
import abipy.core.abinit_units as abu

from abipy.core.testing import AbipyTest

#from abipy.dfpt.phonons import factor_ev2units, unit_tag, dos_label_from_units, wlabel_from_units


class TestUnitTools(AbipyTest):

    def test_units_api(self):
        for units in ["ev", "meV" ,"ha", "cm-1", "cm^-1", "Thz"]:
            assert abu.phfactor_ev2units(units)
            assert abu.phunit_tag(units)
            assert abu.phdos_label_from_units(units)
            assert abu.wlabel_from_units(units)

        for func in [abu.phfactor_ev2units, abu.phunit_tag, abu.phdos_label_from_units]:
            with self.assertRaises(KeyError):
                func("foo")

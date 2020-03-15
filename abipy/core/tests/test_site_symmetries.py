"""Tests for core.site_symmetries module"""
import os
import numpy as np

from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
#from abipy.core.site_symmetries import SiteSymmetries

import abipy.data as abidata


class TestSiteSymmetries(AbipyTest):
    """Unit tests for SiteSymmetries."""

    def test_si(self):
        """Testing wyckoff positions for Si2"""
        si = Structure.from_file(abidata.cif_file("si.cif"))
        ss = si.site_symmetries
        repr(ss); str(ss)
        assert ss.to_string(verbose=2)
        df = ss.get_wyckoff_dataframe(verbose=2)
        self.assert_equal(np.array(df["xfrac"].values, dtype=float), [0, 0.25])
        self.assert_equal(np.array(df["yfrac"].values, dtype=float), [0, 0.25])
        self.assert_equal(np.array(df["zfrac"].values, dtype=float), [0, 0.25])
        #0  -43m (#31) nsym:24                  0                  0                  0
        #1  -43m (#31) nsym:24  0.250000000000000  0.250000000000000  0.250000000000000

        df = ss.get_tensor_rank2_dataframe(verbose=2)
        ref = ["Tzz", "Tzz"]
        self.assert_equal(df["Txx"].values, ref)
        self.assert_equal(df["Tyy"].values, ref)
        ref = ["-Tzz/3", "-Tzz/3"]
        self.assert_equal(df["Txy"].values, ref)
        self.assert_equal(df["Txz"].values, ref)
        self.assert_equal(df["Tyz"].values, ref)

    def test_alpha_sio2(self):
        """Testing wyckoff positions for alpha-SiO2"""
        asi02 = Structure.from_file(os.path.join(abidata.dirpath, "refs", "mp-7000_DDB.bz2"))
        ss = asi02.site_symmetries
        df = ss.get_wyckoff_dataframe(verbose=2)
        df = df[df["element"] == "Si"]
        self.assert_equal(df["xfrac"].values, ["xfrac", "yfrac", "0.0"])
        self.assert_equal(df["yfrac"].values, ["0.0", "yfrac", "yfrac"])
        self.assert_equal(np.array(df["zfrac"].values, dtype=float), [0.833335, 0.5, 0.166665])

        """
                 wyckoff  site_symmetry    Txx  Tyy  Tzz    Txy    Txz  Tyz
        0      3a  2 (#3,nsym:2)           Txx  Tyy  Tzz  Tyy/2  Tyz/2  Tyz
        1      3a  2 (#3,nsym:2)           Tyy  Tyy  Tzz    Txy   -Tyz  Tyz
        2      3a  2 (#3,nsym:2)         2*Txy  Tyy  Tzz    Txy  2*Tyz  Tyz
        """
        df = ss.get_tensor_rank2_dataframe(verbose=2)
        ref = ["Txx", "Tyy", "2*Txy"]
        self.assert_equal(df["Txx"].values, ref)
        ref = ["Tyy", "Tyy", "Tyy"]
        self.assert_equal(df["Tyy"].values, ref)
        ref = ["Tzz", "Tzz", "Tzz"]
        self.assert_equal(df["Tzz"].values, ref)
        ref = ["Tyy/2", "Txy", "Txy"]
        self.assert_equal(df["Txy"].values, ref)
        ref = ["Tyz/2", "-Tyz", "2*Tyz"]
        self.assert_equal(df["Txz"].values, ref)
        ref = ["Tyz", "Tyz", "Tyz"]
        self.assert_equal(df["Tyz"].values, ref)

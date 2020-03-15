"""Tests for core.restapi module"""
import contextlib
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.core import restapi


class TestMpRestApi(AbipyTest):
    """Test interfaces with the Materials Project REST API."""

    def test_mprester(self):
        """Testing MP Rest API wrappers."""

        with restapi.get_mprester() as rest:
            print(rest.Error)
            #pdr = rest.get_phasediagram_results("Li-Fe-O".split("-"))
            pdr = rest.get_phasediagram_results("Ga-Ge-As".split("-"))  # This one is faster.
            repr(pdr); str(pdr)
            pdr.print_dataframes(verbose=2)
            if self.has_matplotlib():
                plotter = pdr.plot(show_unstable=True, show=False)
                assert hasattr(plotter, "show")

        # Test mp_search
        mp = abilab.mp_search("MgB2")
        repr(mp); str(mp)
        assert mp.structures
        assert "mp-763" in mp.ids
        assert len(mp.structures) == len(mp.data)
        assert hasattr(mp.dataframe, "describe")
        mp.print_results(fmt="abivars", verbose=2)
        new = mp.add_entry(mp.structures[-1], "newid")
        assert len(new.ids) == len(mp.ids) + 1
        assert new.ids == mp.ids + ["newid"]

        with contextlib.redirect_stdout(None):
            new.print_results(fmt="cif", verbose=2)

        # Test mp_match_structure
        mp = abilab.mp_match_structure(abidata.cif_file("al.cif"))
        repr(mp); str(mp)
        assert mp.structures and mp
        assert "mp-134" in mp.ids
        assert mp.data is None and mp.dataframe is None
        mp.print_results(fmt="abivars", verbose=2)

        if self.has_nbformat():
            mp.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_cod(self):
        """Testing COD interface."""
        self.skip_if_not_executable("mysql")
        # Test abilab.cod_search
        cod = abilab.cod_search("MgB2", primitive=True)
        repr(cod); str(cod)
        assert cod.structures and cod
        assert 1000026 in cod.ids
        assert cod.data is not None
        assert hasattr(cod.dataframe, "describe")

        with contextlib.redirect_stdout(None):
            cod.print_results(fmt="POSCAR", verbose=2)

"""Tests for pseudos_models module."""

from pydantic import ValidationError
from abipy.core.testing import AbipyTest
from abipy.htc.pseudos_models import PseudoSpecs


class TestPseudosSpecs(AbipyTest):

    def test_api(self):
        """Testing PseudoSpecs API."""
        # FIXME: This requires an installed repo
        repo_name = "ONCVPSP-PBEsol-SR-PDv0.4"
        specs = PseudoSpecs.from_repo_table_name(repo_name, "standard")
        assert specs.repo_name == repo_name
        assert specs.table_name == "standard"
        assert specs.ps_generator == "ONCVPSP"
        assert specs.ps_type == "NC"
        assert specs.xc_name == "PBEsol"
        assert specs.relativity_type == "SR"
        assert specs.project_name == "PD"
        assert specs.version == "0.4"

        data = specs.dict()

        with self.assertRaises(ValidationError):
            d = data.copy()
            d["relativity_type"] = "foobar"
            PseudoSpecs(**d)

        with self.assertRaises(KeyError):
            PseudoSpecs.from_repo_table_name("ONCVPSP-PBEsol-SRfoo-PDv0.4", "standard")

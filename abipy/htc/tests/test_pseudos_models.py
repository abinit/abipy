from pydantic import ValidationError
from abipy.core.testing import AbipyTest
from abipy.htc.pseudos_models import PseudoSpecs


class TestPseudosSpecs(AbipyTest):

    def test_api(self):
        repo_name = "ONCVPSP-PBEsol-SR-PDv0.4"
        specs = PseudoSpecs.from_repo_name(repo_name)
        assert specs.repo_name == repo_name
        assert specs.ps_generator == "ONCVPSP"
        assert specs.xc_name == "PBEsol"
        assert specs.relativity_type == "SR"
        assert specs.project_name == "PD"
        assert specs.version == "0.4"
        assert specs.table_accuracy == "standard"

        data = specs.dict()

        with self.assertRaises(ValidationError):
            d = data.copy()
            d["table_accuracy"] = "foobar"
            PseudoSpecs(**d)

        with self.assertRaises(ValidationError):
            d = data.copy()
            d["relativity_type"] = "foobar"
            PseudoSpecs(**d)

        with self.assertRaises(ValueError):
            PseudoSpecs.from_repo_name("ONCVPSP-PBEsol-SRfoo-PDv0.4")

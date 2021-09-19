from abipy.core.testing import AbipyTest
from abipy.flowtk.psrepos import OncvpspRepo, get_repo_with_name


class TestPsRepos(AbipyTest):

    def test_base_api(self):
        repo = OncvpspRepo(ps_generator="ONCVPSP", xc_name="PBE", relativity_type="SR",
                           project_name="PD", version="0.1", rid=1, url="http://example.org")

        assert repo.xc_name == "PBE" and repo.version == "0.1"
        assert repo.ps_type == "NC"
        assert repr(repo)
        assert str(repo)
        assert repo.isnc and not repo.ispaw
        assert repo.dirname == "ONCVPSP-PBE-SR-PDv0.1"
        assert repo == repo

        repos_root = "/tmp"
        assert not repo.is_installed(repos_root)
        d = repo.to_rowdict(repos_root, verbose=2)
        assert d["installed"] == "False"
        assert d["version"] == "0.1"

        # For the time being, we don't test the installation procedure.
        #repo.install()

        with self.assertRaises(ValueError):
            OncvpspRepo(ps_generator="ONCVPSP", xc_name="PBE", relativity_type="Foo",
                        project_name="PD", version="0.1", rid=1, url="http://example.org")

    def test_oncvpsp_api(self):
        repo_sr = get_repo_with_name("ONCVPSP-PBEsol-SR-PDv0.4")
        repo_fr = get_repo_with_name("ONCVPSP-PBEsol-FR-PDv0.4")
        repos = [repo_sr, repo_fr]
        assert repo_sr != repo_fr
        assert repo_sr.rid != repo_fr.rid
        assert repo_sr.relativity_type == "SR"
        assert repo_fr.relativity_type == "FR"
        assert all(repo.xc_name == "PBEsol" for repo in repos)
        assert all(repo.version == "0.4" for repo in repos)
        assert all(repo.project_name == "PD" for repo in repos)
        assert all(repo.ps_generator == "ONCVPSP" for repo in repos)

        #repos_root = "/tmp"
        #if repo_sr.is_installed(repos_root)
        #    repo_sr.validate_checksums()
        #    pseudos = repo_sr.get_pseudos()

    def test_atompaw_api(self):
        repo_lda = get_repo_with_name("ATOMPAW-LDA-JTHv1.1")
        repo_pbe = get_repo_with_name("ATOMPAW-PBE-JTHv1.1")
        repos = [repo_lda, repo_pbe]
        assert repo_lda != repo_pbe
        assert repo_lda.xc_name == "LDA"
        assert repo_pbe.xc_name == "PBE"
        assert all(repo.version == "1.1" for repo in repos)
        assert all(repo.ps_generator == "ATOMPAW" for repo in repos)
        assert all(repo.project_name == "JTH" for repo in repos)
        assert all(repo.ispaw for repo in repos)
        assert all(not repo.isnc for repo in repos)

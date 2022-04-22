"""
Tests for psrepos.py module
"""
import tempfile
import os

from abipy.core.testing import AbipyTest
from abipy.flowtk.psrepos import (OncvpspRepo, get_repo_from_name, encode_pseudopath, decode_pseudopath,
                                  download_repo_from_url, tabulate_repos, get_all_registered_repos,
                                  get_installed_repos_and_root)


class TestPsRepos(AbipyTest):

    def test_encode_decode_filepath(self):
        repo = get_repo_from_name("ONCVPSP-PBEsol-SR-PDv0.4")
        if not repo.is_installed():
            raise self.SkipTest("ONCVPSP-PBEsol-SR-PDv0.4 should be installed")
        filepath = os.path.join(repo.dirpath, "Ni/Ni-sp.psp8")
        encoded = encode_pseudopath(filepath)
        assert encoded == "@ONCVPSP-PBEsol-SR-PDv0.4/Ni/Ni-sp.psp8"
        assert decode_pseudopath(encoded) == filepath

        #filepath_with_tilde = "~/.abinit/pseudos/ONCVPSP-PBEsol-SR-PDv0.4/Ni/Ni-sp.psp8"
        #encoded = encode_pseudopath(filepath_with_tilde)
        #assert encoded == "@ONCVPSP-PBEsol-SR-PDv0.4/Ni/Ni-sp.psp8"

        # Path of pseudos in directories that are not known to AbiPy should be left unchanged.
        filepath = "/Users/gmatteo/MyTables/ONCVPSP-PBEsol-SR-PDv0.4/Ni/Ni-sp.psp8"
        encoded = encode_pseudopath(filepath)
        assert encoded == filepath
        assert decode_pseudopath(encoded) == filepath

    def test_download_repo(self):
        """Testing download_repo_from_url"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            url = "https://file-examples-com.github.io/uploads/2017/02/zip_2MB.zip"
            download_repo_from_url(url, tmp_dir, chunk_size=2 * 1024**2, verbose=1)

    def test_base_api(self):
        """Testing base API."""
        repo = OncvpspRepo(ps_generator="ONCVPSP", xc_name="PBE", relativity_type="SR",
                           project_name="PD", version="0.1", url="http://example.org")

        assert repo.xc_name == "PBE" and repo.version == "0.1"
        assert repo.ps_type == "NC"
        assert repr(repo)
        assert str(repo)
        assert repo.isnc and not repo.ispaw
        assert repo.name == "ONCVPSP-PBE-SR-PDv0.1"
        assert repo == repo

        assert not repo.is_installed()
        d = repo.to_rowdict(verbose=2)
        assert d["installed"] is False
        assert d["version"] == "0.1"

        # For the time being, we don't test the installation procedure.
        #repo.install()

        with self.assertRaises(ValueError):
            OncvpspRepo(ps_generator="ONCVPSP", xc_name="PBE", relativity_type="Foo",
                        project_name="PD", version="0.1", url="http://example.org")

    def test_high_level_api(self):

        def check_url(url):
            import requests
            response = requests.head(url)
            # 302 corresponds to redirection and it's returned by github.
            if response.status_code not in (200, 302):
                raise RuntimeError(f'url {url} returned {response.status_code}')

        repos = get_all_registered_repos()
        assert repos
        for repo in repos:
            check_url(repo.url)

        installed_repos, root = get_installed_repos_and_root()
        assert os.path.exists(root)

    def test_oncvpsp_api(self):
        """Testing ONCVPSP API."""
        repo_sr = get_repo_from_name("ONCVPSP-PBEsol-SR-PDv0.4")
        repo_fr = get_repo_from_name("ONCVPSP-PBEsol-FR-PDv0.4")
        repos = [repo_sr, repo_fr]
        assert repo_sr != repo_fr
        assert repo_sr.relativity_type == "SR"
        assert repo_fr.relativity_type == "FR"
        assert all(repo.isnc and not repo.ispaw for repo in repos)
        assert all(repo.xc_name == "PBEsol" for repo in repos)
        assert all(repo.version == "0.4" for repo in repos)
        assert all(repo.project_name == "PD" for repo in repos)
        assert all(repo.ps_generator == "ONCVPSP" for repo in repos)

        assert tabulate_repos(repos, with_citations=True, verbose=1)
        assert tabulate_repos(repos, exclude=["installed"], verbose=1)

        if repo_sr.is_installed():
            repo_sr.validate_checksums(verbose=1)
            for table_name in ("standard", "stringent"):
                pseudos = repo_sr.get_pseudos(table_name=table_name)
                # Test memoized method.
                same_pseudos = repo_sr.get_pseudos(table_name=table_name)
                assert same_pseudos is pseudos
                zmax = 56  # Ba
                assert pseudos.is_complete(zmax=zmax)
                assert pseudos.allnc
                assert all(p.has_hints for p in pseudos)
                assert all(not p.supports_soc for p in pseudos)

        if repo_fr.is_installed():
            repo_fr.validate_checksums(verbose=1)
            for table_name in ("standard", "stringent"):
                pseudos = repo_fr.get_pseudos(table_name=table_name)
                # Test memoized method.
                same_pseudos = repo_fr.get_pseudos(table_name=table_name)
                assert same_pseudos is pseudos
                zmax = 56  # Ba
                assert pseudos.is_complete(zmax=zmax)
                assert pseudos.allnc
                assert all(p.has_hints for p in pseudos)
                assert all(p.supports_soc for p in pseudos)

        with self.assertRaises(KeyError):
            get_repo_from_name("ONCVPSP-PBEEsol-SR-PDv0.4")

    def test_atompaw_api(self):
        """Testing ATOMPAW API."""
        repo_lda = get_repo_from_name("ATOMPAW-LDA-JTHv1.1")
        repo_pbe = get_repo_from_name("ATOMPAW-PBE-JTHv1.1")
        repos = [repo_lda, repo_pbe]
        assert repo_lda != repo_pbe
        assert repo_lda.xc_name == "LDA"
        assert repo_pbe.xc_name == "PBE"
        assert all(repo.version == "1.1" for repo in repos)
        assert all(repo.ps_generator == "ATOMPAW" for repo in repos)
        assert all(repo.project_name == "JTH" for repo in repos)
        assert all(repo.ispaw and not repo.isnc for repo in repos)

        # TODO: Add additional tests once the JTH table fulfills the repo protocol.

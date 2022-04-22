"""
This module provides an API to deal with pseudopotential repositories.
A repo is essentially a collection of versioned pseudopotentials files installed within the same root directory.
A repo has a unique name that encodes the XC functional, relativity type, the kind of pseudopotential and the version.
The default root is ~/.abinit/pseudos although it is possible to change it via the ABIPY_PSREPOS_ROOT env variable

Note that all pseudos in a repo share the same XC functional, the type (NC, PAW) and the
treatment of relativistic corrections although one might have multiple pseudos for the same element.

Due to this ambiguity, a repo cannot be directly used for running calculations in an automatic fashion
hence the user is supposed to specify both the `repo_name` and the `table_name` when constructing a `PseudoTable`.
"""
from __future__ import annotations

import os
import time
import abc
import json
import posixpath
import tempfile
import shutil
import hashlib
import requests

from typing import List, Optional, Tuple
from urllib.parse import urlsplit
from monty.termcolor import cprint, colored
from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
from abipy.tools.decorators import memoized_method
from tqdm import tqdm


REPOS_ROOT = os.environ.get("ABIPY_PSREPOS_ROOT",
                            default=os.path.join(os.path.expanduser("~"), ".abinit", "pseudos"))


def get_repos_root() -> str:
    """
    Return the path to the installation directory. Create directory if needed.
    """
    if not os.path.exists(REPOS_ROOT):
        os.mkdir(REPOS_ROOT)

    return REPOS_ROOT


def encode_pseudopath(filepath: str) -> str:
    for repo in _ALL_REPOS:
        if filepath.startswith(repo.dirpath):
            return filepath.replace(repo.dirpath, f"@{repo.name}", 1)
    return filepath


def decode_pseudopath(filepath: str) -> str:
    for repo in _ALL_REPOS:
        start = f"@{repo.name}"
        if filepath.startswith(start):
            return filepath.replace(start, repo.dirpath, 1)
    return filepath


def download_repo_from_url(url: str, save_dirpath: str,
                           chunk_size: int = 2 * 1024**2, verbose: int = 0) -> None:
    """
    Dowload file from url.

    Args:
        url: The url from which the targz is taken.
        save_dirpath: The directory in which the tarball is unpacked.
        chunk_size: Chunk size used for downloading the file.
        verbose: Verbosity level
    """
    path = urlsplit(url).path
    filename = posixpath.basename(path)
    #print(path, filename)

    # stream = True is required by the iter_content below
    with requests.get(url, stream=True) as r:
        #tmp_dir = tempfile.mkdtemp()
        with tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None) as tmp_dir:
            tmp_filepath = os.path.join(tmp_dir, filename)
            if verbose:
                print("Writing temporary file:", tmp_filepath)

            total_size_in_bytes = int(r.headers.get('content-length', 0))
            progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)

            with open(tmp_filepath, 'wb') as fd:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    fd.write(chunk)
                    progress_bar.update(len(chunk))

            progress_bar.close()
            if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
                raise RuntimeError(f"Something went wrong while donwloading url: {url}")

            shutil.unpack_archive(tmp_filepath, extract_dir=tmp_dir)

            dirpaths = [os.path.join(tmp_dir, basen) for basen in os.listdir(tmp_dir) if basen != filename]

            if len(dirpaths) != 1:
                raise RuntimeError(f"Expecting single directory, got {dirpaths}")
            if not os.path.isdir(dirpaths[0]):
                raise RuntimeError(f"Expecting single directory, got {dirpaths}")

            if verbose: print(f"Moving {dirpaths[0]} to {save_dirpath}")
            shutil.move(dirpaths[0], save_dirpath)


def md5_for_filepath(filepath: str) -> str:
    """
    Compute and return the md5 of a file.
    """
    with open(filepath, "rt") as fh:
        text = fh.read()
        m = hashlib.md5(text.encode("utf-8"))
        return m.hexdigest()


def get_repo_from_name(repo_name: str) -> PseudosRepo:
    """
    Return a PseudosRepo from its name ``repo_name``.
    Raises KeyError if ``repo_name`` is not registered.
    """
    for repo in _ALL_REPOS:
        if repo.name == repo_name:
            return repo
    else:
        all_names = [repo.name for repo in _ALL_REPOS]
        raise KeyError(f"Couldn't find {repo_name} in the list of registered repos:\n{all_names}")


def get_installed_repos_and_root(dirpath: Optional[str] = None) -> Tuple[List[PseudosRepo], str]:
    """
    Return (all_repos, dirpath)
    """
    dirpath = REPOS_ROOT if not dirpath else dirpath
    dir_basenames = [name for name in os.listdir(dirpath) if os.path.isdir(os.path.join(dirpath, name))]
    dirname2repo = {repo.name: repo for repo in _ALL_REPOS}
    return [dirname2repo[dirname] for dirname in dir_basenames if dirname in dirname2repo], dirpath


class Citation:

    def __init__(self, title: str, doi: str):
        self.title = title
        self.doi = doi

    def __str__(self):
        return f"{self.title}\ndoi:{self.doi}"

    #def __hash__(self) -> int:
    #    return hash(self.doi)

    #def __eq__(self, other):
    #    return self.doi == other.doi


class PseudosRepo(abc.ABC):
    """
    Base abstract class for a Pseudopotential repository that is essentially a collection.
    of pseudopotential files with some metadata and some external files that allow us to
    construct a PseudoTable object.
    """

    #accuracies = ["standard", "stringent"]

    def __init__(self, ps_generator: str, xc_name: str, relativity_type: str, project_name: str,
                 version: str, url: str):
        """
        Args:
            ps_generator: Name of the pseudopotential generator
            xc_name: XC functional.
            relativity_type: SR for scalar-relativistic or FR for fully relativistic.
            project_name: Name of the project associated to this repository.
            version: Version string.
            url: URL from which the targz will be taken.
        """
        if relativity_type not in {"SR", "FR"}:
            raise ValueError(f"Invalid relativity_type: {relativity_type}")

        self.ps_generator = ps_generator
        self.xc_name = xc_name
        self.version = version
        self.project_name = project_name
        self.relativity_type = relativity_type
        self.url = url

    def to_rowdict(self, verbose: int = 0) -> dict:
        row = dict(
            ps_generator=self.ps_generator,
            ps_type=self.ps_type,
            xc_name=self.xc_name,
            relativity_type=self.relativity_type,
            project_name=self.project_name,
            version=self.version,
            installed=self.is_installed(),
            name=self.name,
        )

        if verbose:
            row.update(url=self.url)

        return row

    def __repr__(self) -> str:
        return self.name

    def __str__(self) -> str:
        lines = [self.name]
        app = lines.append
        return "\n".join(lines)

    def __eq__(self, other):
        return self.name == other.name

    @property
    def isnc(self) -> bool:
        """True if norm-conserving repo."""
        return self.ps_type == "NC"

    @property
    def ispaw(self) -> bool:
        """True if PAW repo."""
        return self.ps_type == "PAW"

    #@property
    #def all_table_names(self) -> List[str]:

    @property
    def dirpath(self) -> str:
        """Absolute path to the directory with the pseudopotentials."""
        return os.path.join(REPOS_ROOT, self.name)

    def is_installed(self) -> bool:
        """True if the repo is already installed in REPOS_ROOT."""
        return os.path.exists(os.path.join(REPOS_ROOT, self.name))

    def install(self, verbose: int = 0) -> None:
        """
        Install the repository in the standard location relative to the `REPOS_ROOT` directory.
        """
        print(f"Downloading repository from: {self.url} ...")
        print(f"Installing {repr(self)} in: {self.dirpath}")
        start = time.time()
        if not os.path.exists(REPOS_ROOT): os.mkdir(REPOS_ROOT)
        download_repo_from_url(self.url, self.dirpath, verbose=verbose)
        self.validate_checksums(verbose)
        print(f"Installation completed successfully in {time.time() - start:.2f} [s]")

    #####################
    # Abstract interface.
    #####################

    @abc.abstractmethod
    def validate_checksums(self, verbose: int) -> None:
        """Validate md5 checksums after download."""

    @property
    @abc.abstractmethod
    def ps_type(self) -> str:
        """Pseudopotential type e.g. NC or PAW"""

    @property
    @abc.abstractmethod
    def name(self):
        """The name of repository built from the metadata. Must be unique"""

    @abc.abstractmethod
    def get_citations(self) -> List[Citation]:
        """Return list of citations for this repository."""

    @abc.abstractmethod
    def get_pseudos(self, table_name: str) -> PseudoTable:
        """
        Build and return the PseudoTable associated to the given table_name
        """


class OncvpspRepo(PseudosRepo):

    @classmethod
    def from_github(cls, xc_name: str, relativity_type: str, version: str) -> OncvpspRepo:
        """
        Build a OncvpsRepo assuming a github repository.
        """
        ps_generator, project_name = "ONCVPSP", "PD"

        if relativity_type == "FR":
            # https://github.com/PseudoDojo/ONCVPSP-PBE-FR-PDv0.4/archive/refs/heads/master.zip
            sub_url = f"{ps_generator}-{xc_name}-FR-{project_name}v{version}"
        elif relativity_type == "SR":
            # https://github.com/PseudoDojo/ONCVPSP-PBE-PDv0.4/archive/refs/heads/master.zip
            sub_url = f"{ps_generator}-{xc_name}-{project_name}v{version}"
        else:
            raise ValueError(f"Invalid relativity_type {relativity_type}")

        url = f"https://github.com/PseudoDojo/{sub_url}/archive/refs/heads/master.zip"
        return cls(ps_generator, xc_name, relativity_type, project_name, version, url)

    @property
    def ps_type(self) -> str:
        return "NC"

    @property
    def name(self) -> str:
        # ONCVPSP-PBEsol-PDv0.4/
        # ONCVPSP-PBE-FR-PDv0.4/
        return f"{self.ps_generator}-{self.xc_name}-{self.relativity_type}-{self.project_name}v{self.version}"

    def get_citations(self) -> List[Citation]:
        return [
            Citation(
                title="The PseudoDojo: Training and grading a 85 element optimized norm-conserving pseudopotential table",
                doi="https://doi.org/10.1016/j.cpc.2018.01.012"),
        ]

    def validate_checksums(self, verbose: int) -> None:
        print(f"\nValidating md5 checksums of {repr(self)}...")
        djson_paths = [os.path.join(self.dirpath, jfile) for jfile in ("standard.djson", "stringent.djson")]

        seen = set()
        errors = []
        for djson_path in djson_paths:
            with open(djson_path, "rt") as fh:
                djson = json.load(fh)

            for symbol, d in djson["pseudos_metadata"].items():
                bname = d["basename"]
                if bname in seen: continue
                seen.add(bname)
                ref_md5 = d["md5"]
                this_path = os.path.join(self.dirpath, symbol, bname)
                this_md5 = md5_for_filepath(this_path)
                if ref_md5 != this_md5:
                    errors.append(f"Different md5 checksums for {this_path}")
                else:
                    if verbose:
                        print(f"MD5 checksum for {this_path} is OK")

        if errors:
            cprint("Checksum test: FAILED", color="red")
            errstr = "\n".join(errors)
            raise ValueError(f"Checksum test failed for the following pseudos:\n{errstr}\n"
                             f"Data is corrupted. Try to download {repr(self)} again")
        else:
            cprint("Checksum test: OK", color="green")

    @memoized_method()
    def get_pseudos(self, table_name: str) -> PseudoTable:
        """
        Build and return the PseudoTable associated to the given table_name.
        Note that we use a per-instance cache to store the results.
        """
        djson_path = os.path.join(self.dirpath, f"{table_name}.djson")
        pseudos = []
        with open(djson_path, "rt") as fh:
            djson = json.load(fh)
            for symbol, d in djson["pseudos_metadata"].items():
                bname = d["basename"]
                pseudo_path = os.path.join(self.dirpath, symbol, bname)
                #print(f"Reading pseudo from {pseudo_path}")
                # FIXME: Bug if ~ in pseudo_path
                pseudo_path = os.path.expanduser(pseudo_path)
                pseudo = Pseudo.from_file(pseudo_path)
                # Attach a fake dojo_report
                # TODO: This part should be rationalized/rewritten
                hints = d["hints"]
                dojo_report = {"hints": hints}
                pseudo.dojo_report = dojo_report
                #print(f"pseudo.filepath after {pseudo.filepath}")
                pseudos.append(pseudo)

        return PseudoTable(pseudos)


class JthRepo(PseudosRepo):

    @classmethod
    def from_abinit_website(cls, xc_name: str, relativity_type: str, version: str) -> JthRepo:
        ps_generator, project_name = "ATOMPAW", "JTH"
        # https://www.abinit.org/ATOMICDATA/JTH-LDA-atomicdata.tar.gz
        # ATOMPAW-LDA-JTHv0.4
        url = f"https://www.abinit.org/ATOMICDATA/JTH-{xc_name}-atomicdata.tar.gz"
        return cls(ps_generator, xc_name, relativity_type, project_name, version, url)

    @property
    def ps_type(self) -> str:
        return "PAW"

    @property
    def name(self) -> str:
        # ATOMPAW-LDA-JTHv0.4
        return f"{self.ps_generator}-{self.xc_name}-{self.project_name}v{self.version}"

    def validate_checksums(self, verbose: int) -> None:
        print(f"\nValidating md5 checksums of {repr(self)} ...")
        cprint("WARNING: JTH-PAW repository does not support md5 checksums!!!", color="red")

    @memoized_method()
    def get_pseudos(self, table_name: str) -> PseudoTable:
        """
        Build and return the PseudoPotential table associated to the given table_name
        """
        if table_name != "standard":
            raise ValueError(f"JTH table does not support table_name: {table_name}")

        # Read the list of pseudopotential paths (relative to dirpath)
        txt_path = os.path.join(self.dirpath, f"{table_name}.txt")
        with open(txt_path, "rt") as fh:
            relpaths = fh.readlines()
            relpaths = [l for l in relpaths if l.strip()]

        pseudos = []
        for rpath in relpaths:
            pseudo = Pseudo.from_file(os.path.join(self.dirpath, rpath))
            # TODO: Get hints
            pseudos.append(pseudo)

        return PseudoTable(pseudos)

    def get_citations(self) -> List[Citation]:
        return [
            Citation(
                title="Generation of Projector Augmented-Wave atomic data: A 71 element validated table in the XML format",
                doi="https://www.sciencedirect.com/science/article/abs/pii/S0010465513004359?via%3Dihub"),
        ]


def repos_from_names(repo_names: List[str]) -> List[PseudosRepo]:
    """
    Return list of PP Repos from a list of repo_names
    """
    id2repo = {repo.name: repo for repo in _ALL_REPOS}
    # This will fail if we have received an invalid repo_name
    return [id2repo[repo_name] for repo_name in repo_names]


def repo_from_name(repo_name: str) -> PseudosRepo:
    """
    Return PseudosRepo from its name ``repo_name``.
    """
    id2repo = {repo.name: repo for repo in _ALL_REPOS}
    # This will fail if we have received an invalid repo_name.
    return id2repo[repo_name]


def tabulate_repos(repos: List[PseudosRepo], exclude: Optional[List[str]] = None,
                   with_citations: bool = False, verbose: int = 0) -> str:

    bool2color = {True: "green", False: "red"}

    rows = []
    for i, repo in enumerate(repos):
        d = repo.to_rowdict(verbose=verbose)
        if exclude:
            d = {k: d[k] for k in d if k not in exclude}

        if i == 0:
            headers = list(d.keys())

        if "installed" in d:
            d["installed"] = colored(d["installed"], bool2color[d["installed"]])

        rows.append(list(d.values()))

    from tabulate import tabulate
    s = tabulate(rows, headers)
    if not with_citations: return s

    lines = [s, ""]

    proj2citations = {}
    for repo in repos:
        proj2citations[repo.project_name] = repo.get_citations()

    for projname, citations in proj2citations.items():
        lines.append(f"References for the {projname} project:")
        for i, citation in enumerate(citations):
            lines.append(f"\t- {citation.doi}")

    s = "\n".join(lines)

    return s


#############################################################
# Here we register the repositories and build _ALL_REPOS list
#############################################################

def get_all_registered_repos():
    return _ALL_REPOS[:]


_mk_onc = OncvpspRepo.from_github

_ONCVPSP_REPOS = [
    _mk_onc(xc_name="PBEsol", relativity_type="SR", version="0.4"),
    _mk_onc(xc_name="PBEsol", relativity_type="FR", version="0.4"),
    _mk_onc(xc_name="PBE", relativity_type="SR", version="0.4"),
    #_mk_onc(xc_name="PBE", relativity_type="FR", version="0.4"),  FIXME: checksum fails
]

_mk_jth = JthRepo.from_abinit_website

_PAW_REPOS = [
    _mk_jth(xc_name="LDA", relativity_type="SR", version="1.1"),
    _mk_jth(xc_name="PBE", relativity_type="SR", version="1.1"),
]

_ALL_REPOS = _ONCVPSP_REPOS + _PAW_REPOS

# Make sure repo name is unique.
_repo_names = [_repo.name for _repo in _ALL_REPOS]
if len(set(_repo_names)) != len(_repo_names):
    raise RuntimeError(f"Found duplicated repo_names in ALL_REPOS:\nids: {_repo_names}")

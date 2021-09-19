"""
Pseudopotential repositories
"""
from __future__ import annotations

import sys
import os
import time
import abc
import json
import posixpath
import tempfile
import shutil
import hashlib

from typing import List
from urllib.parse import urlsplit
from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable


def download_url(url: str, save_dirpath: str, chunk_size: int = 128, verbose: int = 0) -> None:
    """

    Args:
        url:
        save_dirpath:
        chunk_size:
        verbose:
    """

    import requests
    path = urlsplit(url).path
    filename = posixpath.basename(path)
    #print(path, filename)

    # stream = True is required by the iter_content below
    with requests.get(url, stream=True) as r:
        tmp_dir = tempfile.mkdtemp()
        #with tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None) as tmp_dir:
        tmp_filepath = os.path.join(tmp_dir, filename)
        if verbose: print("Writing temporary file:", tmp_filepath)

        with open(tmp_filepath, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=chunk_size):
                fd.write(chunk)

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


def pprint_rows(rows: list, out=sys.stdout, rstrip: bool = False) -> None:
    """
    Prints out a table of data, padded for alignment
    Each row must have the same number of columns.

    Args:
        out: Output stream (file-like object)
        rows: The table to print. A list of lists.
        rstrip: if true, trailing withespaces are removed from the entries.
    """
    def max_width_col(table, col_idx):
        """Get the maximum width of the given column index"""
        return max([len(r[col_idx]) for r in table])

    if rstrip:
        for row_idx, row in enumerate(rows):
            rows[row_idx] = [c.rstrip() for c in row]

    col_paddings = []
    ncols = len(rows[0])
    for i in range(ncols):
        col_paddings.append(max_width_col(rows, i))

    for row in rows:
        # left col
        out.write( row[0].ljust(col_paddings[0] + 1) )
        # rest of the cols
        for i in range(1, len(row)):
            col = row[i].rjust(col_paddings[i] + 2)
            out.write(col)
        out.write("\n")


def get_repo_with_name(repo_name) -> PseudosRepo:
    """
    Return a PseudoRepo from its name ``repo_name``.
    Raises ValueError if ``repo_name`` is not registered.
    """
    for repo in ALL_REPOS:
        if repo.dirname == repo_name:
            return repo
    else:
        raise ValueError(f"Couldn't find {repo_name} in the list of registered repos.")



class Citation:

    def __init__(self, title: str, doi: str):
        self.title = title
        self.doi = doi



class PseudosRepo(abc.ABC):
    """
    Base abstract class for Pseudopotential repositories.
    """

    def __init__(self, ps_generator: str, xc_name: str, relativity_type: str, project_name: str,
                 version: str, rid: int, url: str):
        """
        Args:
            ps_generator:
            xc_name:
            relativity_type:
            project_name:
            version:
            rid:
            url:
        """
        if relativity_type not in {"SR", "FR"}:
            raise ValueError(f"Invalid relativity_type: {relativity_type}")

        self.ps_generator = ps_generator
        self.xc_name = xc_name
        self.version = version
        self.project_name = project_name
        self.relativity_type = relativity_type
        self.rid = rid
        self.url = url

    @classmethod
    def from_dirpath(cls, dirpath: str) -> PseudosRepo:
        dirname = os.path.basename(dirpath)
        return cls.from_dirname(dirname)

    @classmethod
    def from_dirname(cls, dirname: str) -> PseudosRepo:
        """Return a PseudosRepo for the installation directory."""
        for repo in ALL_REPOS:
            if repo.dirname == dirname:
                return repo

        raise ValueError(f"Cannot find `{dirname}` in the list of registered Tables!")

    def to_rowdict(self, repos_root: str, verbose: int = 0) -> dict:
        row = dict(
            ps_generator=self.ps_generator,
            ps_type=self.ps_type,
            xc_name=self.xc_name,
            relativity_type=self.relativity_type,
            project_name=self.project_name,
            version=self.version,
            installed=str(self.is_installed(repos_root)),
            dirname=self.dirname,
            repo_id=str(self.rid),
        )

        if verbose:
            row.update(url=self.url)

        return row

    @property
    def dirname(self) -> str:
        if self.isnc:
            # ONCVPSP-PBEsol-PDv0.4/
            # ONCVPSP-PBE-FR-PDv0.4/
            return f"{self.ps_generator}-{self.xc_name}-{self.relativity_type}-{self.project_name}v{self.version}"
        elif self.ispaw:
            # ATOMPAW-LDA-JTHv0.4
            return f"{self.ps_generator}-{self.xc_name}-{self.project_name}v{self.version}"
        else:
            raise ValueError(f"Invalid ps_generator: {self.ps_generator}")

    def __repr__(self) -> str:
        return self.dirname

    def __str__(self) -> str:
        lines = [self.dirname]
        app = lines.append
        return "\n".join(lines)

    def __eq__(self, other):
        return self.dirname == other.dirname

    @property
    def isnc(self) -> bool:
        """True if norm-conserving repo."""
        return self.ps_type == "NC"

    @property
    def ispaw(self) -> bool:
        """True if PAW repo."""
        return self.ps_type == "PAW"

    def is_installed(self, repos_root: str) -> bool:
        """True if the repo is already installed in repos_root."""
        return os.path.exists(os.path.join(repos_root, self.dirname))

    def install(self, repos_root: str, verbose: int) -> None:
        """
        Install the repo in the standard location relative to the `repos_root` directory.
        """
        save_dirpath = os.path.join(repos_root, self.dirname)
        print(f"Installing {repr(self)} in {save_dirpath} directory")
        print(f"Downloading repository from: {self.url} ...")
        start = time.time()
        download_url(self.url, save_dirpath, verbose=verbose)
        self.validate_checksums(repos_root, verbose)
        print(f"Completed in {time.time() - start:.2f} [s]")

    #####################
    # Abstract interface.
    #####################

    @abc.abstractmethod
    def validate_checksums(self, repos_root: str, verbose: int) -> None:
        """Validate checksums after download."""

    @property
    @abc.abstractmethod
    def ps_type(self) -> str:
        """Pseudopotential type e.g. NC or PAW"""

    @abc.abstractmethod
    def get_citations(self) -> List[Citation]:
        """Return list of citations for this repository."""

    @abc.abstractmethod
    def get_pseudos(self, repos_root: str, table_accuracy: str) -> PseudoTable:
        """Build and return the PseudoPotential table associated to the given table_accuracy"""


class OncvpspRepo(PseudosRepo):

    @classmethod
    def from_github(cls, xc_name: str, relativity_type: str, rid: int, version: str) -> OncvpspRepo:
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
        return cls(ps_generator, xc_name, relativity_type, project_name, version, rid, url)

    @property
    def ps_type(self) -> str:
        return "NC"

    def get_citations(self) -> List[Citation]:
        return [
            Citation(
                title="The PseudoDojo: Training and grading a 85 element optimized norm-conserving pseudopotential table",
                doi="https://doi.org/10.1016/j.cpc.2018.01.012"),
        ]

    def validate_checksums(self, repos_root: str, verbose: int) -> None:
        print(f"\nValidating checksums of {repr(self)} ...")
        dirpath = os.path.join(repos_root, self.dirname)
        djson_paths = [os.path.join(dirpath, jfile) for jfile in ("standard.djson", "stringent.djson")]

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
                this_path = os.path.join(dirpath, symbol, bname)
                this_md5 = md5_for_filepath(this_path)
                if ref_md5 != this_md5:
                    errors.append(f"Different md5 checksums for {this_path}")
                else:
                    if verbose:
                        print(f"MD5 checksum for {this_path} is OK")

        if errors:
            print("Checksum test: FAILED")
            errstr = "\n".join(errors)
            raise ValueError(f"Checksum test failed for the following pseudos:\n{errstr}\n"
                             f"Data is corrupted. Try to download {repr(self)} again")
        else:
            print("Checksum test: OK")

    def get_pseudos(self, repos_root: str, table_accuracy: str) -> PseudoTable:

        dirpath = os.path.join(repos_root, self.dirname)
        djson_path = os.path.join(dirpath, f"{table_accuracy}.djson")
        pseudos = []
        with open(djson_path, "rt") as fh:
            djson = json.load(fh)
            for symbol, d in djson["pseudos_metadata"].items():
                bname = d["basename"]
                pseudo_path = os.path.join(dirpath, symbol, bname)
                pseudo = Pseudo.from_file(pseudo_path)
                # Attach a fake dojo_report
                # TODO: This part should be rationalized
                hints = d["hints"]
                dojo_report = {"hints": hints}
                pseudo.dojo_report = dojo_report
                print(pseudo)
                pseudos.append(pseudo)

        return PseudoTable(pseudos)


class JthRepo(PseudosRepo):

    @classmethod
    def from_abinit_website(cls, xc_name: str, relativity_type: str, rid: int, version: str) -> JthRepo:
        ps_generator, project_name = "ATOMPAW", "JTH"
        # https://www.abinit.org/ATOMICDATA/JTH-LDA-atomicdata.tar.gz
        # ATOMPAW-LDA-JTHv0.4
        url = f"https://www.abinit.org/ATOMICDATA/JTH-{xc_name}-atomicdata.tar.gz"
        return cls(ps_generator, xc_name, relativity_type, project_name, version, rid, url)

    @property
    def ps_type(self) -> str:
        return "PAW"

    def validate_checksums(self, repos_root: str, verbose: int) -> None:
        print(f"\nValidating checksums of {repr(self)} ...")
        print("WARNING: JTH-PAW repository does not support md5 checksums!!!")

    def get_pseudos(self, repos_root: str, table_accuracy: str) -> PseudoTable:
        if table_accuracy != "standard":
            raise ValueError(f"JTH table does not support table_accuracy: {table_accuracy}")
        dirpath = os.path.join(repos_root, self.dirname)
        pseudos = []
        raise NotImplementedError()
        #return PseudoTable(pseudos)

    def get_citations(self) -> List[Citation]:
        return [
            Citation(
                title="Generation of Projector Augmented-Wave atomic data: A 71 element validated table in the XML format",
                doi="https://www.sciencedirect.com/science/article/abs/pii/S0010465513004359?via%3Dihub"),
        ]


def repos_from_id_list(id_list: list) -> List[PseudosRepo]:
    """
    Return list of PP Repos from a list of table identifiers.
    """
    ids = sorted(set([int(i) for i in id_list]))
    id2repo = {repo.rid: repo for repo in ALL_REPOS}
    return [id2repo[rid] for rid in ids]   # This will fail if we have received an invalid id.


def pprint_repos(repos: List[PseudosRepo], repos_root: str, out=sys.stdout,
                 rstrip: bool = False, verbose: int = 0) -> None:

    rows = None
    for i, repo in enumerate(repos):
        d = repo.to_rowdict(repos_root, verbose=verbose)
        if i == 0:
            rows = [list(d.keys())]
        rows.append(list(d.values()))

    pprint_rows(rows, out=out, rstrip=rstrip)


###################################
# Here we register the repositories
# Client code should use ALL_REPOS
###################################

mk_onc = OncvpspRepo.from_github

ONCVPSP_REPOS = [
    mk_onc(xc_name="PBEsol", relativity_type="SR", version="0.4", rid=1),
    mk_onc(xc_name="PBEsol", relativity_type="FR", version="0.4", rid=2),
    mk_onc(xc_name="PBE", relativity_type="SR", version="0.4", rid=3),
    #mk_onc(xc_name="PBE", relativity_type="FR", version="0.4", rid=4),  FIXME: checksum fails
]

mk_jth = JthRepo.from_abinit_website

PAW_REPOS = [
    mk_jth(xc_name="LDA", relativity_type="SR", version="1.1", rid=21),
    mk_jth(xc_name="PBE", relativity_type="SR", version="1.1", rid=22),
]

ALL_REPOS = ONCVPSP_REPOS + PAW_REPOS

# Check for possible duplications in the IDs.
_ids = [_repo.rid for _repo in ALL_REPOS]
if len(set(_ids)) != len(_ids):
    raise RuntimeError(f"Found duplicated ids in ALL_REPOS:\nids: {_ids}")


_repo_names = [_repo.dirname for _repo in ALL_REPOS]
if len(set(_repo_names)) != len(_repo_names):
    raise RuntimeError(f"Found duplicated repo_names in ALL_REPOS:\nids: {_repo_names}")
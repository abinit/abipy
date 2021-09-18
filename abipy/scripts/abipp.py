#!/usr/bin/env python
from __future__ import annotations

import sys
import os
import argparse
import time
import abc
import warnings
import json
import posixpath
import urllib.request
import tempfile
import shutil
import hashlib

from typing import List
from urllib.parse import urlsplit


#__version__ = "0.4.0"

# https://stackoverflow.com/questions/9419162/download-returned-zip-file-from-url
try:
    import requests

    def download_url(url: str, save_dirpath: str, chunk_size: int = 128, verbose: int = 0) -> None:

        path = urlsplit(url).path
        filename = posixpath.basename(path)
        #print(path, filename)

        # stream = true is required by the iter_content below
        with requests.get(url, stream=True) as r:
            tmp_dir = tempfile.mkdtemp()
            #with tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None) as tmp_dir:
            tmp_filepath = os.path.join(tmp_dir, filename)
            if verbose:
                print("Writing temporary file:", tmp_filepath)

            with open(tmp_filepath, 'wb') as fd:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    fd.write(chunk)

            shutil.unpack_archive(tmp_filepath, extract_dir=tmp_dir)

            dirpaths = [os.path.join(tmp_dir, basen) for basen in os.listdir(tmp_dir) if basen != filename]
            if len(dirpaths) != 1:
                raise RuntimeError(f"Expecting single directory, got {dirpaths}")
            if not os.path.isdir(dirpaths[0]):
                raise RuntimeError(f"Expecting single directory, got {dirpaths}")

            if verbose:
                print(f"Moving {dirpaths[0]} to {save_dirpath}")
            shutil.move(dirpaths[0], save_dirpath)


except ImportError:
    warnings.warn("""
Cannot import requests package.
For performance and reliabilty reasons, it is highly recommended to install the package with e.g.:

    pip install requests --user
    
Continuing with urllib.request from python stdlib.
""")

    import urllib.request

    def download_url(url: str, save_dirpath: str, chunk_size: int = 128, verbose: int = 0) -> None:
        with urllib.request.urlopen(url) as dl_file:
            with open(save_dirpath, 'wb') as out_file:
                out_file.write(dl_file.read())


def user_wants_to_abort():
    """Interactive problem, return False if user entered `n` or `no`."""
    try:
        answer = input("\nDo you want to continue [Y/n]")
    except EOFError:
        return False

    return answer.lower().strip() in ["n", "no"]


def md5_for_filepath(filepath: str) -> str:
    with open(filepath, "rt") as fh:
        text = fh.read()
        m = hashlib.md5(text.encode("utf-8"))
        return m.hexdigest()


def pprint_rows(rows: list, out=sys.stdout, rstrip: bool = False) -> None:
    """
    Prints out a table of data, padded for alignment
    Each row must have the same number of columns.

    Args:
        out:
            Output stream (file-like object)
        rows:
            The table to print. A list of lists.
        rstrip:
            if true, trailing withespaces are removed from the entries.
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


class Repo(abc.ABC):

    def __init__(self, pp_generator: str, xc_name: str, relativity_type: str, project_name: str,
                 version: str, rid: int, url: str):

        if relativity_type not in {"SR", "FR"}:
            raise ValueError(f"Invalid relativity_type: {relativity_type}")

        self.pp_generator = pp_generator
        self.xc_name = xc_name
        self.version = version
        self.project_name = project_name
        self.relativity_type = relativity_type
        self.rid = rid
        self.url = url

    def to_rowdict(self, abipp_home: str, verbose: int = 0) -> dict:
        row = dict(
            pp_generator=self.pp_generator,
            pp_type=self.pp_type,
            xc_name=self.xc_name,
            relativity_type=self.relativity_type,
            project_name=self.project_name,
            version=self.version,
            installed=str(self.is_installed(abipp_home)),
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
            return f"{self.pp_generator}-{self.xc_name}-{self.relativity_type}-{self.project_name}v{self.version}"
        elif self.ispaw:
            # ATOMPAW-LDA-JTHv0.4
            return f"{self.pp_generator}-{self.xc_name}-{self.project_name}v{self.version}"
        else:
            raise ValueError(f"Invalid pp_generator: {self.pp_generator}")

    @classmethod
    def from_dirpath(cls, dirpath: str) -> Repo:
        """Return a Repo for the installation directory."""
        dirname = os.path.basename(dirpath)
        for repo in ALL_REPOS:
            if repo.dirname == dirname:
                return repo

        raise ValueError(f"Cannot find `{dirname}` in the list of registered Tables!")

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
        return self.pp_type == "NC"

    @property
    def ispaw(self) -> bool:
        """True if PAW repo."""
        return self.pp_type == "PAW"

    def is_installed(self, abipp_home: str) -> bool:
        """True if the repo is already installed in abipp_home."""
        return os.path.exists(os.path.join(abipp_home, self.dirname))

    def install(self, abipp_home: str, verbose: int) -> None:
        """
        Install the repo in the standard location relative to the `abipp_home` directory.
        """
        save_dirpath = os.path.join(abipp_home, self.dirname)
        print(f"Installing {repr(self)} in {save_dirpath} directory")
        print(f"Downloading repo from: {self.url} ...")
        start = time.time()
        download_url(self.url, save_dirpath, verbose=verbose)
        self.validate_checksums(abipp_home, verbose)
        print(f"Completed in {time.time() - start:.2f} [s]")

    @abc.abstractmethod
    def validate_checksums(self, abipp_home: str, verbose: int) -> None:
        """Validate checksums after download."""

    @property
    @abc.abstractmethod
    def pp_type(self) -> str:
        """Pseudopotential type e.g. NC or PAW"""

    #@abc.abstractmethod
    #def build_pseudotable(self, accuracy: str):



class OncvpspRepo(Repo):

    @classmethod
    def from_github(cls, xc_name: str, relativity_type: str, rid: int, version: str) -> OncvpspRepo:
        pp_generator, project_name = "ONCVPSP", "PD"

        if relativity_type == "FR":
            # https://github.com/PseudoDojo/ONCVPSP-PBE-FR-PDv0.4/archive/refs/heads/master.zip
            sub_url = f"{pp_generator}-{xc_name}-FR-{project_name}v{version}"
        elif relativity_type == "SR":
            # https://github.com/PseudoDojo/ONCVPSP-PBE-PDv0.4/archive/refs/heads/master.zip
            sub_url = f"{pp_generator}-{xc_name}-{project_name}v{version}"
        else:
            raise ValueError(f"Invalid relativity_type {relativity_type}")

        url = f"https://github.com/PseudoDojo/{sub_url}/archive/refs/heads/master.zip"
        return cls(pp_generator, xc_name, relativity_type, project_name, version, rid, url)

    @property
    def pp_type(self) -> str:
        return "NC"

    def validate_checksums(self, abipp_home: str, verbose: int) -> None:
        print(f"\nValidating checksums of {repr(self)} ...")
        dirpath = os.path.join(abipp_home, self.dirname)
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

    def build_pseudotable(self, abipp_home: str, accuracy: str):
        from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
        dirpath = os.path.join(abipp_home, self.dirname)
        djson_path = os.path.join(dirpath, f"{accuracy}.djson")
        with open(djson_path, "rt") as fh:
            djson = json.load(fh)
            for symbol, d in djson["pseudos_metadata"].items():
                bname = d["basename"]
                pseudo_path = os.path.join(dirpath, symbol, bname)
                pseudo = Pseudo.from_file(pseudo_path)

        #return PseudoTable


class JthRepo(Repo):

    @classmethod
    def from_abinit_website(cls, xc_name: str, relativity_type: str, rid: int, version: str) -> JthRepo:
        pp_generator, project_name = "ATOMPAW", "JTH"
        # https://www.abinit.org/ATOMICDATA/JTH-LDA-atomicdata.tar.gz
        # ATOMPAW-LDA-JTHv0.4
        url = f"https://www.abinit.org/ATOMICDATA/JTH-{xc_name}-atomicdata.tar.gz"
        return cls(pp_generator, xc_name, relativity_type, project_name, version, rid, url)

    @property
    def pp_type(self) -> str:
        return "PAW"

    def validate_checksums(self, abipp_home: str, verbose: int) -> None:
        print(f"\nValidating checksums of {repr(self)} ...")
        print("WARNING: JTH-PAW repo does not support md5 checksums!!!!!!!!!!")


def repos_from_id_list(id_list: list) -> List[Repo]:
    """
    Return list of PP Repos from list of table identifiers.
    """
    ids = sorted(set([int(i) for i in id_list]))
    id2repo = {repo.rid: repo for repo in ALL_REPOS}
    return [id2repo[rid] for rid in ids]   # This will fail if we have received an invalid id.


def pprint_repos(repos: List[Repo], abipp_home: str, out=sys.stdout,
                 rstrip: bool = False, verbose: int = 0) -> None:

    rows = None
    for i, repo in enumerate(repos):
        d = repo.to_rowdict(abipp_home, verbose=verbose)
        if i == 0:
            rows = [list(d.keys())]
        rows.append(list(d.values()))

    pprint_rows(rows, out=out, rstrip=rstrip)


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


def get_abipp_home(options) -> str:
    """
    Return the path to the PseudoDojo installation directory.
    Create the directory if needed.
    """
    abipp_home = options.abipp_home
    if not os.path.exists(abipp_home):
        os.mkdir(abipp_home)

    return abipp_home


def abipp_list(options):
    """
    List all installed pseudopotential repos.
    """
    abipp_home = get_abipp_home(options)
    dirpaths = [os.path.join(abipp_home, name) for name in os.listdir(abipp_home) if
                os.path.isdir(os.path.join(abipp_home, name))]

    if not dirpaths:
        print("Could not find any pseudopotential repository installed in:", abipp_home)
        return 0

    print(f"The following repositories are installed in {abipp_home}:\n")
    repos = [Repo.from_dirpath(dirpath) for dirpath in dirpaths]
    # Keep the list sorted by ID.
    repos = sorted(repos, key=lambda repo: repo.rid)
    pprint_repos(repos, abipp_home=abipp_home)

    if not options.checksums:
        return 0

    exc_list = []
    for repo in repos:
        try:
            repo.validate_checksums(abipp_home, options.verbose)
        except Exception as exc:
            exc_list.append(exc)

    if exc_list:
        print("\nList of exceptions raised by validate_checksums:")
        for exc in exc_list:
            print(exc)

    return len(exc_list)


def abipp_avail(options):
    """
    Show available repos.
    """
    print("List of available pseudopotential repositories:\n")
    abipp_home = get_abipp_home(options)
    pprint_repos(ALL_REPOS, abipp_home)


def abipp_nc_get(options):
    """
    Get NC repo. Can choose among three formats: psp8, upf2 and psml.
    By default we fetch all formats.
    """
    abipp_home = get_abipp_home(options)
    repos = [repo for repo in ALL_REPOS if repo.isnc and not repo.is_installed(abipp_home)]
    if not repos:
        print(f"All registered NC repositories are already installed in {abipp_home}. Returning")
        return 0

    print("The following NC repositories will be installed:\n")
    pprint_repos(repos, abipp_home=abipp_home)
    if not options.yes and user_wants_to_abort(): return 2

    print("Fetching NC repositories. It may take some time ...")
    for repo in repos:
        repo.install(abipp_home, options.verbose)

    abipp_list(options)
    return 0


def abipp_paw_get(options):
    """
    Get PAW repositories in PAWXML format.
    """
    abipp_home = get_abipp_home(options)
    repos = [repo for repo in ALL_REPOS if repo.ispaw and not repo.is_installed(abipp_home)]
    if not repos:
        print(f"All registered PAW repositories are already installed in {abipp_home}. Returning")
        return 0

    print("The following PAW repositories will be installed:")
    pprint_repos(repos, abipp_home=abipp_home)
    if not options.yes and user_wants_to_abort(): return 2

    print("Fetching PAW repositories. It may take some time ...")
    for repo in repos:
        repo.install(abipp_home, options.verbose)

    abipp_list(options)
    return 0


def abipp_get_byid(options):
    """
    Get list of repos by their IDs.
    Use the `avail` command to get the repo ID.
    """
    abipp_home = get_abipp_home(options)
    repos = repos_from_id_list(options.id_list)
    repos = [repo for repo in repos if not repo.is_installed(abipp_home)]

    if not repos:
        print("Tables are already installed!")
        abipp_list(options)
        return 1

    print("The following repositories will be installed:")
    pprint_repos(repos, abipp_home=abipp_home)
    if not options.yes and user_wants_to_abort(): return 2

    for repo in repos:
        repo.install(abipp_home, options.verbose)

    abipp_list(options)
    return 0


#def abipp_apropos(options):


def get_epilog():
    return """\

Usage example:

  abipp.py avail                          --> Show registered repositories and the associated IDs.
  abipp.py list                           --> List installed repositories.
  abipp.py get_byid 1 3                   --> Download repositories by ID(s).          
  abipp.py nc_get                         --> Get all NC repositories (most recent version)
  abipp.py nc_get -xc PBE -fr -sr --version 0.4 
  abipp.py paw_get                        --> Get all PAW repositories (most recent version)
"""


def get_parser(with_epilog=False):

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                              help='verbose, can be supplied multiple times to increase verbosity.')

    copts_parser.add_argument('--abipp-home', type=str,
                              default=os.path.expanduser(os.path.join("~", ".abinit", "pseudos")),
                              help='Installation directory. Default: $HOME/.abinit/pseudos')

    copts_parser.add_argument('-y', "--yes", action="store_true", default=False,
                              help="Do not ask for confirmation when installing repositories.")

    copts_parser.add_argument("-c", "--checksums", action="store_true", default=False,
                              help="Validate checksums")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version')  #, version=__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for list command.
    p_list = subparsers.add_parser("list", parents=[copts_parser], help=abipp_list.__doc__)

    # Subparser for avail command.
    subparsers.add_parser("avail", parents=[copts_parser], help=abipp_avail.__doc__)

    # Subparser for nc_get command.
    p_nc_get = subparsers.add_parser("nc_get", parents=[copts_parser], help=abipp_nc_get.__doc__)

    # Subparser for paw_get command.
    p_paw_get = subparsers.add_parser("paw_get", parents=[copts_parser], help=abipp_paw_get.__doc__)

    # Subparser for get_byid command.
    p_get_byid = subparsers.add_parser("get_byid", parents=[copts_parser], help=abipp_get_byid.__doc__)
    p_get_byid.add_argument("id_list", type=int, nargs="+", help="List of PseudoPotential Repo IDs to download.")

    return parser


def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    return globals()[f"abipp_{options.command}"](options)


if __name__ == "__main__":
    sys.exit(main())

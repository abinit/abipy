#!/usr/bin/env python
"""
Script to download and install pseudopotential tables from the web
"""
from __future__ import annotations

import sys
#import os
import argparse

from abipy.core.release import __version__
from abipy.flowtk.psrepos import (tabulate_repos, repos_from_names,
                                  get_all_registered_repos, get_installed_repos_and_root)


def user_wants_to_abort():
    """Interactive problem, return False if user entered `n` or `no`."""
    try:
        answer = input("\nDo you want to continue [Y/n]")
    except EOFError:
        return False

    return answer.lower().strip() in ["n", "no"]


def abips_list(options):
    """
    List all installed pseudopotential repos.
    """
    repos, repos_root = get_installed_repos_and_root()

    if not repos:
        print("Could not find any pseudopotential repository installed in:", repos_root)
        return 0

    print(f"The following pseudopotential repositories are installed in {repos_root}:\n")
    print(tabulate_repos(repos, exclude=["installed"], verbose=options.verbose))

    #if options.verbose:
    #    for repo in repos:
    #        if repo.ispaw: continue
    #        pseudos = repo.get_pseudos(table_name="standard")
    #        print(pseudos)
    #else:
    #    print("\nUse -v to print the pseudos")

    if not options.checksums:
        print("\nUse -c to validate the md5 checksum")
        return 0

    exc_list = []
    for repo in repos:
        try:
            repo.validate_checksums(options.verbose)
        except Exception as exc:
            exc_list.append(exc)

    if exc_list:
        print("\nList of exceptions raised by validate_checksums:")
        for exc in exc_list:
            print(exc)

    return len(exc_list)


def abips_avail(options):
    """
    Show available repos.
    """
    print("List of available pseudopotential repositories:\n")
    all_repos = get_all_registered_repos()
    print(tabulate_repos(all_repos, with_citations=True, verbose=options.verbose))


#def abips_nc_install(options):
#    """
#    Get NC repo. Can choose among three formats: psp8, upf2 and psml.
#    By default we fetch all formats.
#    """
#    repos_root = get_repos_root(options)
#    all_repos = get_all_registered_repos()
#    repos = [repo for repo in ALL_REPOS if repo.isnc and not repo.is_installed(repos_root)]
#    if not repos:
#        print(f"All registered NC repositories are already installed in {repos_root}. Returning")
#        return 0
#
#    print("The following NC repositories will be installed:\n")
#    pprint_repos(repos, repos_root=repos_root)
#    if not options.yes and user_wants_to_abort():
#        return 2
#
#    print("Fetching NC repositories. It may take some time ...")
#    for repo in repos:
#        repo.install(repos_root, options.verbose)
#
#    if options.verbose: abips_list(options)
#    return 0
#
#
#def abips_paw_install(options):
#    """
#    Get PAW repositories in PAWXML format.
#    """
#    repos_root = get_repos_root(options)
#    all_repos = get_all_registered_repos()
#    repos = [repo for repo in ALL_REPOS if repo.ispaw and not repo.is_installed(repos_root)]
#    if not repos:
#        print(f"All registered PAW repositories are already installed in {repos_root}. Returning")
#        return 0
#
#    print("The following PAW repositories will be installed:")
#    pprint_repos(repos, repos_root=repos_root)
#    if not options.yes and user_wants_to_abort():
#        return 2
#
#    print("Fetching PAW repositories. It may take some time ...")
#    for repo in repos:
#        repo.install(repos_root, options.verbose)
#
#    if options.verbose: abips_list(options)
#    return 0


def abips_install(options):
    """
    Install pseudopotential repositories by name(s).
    Use the `avail` command to get the repo name.
    """
    repos = repos_from_names(options.repo_names)
    repos = [repo for repo in repos if not repo.is_installed()]

    if not repos:
        print("Table(s) are already installed! Nothing to do. Returning")
        return 0

    print("The following pseudopotential repositories will be installed:")
    print(tabulate_repos(repos, verbose=options.verbose))
    #if not options.yes and user_wants_to_abort(): return 2

    for repo in repos:
        repo.install(verbose=options.verbose)

    if options.verbose: abips_list(options)
    return 0


def abips_show(options):
    """
    Show pseudopotential tables
    """
    repos = repos_from_names(options.repo_names)
    repos = [repo for repo in repos if repo.is_installed()]

    if not repos:
        print(f"There's no installed repository with name in: {options.repo_names}")
        return 1

    for repo in repos:
        print(repo)
        table_names = ["standard", "stringent"]
        for table_name in table_names:
            pseudos = repo.get_pseudos(table_name=table_name)
            #for pseudo in pseudos:
            #    print(pseudo.filepath)
            #    print(pseudo.as_dict()["filepath"])
            print(f"For accuracy: {table_name}:")
            #print(pseudos.summarize())
            print(pseudos)

    return 0


def abips_mkff(options):
    """
    Call Abinit to compute PSPS.nc file and show results
    """
    from abipy.electrons.psps import PspsFile
    with PspsFile.from_abinit_run(options.pseudo_path) as abifile:
        print(abifile)
        abifile.expose(use_web=True,
                       #slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                       #use_web=options.expose_web, verbose=options.verbose
                       )


def get_epilog():
    return """\

Usage example:

  abips.py avail                             --> Show all registered pseudopotential repositories
  abips.py list                              --> List repositories installed on this machine
  abips.py install ONCVPSP-PBEsol-SR-PDv0.4  --> Install repository by name (requires internet connection).
"""

  #abips.py show ONCVPSP-PBEsol-SR-PDv0.4  --> Install repository.
  #abips.py nc_install                        --> Get all NC repositories (most recent version)
  #abips.py nc_install -xc PBE -fr -sr --version 0.4
  #abips.py paw_install                        --> Get all PAW repositories (most recent version)


def get_parser(with_epilog=False):

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                              help='verbose, can be supplied multiple times to increase verbosity.')

    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                              help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    #copts_parser.add_argument('--repos-root', "-r", type=str,
    #                          default=os.path.expanduser(os.path.join("~", ".abinit", "pseudos")),
    #                          help='Installation directory. Default: $HOME/.abinit/pseudos')

    #copts_parser.add_argument('-y', "--yes", action="store_true", default=False,
    #                          help="Do not ask for confirmation when installing repositories.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for avail command.
    subparsers.add_parser("avail", parents=[copts_parser], help=abips_avail.__doc__)

    # Subparser for list command.
    p_list = subparsers.add_parser("list", parents=[copts_parser], help=abips_list.__doc__)
    p_list.add_argument("-c", "--checksums", action="store_true", default=False,
                        help="Validate md5 checksums")

    # Subparser for install command.
    p_install = subparsers.add_parser("install", parents=[copts_parser], help=abips_install.__doc__)
    p_install.add_argument("repo_names", type=str, nargs="+", help="List of repositories to download.")
    p_install.add_argument("-c", "--checksums", action="store_true", default=False,
                           help="Validate md5 checksums")

    # Subparser for nc_install command.
    #p_nc_install = subparsers.add_parser("nc_install", parents=[copts_parser], help=abips_nc_install.__doc__)

    # Subparser for paw_install command.
    #p_paw_install = subparsers.add_parser("paw_install", parents=[copts_parser], help=abips_paw_install.__doc__)

    # Subparser for show command.
    p_show = subparsers.add_parser("show", parents=[copts_parser], help=abips_show.__doc__)
    p_show.add_argument("repo_names", type=str, nargs="+", help="List of Repo names.")

    # Subparser for show command.
    p_mkff = subparsers.add_parser("mkff", parents=[copts_parser], help=abips_show.__doc__)
    p_mkff.add_argument("pseudo_path", type=str, help="Pseudopotential path")

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
        print(exc)
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    #if options.repos_root:

    return globals()[f"abips_{options.command}"](options)


if __name__ == "__main__":
    sys.exit(main())

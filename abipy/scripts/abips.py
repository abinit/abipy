#!/usr/bin/env python
"""
Script to download and install pseudopotential tables from the web.
"""
from __future__ import annotations

import sys
#import os
import argparse
import abipy.tools.cli_parsers as cli

from monty.termcolor import cprint
from pymatgen.core.periodic_table import Element
from abipy.core.release import __version__
from abipy.tools import duck
from abipy.flowtk.pseudos import PseudoTable
from abipy.flowtk.psrepos import (tabulate_repos, repos_from_names,
                                  get_all_registered_repos, get_installed_repos_and_root)


def abips_list(options) -> list:
    """
    List installed pseudopotential repos.
    """
    repos, repos_root = get_installed_repos_and_root()
    if not repos:
        print("Could not find any pseudopotential repository installed in:", repos_root)
        return 0

    print(f"The following pseudopotential repositories are installed in {repos_root}:\n")
    print(tabulate_repos(repos, exclude=["installed"], verbose=options.verbose))

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


def abips_avail(options) -> int:
    """
    Show available pseudopotential repos.
    """
    print("List of available pseudopotential repositories:\n")
    all_repos = get_all_registered_repos()
    print(tabulate_repos(all_repos, with_citations=True, verbose=options.verbose))
    return 0


#def abips_nc_install(options) -> int:
#    """
#    Get all NC repos for a given version,
#    Can choose among three formats: psp8, upf2 and psml. By default we fetch all formats.
#    """
#    repos_root = get_repos_root(options)
#    all_repos = get_all_registered_repos()
#    repos = [repo for repo in all_repos if repo.isnc and not repo.is_installed(repos_root)]
#    if not repos:
#        print(f"All registered NC repositories are already installed in {repos_root}. Returning")
#        return 0
#
#    print("The following NC repositories will be installed:\n")
#    pprint_repos(repos, repos_root=repos_root)
#    if not options.yes and cli.user_wants_to_abort():
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
#    Get all JTH PAW repositories in PAWXML format for the given version.
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
#    if not options.yes and cli.user_wants_to_abort():
#        return 2
#
#    print("Fetching PAW repositories. It may take some time ...")
#    for repo in repos:
#        repo.install(repos_root, options.verbose)
#
#    if options.verbose: abips_list(options)
#    return 0


def abips_install(options) -> int:
    """
    Install pseudopotential repositories by name(s).
    Use `avail` command to get repo names.
    """
    repos = repos_from_names(options.repo_names)
    repos = [repo for repo in repos if not repo.is_installed()]

    if not repos:
        print("Table(s) are already installed! Nothing to do. Returning")
        return 0

    print("The following pseudopotential repositories will be installed:")
    print(tabulate_repos(repos, verbose=options.verbose))
    #if not options.yes and cli.user_wants_to_abort(): return 2

    for repo in repos:
        repo.install(verbose=options.verbose)

    if options.verbose:
        abips_list(options)

    return 0


def abips_show(options) -> int:
    """
    Show info on pseudopotential table(s).
    """
    repos = repos_from_names(options.repo_names)
    repos = [repo for repo in repos if repo.is_installed()]

    if not repos:
        print(f"There's no installed repository with name in: {options.repo_names}")
        return 1

    for repo in repos:
        print(repo)
        for table_name in repo.table_names:
            pseudos = repo.get_pseudos(table_name=table_name)
            print(f"For table_name: {table_name}:\n")
            if options.symbol is not None:
                print("Selecting pseudos with symbol:", options.symbol)
                pseudos = pseudos.pseudo_with_symbol(options.symbol, allow_multi=False)
            print(pseudos)

    return 0


def abips_element(options) -> int:
    """
    Find all pseudos in the installed tables for the given element (symbol or znucl).
    """
    # Accept symbol string or Z
    symbol = options.element
    if symbol.isnumeric():
        symbol = Element.from_Z(int(symbol)).symbol

    repos, repos_root = get_installed_repos_and_root()
    if not repos:
        print("Could not find any pseudopotential repository installed in:", repos_root)
        return 1

    pseudo_list = []
    for repo in repos:
        for table_name in repo.table_names:
            try:
                pseudos = repo.get_pseudos(table_name=table_name)
                pseudo = pseudos.pseudo_with_symbol(symbol, allow_multi=False)
            except Exception as exc:
                cprint(f"{str(exc)}", "red")
                continue

            if pseudo not in pseudo_list:
                pseudo_list.append(pseudo)

    pseudos = PseudoTable(pseudo_list)
    pseudos.print_table()

    return 0


def abips_mkff(options) -> int:
    """
    Call Abinit to compute PSPS.nc files from a list of pseudos and show results.
    """
    from abipy.electrons.psps import PspsFile, PspsRobot
    ecut = options.ecut

    if len(options.pseudo_paths) == 1:
        if options.vloc_rcut_list is None:
            with PspsFile.from_abinit_run(options.pseudo_paths[0], ecut) as abifile:
                if options.verbose: print(abifile)
                abifile.expose(use_web=options.expose_web,
                               slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                               verbose=options.verbose
                               )

        else:
            robot = PspsRobot.from_vloc_rcut_list(options.pseudo_paths[0], options.vloc_rcut_list, ecut)
            if options.verbose: print(robot)
            robot.expose(use_web=options.expose_web,
                         slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                         verbose=options.verbose
                         )

    else:
        if options.vloc_rcut_list is not None:
            raise ValueError("vloc_rcut_list does not support more than one pseudo!")
        robot = PspsRobot.from_abinit_run(options.pseudo_paths, ecut)

        if options.verbose: print(robot)
        robot.expose(use_web=options.expose_web,
                     slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                     verbose=options.verbose
                     )

    return 0


def get_epilog() -> str:
    return """\

Usage example:

  abips.py avail                             --> Show all registered pseudopotential repositories.
  abips.py list                              --> List repositories installed on this machine.
  abips.py install ONCVPSP-PBEsol-SR-PDv0.4  --> Install repository by name (requires internet connection).
  abips.py show ONCVPSP-PBEsol-SR-PDv0.4     --> Show info on pseudos in repository.
  abips.py element O                         --> Show all installed pseudos for element.
  abips.py mkff PSEUDO1 [PSEUDO2 ...]        --> Compute form factors for pseudos and show them.
"""
  #abips.py onc_install                        --> Get all NC repositories (most recent version)
  #abips.py onc_install -xc PBE -fr -sr -v 0.4
  #abips.py jth_install                        --> Get all JTH PAW repositories (most recent version)


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
                        help="Validate md5 checksums.")

    # Subparser for install command.
    p_install = subparsers.add_parser("install", parents=[copts_parser], help=abips_install.__doc__)
    p_install.add_argument("repo_names", type=str, nargs="+", help="List of repositories to download.")
    p_install.add_argument("-c", "--checksums", action="store_true", default=False,
                           help="Validate md5 checksums.")

    # Subparser for onc_install command.
    #p_onc_install = subparsers.add_parser("onc_install", parents=[copts_parser], help=abips_onc_install.__doc__)
    #p_onc_install.add_argument("-v", type=str, default=None, help="Table version. Default: latest one ")

    # Subparser for jth_install command.
    #p_jth_install = subparsers.add_parser("jth_install", parents=[copts_parser], help=abips_paw_install.__doc__)
    #p_jth_install.add_argument("-v", type=str, default=None, help="Table version. Default: latest one ")

    # Subparser for show command.
    p_show = subparsers.add_parser("show", parents=[copts_parser], help=abips_show.__doc__)
    p_show.add_argument("repo_names", type=str, nargs="+", help="List of repo names.")
    p_show.add_argument("-s", "--symbol", type=str, default=None, help="Select pseudo by element symbol.")

    # Subparser for list command.
    p_element = subparsers.add_parser("element", parents=[copts_parser], help=abips_element.__doc__)
    p_element.add_argument("element", type=str, help="Element symbol or atomic number.")

    # Subparser for mkff command.
    p_mkff = subparsers.add_parser("mkff", parents=[copts_parser], help=abips_mkff.__doc__)
    p_mkff.add_argument("pseudo_paths", nargs="+", type=str, help="Pseudopotential path.")
    p_mkff.add_argument("--ecut", type=float, required=True, help="Cutoff energy in Ha.")
    p_mkff.add_argument("-rc", "--vloc-rcut-list", nargs="+", default=None, type=float,
                        help="List of cutoff radii for vloc in Bohr.")
    cli.add_expose_options_to_parser(p_mkff)

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

    #cli.set_loglevel(options.loglevel)

    # Use seaborn settings.
    if hasattr(options, "seaborn") and options.seaborn:
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)

    return globals()[f"abips_{options.command}"](options)


if __name__ == "__main__":
    sys.exit(main())

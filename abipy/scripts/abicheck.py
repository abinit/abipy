#!/usr/bin/env python
"""
Diagnostic tool to verify the AbiPy environment and configuration.

This script checks that the required Python packages are installed, the 
Abinit executable is reachable, and the configuration files (manager.yml, 
scheduler.yml) are correctly set up.

It can also be used to automatically install template configuration files 
and run a small test flow to verify the entire toolchain.

Examples:
    Run basic consistency checks:
        $ abicheck.py

    Run a full test including the execution of a small Abinit flow:
        $ abicheck.py --with-flow

    Install template configuration files in ~/.abinit/abipy/:
        $ abicheck.py -c

    Show available manager templates provided by AbiPy:
        $ abicheck.py -m
"""

from __future__ import annotations

import argparse
import os
import sys

from monty import termcolor
from monty.termcolor import cprint

import abipy.data as abidata
import abipy.tools.cli_parsers as cli
from abipy import abilab, flowtk


def show_managers(options):
    """
    Print a table with manager configuration files provided by AbiPy.

    Args:
        options: Namespace object containing command-line options.

    Returns:
        int: System exit code.
    """
    from tabulate import tabulate

    table = []
    root = os.path.join(abidata.dirpath, "managers")
    yaml_paths = [os.path.join(root, f) for f in os.listdir(root) if f.endswith(".yml") and "_manager" in f]
    for path in yaml_paths:
        manager = flowtk.TaskManager.from_file(path)
        hostname = os.path.basename(path).split("_")[0]
        table.append([hostname, manager.qadapter.QTYPE, path])
        if options.verbose > 1:
            print(manager)
    print(tabulate(table, headers=["hostname", "queue-type", "filepath"], tablefmt="rst"))
    return 0


def get_epilog() -> str:
    """
    Return the epilog string for the command-line parser.
    """
    return """\
Usage example:
    abicheck.py                ==> Test abipy installation and requirements.
    abicheck.py -m             ==> Print table with manager files provided by AbiPy.
    abicheck.py -c             ==> Install Yaml configuration files for manager and scheduler in ~/.abinit/abipy dir.
    abicheck.py --with-flow    ==> Consistency check followed by execution of AbiPy flow.
"""


def get_parser(with_epilog=False):
    """
    Return the ArgumentParser object for the script.

    Args:
        with_epilog: If True, include the epilog in the parser.
    """
    parser = argparse.ArgumentParser(
        epilog=get_epilog() if with_epilog else "", formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--loglevel",
        default="ERROR",
        type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG",
    )
    parser.add_argument("-V", "--version", action="version", version=abilab.__version__)
    parser.add_argument(
        "-v",
        "--verbose",
        default=0,
        action="count",  # -vv --> verbose=2
        help="verbose, can be supplied multiple times to increase verbosity.",
    )
    parser.add_argument("--no-colors", default=False, action="store_true", help="Disable ASCII colors.")
    parser.add_argument(
        "--with-flow", default=False, action="store_true", help="Build and run small abipy flow for testing."
    )
    parser.add_argument(
        "-d",
        "--flow-dir",
        type=str,
        default=None,
        help="Create AbiPy flow in this directory. If None, a default directory is used,",
    )
    parser.add_argument(
        "-m",
        "--show-managers",
        default=False,
        action="store_true",
        help="Print table with manager files provided by AbiPy.",
    )

    parser.add_argument(
        "-c",
        "--create-config",
        default=False,
        action="store_true",
        help="Create yaml configuration files in ~/abinit/.abipy with predefined settings.",
    )
    parser.add_argument(
        "-f",
        "--force-reinstall",
        default=False,
        action="store_true",
        help="Overwrite yaml configuration files if --create-config and files already exist.",
    )
    return parser


@cli.prof_main
def main():
    """
    Main entry point for the script.
    """

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    cli.set_loglevel(options.loglevel)

    if options.no_colors:
        # Disable colors
        termcolor.enable(False)

    if options.show_managers:
        # Show table with manager configuration files.
        return show_managers(options)

    if options.create_config:
        # Install Yaml configuration files for manager and scheduler.
        abilab.install_config_files(workdir=None, force_reinstall=options.force_reinstall)

    errmsg = abilab.abicheck(verbose=options.verbose)
    if errmsg:
        cprint(errmsg, "red")
        cprint(
            "TIP: Use `--show-managers` to print the manager files provided by AbiPy.\n"
            "If abicheck.py is failing because it cannot find the manager.yml configuration file",
            "yellow",
        )
        return 2
    cprint("\nAbipy requirements are properly configured\n", "green")

    if not options.with_flow:
        return 0

    retcode = run_flow(options)
    if retcode == 0:
        cprint("\nTest flow completed successfully\n", "green")

    return retcode


def make_scf_nscf_inputs(paral_kgb=0):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos("14si.pspnc")
    # pseudos = data.pseudos("Si.GGA_PBE-JTH-paw.xml")

    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"), pseudos=pseudos, ndtset=2)

    # Global variables
    ecut = 6
    global_vars = dict(
        ecut=ecut,
        nband=8,
        timopt=-1,
        istwfk="*1",
        nstep=15,
        paral_kgb=paral_kgb,
        iomode=3,
    )

    if multi.ispaw:
        global_vars.update(pawecutdg=2 * ecut)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(ngkpt=[8, 8, 8], shiftk=[0, 0, 0])
    multi[0].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0],  # L point
        [0.0, 0.0, 0.0],  # Gamma point
        [0.0, 0.5, 0.5],  # X point
    ]

    multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
    multi[1].set_vars(tolwfr=1e-12)

    # Generate two input files for the GS and the NSCF run
    scf_input, nscf_input = multi.split_datasets()
    return scf_input, nscf_input


def run_flow(options):
    """
    Build and run a small Abinit flow to test the environment.

    Args:
        options: Namespace object containing command-line options.

    Returns:
        int: Exit code of the scheduler.
    """
    import tempfile

    workdir = tempfile.mkdtemp(dir=options.flow_dir)
    cprint("Running small flow in workdir: %s" % workdir, "yellow")
    print()

    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs()

    # Build the flow.
    flow = flowtk.bandstructure_flow(workdir, scf_input, nscf_input, manager=None)
    flow.build_and_pickle_dump()
    scheduler = flow.make_scheduler()
    retcode = scheduler.start()
    if retcode != 0:
        cprint("Scheduler returned retcode %s" % retcode, "red")
        return retcode

    flow.show_status()
    if not flow.all_ok:
        cprint("Not all tasks in flow reached all_ok", "red")
        retcode = 1

    return retcode


if __name__ == "__main__":
    sys.exit(main())

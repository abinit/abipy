#!/usr/bin/env python
"""
Interface to the database of ABINIT input variables
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import abipy.flowtk as flowtk

from pprint import pprint
from monty.functools import prof_main
from monty.termcolor import cprint
from abipy.core.release import __version__
from abipy import abilab


def print_vlist(vlist, options):
    for v in vlist:
        print(repr(v))

    if options.verbose:
        for v in vlist: abilab.abinit_help(v)
    else:
        print("\nUse -v for more info")


def get_epilog():
    return """\
Usage example:

    abidoc.py man ecut        --> Show documentation for ecut input variable.
    abidoc.py graphviz acell  --> Draw parents and children of variable with graphviz package.
    abidoc.py browse acell    --> Open url in external browser.
    abidoc.py apropos ecut    --> To search in the database for the variables related to ecut.
    abidoc.py find paw        --> To search in the database for the variables whose name contains paw.
    abidoc.py list            --> Print full list of variables.
    abidoc.py withdim natom   --> Print arrays depending on natom.
    abidoc.py scheduler       --> Document scheduler options.
    abidoc.py manager         --> Document manager options.

Use `abidoc.py --help` for help and `abidoc.py COMMAND --help` to get the documentation for `COMMAND`.
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).
"""


def get_parser(with_epilog=False):

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                              help='verbose, can be supplied multiple times to increase verbosity')
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                              help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    copts_parser.add_argument("-c", '--codename', type=str, default="abinit",
                              help="Code name e.g. anaddb, optic... Default: abinit")

    var_parser = argparse.ArgumentParser(add_help=False)
    var_parser.add_argument('varname', help="ABINIT variable")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for man.
    p_man = subparsers.add_parser('man', parents=[copts_parser, var_parser], help="Show documentation for varname.")

    # Subparser for graphviz.
    p_graphviz = subparsers.add_parser('graphviz', parents=[copts_parser, var_parser],
        help=("Draw variable dependencies with graphviz package."
              "See https://graphviz.readthedocs.io/."))
    p_graphviz.add_argument("-e", "--engine", type=str, default="automatic",
        help=("graphviz engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']. "
            "Default: automatic i.e. the engine is automatically selected. See http://www.graphviz.org/pdf/dot.1.pdf "
            "Use `conda install python-graphviz` or `pip install graphviz` to install the python package"))

    # Subparser for browse.
    p_browse = subparsers.add_parser('browse', parents=[copts_parser, var_parser], help="Open documentation in browser.")

    # Subparser for apropos.
    p_apropos = subparsers.add_parser('apropos', parents=[copts_parser, var_parser],
                                      help="Find variables related to varname.")
    # Subparser for find.
    p_find = subparsers.add_parser('find', parents=[copts_parser, var_parser],
                                    help="Find all variables whose name contains varname.")
    # Subparser for require.
    #p_require = subparsers.add_parser('require', parents=[copts_parser], help="Find all variables required by varname.")

    # Subparser for withdim.
    p_withdim = subparsers.add_parser('withdim', parents=[copts_parser],
                                      help="Find all arrays depending on the given dimension.")
    p_withdim.add_argument("dimname", help="Dimension name")

    # Subparser for list.
    p_list = subparsers.add_parser('list', parents=[copts_parser], help="List all variables.")
    p_list.add_argument('--mode', default="a",
                        help="Sort mode, `a` for alphabethical, `s` for varset, `c` for characteristics.")

    # Subparser for manager.
    p_manager = subparsers.add_parser('manager', parents=[copts_parser], help="Document the TaskManager options.")
    p_manager.add_argument("qtype", nargs="?", default=None, help=("Write job script to terminal if qtype='script' else "
        "document the qparams for the given QueueAdapter qtype e.g. slurm."))

    # Subparser for scheduler
    p_docsched = subparsers.add_parser('scheduler', parents=[copts_parser],
        help="Document the options available in scheduler.yml.")

    # Subparser for abibuild
    p_abibuild = subparsers.add_parser('abibuild', parents=[copts_parser],
        help="Show ABINIT build information and exit.")

    return parser


@prof_main
def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    # Get the dabase of variables for codename.
    from abipy.abio.abivar_database.variables import get_codevars
    vdb = get_codevars()[options.codename]

    if options.command == "man":
        abilab.abinit_help(options.varname)

    elif options.command == "browse":
        return vdb[options.varname].browse()

    elif options.command == "graphviz":
        if options.varname in vdb.my_varset_list:
            graph = vdb.get_graphviz(varset=options.varname, vartype=None, engine=options.engine)
        else:
            graph = vdb.get_graphviz_varname(varname=options.varname, engine=options.engine)

        import tempfile
        directory = tempfile.mkdtemp()
        print("Producing source files in:", directory)
        graph.view(directory=directory, cleanup=False)

    elif options.command == "apropos":
        vlist = vdb.apropos(options.varname)
        print_vlist(vlist, options)

    elif options.command == "find":
        vlist = [v for v in vdb.values() if options.varname in v.name]
        print("Find results:\n")
        print_vlist(vlist, options)

    elif options.command == "list":
        if options.mode == "a":
            # Alphabetical
            for i, var in enumerate(vdb.values()):
                print(i, repr(var))

        elif options.mode == "s":
            # Grouped by varset
            for section in vdb.my_varset_list:
                print(30*"#" +  " Section: " + section + " " + 30*"#")
                print_vlist(vdb.vars_with_section(section), options)

        elif options.mode == "c":
            # Grouped by characteristics.
            for char in vdb.my_characteristics:
                print(30*"#" +  " Characteristic: " + char + 30*"#")
                print_vlist(vdb.vars_with_char(char), options)

        else:
            raise ValueError("Wrong mode %s" % options.mode)

    elif options.command == "withdim":
        for var in vdb.values():
            if var.depends_on_dimension(options.dimname):
                cprint(repr(var), "yellow")
                print("dimensions:", str(var.dimensions), "\n")

    elif options.command == "scheduler":
        print("Options that can be specified in scheduler.yml:")
        print(flowtk.PyFlowScheduler.autodoc())
        return 0

    elif options.command == "manager":
        # Document TaskManager options and qparams.
        qtype = options.qtype

        if qtype == "script":
            manager = flowtk.TaskManager.from_user_config()
            script = manager.qadapter.get_script_str(
                job_name="job_name", launch_dir="workdir", executable="executable",
                qout_path="qout_file.path", qerr_path="qerr_file.path",
                stdin="stdin", stdout="stdout", stderr="stderr")
            print(script)

        else:
            print(flowtk.TaskManager.autodoc())
            print("qtype supported: %s" % flowtk.all_qtypes())
            print("Use `abidoc.py manager slurm` to have the list of qparams for slurm.\n")

            if qtype is not None:
                print("QPARAMS for %s" % qtype)
                flowtk.show_qparams(qtype)

    elif options.command == "abibuild":
        abinit_build = flowtk.AbinitBuild()
        print()
        print(abinit_build)
        print()
        if not options.verbose:
            print("Use --verbose for additional info")
        else:
            print(abinit_build.info)

    else:
        raise ValueError("Don't know how to handle command %s" % options.command)

    return 0


if __name__ == "__main__":
    sys.exit(main())

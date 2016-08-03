#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse

from abipy import abilab


def abicomp_struct(options):
    """
    Compare crystalline structures.
    """
    paths = options.paths
    index = [os.path.relpath(p) for p in paths]
    frame = abilab.frame_from_structures(paths, index=None)
    print("File list:")
    for i, p in enumerate(paths):
        print("%d %s" % (i, p))
    print()
    abilab.print_frame(frame)
    return 0


def abicomp_ebands(options):
    """
    Plot electron bands on a grid.
    """
    paths = options.paths
    eb_objects = paths
    titles = paths
    e0 = "fermie"
    mode = "single"
    mode = "grid"

    #if mode == "grid":
    #    abilab.ebands_gridplot(eb_objects, titles=titles, edos_objects=None, edos_kwargs=None)
    #elif mode == "single":
    #    abilab.ebands_singleplot(eb_objects, edos_objects=None, edos_kwargs=None, e0=e0)
    #elif mode == "animate":
    #    abilab.ebands_animate(eb_objects, edos_objects=None, edos_kwargs=None, e0=e0,
    #                          interval=250, savefile=None)

    return 0


def abicomp_phbands(options):
    """
    Plot phonon bands on a grid.
    """
    paths = options.paths
    phb_objects = paths
    titles = paths
    abilab.phbands_gridplot(phb_objects, titles=titles, phdos_objects=None, phdos_kwargs=None)
    return 0


#def abicomp_pseudos(options):
#    paths = options.paths
#    index = [os.path.relpath(p) for p in paths]
#    frame = abilab.frame_from_pseudos(paths, index=None)
#    print("File list:")
#    for i, p in enumerate(paths):
#        print("%d %s" % (i, p))
#    print()
#    abilab.print_frame(frame)
#    return 0


def abicomp_robot(options):
    """
    Analyze multiple files with a robot. The command has two different variants.

    [1] The files can be listed explicitly as in:

            abicomp.py robot out1_GSR.nc out1_GSR.nc

        or, alternatively, with the shell syntax:

            abicomp.py robot *_GSR.nc

    [2] It's possible to use a directory as the first (and only) argument
        of the robot command followed by the Abinit file extension as in:

            abicomp.py robot . GSR.nc

    By default, the script with call `robot.get_dataframe()` to create an print a table
    with the most important results. For finer control, use --ipython to start an ipython
    console to interact with the robot directly.
    """
    # To define an Help action
    #http://stackoverflow.com/questions/20094215/argparse-subparser-monolithic-help-output?rq=1
    paths = options.paths

    # Temp code
    robot = abilab.abirobot(".", "GSR")
    print(robot)
    abilab.print_frame(robot.get_dataframe())
    return 0

    if os.path.isdir(paths[0]):
        # Directory + extension
        top, ext = paths[:2]
        cls = Robot.class_for_ext(ext)
        robot = cls.from_dir(top, walk=True)
    else:
        # File list.
        #cls = Robot.class_for_ext(ext)
        robot = Robot.from_files(paths)

    if options.ipython:
        import IPython
        IPython.embed(header=str(robot) + "\nType `robot` in the terminal and use <TAB> to list its methods",  robot=robot)
    elif options.notebook:
        robot.make_and_open_notebook(nbpath=None, daemonize=True)
    else:
        print(robot)
        abilab.print_frame(robot.get_dataframe())

    return 0

def abicomp_gs_scf(options):
    """
    Compare ground-state SCF cycles.
    """
    paths = options.paths
    f0 = abilab.AbinitOutputFile(paths[0])
    f0.compare_gs_scf_cycles(paths[1:])
    return 0


def abicomp_dfpt2_scf(options):
    """
    Compare DFPT SCF cycles.
    """
    paths = options.paths
    f0 = abilab.AbinitOutputFile(paths[0])
    f0.compare_d2de_scf_cycles(paths[1:])
    return 0


def main():
    def str_examples():
        return """\
Usage example:
  abidiff.py struct */*/outdata/out_GSR.nc        => Compare structures in multiple files.
  abidiff.py ebands out1_GSR.nc out2_GSR.nc       => Plot electron bands on a grid.
  abidiff.py phbands out1_PHBST.nc out2_PHBST.nc  => Plot electron bands on a grid.
  abidiff.py gs_scf run1.abo run2.abo             => Compare the SCF cycles in two output files.
  abidiff.py dfpt2_scf                            => Compare the DFPT SCF cycles in two output files.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('paths', nargs="+", help="List of files to compare")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='Verbose, can be supplied multiple times to increase verbosity')
    copts_parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + abilab.__version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for struct command.
    p_struct = subparsers.add_parser('struct', parents=[copts_parser], help=abicomp_struct.__doc__)

    # Subparser for ebands command.
    p_ebands = subparsers.add_parser('ebands', parents=[copts_parser], help=abicomp_ebands.__doc__)

    # Subparser for phbands command.
    p_phbands = subparsers.add_parser('phbands', parents=[copts_parser], help=abicomp_phbands.__doc__)

    # Subparser for pseudos command.
    #p_pseudos = subparsers.add_parser('pseudos', parents=[copts_parser], help=abicomp_pseudos.__doc__)

    # Subparser for robot command.
    p_robot = subparsers.add_parser('robot', parents=[copts_parser], help=abicomp_robot.__doc__)

    # Subparser for gs_scf command.
    p_gs_scf = subparsers.add_parser('gs_scf', parents=[copts_parser], help=abicomp_gs_scf.__doc__)

    # Subparser for dfpt2_scf command.
    p_dftp2_scf = subparsers.add_parser('dfpt2_scf', parents=[copts_parser], help=abicomp_dfpt2_scf.__doc__)

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.seaborn:
        import seaborn as sns
        #sns.set(style="dark", palette="Set2")
        #sns.set(style='ticks', palette='Set2')
        #sns.despine()

    # Dispatch
    return globals()["abicomp_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())

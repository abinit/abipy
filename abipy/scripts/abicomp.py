#!/usr/bin/env python
"""
Script to analyze/compare results stored in multiple netcdf files.
By default the script displays the results/plots in the shell.
Use --ipython to start an ipython terminal or -nb to generate an ipython notebook.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse

from monty.functools import prof_main
from monty.termcolor import cprint
from abipy import abilab


def abicomp_structure(options):
    """
    Compare crystalline structures.
    """
    paths = options.paths
    index = [os.path.relpath(p) for p in paths]

    if options.notebook:
        import nbformat
        nbv = nbformat.v4
        nb = nbv.new_notebook()

        nb.cells.extend([
            nbv.new_code_cell("""\
from __future__ import print_function, division, unicode_literals, absolute_import
import sys
import os

%matplotlib notebook
from IPython.display import display
#import seaborn as sns

from abipy import abilab"""),
            nbv.new_code_cell("dfs = abilab.frames_from_structures(%s, index=%s)" % (paths, index)),
            # Analyze dataframes.
            nbv.new_code_cell("dfs.lattice"),
            nbv.new_code_cell("dfs.coords"),
            nbv.new_code_cell("# for structure in dfs.structures: display(structure)"),
        ])

        import io, tempfile # os,
        _, nbpath = tempfile.mkstemp(prefix="abinb_", suffix='.ipynb', dir=os.getcwd(), text=True)

        # Write notebook
        import nbformat
        with io.open(nbpath, 'wt', encoding="utf8") as fh:
            nbformat.write(nb, fh)

        cmd = "jupyter notebook %s" % nbpath
        return os.system(cmd)

    dfs = abilab.frames_from_structures(paths, index=index)

    if options.ipython:
        import IPython
        IPython.embed(header="Type `dfs` in the terminal and use <TAB> to list its methods", dfs=dfs)
    else:
        print("File list:")
        for i, p in enumerate(paths):
            print("%d %s" % (i, p))
        print()
        abilab.print_frame(dfs.lattice, title="Lattice parameters:")
        abilab.print_frame(dfs.coords, title="Atomic positions (columns give the site index):")

    return 0


def abicomp_ebands(options):
    """
    Plot electron bands on a grid.
    """
    paths, e0 = options.paths, options.e0
    plotter = abilab.ElectronBandsPlotter(key_ebands=[(os.path.relpath(p), p) for p in paths])

    if options.ipython:
        import IPython
        IPython.embed(header=str(plotter) + "\nType `plotter` in the terminal and use <TAB> to list its methods",
                      plotter=plotter)

    elif options.notebook:
        plotter.make_and_open_notebook(foreground=options.foreground)

    else:
        # Print pandas Dataframe.
        frame = plotter.get_ebands_frame()
        abilab.print_frame(frame)

        # Optionally, print info on gaps and their location
        if not options.verbose:
            print("\nUse --verbose for more information")
        else:
            for ebands in plotter.ebands_list:
                print(ebands)

        # Here I select the plot method to call.
        if options.plot_mode != "None":
            plotfunc = getattr(plotter, options.plot_mode, None)
            if plotfunc is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % options.plot_mode)
            plotfunc(e0=e0)

    return 0


def abicomp_edos(options):
    """
    Plot electron DOSes on a grid.
    """
    paths, e0 = options.paths, options.e0
    plotter = abilab.ElectronDosPlotter(key_edos=[(os.path.relpath(p), p) for p in paths])

    if options.ipython:
        import IPython
        IPython.embed(header=str(plotter) + "\nType `plotter` in the terminal and use <TAB> to list its methods",
                      plotter=plotter)

    elif options.notebook:
        plotter.make_and_open_notebook(foreground=options.foreground)

    else:
        # Print pandas Dataframe.
        #frame = plotter.get_ebands_frame()
        #abilab.print_frame(frame)

        # Optionally, print info on gaps and their location
        if not options.verbose:
            print("\nUse --verbose for more information")
        else:
            for edos in plotter.edos_list:
                print(edos)

        # Here I select the plot method to call.
        if options.plot_mode != "None":
            plotfunc = getattr(plotter, options.plot_mode, None)
            if plotfunc is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % options.plot_mode)
            plotfunc(e0=e0)

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
    with the most important results. For finer control, use --ipy to start an ipython
    console to interact with the robot directly or --nb to generate a jupyter notebook.
    """
    # To define an Help action
    #http://stackoverflow.com/questions/20094215/argparse-subparser-monolithic-help-output?rq=1
    paths = options.paths

    # Temporary code
    #robot = abilab.abirobot(".", "GSR")

    if os.path.isdir(paths[0]):
        # Assume directory + extension
        top, ext = paths[:2]
        cls = abilab.Robot.class_for_ext(ext)
        robot = cls.from_dir(top, walk=True)
    else:
        # Assume file list.
        ext = "GSR"
        cls = abilab.Robot.class_for_ext(ext)
        robot = cls.from_files(paths)

    if options.ipython:
        import IPython
        IPython.embed(header=str(robot) + "\nType `robot` in the terminal and use <TAB> to list its methods",
                     robot=robot)
    elif options.notebook:
        robot.make_and_open_notebook(foreground=options.foreground)
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
    figures = f0.compare_gs_scf_cycles(paths[1:])
    if not figures:
        cprint("Cannot find GS-SCF sections in output files.", "yellow")
    return 0


def abicomp_dfpt2_scf(options):
    """
    Compare DFPT SCF cycles.
    """
    paths = options.paths
    f0 = abilab.AbinitOutputFile(paths[0])
    figures = f0.compare_d2de_scf_cycles(paths[1:])
    if not figures:
        cprint("Cannot find DFPT-SCF sections in output files.", "yellow")
    return 0


def abicomp_time(options):
    """
    Analyze/plot the timing data of single or multiple runs.
    """
    paths = options.paths

    if len(options.paths) == 1 and os.path.isdir(paths[0]):
        top = options.paths[0]
        print("Walking directory tree from top:", top, "Looking for file extension:", options.ext)
        parser, paths_found, okfiles = abilab.AbinitTimerParser.walk(top=top, ext=options.ext)

        if not paths_found:
            cprint("Empty file list!", color="magenta")
            return 1

        print("Found %d files\n" % len(paths_found))
        if okfiles != paths_found:
            badfiles = [f for f in paths_found if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)

    else:
        parser = abilab.AbinitTimerParser()
        okfiles = parser.parse(options.paths)

        if okfiles != options.paths:
            badfiles = [f for f in options.paths if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)

    if parser is None:
        cprint("Cannot analyze timing data. parser is None", color="magenta")
        return 1

    print(parser.summarize())

    if options.verbose:
        for timer in parser.timers():
            print(timer.get_dataframe())

    if options.ipython:
        cprint("Invoking ipython shell. Use parser to access the object inside ipython", color="blue")
        import IPython
        IPython.start_ipython(argv=[], user_ns={"parser": parser})
    elif options.notebook:
        parser.make_and_open_notebook(foreground=options.foreground)
    else:
        parser.plot_all()

    return 0


@prof_main
def main():
    def str_examples():
        return """\
Usage example:

  abicomp.py structure */*/outdata/out_GSR.nc         => Compare structures in multiple files.
  abicomp.py ebands out1_GSR.nc out2_GSR.nc           => Plot electron bands on a grid (Use `-p` to change plot mode)
  abicomp.py ebands *_GSR.nc -ipy                     => Build plotter object and start ipython console.
  abicomp.py ebands *_GSR.nc -nb                      => Interact with the plotter via the jupyter notebook.
  abicomp.py edos *_WFK.nc -nb                        => Compare electron DOS in the jupyter notebook.
  abicomp.py gs_scf run1.abo run2.abo                 => Compare the SCF cycles in two output files.
  abicomp.py dfpt2_scf run1.abo run2.abo              => Compare the DFPT SCF cycles in two output files.
  abicomp.py.py time [OUT_FILES]                      => Parse timing data in files and plot results
  abicomp.py.py time . --ext=abo                      => Scan directory tree from `.`, look for files with extension `abo`
                                                         parse timing data and plot results.
"""
  #abicomp.py phbands out1_PHBST.nc out2_PHBST.nc      => Plot electron bands on a grid.

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

    # Parent parser for commands support (ipython/jupyter)
    ipy_parser = argparse.ArgumentParser(add_help=False)
    ipy_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    ipy_parser.add_argument('--foreground', action='store_true', default=False,
                            help="Run jupyter notebook in the foreground.")
    ipy_parser.add_argument('-ipy', '--ipython', default=False, action="store_true", help='Invoke ipython terminal.')

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + abilab.__version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for structure command.
    p_struct = subparsers.add_parser('structure', parents=[copts_parser, ipy_parser], help=abicomp_structure.__doc__)

    # Subparser for ebands command.
    p_ebands = subparsers.add_parser('ebands', parents=[copts_parser, ipy_parser], help=abicomp_ebands.__doc__)
    p_ebands.add_argument("-p", "--plot-mode", default="gridplot",
                          choices=["gridplot", "combiplot", "boxplot", "combiboxplot", "animate", "None"],
                          help="Plot mode e.g. `-p combiplot` to plot bands on the same figure. Default is `gridplot`")
    p_ebands.add_argument("-e0", default="fermie", choices=["fermie", "None"],
                          help="Option used to define the zero of energy in the band structure plot. Default is `fermie`")

    # Subparser for edos command.
    p_edos = subparsers.add_parser('edos', parents=[copts_parser, ipy_parser], help=abicomp_edos.__doc__)
    p_edos.add_argument("-p", "--plot-mode", default="gridplot",
                          choices=["gridplot", "combiplot", "None"],
                          help="Plot mode e.g. `-p combiplot` to plot DOSes on the same figure. Default is `gridplot`")
    p_edos.add_argument("-e0", default="fermie", choices=["fermie", "None"],
                          help="Option used to define the zero of energy in the DOS plot. Default is `fermie`")

    # Subparser for phbands command.
    #p_phbands = subparsers.add_parser('phbands', parents=[copts_parser, ipy_parser], help=abicomp_phbands.__doc__)

    # Subparser for pseudos command.
    #p_pseudos = subparsers.add_parser('pseudos', parents=[copts_parser, ipy_parser], help=abicomp_pseudos.__doc__)

    # Subparser for robot command.
    # TODO
    #p_robot = subparsers.add_parser('robot', parents=[copts_parser, ipy_parser], help=abicomp_robot.__doc__)

    # Subparser for time command.
    p_time = subparsers.add_parser('time', parents=[copts_parser, ipy_parser], help=abicomp_time.__doc__)
    p_time.add_argument("-e", "--ext", default=".abo", help=("File extension for Abinit output files. "
                        "Used when the first argument is a directory. Default is `.abo`"))

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

    # Dispatch
    return globals()["abicomp_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())

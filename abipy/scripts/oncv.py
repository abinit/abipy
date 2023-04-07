#!/usr/bin/env python
"""
Script to generate/analyze/plot ONCVPSP pseudopotentials.
"""
import sys
import os
import argparse
import json
import shutil
import abipy.tools.cli_parsers as cli

from pprint import pformat
from monty.termcolor import cprint
from abipy.flowtk.pseudos import Pseudo
from abipy.ppcodes.ppgen import OncvGenerator
from abipy.ppcodes.oncv_parser import OncvParser
from abipy.ppcodes.oncv_plotter import OncvPlotter, oncv_make_open_notebook, MultiOncvPlotter


def _find_oncv_output(path):
    """
    Fix possible error in the specification of filepath when we want a `.out` file.
    Return output path.
    """
    if path.endswith(".out"): return path
    root, _ = os.path.splitext(path)
    new_path = root + ".out"
    if not os.path.exists(new_path):
        raise ValueError("Cannot find neither %s nor %s" % (path, new_path))
    cprint("Maybe you meant %s" % new_path, "yellow")
    return new_path


def oncv_nbplot(options):
    """
    Generate jupyter notebook to plot data. Requires oncvpsp output file.
    """
    out_path = _find_oncv_output(options.filepath)
    return oncv_make_open_notebook(out_path, foreground=options.foreground, classic_notebook=options.classic_notebook,
                                   no_browser=options.no_browser)


def oncv_gnuplot(options):
    """
    Plot data with gnuplot.
    """
    out_path = _find_oncv_output(options.filepath)

    # Parse output file.
    onc_parser = OncvParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not completed. Exiting", "red")
        return 1

    onc_parser.gnuplot()
    return 0


def oncv_print(options):
    """
    Print result to terminal.
    """
    p = OncvParser(options.filepath)
    p.scan()
    if not p.run_completed:
        raise RuntimeError("ocnvpsp output file is not completed")

    print(p)

    #results = p.get_results()
    #pprint(results)
    #for l in range(p.lmax + 1):
    #    # Get AE/PS logder(l) as a function of energy in Ha.
    #    f1, f2 = p.atan_logders.ae[l], p.atan_logders.ps[l]
    #    print(f1)
    #    print("f1", f1.energies.shape, f1.values.shape)


def oncv_plot(options):
    """
    Plot data with matplotlib. Requires oncvpsp output file.
    """
    cli.customize_mpl(options)

    out_path = _find_oncv_output(options.filepath)

    plotter = OncvPlotter.from_file(out_path)
    #plotter.expose(use_web=True)
    plotter.expose(use_web=False)

    return 0


def oncv_compare(options):
    """
    Compare multiple oncvpsp output fies.
    """
    cli.customize_mpl(options)

    plotter = MultiOncvPlotter.from_files(options.filepaths)

    plotter.plot_atan_logders(show=True, fontsize=12)
    return

    # Plot data
    use_web = False
    #use_web = True
    #import matplotlib.pyplot as plt
    #with plt.xkcd():
    plotter.expose(use_web=use_web)

    return 0


def oncv_json(options):
    """
    Produce a string with the results in a JSON dictionary and exit
    Requires oncvpsp output file.
    """
    out_path = _find_oncv_output(options.filepath)
    onc_parser = OncvParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    # Generate json files with oncvpsp results.
    print(json.dumps(onc_parser.to_dict, indent=-1))
    return 0


def oncv_run(options):
    """
    Run oncvpsp, generate djrepo file, plot results. Requires input file.
    """
    # Select calc_type
    calc_type = dict(nor="non-relativistic",
                     sr="scalar-relativistic",
                     fr="fully-relativistic")[options.rel]

    # Build names of psp8 and djson files from input and relativistic mode.
    in_path = options.filepath
    root, _ = os.path.splitext(in_path)

    # Enforce convention on output files.
    if options.rel == "nor":
        if not root.endswith("_nor"): root += "_nor"
    elif options.rel == "fr":
        if not root.endswith("_r"):
            root += "_r"
            cprint("FR calculation with input file without `_r` suffix. Will add `_r` to output files", "yellow")

    # Build names of output files.
    psp8_path = root + ".psp8"
    djrepo_path = root + ".djrepo"
    out_path = root + ".out"
    if os.path.exists(psp8_path):
        cprint("%s already exists and will be overwritten" % os.path.relpath(psp8_path), "yellow")
    if os.path.exists(djrepo_path):
        cprint("%s already exists and will be overwritten" % os.path.relpath(djrepo_path), "yellow")
    if os.path.exists(out_path):
        cprint("%s already exists and will be overwritten" % os.path.relpath(out_path), "yellow")

    # Build Generator and start generation.
    psgen = OncvGenerator.from_file(in_path, calc_type, workdir=None)
    #print(psgen)
    print(psgen.input_str)
    print("Using oncvpsp executable:\n\t", psgen.executable)
    print(f"Output files produced in directory:\n\t{psgen.workdir}")

    psgen.start()
    retcode = psgen.wait()

    if retcode != 0:
        cprint("oncvpsp returned %s. Exiting" % retcode, "red")
        return retcode

    if psgen.status != psgen.S_OK:
        cprint(f"psgen.status = {psgen.status} != psgen.S_OK", "red")
        if psgen.parser.warnings:
            print(2 * "\n")
            print("List of WARNINGS:")
            for w in psgen.parser.warnings:
                print(w)

        if psgen.parser.errors:
            print(2 * "\n")
            print("List of ERRORS:")
            for e in psgen.parser.errors:
                print(e)

        #return 1

    # Transfer final output file.
    shutil.copy(psgen.stdout_path, out_path)

    # Parse the output file
    onc_parser = OncvParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    # Extract psp8 files from the oncvpsp output and write it to file.
    psp8_str = onc_parser.get_psp8_str()
    with open(psp8_path, "wt") as fh:
        fh.write(psp8_str)

    # Write upf if available.
    upf_str = onc_parser.get_upf_str()
    if upf_str is not None:
        with open(psp8_path.replace(".psp8", ".upf"), "wt") as fh:
            fh.write(upf_str)

    pseudo = Pseudo.from_file(psp8_path)
    if pseudo is None:
        cprint("Cannot parse psp8 file: %s" % psp8_path, "red")
        return 1

    # Initialize and write djson file.
    # FIXME: This part can be removed when we migrate to the new lightweight testing
    # infrastructure with small json files.
    from pseudo_dojo.core.dojoreport import DojoReport
    report = DojoReport.empty_from_pseudo(pseudo, onc_parser.hints, devel=False)
    report.json_write()

    cli.customize_mpl(options)

    # Build the plotter
    plotter = onc_parser.get_plotter()

    plotter.expose(use_web=True)

    return 0


def oncv_gui(options):
    """
    Start a panel web app to generate pseudopotentials.
    """
    import panel as pn
    from abipy.panels.oncvpsp_gui import OncvGui
    from abipy.panels.core import abipanel, get_abinit_template_cls_kwds, AbipyParameterized

    # Load abipy/panel extensions and set the default template
    #tmpl_kwds.update(dict(
    #    sidebar_width=240,
    #    #sidebar_width=280,
    #    #background_color="yellow",
    #))
    abipanel(panel_template=options.panel_template)

    if options.has_remote_server:
        print("Enforcing limitations on what the user can do on the abinit server")
        AbipyParameterized.has_remote_server = options.has_remote_server

    tmpl_cls, tmpl_kwds = get_abinit_template_cls_kwds()

    if options.verbose:
        print("Using default template:", tmpl_cls, "with kwds:\n", pformat(tmpl_kwds), "\n")

    #if options.seaborn:
    # Use seaborn settings.
    import seaborn as sns
    context = "paper"
    context = "notebook"
    sns.set(context=context, style='darkgrid', palette='deep',
            font='sans-serif', font_scale=1, color_codes=False, rc=None)

    def build():
        gui = OncvGui.from_file(os.path.abspath(options.filepath))
        template = tmpl_cls(main=gui.get_panel(), title="Oncvpsp GUI", **tmpl_kwds)
        return template

    # Call pn.serve to serve the multipage app.
    serve_kwargs = cli.get_pn_serve_kwargs(options)

    return pn.serve(build, **serve_kwargs)


def main():

    def str_examples():
        return """\
Usage example:

    oncv.py run H.in         ==> Run oncvpsp input file (scalar relativistic mode).
    oncv.py plot H.out       ==> Use matplotlib to plot oncvpsp results for pseudo H.psp8.
    oncv.py print H.out      ==> Use matplotlib to plot oncvpsp results for pseudo H.psp8.
    oncv.py gui H.out        ==> Run oncvpsp input file (scalar relativistic mode).
    oncv.py gnuplot H.out    ==> Use gnuplot to plot oncvpsp results for pseudo H.psp8.
    oncv.py nbplot H.out     ==> Generate jupyter notebook to plot oncvpsp results.
    oncv.py compare H.out other_H.out  ==> Compare multiple output files.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    def get_copts_parser(multi=False):

        # Parent parser implementing common options.
        p = argparse.ArgumentParser(add_help=False)
        p.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                       help='Verbose, can be supplied multiple times to increase verbosity')

        p.add_argument('--loglevel', default="ERROR", type=str,
                            help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

        if multi:
            p.add_argument('filepaths', nargs="+", help="List of files to compare.")
        else:
            p.add_argument('filepath', default="", help="Path to the input/output file")

        from abipy.core.release import __version__
        p.add_argument('-V', '--version', action='version', version=__version__)

        return p

    copts_parser = get_copts_parser(multi=False)

    # Parent parser for commands supporting MplExpose.
    plot_parser = argparse.ArgumentParser(add_help=False)

    plot_parser.add_argument("-s", "--slide-mode", default=False, action="store_true",
            help="Iterate over figures. Expose all figures at once if not given on the CLI.")
    plot_parser.add_argument("-t", "--slide-timeout", type=int, default=None,
            help="Close figure after slide-timeout seconds (only if slide-mode). Block if not specified.")
    plot_parser.add_argument('-sns', "--seaborn", const="paper", default=None, action='store', nargs='?', type=str,
        help='Use seaborn settings. Accept value defining context in ("paper", "notebook", "talk", "poster"). Default: paper')
    plot_parser.add_argument('-mpl', "--mpl-backend", default=None,
        help=("Set matplotlib interactive backend. "
              "Possible values: GTKAgg, GTK3Agg, GTK, GTKCairo, GTK3Cairo, WXAgg, WX, TkAgg, Qt4Agg, Qt5Agg, macosx."
              "See also: https://matplotlib.org/faq/usage_faq.html#what-is-a-backend."))

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Create the parsers for the sub-commands
    p_run = subparsers.add_parser('run', parents=[copts_parser, plot_parser], help=oncv_run.__doc__)
    p_run.add_argument("--rel", default="sr", help=("Relativistic treatment: `nor` for non-relativistic, "
                       "`sr` for scalar-relativistic, `fr` for fully-relativistic."))

    p_print = subparsers.add_parser('print', parents=[copts_parser], help=oncv_print.__doc__)

    # Create the parsers for the sub-commands
    p_plot = subparsers.add_parser('plot', parents=[copts_parser, plot_parser],
                                   help=oncv_plot.__doc__)

    copts_parser_multi = get_copts_parser(multi=True)
    p_compare = subparsers.add_parser('compare', parents=[copts_parser_multi, plot_parser],
                                      help=oncv_compare.__doc__)

    p_nbplot = subparsers.add_parser('nbplot', parents=[copts_parser], help=oncv_nbplot.__doc__)
    # notebook options.
    p_nbplot.add_argument('-nb', '--notebook', action='store_true', default=False, help="Open file in jupyter notebook")
    p_nbplot.add_argument('--classic-notebook', "-cnb", action='store_true', default=False,
                        help="Use classic jupyter notebook instead of jupyterlab.")
    p_nbplot.add_argument('--no-browser', action='store_true', default=False,
                        help=("Start the jupyter server to serve the notebook "
                              "but don't open the notebook in the browser.\n"
                              "Use this option to connect remotely from localhost to the machine running the kernel"))
    p_nbplot.add_argument('--foreground', action='store_true', default=False,
                        help="Run jupyter notebook in the foreground.")

    parents = [copts_parser, cli.pn_serve_parser(), plot_parser]
    p_gui = subparsers.add_parser('gui', parents=parents, help=oncv_gui.__doc__)

    p_gnuplot = subparsers.add_parser('gnuplot', parents=[copts_parser], help=oncv_gnuplot.__doc__)

    p_json = subparsers.add_parser('json', parents=[copts_parser], help=oncv_json.__doc__)

    # Parse command line.
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

    # Dispatch
    return globals()["oncv_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())

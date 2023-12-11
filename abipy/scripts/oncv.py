#!/usr/bin/env python
"""
Script to generate/analyze/plot ONCVPSP pseudopotentials.
"""
from __future__ import annotations

import sys
import os
import argparse
import shutil
import abipy.tools.cli_parsers as cli

from pprint import pformat
from monty.termcolor import cprint
from abipy.flowtk.pseudos import Pseudo
from abipy.ppcodes.ppgen import OncvGenerator
from abipy.ppcodes.oncv_parser import OncvParser
from abipy.ppcodes.oncv_plotter import OncvPlotter, oncv_make_open_notebook, MultiOncvPlotter


def _find_oncv_output(path: str) -> str:
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


def oncv_notebook(options):
    """
    Generate jupyter notebook to plot data. Requires oncvpsp output file.
    """
    out_path = _find_oncv_output(options.filepath)
    return oncv_make_open_notebook(out_path, foreground=options.foreground,
                                   classic_notebook=options.classic_notebook,
                                   no_browser=options.no_browser)


def oncv_gnuplot(options):
    """
    Plot data with gnuplot.
    """
    out_path = _find_oncv_output(options.filepath)

    # Parse output file.
    onc_parser = OncvParser(out_path).scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not completed. Exiting", "red")
        return 1

    onc_parser.gnuplot()
    return 0


def oncv_print(options) -> int:
    """
    Parse oncvps output and print results to terminal.
    """
    out_path = _find_oncv_output(options.filepath)
    p = OncvParser(out_path).scan()
    if not p.run_completed:
        raise RuntimeError("ocnvpsp output file is not completed")

    print(p.to_string(verbose=options.verbose))
    return 0


def oncv_plot(options) -> int:
    """
    Plot data with matplotlib. Requires oncvpsp output file.
    """
    cli.customize_mpl(options)

    out_path = _find_oncv_output(options.filepath)
    plotter = OncvPlotter.from_file(out_path)
    plotter.plotly_atan_logders().show()
    return 0

    plotter.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                   use_web=options.expose_web, verbose=options.verbose)
    return 0


def oncv_plot_pseudo(options) -> int:
    """
    Plot data with matplotlib. Requires pseudopotential file (UPF2 or pawxml).
    """
    cli.customize_mpl(options)

    pseudo = Pseudo.from_file(options.filepath)
    print(pseudo)

    exposer = "panel" if options.expose_web else "mpl"
    from abipy.tools.plotting import Exposer
    with Exposer.as_exposer(exposer) as e:
        e.add_obj_with_yield_figs(pseudo)

    return 0


def oncv_compare(options) -> int:
    """
    Compare multiple oncvpsp output files.
    """
    cli.customize_mpl(options)

    out_paths = [_find_oncv_output(p) for p in options.filepaths]
    plotter = MultiOncvPlotter.from_files(out_paths)

    # Plot data
    plotter.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                   use_web=options.expose_web, verbose=options.verbose)

    return 0


#@flowtk.flow_main
#def main(options):
#    return build_flow(options)

#def oncv_hints(options):
#    """
#    """
#    from abipy import flowtk
#    from abipy.flowtk.pseudo_works import GsEcutConvWork
#    flow = flowtk.Flow(workdir=options.workdir)
#
#    ecut_list = [35, 40, 45, 50, 55, 60, 65]
#    #ecut_list = [35, 40, 45]
#    work = GsEcutConvWork.from_pseudo(pseudo, ecut_list)
#    flow.register_work(work)
#    return flow


def oncv_run(options):
    """
    Run oncvpsp, generate djrepo file, plot results. Requires oncvps input file.
    """
    # Build names of psp8 and djson files from input and relativistic mode.
    in_path = options.filepath
    root, _ = os.path.splitext(in_path)

    # Enforce convention on output files.
    calc_type = None
    if options.rel == "nor":
        if not root.endswith("_nor"): root += "_nor"

    elif options.rel == "fr":
        if not root.endswith("_r"):
            root += "_r"
            cprint("FR calculation with input file without `_r` suffix. Will add `_r` to output files", "yellow")

    elif options.rel == "from_file":
        calc_type  = "scalar-relativistic"
        if root.endswith("_r"): calc_type = "fully-relativistic"
        if root.endswith("_nor"): calc_type = "non-relativistic"

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

    # Select calc_type
    if calc_type is None:
        calc_type = dict(nor="non-relativistic",
                         sr="scalar-relativistic",
                         fr="fully-relativistic")[options.rel]

    # Build Generator and start generation.
    psgen = OncvGenerator.from_file(in_path, calc_type, workdir=None)

    print(psgen.input_str)
    print("Using executable:\n\t", psgen.executable)
    print(f"Output files produced in directory:\n\t{psgen.workdir}")

    if retcode := psgen.start_and_wait() != 0:
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

    # Transfer final output file.
    shutil.copy(psgen.stdout_path, out_path)

    # Parse the output file
    onc_parser = OncvParser(out_path).scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not completed. Exiting", "red")
        return 1

    # Extract psp8 files from the oncvpsp output and write it to file.
    with open(psp8_path, "wt") as fh:
        fh.write(onc_parser.get_psp8_str())

    # Write UPF2 file if available.
    upf_str = onc_parser.get_upf_str()
    if upf_str is not None:
        with open(psp8_path.replace(".psp8", ".upf"), "wt") as fh:
            fh.write(upf_str)
    else:
        cprint("UPF2 file has not been produced. Use `both` in input file!", "red")

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
    plotter.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                   use_web=options.expose_web, verbose=options.verbose)
    return 0


def oncv_gui(options):
    """
    Start a panel web app to generate pseudopotentials.
    """
    import panel as pn
    from abipy.panels.oncvpsp_gui import OncvGui
    from abipy.panels.core import abipanel, get_abinit_template_cls_kwds, AbipyParameterized

    # Load abipy/panel extensions and set the default template
    abipanel(panel_template=options.panel_template)

    if options.has_remote_server:
        print("Enforcing limitations on what the user can do on the abinit server")
        AbipyParameterized.has_remote_server = options.has_remote_server

    tmpl_cls, tmpl_kwds = get_abinit_template_cls_kwds()

    if options.verbose:
        print("Using default template:", tmpl_cls, "with kwds:\n", pformat(tmpl_kwds), "\n")

    def build():
        gui = OncvGui.from_file(os.path.abspath(options.filepath), plotlyFlag=options.plotly)
        return tmpl_cls(main=gui.get_panel(), title="Oncvpsp GUI", **tmpl_kwds)

    # Call pn.serve to serve the multipage app.
    serve_kwargs = cli.get_pn_serve_kwargs(options)

    return pn.serve(build, **serve_kwargs)


def main():

    def str_examples():
        return """\
Usage example:

    oncv.py run H.in                   ==> Run oncvpsp input file (scalar relativistic mode).
    oncv.py run H_r.in                 ==> Run oncvpsp input file (relativistic mode).
    oncv.py plot H.out                 ==> Use matplotlib to plot oncvpsp results for pseudo H.psp8.
    oncv.py plot H.out -ew             ==> Show matplotlib figures in browser.
    oncv.py gui H.in                   ==> Start a panel web app to generate pseudopotential.
    oncv.py print H.out                ==> Parse oncvps output and print results to terminal.
    oncv.py compare H.out other_H.out  ==> Compare multiple output files (supports -ew option as well)
    oncv.py gnuplot H.out              ==> Use gnuplot to plot oncvpsp results for pseudo H.psp8.
    oncv.py notebook H.out             ==> Generate jupyter notebook to plot oncvpsp results.
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

        return p

    copts_parser = get_copts_parser(multi=False)

    # Parent parser for commands supporting MplExposer.
    plot_parser = argparse.ArgumentParser(add_help=False)
    cli.add_expose_options_to_parser(plot_parser)

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    from abipy.core.release import __version__
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for run command.
    p_run = subparsers.add_parser('run', parents=[copts_parser, plot_parser], help=oncv_run.__doc__)
    p_run.add_argument("--rel", default="from_file", help=("Relativistic treatment: `nor` for non-relativistic, "
        "`sr` for scalar-relativistic, `fr` for fully-relativistic. Default: `from_file` i.e. detected from file"))

    # Subparser for print command.
    p_print = subparsers.add_parser('print', parents=[copts_parser], help=oncv_print.__doc__)

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[copts_parser, plot_parser], help=oncv_plot.__doc__)

    # Subparser for plot command.
    p_plot_pseudo = subparsers.add_parser('plot_pseudo', parents=[copts_parser, plot_parser],
                                          help=oncv_plot_pseudo.__doc__)

    # Subparser for compare command.
    copts_parser_multi = get_copts_parser(multi=True)
    p_compare = subparsers.add_parser('compare', parents=[copts_parser_multi, plot_parser],
                                      help=oncv_compare.__doc__)

    # notebook options.
    p_nb = subparsers.add_parser('notebook', parents=[copts_parser], help=oncv_notebook.__doc__)
    p_nb.add_argument('-nb', '--notebook', action='store_true', default=False, help="Open file in jupyter notebook")
    p_nb.add_argument('--classic-notebook', "-cnb", action='store_true', default=False,
                          help="Use classic jupyter notebook instead of jupyterlab.")
    p_nb.add_argument('--no-browser', action='store_true', default=False,
                          help=("Start the jupyter server to serve the notebook "
                                "but don't open the notebook in the browser.\n"
                                "Use this option to connect remotely from localhost to the machine running the kernel"))
    p_nb.add_argument('--foreground', action='store_true', default=False,
                          help="Run jupyter notebook in the foreground.")

    parents = [copts_parser, cli.pn_serve_parser(), plot_parser]

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', parents=parents, help=oncv_gui.__doc__)

    # Subparser for gnuplot command.
    p_gnuplot = subparsers.add_parser('gnuplot', parents=[copts_parser], help=oncv_gnuplot.__doc__)

    # Subparser for hints command.
    #p_hints = subparsers.add_parser('hints', parents=[copts_parser], help=oncv_hints.__doc__)
    #p_hints.add_argument("pseudo_paths", nargs="+", type=str, help="Pseudopotential path.")
    #p_hints.add_argument("--ecut", type=float, required=True, help="Cutoff energy in Ha.")
    #p_hints.add_argument("-rc", "--vloc-rcut-list", nargs="+", default=None, type=float,
    #                    help="List of cutoff radii for vloc in Bohr.")
    #cli.add_expose_options_to_parser(p_hints)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    cli.set_loglevel(options.loglevel)

    # Use seaborn settings.
    if getattr(options, "seaborn", None):
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)

    # Dispatch
    return globals()["oncv_" + options.command](options)

if __name__ == "__main__":
    sys.exit(main())

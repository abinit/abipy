#!/usr/bin/env python
"""
This script opens one of the output files produced by Abinit (usually in netcdf format but
other files are supported as well). By default the script starts an interactive ipython
session so that one can interact with the file and call its methods.
Alternatively, it is possible to generate automatically a jupyter notebook to execute code.
"""
import sys
import os
import argparse
import subprocess

from pprint import pprint
from monty.os.path import which
from monty.termcolor import cprint
from monty.functools import prof_main
from abipy import abilab


def make_and_open_notebook(options):
    """
    Generate a jupyter notebook and open it in the browser.
    Return system exit code.

    Raise:
        RuntimeError if jupyter is not in $PATH
    """
    import os
    import nbformat
    nbf = nbformat.v4
    nb = nbf.new_notebook()

    nb.cells.extend([
        nbf.new_markdown_cell("## This is an auto-generated notebook for %s" % os.path.relpath(options.filepath)),
        nbf.new_code_cell("""\

%matplotlib notebook
import numpy as np
#import seaborn as sns
#sns.set(context='notebook', style='darkgrid', palette='deep',
#        font='sans-serif', font_scale=1, color_codes=False, rc=None)
from abipy import abilab

"""),
        nbf.new_code_cell("abifile = abilab.abiopen('%s')" % options.filepath)
    ])

    import io, tempfile
    _, nbpath = tempfile.mkstemp(prefix="abinb_", suffix='.ipynb', dir=os.getcwd(), text=True)

    with io.open(nbpath, 'wt', encoding="utf8") as f:
        nbformat.write(nb, f)

    if which("jupyter") is None:
        raise RuntimeError("Cannot find jupyter in PATH. Install it with `pip install`")

    if not options.classic_notebook:
        # Use jupyter-lab instead of classic notebook
        has_jupyterlab = which("jupyter-lab") is not None
        appname = "jupyter-lab" if has_jupyterlab else "jupyter notebook"
    else:
        appname = "jupyter notebook"

    if options.foreground:
        return os.system("%s %s" % (appname, nbpath))
    else:
        fd, tmpname = tempfile.mkstemp(text=True)
        print(tmpname)
        cmd = "%s %s" % (appname, nbpath)
        print("Executing:", cmd, "\nstdout and stderr redirected to %s" % tmpname)
        process = subprocess.Popen(cmd.split(), shell=False, stdout=fd, stderr=fd)
        cprint("pid: %s" % str(process.pid), "yellow")
        return 0


def get_epilog():
    s = """\
======================================================================================================
Usage example:

    abiopen.py FILE          => Open file in ipython shell.
    abiopen.py FILE -p       => Print info on object to terminal.
    abiopen.py FILE -e       => Generate matplotlib figures automatically.
                                Use -sns to activate seaborn settings.
    abiopen.py FILE -eweb    => Generate matplotlib figures, show them in the $BROWSER.
    abiopen.py FILE -ply     => Generate plotly figures automatically. Show them in the $BROWSER.
                                Note that not all FILEs support plotly.
    abiopen.py FILE -pn      => Generate GUI in web BROWSER to interact with FILE
                                Requires panel package (WARNING: still under development!)
    abiopen.py FILE -nb      => Generate jupyter-lab notebook.
    abiopen.py FILE -cnb     => Generate classic jupyter notebook.

where `FILE` is any file supported by abipy/pymatgen e.g. Netcdf files, Abinit input, POSCAR, xsf.
File extensions supported (including zipped files with extension in ".bz2", ".gz", ".z"):
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).

JSON file are supported as well. In this case, abiopen.py tries to reconstruct python objects
assuming JSON document in MSONable format and then invokes ipython with the `data` object.
Use `-e` or `--notebook` or `--panel` to print the JSON dictionary without reconstructing python objects.
======================================================================================================

Table mapping file extension to AbiPy object:

"""
    return s + abilab.abiopen_ext2class_table()


def get_parser(with_epilog=False):
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument("filepath", help="File to open. See table below for the list of supported extensions.")

    # notebook options.
    parser.add_argument('-nb', '--notebook', action='store_true', default=False, help="Open file in jupyter notebook")
    parser.add_argument('--classic-notebook', "-cnb", action='store_true', default=False,
                        help="Use classic jupyter notebook instead of jupyterlab.")
    parser.add_argument('--no-browser', action='store_true', default=False,
                        help=("Start the jupyter server to serve the notebook "
                              "but don't open the notebook in the browser.\n"
                              "Use this option to connect remotely from localhost to the machine running the kernel"))
    parser.add_argument('--foreground', action='store_true', default=False,
                        help="Run jupyter notebook in the foreground.")

    # print option
    parser.add_argument('-p', '--print', action='store_true', default=False, help="Print python object and return.")

    # panel option
    parser.add_argument("-pn", '--panel', action='store_true', default=False,
                        help="Open Dashboard in web browser, requires panel package.")

    parser.add_argument("-pnt", "--panel-template", default="FastList", type=str,
                        help="Specify template for panel dasboard." +
                             "Possible values are: FastList, FastGrid, Golden, Bootstrap, Material, React, Vanilla." +
                             "Default: FastList"
                        )

    # Expose option.
    parser.add_argument('-e', '--expose', action='store_true', default=False,
        help="Open file and generate matplotlib figures automatically by calling expose method.")
    parser.add_argument("-s", "--slide-mode", default=False, action="store_true",
        help="Iterate over figures. Expose all figures at once if not given on the CLI.")
    parser.add_argument("-t", "--slide-timeout", type=int, default=None,
        help="Close figure after slide-timeout seconds (only if slide-mode). Block if not specified.")
    parser.add_argument('-sns', "--seaborn", const="paper", default=None, action='store', nargs='?', type=str,
        help='Use seaborn settings. Accept value defining context in ("paper", "notebook", "talk", "poster"). Default: paper')
    parser.add_argument('-mpl', "--mpl-backend", default=None,
        help=("Set matplotlib interactive backend. "
              "Possible values: GTKAgg, GTK3Agg, GTK, GTKCairo, GTK3Cairo, WXAgg, WX, TkAgg, Qt4Agg, Qt5Agg, macosx."
              "See also: https://matplotlib.org/faq/usage_faq.html#what-is-a-backend."))
    parser.add_argument("-ew", "--expose-web", default=False, action="store_true",
            help='Generate matplotlib plots in $BROWSER instead of X-server. WARNING: Not all the features are supported.')
    parser.add_argument("-ply", "--plotly", default=False, action="store_true",
            help='Generate plotly plots in $BROWSER instead of matplotlib. WARNING: Not all the features are supported.')

    return parser


@prof_main
def main():
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

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    ##############################################################################################
    # Handle meta options i.e. options that set other options.
    # OK, it's not very clean but I haven't find any parse API to express this kind of dependency.
    ##############################################################################################
    if options.plotly: options.expose = True
    if options.expose_web: options.expose = True
    if options.classic_notebook: options.notebook = True

    if options.verbose > 2: print(options)

    if options.mpl_backend is not None:
        # Set matplotlib backend
        import matplotlib
        matplotlib.use(options.mpl_backend)

    if options.seaborn:
        # Use seaborn settings.
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)

    if not os.path.exists(options.filepath):
        raise RuntimeError("%s: no such file" % options.filepath)

    if options.filepath.endswith(".json"):
        return handle_json(options)

    if not options.notebook:
        abifile = abilab.abiopen(options.filepath)

        if options.print:
            # Print object to terminal.
            if hasattr(abifile, "to_string"):
                print(abifile.to_string(verbose=options.verbose))
            else:
                print(abifile)
            return 0

        elif options.expose:
            # Print info to terminal
            if hasattr(abifile, "to_string"):
                print(abifile.to_string(verbose=options.verbose))
            else:
                print(abifile)

            # Generate plots automatically.
            if options.plotly:
                # plotly version
                if hasattr(abifile, "plotly_expose"):
                    abifile.plotly_expose(verbose=options.verbose)
                else:
                    cprint("`%s` does not implement plotly_expose method" % type(abifile), "red")

            elif hasattr(abifile, "expose"):
                # matplotlib version
                abifile.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                               use_web=options.expose_web, verbose=options.verbose)
            else:
                if not hasattr(abifile, "yield_figs"):
                    raise TypeError("Object of type `%s` does not implement (expose or yield_figs methods" % type(abifile))
                from abipy.tools.plotting import MplExpose
                with MplExpose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                               verbose=options.verbose) as e:
                    e(abifile.yield_figs())

            return 0

        elif options.panel:
            import matplotlib
            matplotlib.use("Agg")
            abilab.abipanel()

            if not hasattr(abifile, "get_panel"):
                raise TypeError("Object of type `%s` does not implement get_panel method" % type(abifile))

            app = abifile.get_panel(template=options.panel_template)
            app.show(debug=options.verbose > 0)
            return 0

        # Start ipython shell with namespace
        # Use embed because I don't know how to show a header with start_ipython.
        import IPython
        IPython.embed(header="""
The Abinit file object is associated to the `abifile` python variable.
Use `abifile.<TAB>` to list available methods.
Use e.g. `abifile.plot?` to access docstring and `abifile.plot??` to visualize source.
Use `print(abifile)` to print the object.
""")

    else:
        # Call specialized method if the object is a NotebookWriter
        # else generate simple notebook by calling `make_and_open_notebook`
        cls = abilab.abifile_subclass_from_filename(options.filepath)
        if hasattr(cls, "make_and_open_notebook"):
            if hasattr(cls, "__exit__"):
                with abilab.abiopen(options.filepath) as abifile:
                    return abifile.make_and_open_notebook(foreground=options.foreground,
                                                          classic_notebook=options.classic_notebook,
                                                          no_browser=options.no_browser)
            else:
                abifile = abilab.abiopen(options.filepath)
                return abifile.make_and_open_notebook(foreground=options.foreground,
                                                      classic_notebook=options.classic_notebook,
                                                      no_browser=options.no_browser)
        else:
            return make_and_open_notebook(options)

    return 0


def handle_json(options):
    """Handle JSON file."""

    if options.notebook:
        # Visualize JSON document in jupyter
        cmd = "jupyter-lab %s" % options.filepath
        print("Executing:", cmd)
        process = subprocess.Popen(cmd.split(), shell=False) #, stdout=fd, stderr=fd)
        cprint("pid: %s" % str(process.pid), "yellow")
        return 0

    elif options.panel:
        # Visualize JSON document in panel dashboard
        import json
        import panel as pn
        with open(options.filepath, "rt") as fh:
            d = json.load(fh)
        json_pane = pn.pane.JSON(d, name='JSON', height=300, width=500)
        app = pn.Row(json_pane.controls(jslink=True), json_pane)
        app.show()
        return 0

    else:
        if options.print:
            # Print python object to terminal.
            data = abilab.mjson_load(options.filepath)
            pprint(data, indent=4)
            return 0
        elif options.expose:
            # Pretty-print dict to terminal.
            import json
            with open(options.filepath, "rt") as fh:
                data = json.load(fh)
            pprint(data, indent=4)
            return 0

        data = abilab.mjson_load(options.filepath)
        # Start ipython shell with namespace
        # Use embed because I don't know how to show a header with start_ipython.
        import IPython
        IPython.embed(header="""
The object initialized from JSON (MSONable) is associated to the `data` python variable.
""")


if __name__ == "__main__":
    sys.exit(main())

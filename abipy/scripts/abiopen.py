#!/usr/bin/env python
"""
This script opens an output file produced by Abinit (usually in netcdf format but
other files are supported as well). By default the script starts an interactive ipython
session so that one can interact with the file and call its methods.
Alternatively, it is possible to generate automatically a jupyter notebook to execute code.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import io
import argparse
import tempfile

from monty.os.path import which
from monty.termcolor import cprint
from monty.functools import prof_main
from abipy import abilab


def make_and_open_notebook(options):
    """
    Generate an jupyter notebook and open it in the browser.
    Return system exit code.

    Raise:
        RuntimeError if jupyther is not in $PATH
    """
    import os
    import nbformat
    nbf = nbformat.v4
    nb = nbf.new_notebook()

    nb.cells.extend([
        nbf.new_markdown_cell("## This is an auto-generated notebook for %s" % os.path.relpath(options.filepath)),
        nbf.new_code_cell("""\
from __future__ import print_function, division, unicode_literals, absolute_import
%matplotlib notebook
import numpy as np
#import seaborn as sns
#sns.set(context='notebook', style='darkgrid', palette='deep',
#        font='sans-serif', font_scale=1, color_codes=False, rc=None)
from abipy import abilab\
"""),

    nbf.new_code_cell("abifile = abilab.abiopen('%s')" % options.filepath)
    ])

    import io, tempfile
    _, nbpath = tempfile.mkstemp(prefix="abinb_", suffix='.ipynb', dir=os.getcwd(), text=True)

    with io.open(nbpath, 'wt', encoding="utf8") as f:
        nbformat.write(nb, f)

    if which("jupyter") is None:
        raise RuntimeError("Cannot find jupyter in PATH. Install it with `pip install`")

    if options.foreground:
        return os.system("jupyter notebook %s" % nbpath)
    else:
        fd, tmpname = tempfile.mkstemp(text=True)
        print(tmpname)
        cmd = "jupyter notebook %s" % nbpath
        print("Executing:", cmd)
        print("stdout and stderr redirected to %s" % tmpname)
        import subprocess
        process = subprocess.Popen(cmd.split(), shell=False, stdout=fd, stderr=fd)
        cprint("pid: %s" % str(process.pid), "yellow")


def get_epilog():
    s = """\
Usage example:

    abiopen.py FILE        => Open file in ipython shell.
    abiopen.py FILE -nb    => Generate jupyter notebook.
    abiopen.py FILE -p     => Print info on object to terminal.
    abiopen.py FILE -e     => Generate matplotlib figures automatically.
                              Use -sns to activate seaborn settings.

`FILE` is any file supported by abipy/pymatgen e.g Netcdf files, Abinit input, POSCAR, xsf ...
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).

File extensions supported:
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

    # notebook option
    parser.add_argument('-nb', '--notebook', action='store_true', default=False, help="Open file in jupyter notebook")
    parser.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")

    # print option
    parser.add_argument('-p', '--print', action='store_true', default=False, help="Print python object and return.")

    # expose option.
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

    if options.verbose > 2:
        print(options)

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

    #if options.filepath.endswith(".ipynb"):
    #    import papermill as pm
    #    nb = pm.read_notebook('notebook.ipynb')
    #    import IPython
    #    IPython.embed(header="The Abinit file is bound to the `nb` variable.\n")
    #    return 0

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
            # Generate matplotlib plots automatically.
            if hasattr(abifile, "to_string"):
                print(abifile.to_string(verbose=options.verbose))
            else:
                print(abifile)

            if hasattr(abifile, "expose"):

                abifile.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                               verbose=options.verbose)
            else:
                if not hasattr(abifile, "yield_figs"):
                    raise TypeError("Object of type `%s` does not implement (expose or yield_figs methods" % type(abifile))
                from abipy.tools.plotting import MplExpose
                with MplExpose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                               verbose=options.verbose) as e:
                    e(abifile.yield_figs())

            return 0

        # Start ipython shell with namespace
        # Use embed because I don't know how to show a header with start_ipython.
        import IPython
        IPython.embed(header="""
The Abinit file is bound to the `abifile` variable.
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
                    return abifile.make_and_open_notebook(foreground=options.foreground)
            else:
                abifile = abilab.abiopen(options.filepath)
                return abifile.make_and_open_notebook(foreground=options.foreground)
        else:
            return make_and_open_notebook(options)

    return 0


if __name__ == "__main__":
    sys.exit(main())

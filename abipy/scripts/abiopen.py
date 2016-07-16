#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import io
import argparse
import tempfile

from monty.os.path import which
from abipy import abilab


def make_open_notebook(options):
    """
    Generate an ipython notebook and open it in the browser.
    Return system exit code.

    Raise:
        RuntimeError if jupyther is not in $PATH
    """
    import nbformat
    nbf = nbformat.v4
    nb = nbf.new_notebook()

    nb.cells.extend([
        nbf.new_markdown_cell("# This is an auto-generated notebook for %s" % os.path.relpath(filepath)),
        nbf.new_code_cell("""\
from __future__ import print_function, division, unicode_literals, absolute_import
%matplotlib notebook
#import numpy as np
#import seaborn as sns
from abipy import abilab\
"""),

    nbf.new_code_cell("abifile = abilab.abiopen('%s')" % options.filepath)
    ])

    _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)

    with io.open(nbpath, 'wt', encoding="utf8") as f:
        nbformat.write(nb, f)

    if which("jupyter") is None:
        raise RuntimeError("Cannot find jupyter in PATH. Install it with `pip install`")
    return os.system("jupyter notebook %s" % nbpath)


#from monty.functools import prof_main
#@prof_main
def main():

    def str_examples():
        s = """\
Usage example:
    abiopen.py out_GSR
    abiopen.py out_DDB

File extensions supported:
"""
        return s + abilab.abiopen_ext2class_table()

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + abilab.__version__)

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                     help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('-nb', '--notebook', action='store_true', default=False, help="Open file in jupyter notebook")
    parser.add_argument("filepath", help="File to open. See table below for the list of supported extensions.")

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

    options.filepath = os.path.abspath(options.filepath)
    if not os.path.exists(options.filepath):
        raise RuntimeError("%s: no such file" % options.filepath)

    if not options.notebook:
        # Start ipython shell with namespace
        abifile = abilab.abiopen(options.filepath)
        import IPython
        # Use embed because I don't know how to show a header with start_ipython.
        IPython.embed(header="The Abinit file is bound to the `abifile` variable.\nTry `print(abifile)`")
        #IPython.start_ipython(argv=options.argv,
        #                      user_ns={"abifile": abifile},
        #                      banner="hello",
        #                      banner1="hello1",
        #                      header="hello_header",
        #                      )
        #
    else:
        import daemon
        with daemon.DaemonContext(detach_process=True):
            return make_open_notebook(options)

    return 0


if __name__ == "__main__":
    sys.exit(main())

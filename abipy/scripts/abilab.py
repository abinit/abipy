#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import sys

from IPython.config.loader import Config
# First import the embeddable shell class
from IPython.frontend.terminal.embed import InteractiveShellEmbed

import argparse
import pymatgen as pymatgen
import abipy.data as abidata

from abipy.abilab import *
from abipy.core.release import __version__

def main():
    # Build the main parser.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

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

    try:
        get_ipython
    except NameError:
        nested = 0
        cfg = Config()
        prompt_config = cfg.PromptManager
        prompt_config.in_template = 'In [\\#]: '
        prompt_config.in2_template = '   .\\D.: '
        prompt_config.out_template = 'Out[\\#]: '
    else:
        print("Running nested copies of IPython.")
        print("The prompts for the nested copy have been modified")
        cfg = Config()
        nested = 1

    # Now create an instance of the embeddable shell. The first argument is a
    # string with options exactly as you would type them if you were starting
    # IPython at the system command line. Any parameters you want to define for
    # configuration can thus be specified here.

    from abipy import abilab
    _abi_builtins_ = [n for n in dir(abilab) if not n.startswith("_")]
                     #[n for n in dir(pymatgen) if not n.startswith("_")] + \
    del abilab #, pymatgen

    _abi_builtins_ = sorted(set(_abi_builtins_))

    import textwrap
    banner = textwrap.fill(str(_abi_builtins_), width=70)
    del textwrap

    banner = ("Custom ipython environment for abipy. Useful aliases such as abiopen:\n" + 
              #banner + "\n" + 
              "have been loaded.\n" + 
              "Type _abi_builtins_ to get the complete list."
              )

    ipshell = InteractiveShellEmbed(
        config=cfg,
        banner1=banner,
        exit_msg='Leaving abipy interpreter, back to program.')

    ipshell()
    return 0


if __name__ == "__main__":
    sys.exit(main())



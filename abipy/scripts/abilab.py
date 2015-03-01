#!/usr/bin/env python
from IPython.config.loader import Config

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

# First import the embeddable shell class
from IPython.frontend.terminal.embed import InteractiveShellEmbed

# Now create an instance of the embeddable shell. The first argument is a
# string with options exactly as you would type them if you were starting
# IPython at the system command line. Any parameters you want to define for
# configuration can thus be specified here.
import pymatgen as pymatgen
import abipy.abilab as abilab

#from pymatgen import *
from abipy.abilab import *

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

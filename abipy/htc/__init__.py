"""
Abinit files handling library
"""
from __future__ import print_function, division, unicode_literals

__all__ = []

_mods = [
    'input',
    'abinitfiles',
    'filesfile',
    'variable',
    'inputfile',
    'jobfile',
    'abinitinput',
    'launcher',
]

for _mod in _mods:
    exec('import ' + _mod)
    exec('from ' + _mod + ' import *')
    exec('__all__.extend(' + _mod + '.__all__)')
    exec('del ' + _mod)

#from .input import *
#from .abinitfiles import *
#from .filesfile import *
#from .variable import *
#from .inputfile import *
#from .jobfile import *
#from .abinitinput import *
#from .launcher import *

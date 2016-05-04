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

#for _mod in _mods:
#    exec('from .' + _mod + ' import *')
#    exec('__all__.extend(' + _mod + '.__all__)')
#    exec('del ' + _mod)
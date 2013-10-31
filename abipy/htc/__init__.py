"""
Abinit files handling library
"""

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

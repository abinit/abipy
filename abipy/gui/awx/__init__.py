#import os
__all__ = []

_MODS = [
    "core",
    "apps",
    "grids",
    "dialogs",
    "panels",
    "eggs",
    "func1dframe",
    "threads",
]

#_MODS = [os.path.join(os.path.dirname(__file__), f) for f in _MODS]

for _mod in _MODS:
    exec('import ' + _mod)
    exec('from ' + _mod + ' import *')
    exec('__all__.extend(' + _mod + '.__all__)')
    exec('del ' + _mod)
    #_mod = __import__(_mod)
    #from _mod import *
    #__all__.extend(_mod.__all__)
    #del _mod
   
#del os

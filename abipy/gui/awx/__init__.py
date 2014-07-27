#import os

__all__ = []

_MODS = [
    "core",
    "apps",
    "frames",
    "grids",
    "dialogs",
    "panels",
    "eggs",
    "func1dframe",
    "threads",
    "utils",
]

#_MODS = [os.path.join(os.path.dirname(__file__), f) for f in _MODS]

for _mod in _MODS:
    try:
        exec('import ' + _mod)
        exec('from ' + _mod + ' import *')
        exec('__all__.extend(' + _mod + '.__all__)')
        exec('del ' + _mod)

    except ImportError:
        pass

    #dir = os.path.dirname(os.path.abspath(__file__))
    #self.scripts = [f.replace(".py", "") for f in wildcard.filter(os.listdir(dir))]
    #_mod = __import__(_mod)
    #from _mod import *
    #__all__.extend(_mod.__all__)
    #del _mod
#del os

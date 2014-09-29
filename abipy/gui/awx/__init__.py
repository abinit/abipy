#import os

__all__ = []

_MODS = [
    "buttons.py",
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
        exec('from .' + _mod + ' import *')
        exec('__all__.extend(' + _mod + '.__all__)')
        exec('del ' + _mod)
    except ImportError:
        pass

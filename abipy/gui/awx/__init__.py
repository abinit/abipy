__all__ = []

_MODS = [
    "apps",
    "frames",
    "tools",
    "dialogs",
    "panels",
    "hg",
]

for _mod in _MODS:
    exec('import ' + _mod)
    exec('from ' + _mod + ' import *')
    exec('__all__.extend(' + _mod + '.__all__)')
    exec('del ' + _mod)

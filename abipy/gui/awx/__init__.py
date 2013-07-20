__all__ = []

_MODS = [
    "core",
    "apps",
    "grids",
    "dialogs",
    "panels",
    "egs",
    "func1dframe",
]

for _mod in _MODS:
    exec('import ' + _mod)
    exec('from ' + _mod + ' import *')
    exec('__all__.extend(' + _mod + '.__all__)')
    exec('del ' + _mod)


#import pkgutil
#import inspect
#
#for loader, name, is_pkg in pkgutil.walk_packages(__path__):
#    module = loader.find_module(name).load_module(name)
#
#    for name, value in inspect.getmembers(module):
#        if name.startswith('_'): continue
#
#        globals()[name] = value
#        __all__.append(name)


# flake8: noqa
from pymatgen.core import __version__ as pmg_version
from pymatgen.io.abinit.abitimer import AbinitTimerParser, AbinitTimerSection

<<<<<<< HEAD
=======

>>>>>>> 7d0c08a545701bb30cd2a3ea6b32ccb16f001250
#print(f"{pmg_version=}")
#if pmg_version < '2023.7.10':

try:
    from pymatgen.io.abinit.abitimer import AbinitTimerParserError as AbinitTimerParseError
except ImportError:
    from pymatgen.io.abinit.abitimer import AbinitTimerParseError

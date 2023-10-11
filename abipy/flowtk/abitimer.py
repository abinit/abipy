# flake8: noqa
from pymatgen.core import __version__ as pmg_version
from pymatgen.io.abinit.abitimer import AbinitTimerParser, AbinitTimerSection


if pmg_version < '2023.7.10':
    from pymatgen.io.abinit.abitimer import AbinitTimerParserError as AbinitTimerParseError
else:
    from pymatgen.io.abinit.abitimer import AbinitTimerParseError

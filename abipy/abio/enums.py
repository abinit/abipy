"""
This module defines enumerators associated to important Abinit input variables
"""
from __future__ import annotations

from enum import Enum, IntEnum #, StrEnum


class ChangeEnumStr:
    def __str__(self):
        return str(self.value)


class RUNL(ChangeEnumStr, IntEnum):
    """
    Values of optdriver corresponding to the different run-levels. See defs_basis.F90
    """
    GSTATE     = 0
    RESPFN     = 1
    SCREENING  = 3
    SIGMA      = 4
    NONLINEAR  = 5
    GWR        = 6
    EPH        = 7
    WFK        = 8
    RTTDDFT    = 9
    GWLS       = 66
    BSE        = 99
    LONGWAVE   = 10


class WFK_TASK(ChangeEnumStr, Enum):
    """
    Integer flags defining the task to be performed in wfk_analyze. See defs_basis.F90
    """
    NONE      = 0
    FULLBZ    = 1
    CLASSIFY  = 2
    PAW_AEPSI = 3
    EINTERP   = 4
    DDK       = 5
    DDK_DIAGO = 6
    OPTICS_FULLBZ = 7
    KPTS_ERANGE= 8
    CHECK_SYMTAB = 9


class GWR_TASK(ChangeEnumStr, Enum):  # StrEnum added in 3.11
    """
    String flags defining the task to be performed in the GWR code.
    """
    HDIAGO = "HDIAGO"
    HDIAGO_FULL = "HDIAGO_FULL"
    CC4S = "CC4S"
    CC4S_FULL = "CC4S_FULL"
    G0W0 = "G0W0"
    G0V = "G0V"
    EGEW = "EGEW"
    EGW0 = "EGW0"
    G0EW = "G0EW"
    RPA_ENERGY = "RPA_ENERGY"
    GAMMA_GW = "GAMMA_GW"
    CHI0 = "CHI0"

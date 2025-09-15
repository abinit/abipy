"""
This module defines enumerators associated to important Abinit input variables
"""
from __future__ import annotations

import enum


class EnumMixin:
    """Mixin for enums providing extra capabilities."""

    @classmethod
    def validate(cls, value) -> None:
        """Validate value"""
        values = [member.value for member in cls]
        if value not in values:
            raise ValueError(f"{value=} is not valid. Must be among: {values=}")


class StrEnum(str, enum.Enum):
    """StrEnum were added in version 3.11"""

    def __new__(cls, *args):
        for arg in args:
            if not isinstance(arg, (str, enum.auto)):
                raise TypeError(
                    "Values of StrEnums must be strings: {} is a {}".format(
                        repr(arg), type(arg)
                    )
                )
        return super().__new__(cls, *args)

    def __str__(self):
        return self.value

    # The first argument to this function is documented to be the name of the enum member, not `self`:
    # https://docs.python.org/3.6/library/enum.html#using-automatic-values
    def _generate_next_value_(name, *_):
        return name


class RUNL(EnumMixin, enum.IntEnum):
    """
    Values of optdriver corresponding to the different run-levels. See defs_basis.F90
    """
    GSTATE = 0
    RESPFN = 1
    SCREENING = 3
    SIGMA = 4
    NONLINEAR = 5
    GWR = 6
    EPH = 7
    WFK = 8
    RTTDDFT = 9
    GWLS = 66
    BSE = 99
    LONGWAVE = 10

    def __str__(self):
        return str(self.value)


class WFK_TASK(EnumMixin, enum.IntEnum):
    """
    Integer flags defining the task to be performed in wfk_analyze. See defs_basis.F90
    """
    NONE = 0
    FULLBZ = 1
    CLASSIFY = 2
    PAW_AEPSI = 3
    EINTERP = 4
    DDK = 5
    DDK_DIAGO = 6
    OPTICS_FULLBZ = 7
    KPTS_ERANGE = 8
    CHECK_SYMTAB = 9

    def __str__(self):
        return str(self.value)


class GWR_TASK(EnumMixin, StrEnum):  # StrEnum added in 3.11
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

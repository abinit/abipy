"""
Structures used in the VASP GWR paper
Table I of https://journals.aps.org/prb/pdf/10.1103/PhysRevB.75.235102
"""
from __future__ import annotations

from abipy.core.structure import Structure

# Chemical formula to (crystal_type, species, a_ang)
_SYMB_TO_DATA = {
    "Si":   ("zinc-blende", ("Si", "Si"), 5.430),
    "GaAs": ("zinc-blende", ("Ga", "As"), 5.648),
    "SiC":  ("zinc-blende", ("Si", "C"),  4.350),
    "ZnO":  ("zinc-blende", ("Zn", "O"),  4.580),
    "C":    ("zinc-blende", ("C", "C"),   3.567),
    "BN":   ("zinc-blende", ("B", "N"),   3.615),
    "MgO":  ("rock-salt",   ("Mg", "O"),  4.213),
    "LiF":  ("rock-salt",   ("Li", "F"),  4.010),
}


def get_gwr_structure(symbol: str) -> Structure:
    """Return Structure from symbol"""
    stype, species, a = _SYMB_TO_DATA[symbol]
    if stype == "zinc-blende":
        return Structure.zincblende(a, species)
    if stype == "rock-salt":
        return Structure.rocksalt(a, species)

    raise ValueError(f"Invalid {stype=}")


if __name__ == "__main__":
    import numpy as np
    for symbol in _SYMB_TO_DATA:
        structure = get_gwr_structure(symbol)
        print("formula:", structure.formula, "a:", structure.lattice.abc[0] * np.sqrt(2))
"""
This module defines convenience types for type hinting purposes.
It extends the types provided by pymatgen with Abipy-specific ones.
"""
from __future__ import annotations

import numpy as np

from typing import TYPE_CHECKING, Any, Union, Sequence
from pymatgen.util.typing import PathLike


if TYPE_CHECKING:  # needed to avoid circular imports
    from matplotlib.pyplot import Axes
    from matplotlib.figure import Figure
    from abipy.core.kpoints import Kpoint
else:
    Axes = Any
    Figure = Any
    Kpoint = Any

VectorLike = Union[Sequence[float], np.ndarray]

# matplotlib objects
#AxList = list[Axes]

# Abipy objects
KptLike = Union["Kpoint", VectorLike]
KptSelect = Union[int, "Kpoint", "VectorLike"]


"""
This module defines convenience types for type hinting purposes.
It extends the types provided by pymatgen with Abipy-specific ones.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union

import numpy as np

if TYPE_CHECKING:
    # needed to avoid circular imports
    from matplotlib.figure import Figure
    from matplotlib.pyplot import Axes

    from abipy.core.kpoints import Kpoint
else:
    Axes = Any
    Figure = Any
    Kpoint = Any

VectorLike = Union[Sequence[float], np.ndarray]
IVectorLike = Union[Sequence[int], np.ndarray]

# matplotlib objects
AxList = list[Axes]

# Abipy objects
KptLike = Union["Kpoint", VectorLike]

KptSelect = Union[int, "Kpoint", "VectorLike"]

GvecSelect = Union[int, IVectorLike]

PathLike = Union[str, Path]

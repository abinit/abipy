"""
This module defines convenience types for type hinting purposes.
It extends the types provided by pymatgen with Abipy-specific ones.
"""
from pymatgen.util.typing import *

if TYPE_CHECKING:  # needed to avoid circular imports
    from matplotlib.pyplot import Axes
    from matplotlib.figure import Figure
    from abipy.core.kpoints import Kpoint
else:
    Axes = Any
    Figure = Any
    Kpoint = Any

# matplotlib objects
AxList = list[Axes]

# Abipy objects
KptLike = Union["Kpoint", VectorLike]
KptSelect = Union[int, "Kpoint", "VectorLike"]


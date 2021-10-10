"""
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, Any
from pydantic import Field
from .base_models import AbipyModel
from .structure_models import StructureData
from .pseudos_models import PseudoSpecs


class BaseGsProcotol(AbipyModel, ABC):
    """
    A protocol contains metaparameters used to generate input files
    """
    in_structure_data: StructureData

    pseudo_specs: PseudoSpecs

    ##################
    # Input parameters
    ##################
    accuracy: str = "normal"

    kppa: int = Field(1000, description="Defines the sampling used for the SCF run. Defaults to 1000 if not given.")

    spin_mode: str = Field("unpolarized", description="Spin polarization")

    charge: float = Field(0.0, description="Electronic charge added to the unit cell.")

    smearing: str = Field("fermi_dirac:0.1 eV", description="Smearing technique.")

    extra_abivars: Dict[str, Any] = Field(None, description="")

    # validators
    #_normalize_name = validator('name', allow_reuse=True)(normalize)
    #_normalize_name = reuse_validator('name')(normalize)

    #@root_validator
    #def check_pseudos(cls, values):
    #    extra_abivars = values.get("extra_abivars")
    #    return values

    #@root_validator
    #def check_pseudos(cls, values):
    #    extra_abivars = values.get("extra_abivars")
    #    return values

    #@root_validator
    #def check_extra_abivars(cls, values):
    #    extra_abivars = values.get("extra_abivars")
    #    return values

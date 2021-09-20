from __future__ import annotations

#import os

#from typing import Dict, Tuple
from abc import ABC, abstractmethod
from pydantic import BaseModel, Field, validator
from pymatgen.io.abinit.pseudos import PseudoTable
from abipy.flowtk.psrepos import get_repo_from_name


class _PseudosProvider(BaseModel, ABC):

    @abstractmethod
    def get_pseudos(self) -> PseudoTable:
        """Return PseudoPotential Table."""


class PseudoSpecs(_PseudosProvider):
    """
    WIP: Model with the parameters needed to retrieve a PseudoDojo Table
    """

    repo_name: str = Field(..., description="Name of the pseudopotential repository.")
    table_accuracy: str = Field(..., description="Accuracy of the table: normal or stringent")

    ps_generator: str = Field(..., description="Pseudopotential type e.g NC or PAW.")
    xc_name: str = Field(..., description="XC name e.g. LDA, PBE, PBEsol...")
    relativity_type: str = Field(..., description="Scalar-relativistic (SR) or fully-relativistic (FR)")
    project_name: str = Field(..., description="Pseudopotential provider e.g. PD for PseudoDojo.")
    version: str = Field(..., description="Version of the table.")

    @classmethod
    def from_repo_name(cls, repo_name: str, table_accuracy: str = "standard") -> PseudoSpecs:
        """
        Build the object from a repository name and the table_accuracy.
        Note that a PseudoPotential repository may provide different set of pseudos for a given ``table_accuracy``.
        """
        repo = get_repo_from_name(repo_name)
        data = dict(repo_name=repo_name, table_accuracy=table_accuracy)
        attr_list = ["ps_generator", "xc_name", "relativity_type", "project_name", "version"]
        data.update({aname: getattr(repo, aname) for aname in attr_list})

        return cls(**data)

    @validator('table_accuracy')
    def validate_table_accuracy(cls, value):
        """Validate table_accuracy"""
        allowed_values = ("standard", "stringent")
        if value not in allowed_values:
            raise ValueError(f"{value} not in {allowed_values}")

        return value

    @validator('relativity_type')
    def validate_relativity_type(cls, value):
        """Validate relativity_type"""
        allowed_values = ("SR", "FR")
        if value not in allowed_values:
            raise ValueError(f"{value} not in {allowed_values}")

        return value

    def get_pseudos(self) -> PseudoTable:
        """
        Return the PseudoTable associated to the initial specs.
        """
        repo = get_repo_from_name(self.repo_name)
        if not repo.is_installed():
            print(f"Warning: Could not find {self.repo_name} installed in {repo.dirpath}.\n"
                  f"Will try to install it at runtime...")
            repo.install(verbose=1)

        # Note that the repo instance keeps an internal cache of tables.
        return repo.get_pseudos(table_accuracy=self.table_accuracy)


#class PseudoSpecs(_PseudosProvider):

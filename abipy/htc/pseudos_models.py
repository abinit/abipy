from __future__ import annotations

from typing import Dict, Tuple

from abc import ABC, abstractmethod
from pydantic import BaseModel, Field, validator
from pymatgen.io.abinit.pseudos import PseudoTable
from abipy.flowtk.psrepos import get_repo_with_name


class _PseudosProvider(BaseModel, ABC):

    @abstractmethod
    def get_pseudos(self) -> PseudoTable:
        """Return PseudoPotential Table."""


_PSEUDOTABLES_CACHE: Dict[Tuple[str, str], PseudoTable] = {}


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

        repo = get_repo_with_name(repo_name)
        data = dict(repo_name=repo_name, table_accuracy=table_accuracy)
        attr_list = ["ps_generator", "xc_name", "relativity_type", "project_name", "version"]
        data.update({aname: getattr(repo, aname) for aname in attr_list})

        return cls(**data)

    @validator('table_accuracy')
    def validate_table_accuracy(cls, value):
        known_accuracies = ("standard", "stringent")
        if value not in known_accuracies:
            raise ValueError(f"{value} not in {known_accuracies}")

        return value

    def get_pseudos(self) -> PseudoTable:
        key = (self.repo_name, self.table_accuracy)
        if key in _PSEUDOTABLES_CACHE:
            return _PSEUDOTABLES_CACHE[key]

        repo = get_repo_with_name(self.repo_name)
        repos_root = "~/.abinit/pseudos"  # FIXME
        if not repo.is_installed(repos_root):
            repo.install(repos_root, verbose=0)

        # Build PseudoTable and cache it.
        pseudos = repo.get_pseudos(repos_root, table_accuracy=self.table_accuracy)
        _PSEUDOTABLES_CACHE[key] = pseudos
        return pseudos

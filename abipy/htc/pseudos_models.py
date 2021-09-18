from __future__ import annotations

#from typing import List

from abc import ABC, abstractmethod
from pydantic import BaseModel, Field, validator
from pymatgen.io.abinit.pseudos import PseudoTable


class _PseudosProvider(BaseModel, ABC):

    @abstractmethod
    def get_pseudos(self):
        """Return PseudoPotential Table."""


#_PSEUDOTABLES_CACHE = {}


class PseudoDojoSpecs(_PseudosProvider):
    """
    WIP: Model with the parameters needed to retrieve a PseudoDojo Table
    """

    #name: str = Field(..., description="Name of the table")

    #abips_dirname = Field(..., description="Pseudopotential type e.g NC or PAW.")

    #pp_type: str = Field(..., description="Pseudopotential type e.g NC or PAW.")
    #project_name: str = Field(..., description="Pseudopotential provider e.g. PD for PseudoDojo.")
    #xc_name: str = Field(..., description="XC name e.g. LDA, PBE, PBEsol...")
    #relativity_type: str = Field(..., description="Scalar-relativistic (SR) or fully-relativistic (FR)")
    #accuracy: str = Field(..., description="Accuracy of the table: normal or stringent")
    #version: str = Field(..., description="Version of the table.")

    @classmethod
    def from_abips_dirname(cls, abips_dirname: str, accuracy: str = "standard") -> PseudoDojoSpecs:

        #from repo in ALL_REPOS:
        #    if repo.dirname == abips_dirname: break
        #else:
        #    raise ValueError(f"Couldn't find {abips_dirname} in the list of registered repos")

        data = dict()
        return cls(**data)

    #@validator('accuracy')
    #def validate_accuracy(cls, value):
    #    known_accuracies = ("standard", "stringent")
    #    if value not in known_accuracies:
    #        raise ValueError(f"{value} not in {known_accuracies}")
    #    return value

    def get_pseudos(self):
        raise NotImplementedError("get_pseudos")
        #_PSEUDOTABLES_CACHE = {}
        #key = (self.abips_name, self.accuracy)
        #if  key in _PSEUDOTABLES_CACHE:
        #    return _PSEUDOTABLES_CACHE[key]

        #repo = Repo.from_dirname(self.abips_dirname)
        #if not repo.is_installed(abips_home):
        #    repo.install(self.abips_home, verbose=0)

        #pseudos = repo.build_pseudotable(accuracy=self.accuracy)

        # Build PseudoTable and cache it.
        #PseudoTable.from_abips_dirname(self.abips_dirname)
        #_PSEUDOTABLES_CACHE[key] = pseudos
        #return pseudos


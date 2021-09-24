from __future__ import annotations

import os

from typing import List  # Dict, Tuple
from abc import ABC, abstractmethod
from pydantic import BaseModel, Field, validator
#from bson.objectid import ObjectId
from pymatgen.io.abinit.pseudos import PseudoTable
from abipy.flowtk.psrepos import get_repo_from_name
from .base_models import AbipyModel  #, MongoConnector


#class OncvpspModel(AbipyModel):
#
#    psp8_md5: str = Field(..., description="")
#    psp8_string: str = Field(..., description="")
#    upf_md5: str = Field(..., description="")
#    upf_string: str = Field(..., description="")
#    psml_md5: str = Field(..., description="")
#    psml_string: str = Field(..., description="")
#
#    psgen_input_string: str = Field(..., description="")
#    psgen_output_string: str = Field(..., description="")
#
#    @classmethod
#    def from_dirpath(cls, dirpath: str) -> OncvpspModel:
#        from abipy.flowtk.psrepos import md5_for_filepath
#        dirpath = os.path.abspath(dirpath)
#        paths = [os.path.join(dirpath, b) for b in os.listdir(dirpath)]
#        data = {}
#
#        def find_abspath_byext(ext: str):
#            cnt, ret_apath = 0, None
#            for p in paths:
#                if p.endswith(ext):
#                    ret_apath = p
#                    cnt += 1
#
#            if cnt == 1:
#                return ret_apath
#            if cnt > 1:
#                raise RuntimeError(f"Found multiple files with extension {ext} in {paths}")
#
#            raise RuntimeError(f"Cannot find file with extension {ext} in {paths}")
#
#        apath = find_abspath_byext(".in")
#        data["psgen_input_string"] = open(apath, "rt").read()
#        apath = find_abspath_byext(".out")
#        data["psgen_output_string"] = open(apath, "rt").read()
#
#        for fmt in ["psp8", "upf", "psml"]:
#            apath = find_abspath_byext("." + fmt)
#            data[f"{fmt}_md5"] = md5_for_filepath(apath)
#            with open(apath, "rt") as fh:
#                data[f"{fmt}_input_string"] = fh.read()
#
#        return cls(**data)


class _PseudosProvider(BaseModel, ABC):

    @abstractmethod
    def get_pseudos(self) -> PseudoTable:
        """Build and return the PseudoTable associated the model."""


#class PseudoSpecs(_PseudosProvider):
#
#    mongo_connector: MongoConnector = Field(..., description="")
#    oid_list: List[ObjectId] = Field(..., description="")
#
#    def get_pseudos(self) -> PseudoTable:
#        """Build and return the PseudoTable associated the model."""
#        collection = self.mongo_connector.get_collection()
#        query = {"_id": {"$in": self.oid_list}}
#        docs = collection.find(query)
#        if not docs or len(docs) != len(self.oid_list):
#            raise RuntimeError("")
#
#        pseudos = []
#        for doc in docs:
#            s = doc["psp8_string"]
#            filepath = "foo.psp8"
#            with open(filepath, "wt") as fh:
#                fh.write(s)
#            pseudos.append(Pseudo.from_file(filepath))
#        return PseudoTable(pseudos)


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
        """
        Validate table_accuracy
        """
        allowed_values = ("standard", "stringent")
        if value not in allowed_values:
            raise ValueError(f"{value} not in {allowed_values}")

        return value

    @validator('relativity_type')
    def validate_relativity_type(cls, value):
        """
        Validate relativity_type
        """
        allowed_values = ("SR", "FR")
        if value not in allowed_values:
            raise ValueError(f"{value} not in {allowed_values}")

        return value

    def get_pseudos(self) -> PseudoTable:
        """
        Build and return the PseudoTable associated the model.
        """
        repo = get_repo_from_name(self.repo_name)
        if not repo.is_installed():
            print(f"Warning: Could not find {self.repo_name} installed in {repo.dirpath}.\n"
                  f"Will try to install it at runtime...")
            repo.install(verbose=1)

        # Note that the repo instance keeps an internal cache of tables.
        return repo.get_pseudos(table_accuracy=self.table_accuracy)

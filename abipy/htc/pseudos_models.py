from __future__ import annotations

#from typing import List

from abc import ABC, abstractmethod
from pydantic import BaseModel, Field

from pymatgen.io.abinit.pseudos import Pseudo


class _PseudosProvider(BaseModel, ABC):

    @abstractmethod
    def get_pseudos(self):
        """Return PseudoPotential Table."""


#class PseudosFromMongoCollection(_PseudosProvider):
#
#    #collection_name: str = Field(..., description="Name of the MongoDB collection")
#    mongo_connector: MongoConnector = Field(..., description="Name of the MongoDB collection")
#
#    md5_list = List[str]
#
#    def get_pseudos(self):
#        collection = self.mongo_connector.get_collection()
#

class PseudoDojoSpecs(_PseudosProvider):
    """
    WIP: Model with the parameters needed to retrieve a PseudoDojo Table
    """

    name: str = Field(..., description="Name of the table")

    #version: str = Field(..., description="Version of the table")
    #soc_type: str = Field(..., description="Scalar-relativistic or fully-relativistic.")
    #xc_type: str

    @classmethod
    def from_table_name(cls, table_name: str) -> PseudoDojoSpecs:
        return cls(name=table_name)

    #from pydantic import ValidationError, validator
    #@validator('name')
    #def name_must_contain_space(cls, v):
    #    if ' ' not in v:
    #        raise ValueError('must contain a space')
    #    return v.title()

    def get_pseudos(self):
        raise NotImplementedError("get_pseudos")
        #from pseudodojo import ...
        #return PseudoTable
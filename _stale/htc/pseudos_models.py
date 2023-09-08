"""
"""
from __future__ import annotations

#from typing import List  # Dict, Tuple
from abc import ABC, abstractmethod
from pprint import pformat
from pydantic import BaseModel, Field, validator
from bson.objectid import ObjectId
from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
from abipy.core.structure import Structure
from abipy.flowtk.psrepos import get_repo_from_name
from .base_models import AbipyModel, GfsFileDesc, GfsDesc, MongoConnector


#class PseudoSpecs(_PseudosProvider):
class _PseudosProvider(BaseModel, ABC):
    """
    Abstract base class enforcing the implementation of ``get_pseudos``
    """

    ps_type: str = Field(..., description="Pseudopotential type e.g. NC or PAW.")

    ps_generator: str = Field(..., description="Name of the Pseudopotential generator.")

    xc_name: str = Field(..., description="Name of the XC functional e.g. LDA, PBE, PBEsol.")

    relativity_type: str = Field(..., description="Scalar-relativistic (SR) or fully-relativistic (FR)")

    @abstractmethod
    def get_pseudos(self) -> PseudoTable:
        """Build and return the PseudoTable associated the model."""


#class PseudoSpecsWithRepo(_PseudosProvider):
class PseudoSpecs(_PseudosProvider):
    """
    This model stores the parameters needed to retrieve the PseudoPotential Table from a repository.
    The `repo_name` defines the full set of versioned pseudos while `table_name` defines the list of pseudos.
    included in the final Table.
    Note that the PseudoDojo/JTC project may provide multiple pseudos for the same atom.
    This is the reason why the final table with one pseudo per element is uniquely defined by
    the tuple (repo_name, table_name).

    .. code-block::

        specs = PseudoSpecs.from_repo_table_name("ONCVPSP-PBEsol-SR-PDv0.4", "standard")
    """

    repo_name: str = Field(..., description="Name of the pseudopotential repository.")

    table_name: str = Field(..., description="Name of the table e.g. normal, stringent, MD ...")

    project_name: str = Field(..., description="PD for PseudoDojo, JTH for JTH.")

    version: str = Field(..., description="Version of the repository.")

    @classmethod
    def from_repo_table_name(cls, repo_name: str, table_name: str) -> PseudoSpecs:
        """
        Build the object from a repository name and ``table_name``.
        """
        repo = get_repo_from_name(repo_name)
        data = dict(repo_name=repo_name, table_name=table_name)
        attr_list = ["ps_generator", "ps_type", "xc_name", "relativity_type", "project_name", "version"]
        data.update({aname: getattr(repo, aname) for aname in attr_list})

        return cls(**data)

    def __str__(self) -> str:
        return pformat(self.dict())
        #lines = []
        #return "\n".join(lines)

    @validator("relativity_type")
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
        Build and return the PseudoTable associated to the specs.
        """
        repo = get_repo_from_name(self.repo_name)
        if not repo.is_installed():
            print(f"Warning: Could not find {self.repo_name} installed in {repo.dirpath}.\n"
                  f"Will try to install it at runtime...")
            repo.install(verbose=1)

        # Note that the repo instance keeps an internal cache of tables.
        return repo.get_pseudos(table_name=self.table_name)


class MongoDbPseudos(_PseudosProvider):
    """
    A rather specialized provider that takes pseudopotentials by IDs from a MongoDB collection.
    This approach is used for testing/validanting "unofficial" pseudos.
    """

    oid_list: List[ObjectId] = Field(..., description="List of ids in collection")

    mng_connector: MongoConnector = Field(..., description="")

    #pseudo_model_cls =

    def get_pseudos(self) -> PseudoTable:
        """
        Build and return the PseudoTable associated to the specs.
        """
        collection = self.mng_connector.get_collection()
        #pseudo_model = pseudo_model_cls(**d)
        pseudos = []
        for oid in self.oid_list:
            d = collection.find_one({"_id": oid}, {"_id": 0})
            if d is None:
                raise ValueError(f"Cannot find pseudo with id {oid} in collection: {collection}")
            #pseudo_model = pseudo_model_cls(**d)
            p = pseudo_model.get_pseudo_from_db(self.mng_connector)
            pseudos.append(p)

        return PseudoTable(pseudos)



class BasePseudoModel(AbipyModel, ABC):
    """
    Abstract base class for a model associated to a single pseudopotential.
    Used for fetching the pseudopotential file from a MongoDB collection.
    """

    atomic_symbol: str = Field(..., description="")

    ps_type: str = Field(..., description="Pseudopotential type e.g. NC or PAW.")

    ps_generator: str = Field(..., description="Name of the Pseudopotential generator.")

    xc_name: str = Field(..., description="Name of the XC functional e.g. LDA, PBE, PBEsol.")

    relativity_type: str = Field(..., description="Scalar-relativistic (SR) or fully-relativistic (FR)")

    #psgen_input_string: str = Field(..., description="")

    #psgen_output_string: str = Field(..., description="")

    #output_gsfd: GfsFileDesc = Field(..., description="Pointer to the main output file in GridFs")

    @abstractmethod
    def get_pseudo_from_db(self, mng_connector: MongoConnector) -> Pseudo:
        """Build and return a Pseudo object from the model."""



class OncvpspModel(BasePseudoModel):
    """
    This model stores the pseudopotentials generated by oncvpsp with the corresponding
    input and output files.
    """

    psp8_gfsd: GfsFileDesc = Field(..., description="Pointer to the psp8 file in GridFs")
    psp8_md5: str = Field(..., description="md5 checksum for the psp8 file")

    upf_md5: str = Field(..., description="md5 checksum for the upf file")
    upf_gsfd: GfsFileDesc = Field(..., description="Pointer to the UPF file in GridFs")

    psml_md5: str = Field(..., description="md5 checksum for the psml file")
    psml_gsfd: GfsFileDesc = Field(..., description="Pointer to the PSML file in GridFs")

    #@validator("atomic_symbol")
    #def validate_atomic_symbol(cls, value):

    @classmethod
    def from_dirpath(cls, dirpath: str, mng_connector: MongoConnector) -> OncvpspModel:
        """
        Build a the model from a directory containing the output files produced by oncvpsp.
        """
        from abipy.flowtk.psrepos import md5_for_filepath
        dirpath = os.path.abspath(dirpath)
        paths = [os.path.join(dirpath, b) for b in os.listdir(dirpath)]

        def find_abspath_byext(ext: str):
            cnt, ret_apath = 0, None
            for p in paths:
                if p.endswith(ext):
                    ret_apath = p
                    cnt += 1

            if cnt == 1:
                return ret_apath
            if cnt > 1:
                raise RuntimeError(f"Found multiple files with extension {ext} in {paths}")

            raise RuntimeError(f"Cannot find file with extension {ext} in {paths}")

        data = {}
        apath = find_abspath_byext(".in")
        input_string = open(apath, "rt").read()
        #data["psgen_input_string"] = input_string
        apath = find_abspath_byext(".out")
        #data["psgen_output_string"] = open(apath, "rt").read()
        #data[output_gfsd"] = mng_connector.gfs_put_filepath(apath)

        data["atomic_symbol"] = "foo"

        for fmt in ["psp8", "upf", "psml"]:
            apath = find_abspath_byext("." + fmt)
            data[f"{fmt}_md5"] = md5_for_filepath(apath)
            #with open(apath, "rt") as fh:
            #    data[f"{fmt}_input_string"] = fh.read()

            data[f"{fmt}_gfsd"] = mng_connector.gfs_put_filepath(apath)

        return cls(**data)

    def get_pseudo_from_db(self, mng_connector: MongoConnector) -> Pseudo:
        """
        Build and return a Pseudo instance from the model.
        """
        filepath = mng_connector.mktmp_filepath(self.psp8_gfsd)
        return Pseudo.from_file(filepath)


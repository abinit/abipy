"""
Pydantic Models.
"""
import json

from pprint import pprint
#from typing import List, Optional
from monty import design_patterns
from pydantic import BaseModel, Field, BaseConfig
from bson.objectid import ObjectId

from monty.json import MontyEncoder, MontyDecoder
from abipy.core.structure import Structure
from abipy.electrons.ebands import ElectronBands


def monty_json_dumps(obj, **kwargs):
    #print("in monty json dumps with type:", type(obj))
    return json.dumps(obj, cls=MontyEncoder, **kwargs)


def monty_json_loads(string, **kwargs):
    #print("in monty json loads with string", string)
    return json.loads(string, cls=MontyDecoder, **kwargs)


def monty_trick(obj):
    return json.loads(monty_json_dumps(obj))



class MongoConnector(BaseModel):
    """
    Model with the parameters used to connect to the MongoDB server and the name of
    the (default) collection used to perform CRUD operations.
    """

    host: str = Field(..., description="Host address e.g. 'localhost'")

    port: int = Field(..., description="MongoDB server port. e.g. 27017")

    db_name: str = Field("abipy", description="Name of the MongoDB database")

    collection_name: str = Field(..., description="Name of the collection")

    user: str = Field(None, description="User name. Default: None")

    password: str = Field(None, description="password for authentication. Default: None")

    @classmethod
    def for_localhost(cls, collection_name: str, port: int = 27017):
        return cls(host="localhost", port=port, collection_name=collection_name)

    def get_client(self):
        """
        Establish a connection. Return MongoClient
        """
        from pymongo import MongoClient
        from pymongo.errors import ConnectionFailure

        # TODO Handle reconnection, cache client ...
        #self._client =
        client = MongoClient(host=self.host, port=self.port)

        try:
            # The ping command is cheap and does not require auth.
            client.admin.command('ping')
        except ConnectionFailure:
            print("Server not available")

        return client

    def get_collection(self, collection_name=None):
        """
        Returns MongoDb collection
        """
        client = self.get_client()
        db = client[self.db_name]

        # Authenticate if needed
        if self.user and self.password:
            db.autenticate(self.user, password=self.password)

        collection_name = collection_name or self.collection_name
        return db[collection_name]


class MongoModel(BaseModel):
    """
    Base class providing tools to serialize Abipy/Pymatgen objects
    supporting the MSONable protocol.
    It also provide tiny wrappers around the MongoDB API.
    """

    class Config:

        # Here we register specialized encoders
        # to convert MSONable objects to JSON string.
        # subclasses will inherit these encoders.
        json_encoders = {
            Structure: lambda o: monty_trick(o),
            ElectronBands: lambda o: monty_trick(o),
        }

        # This is needed to be able to use MSONable AbiPy objects
        # such as ElectronBands that do not inherit from BaseModel
        arbitrary_types_allowed = True

        #json_loads = monty_json_loads
        #json_dumps = monty_json_dumps

    @classmethod
    def from_mongo_oid(cls, oid, collection):
        """Return a model instance for the ObjectId oid and the collection."""
        oid = ObjectId(oid)
        data = collection.find_one({'_id': oid})
        if data is None: return data
        data.pop("_id")
        data = MontyDecoder().process_decoded(data)

        return cls(**data)
        #return cls(**dict(data, id=id))

    def mongo_insert(self, collection):
        """Insert the model in collection. Return ObjectId."""
        doc = json.loads(self.json())
        #print("inserting doc:", doc)
        return collection.insert_one(doc).inserted_id

    def mongo_full_update_oid(self, oid, collection):
        """Perform a full update of the model given the ObjectId in the collection."""
        oid = ObjectId(oid)
        old_doc = collection.find_one({'_id': oid})
        if old_doc is None:
            raise RuntimeError(f"Cannot find document with oid: {oid}")

        new_doc = json.loads(self.json())
        old_doc.update(new_doc)

        collection.replace_one({"_id": oid}, new_doc, upsert=False)
        #collection.update_one({'_id': oid}, new_doc)


class PseudoDojoSpecs(BaseModel):
    """
    WIP: Model with the parameters needed to retrieve a PseudoDojo Table
    """

    name: str = Field(..., description="Name of the table")

    #version: str = Field(..., description="Version of the table")

    #soc_type: str = Field(..., description="Scalar-relativistic or fully-relativistic.")

    #xc_type: str

    @classmethod
    def from_table_name(cls, table_name: str):
        return cls(name=table_name)

    #def get_pseudos(self, db):
    #    from pseudodojo import ...
    #    return PseudoTable


class FlowModel(MongoModel):
    """
    This model
    """

    flow_status: int = 0

    workdir: str = None

    pseudos_specs: PseudoDojoSpecs = Field(..., description="The input structure.")

    #flow_params: dict

    #flow_created_at: datetime = Field(
    #    description="Timestamp for when this material document was first created",
    #    default_factory=datetime.utcnow,
    #)
    #
    #flow_last_updated: datetime = Field(
    #    description="Timestamp for the most recent calculation update for this property",
    #    default_factory=datetime.utcnow,
    #)

    @classmethod
    def find_runnable_oid_models(cls, collection, limit=0):
        """Return list of models that are ready to run."""
        cursor = collection.find({"flow_status": 0}, limit=limit)
        if cursor is None: return None
        items = []
        for doc in cursor:
            oid = doc.pop("_id")
            doc = MontyDecoder().process_decoded(doc)
            items.append((oid, cls(**doc)))
        return items

    def build_flow(self, workdir, manager):
        """
        API used by the AbiPy Worker to build a Flow from the model.
        Wraps _build_flow implented in the subclass.
        """
        flow = self._build_flow(workdir, manager)
        # Set the workdir of the Flow in the model.
        self.workdir = workdir
        return flow

    def postprocess_flow(self, flow):
        """
        API used by the AbiPy Worker to postprocess a Flow from the model.
        Wraps _postprocess_flow implented in the subclass.
        """
        self._postprocess_flow(flow)
        #self.flow_completed_at



class GsData(MongoModel):
    """GS SCF results: energy, forces, stresses fermi level, gaps"""

    ebands: ElectronBands = Field(
        ..., description="GS SCF results: energy, forces, stresses."
    )

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str):
        """File the model from the GSR filepath."""
        from abipy.electrons.gsr import GsrFile
        with GsrFile(gsr_filepath) as gsr:
            return cls.from_gsr(gsr)

    @classmethod
    def from_gsr(cls, gsr):
        """Fill the model from a |GsrFile|"""
        #print(gsr.ebands.structure)
        # TODO: Use "cartesian_forces_eV/Ang"
        #gsr.ebands.structure.remove_site_property("cartesian_forces")
        return cls(ebands=gsr.ebands)


class BandStructureFlowModel(FlowModel):
    """
    This model defines the input arguments used to build a Flow for band structure calculations
    as well as the submodels used to store the final results.

    Users are supposed to use this model to initialize a MongoDB collection with all
    the input arguments that will be used to generate the flow and provide a concrete
    implementation of:

        - _build_flow.
        - _postprocess_flow

    The first method receinves the input arguments from the MongoDB database
    and use these values to build a flow.

    The second method is invoked by the AbiPy worker when the calculation is completed.
    The function uses the Flow API to fill the ouput part of the model that
    will be then stored in the database collection.
    """

    ########
    # Input
    ########

    structure: Structure = Field(..., description="The input structure.")

    ########
    # Output
    ########

    scf_data: GsData = Field(None, description="Results produced by the GS SCF run.")

    nscf_data: GsData = Field(None, description="Results produced by the GS NSCF run.")

    def _build_flow(self, workdir, manager):
        """
        Build an AbiPy Flow using the input data available in the model and return it.

        Args:
            workdir: Working directory provided by the caller.
            manager: |TaskManager| object

        Return: |Flow| object.
        """
        from abipy.abio.factories import ebands_input
        pseudos = self.pseudo_specs.get_pseudos()  # FIXME
        multi = ebands_input(self.structure, pseudos)
                             #kppa=None, nscf_nband=None, ndivsm=15,
                             #ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
                             #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, dos_kppa=None)

        scf_input, nscf_input = multi.split_datasets()
        from .flows import bandstructure_flow
        return bandstructure_flow(workdir, scf_input, nscf_input, manager=manager)

    def _postprocess_flow(self, flow):
        """
        Analyze the flow and fills the model with output results.
        This function is called by the worker if the flow completed succesfully.
        """
        with flow[0][0].open_gsr() as gsr:
            self.scf_data = GsData.from_gsr(gsr)

        with flow[0][1].open_gsr() as gsr:
            self.nscf_data = GsData.from_gsr(gsr)



if __name__ == "__main__":

    import sys
    import abipy.data as abidata
    si = Structure.from_file(abidata.ref_file("refs/si_ebands/run.abi"))

    from abipy.flowtk.models import MongoConnector
    mongo_connector = MongoConnector.for_localhost(collection_name="ebands")
    collection = mongo_connector.get_collection()
    pseudos_specs = PseudoDojoSpecs.from_table_name("NCPW")

    #items = BandStructureFlowModel.find_runnable_oid_models(collection, limit=1)
    #if items:
    #    for item in items:
    #        print(item)

    #for model_cls in [MongoConnector, PseudoDojoSpecs,]: # BandStructureFlowModel.
    #    print(model_cls.schema_json(indent=2))
    #sys.exit(1)

    structures = [si]
    for structure in structures:
        model = BandStructureFlowModel(structure=structure, pseudos_specs=pseudos_specs)
        oid = model.mongo_insert(collection)

        assert model.flow_status == 0

        model.scf_data = GsData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"))
        model.nscf_data = GsData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"))
        #d = model.dict()
        #for k, v in d.items(): print(f"{type(k)} --> {type(v)}")

        oid = model.mongo_insert(collection)
        model.mongo_full_update_oid(oid, collection)

        same_model = BandStructureFlowModel.from_mongo_oid(str(oid), collection)
        print("same_model", type(same_model)) #, model)
        print("same_model.structure", type(same_model.structure)) #, model)
        assert same_model.flow_status == 0
        assert same_model.structure == structure
        assert isinstance(same_model.scf_data.ebands, ElectronBands)
        assert isinstance(same_model.nscf_data.ebands, ElectronBands)
        print(same_model.scf_data.ebands.structure)

    #from abipy.flowtk.worker import AbipyWorker
    #worker = AbipyWorker.new_with_name("ebands_worker", scratch_dir="/tmp/",
    #                                   manager=manager.
    #                                   mongo_connector=mongo_connector,
    #                                   flow_model=BandStructureFlowModel)

    #worker.serve()
    print("OK")

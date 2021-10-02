"""
Base pydantic models used by Abipy to implement the HTC version of the Flow connected to a MongoDB database.
"""
from __future__ import annotations

import json
import os
import inspect
import time

#import numpy as np

from typing import List, Any, Type, TypeVar, Optional, Iterable
#from uuid import UUID
#from pprint import pprint
from pydantic import BaseModel, Field, PrivateAttr, root_validator
from pydantic.main import ModelMetaclass
from bson.objectid import ObjectId
from pymongo import MongoClient
from pymongo.database import Database
from pymongo.collection import Collection
from gridfs import GridFS
from monty.json import MontyEncoder, MontyDecoder, MSONable
from pymatgen.util.serialization import pmg_serialize
from pymatgen.core import __version__ as pmg_version
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure as pmg_Structure
from pymatgen.core.periodic_table import Element
from abipy.core.release import __version__ as abipy_version
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput
from abipy.electrons.ebands import ElectronBands
from abipy.flowtk.events import EventReport


class AbipyDecoder(MontyDecoder):
    """
    Extends MontyDecoder adding support for decoding classes besides instances..
    """

    def process_decoded(self, d):

        if isinstance(d, dict) and "@module" in d and "@qualname" in d:
            modname = d["@module"]
            qualname = d["@qualname"]
            mod = __import__(modname, None, None, [qualname], 0)
            return getattr(mod, qualname)

        return super().process_decoded(d)


def cls2dict(cls) -> dict:
    """Encodes python class using __qualname__"""
    return {"@qualname": cls.__qualname__, "@module": cls.__module__}


class AbipyEncoder(MontyEncoder):
    """
    Extends MontyEncoder adding support for encoding classes besides instances..
    """

    def default(self, o) -> dict:  # pylint: disable=E0202
        if inspect.isclass(o):
            return cls2dict(o)

        return super().default(o)


def monty_json_dumps(obj: Any, **kwargs) -> str:
    #print("in monty json dumps with type:", type(obj))
    return json.dumps(obj, cls=AbipyEncoder, **kwargs)


def monty_json_loads(string: str, **kwargs) -> Any:
    #print("in monty json loads with string", string)
    return json.loads(string, cls=AbipyDecoder, **kwargs)


def monty_load(obj: Any, **kwargs) -> dict:
    s = monty_json_dumps(obj)
    #print("in monty load:\n", s)
    d = monty_json_loads(s)
    #print(d)
    return d


def monty_trick(obj: Any) -> Any:
    #print(f"in monty trick with {type(obj)}")
    return json.loads(monty_json_dumps(obj))


M = TypeVar('M', bound='AbipyModel')


class AbipyModel(BaseModel, MSONable):
    """
    Base class providing tools to serialize Abipy/Pymatgen objects
    implementing the MSONable protocol.
    Most of the AbiPy models should inherith from this class so that we
    can take advantage of pydantic Config to encode/decode objects that are not pydantci models.
    """

    class Config:

        # This is needed to be able to use MSONable AbiPy objects
        # such as ElectronBands that do not inherit from BaseModel
        arbitrary_types_allowed = True

        # Here we register specialized encoders to convert MSONable objects to JSON strings.
        # Subclasses will inherit these json_encoders.
        json_encoders = {
            ModelMetaclass: lambda cls: cls2dict(cls),
            #ObjectId: lambda oid: str(oid),
            ObjectId: lambda oid: monty_trick(oid),
            #np.ndarray: lambda arr: arr.tolist(),
            #UUID: lambda uuid: str(uuid),
            Element: lambda o: monty_trick(o),
            Composition: lambda o: monty_trick(o),
            pmg_Structure: lambda o: monty_trick(o),
            Structure: lambda o: monty_trick(o),
            AbinitInput: lambda o: monty_trick(o),
            ElectronBands: lambda o: monty_trick(o),
            EventReport: lambda o: monty_trick(o),
        }

        # whether to populate models with the value property of enums, rather than the raw enum.
        #use_enum_values = True

        #json_loads = monty_json_loads
        #json_dumps = monty_json_dumps


    @classmethod
    def from_dict(cls: Type[M], d: dict) -> M:
        """
        Reconstruct the model from a dictionary taking into account the MSONAble protocol
        """
        #print("d1 in from dict:")
        #pprint(d)
        d = monty_load(d)
        if isinstance(d, cls):
            return d
        #print("d2 after monty_load:")
        #pprint(d)
        return cls(**d)

    @pmg_serialize
    def as_dict(self) -> dict:
        """
        Convert the model to dictionary taking into account the MSONable protocol.
        """
        #return monty_trick(self)
        d = self.dict()
        #print("d in as dict:")
        #pprint(d)
        return d

    #@staticmethod
    #def serialize(self, data):

    def to_json(self, **kwargs) -> str:
        """
        Returns a json string representation of the MSONable object.
        """
        # NB: we pass kwargs unlike in monty to_json
        return json.dumps(self.as_dict(), cls=AbipyEncoder, **kwargs)

    @classmethod
    def from_json_file(cls: Type[M], filepath: str) -> M:
        """
        Helper function to reconstruct the model from a json file.
        """
        with open(filepath, "rt") as fp:
            return cls.from_json(fp.read())

    @classmethod
    def from_json(cls: Type[M], json_string: str) -> M:
        """
        Helper function to reconstruct the model from a json string.
        """
        new = monty_json_loads(json_string)
        if isinstance(new, cls): return new
        return cls(**new)

    def json_write(self, filepath: str, **kwargs) -> None:
        """
        Helper function to write a json string with the model to file.
        """
        with open(filepath, "wt") as fp:
            #fp.write(self.json(**kwargs))
            fp.write(self.to_json(**kwargs))

    #@abstractmethod
    #def get_panel_view(self):
    #    """Return panel object with a view of the model"""

    def yaml_dump(self):
        import ruamel.yaml as yaml
        json_string = self.to_json()
        return yaml.safe_dump(yaml.safe_load(json_string), default_flow_style=False)


class GridFsDesc(AbipyModel):
    """
    This submodel lives in a standard MongoDB collection and stores all the information
    needed to retrieve a file from a GridFs collection.
    Users do not need to instanciate this object directly.
    """

    filepath: str = Field(..., description="Absolute path to the file (depends on the host where the file has been produced")

    oid: ObjectId = Field(..., description="ID of the file in the GridFS collection.")

    collection_name: str = Field(..., description="Name of the GridFS collection in which the file is stored.")


#class GridFsDescList(AbipyModel):
#    """
#    This model should be used for defining a list of GridFs descriptors.
#    This is required as in ... we use introspection to find all the GridFSDesc submodels
#    included in the FlowModel and we only look for GridFSDesc or GridFsDescList.
#    """
#
#    __root__: List[GridFSDesc] = Field(..., description="List of GridFsDesc entries")
#
#    def __iter__(self):
#        return iter(self.__root__)
#
#    def __getitem__(self, item):
#        return self.__root__[item]



class MongoConnector(AbipyModel):
    """
    Stores the parameters used to connect to the MongoDB server and the name of
    the collection used to perform CRUD operations.
    This object is usually instanciated from the AbiPy configuration file that defines the host,
    port, dbname, username and password.
    Note, however, that we don't store the name of the collection in the configuration file.
    so the user has to specify it explictly when calling ``from_abipy_config``.
    """

    host: str = Field(..., description="Host address e.g. 'localhost'")

    port: int = Field(27017, description="MongoDB server port")

    dbname: str = Field("abipy", description="MongoDB database")

    collection_name: str = Field(..., description="MongoDB collection")

    username: str = Field(None, description="User name for MongDB authentication. Implies password.")

    password: str = Field(None, description="Password for MongDB authentication")

    # Private attributes
    _client: MongoClient = PrivateAttr(None)

    @classmethod
    def for_localhost(cls, collection_name: str, port: int = 27017, **kwargs) -> MongoConnector:
        """
        Build an instance assuming a MongoDB server running on localhost listening on the default port
        """
        d = dict(host="localhost", port=port, collection_name=collection_name)
        d.update(kwargs)
        return cls(**d)

    @classmethod
    def from_abipy_config(cls, collection_name: str, **kwargs) -> MongoConnector:
        """
        Build an instance from the configuration options stored in the AbiPy config.yml file.
        """
        from abipy.core.config import get_config
        config = get_config()
        d = dict(host=config.mongo_host, port=config.mongo_port,
                 username=config.mongo_username, password=config.mongo_password,
                 dbname=config.mongo_dbname,
                 collection_name=collection_name,
                 )
        d.update(kwargs)

        return cls(**d)

    @root_validator
    def check_username_and_password(cls, values):
        """Both Username and password are needed"""
        username, password = values.get("username"), values.get("password")
        if username is not None and password is None:
            raise ValueError('password must be provided when username is not None')
        return values

    def __str__(self) -> str:
        return self._repr_markdown_()

    def _repr_markdown_(self) -> str:
        """Markdown representation."""
        return f"""

## MongoConnector

- Host: {self.host}
- Port: {self.port}
- Database: {self.dbname}
- Collection: {self.collection_name}
- User: {self.username}
"""

    def get_client(self) -> MongoClient:
        """
        Establish a connection with the MongoDB server. Return MongoClient
        """
        from pymongo.errors import ConnectionFailure

        def client_is_connected(client: MongoClient) -> bool:
            try:
                # The ping command is cheap and does not require auth.
                client.admin.command("ping")
                return True
            except ConnectionFailure:
                print("Server not available. Trying to reconnect...")
                return False

        if self._client is not None and client_is_connected(self._client):
            return self._client

        # Reconnect and cache new client.
        if self.username and self.password:
            # Authenticate if needed
            self._client = MongoClient(host=self.host,
                                       port=self.port,
                                       username=self.username,
                                       password=self.password,
                                       authSource=self.dbname,
                                       #authMechanism='SCRAM-SHA-256'
                                       )

        else:
            self._client = MongoClient(host=self.host, port=self.port)

        return self._client

    def get_db(self) -> Database:
        """Return MongoDB database."""
        client = self.get_client()
        return client[self.dbname]

    def get_collection(self, collection_name: Optional[str] = None) -> Collection:
        """
        Returns MongoDB collection from its name. Use default collection if ``collection_name`` is None.
        """
        db = self.get_db()
        return db[collection_name or self.collection_name]

    def get_gridfs_and_name(self, collection_name: Optional[str] = None, ret_collname=False, **kwargs) -> Tuple[GridFS, str]:
        """
        Returns GridFS collection. Use automatically generated name.
        """
        db = self.get_db()
        coll_name = f"gridfs_{self.collection_name}" if collection_name is None else collection_name
        fs = GridFS(db, collection=coll_name, **kwargs)
        return fs, coll_name

    #def insert_models(models: List[MongoModel], collection_name: Optional[str] = None)) -> List[ObjectId]:
    #    """
    #    Insert list of models in a collection. Return list of objectid
    #    Use default collection if collection_name is None.
    #    """
    #    collection = self.get_collection(collection_name=collection_name)
    #    return mongo_insert_models(models, collection)

    def list_collection_names(self) -> List[str]:
        """
        "Return list of strings with all collection names in the database.
        """
        db = self.get_db()
        _filter = {"name": {"$regex": r"^(?!system\.)"}}
        return db.list_collection_names(filter=_filter)

    def list_collections(self) -> Iterable[Collection]:
        """
        Return iterable with collections available in the MongoDB database.
        """
        db = self.get_db()
        return db.list_collections()

    def gridfs_put_filepath(self, filepath: str) -> GridFsDesc:
        """
        Insert a file in the default GridFs collection.
        Build and return a GridFsDesc instance with metadata and the OID.
        """

        fs, fs_collname = self.get_gridfs_and_name()

        metadata = dict(
            filename=os.path.basename(filepath),
            filepath=filepath,
            parent_collection = self.collection_name,
        )

        with open(filepath, "rb") as fh:
           oid = fs.put(fh, **metadata)

        return GridFsDesc(filepath=filepath, oid=oid, collection_name=fs_collname)

    def abiopen_gfsd(self, gfsd: GridFSDesc):
        """
        Get binary data from the GridFS collection, write buffer to temporary file
        and use ``abiopen`` to open the file.
        """
        if gfsd.oid is None:
            raise RuntimeError("gridfs_ois is None, this means that the file has not been inserted in GridFS")

        fs, fs_collname = self.get_gridfs_and_name()
        grid_out = fs.get(gfsd.oid)

        from tempfile import mkstemp #, TemporaryFile, NamedTemporaryFile
        _, filepath = mkstemp(suffix=grid_out.filename)
        with open(filepath, "wb") as fh:
            fh.write(grid_out.read())

        from abipy.abilab import abiopen
        return abiopen(filepath)

    def open_mongoflow_gui(self, **serve_kwargs):
        """
        Start panel GUI that allows the user to inspect the results stored in the collection.
        """
        collection = self.get_collection()
        from .flow_models import FlowModel
        flow_model_cls = FlowModel.get_subclass_from_collection(collection)
        from abipy.panels.core import abipanel
        from abipy.panels.mongo_gui import MongoGui
        abipanel()
        app = MongoGui(self, flow_model_cls).get_app()
        import panel as pn
        pn.serve(app, **serve_kwargs)

    #def get_panel_view(self):
    #    import panel as pn
    #    return pn.pane.HTML(self._repr_markdown_())


class MockedMongoConnector(MongoConnector):
    """
    Mock a MongoConnector using mongomock. Mainly used for unit tests.
    """

    def get_client(self) -> MongoClient:
        # Cache the client also in the Mock so that we use the same mocked instances.
        if self._client is not None:
            return self._client

        import mongomock
        self._client = mongomock.MongoClient()
        return self._client

    #def get_collection(self, collection_name: Optional[str] = None) -> Collection:
    #    return mongomock.MongoClient().db.collection

    def get_gridfs_and_name(self, collection_name: Optional[str] = None, **kwargs) -> Tuple[GridFS, str]:
        from mongomock.gridfs import enable_gridfs_integration
        enable_gridfs_integration()
        return super().get_gridfs_and_name(collection_name=collection_name, **kwargs)


#class TopLevelModel(AbipyModel):
class MongoModel(AbipyModel):
    """
    Provides tiny wrappers around the MongoDB API.
    """

    abipy_version: str = Field(
        abipy_version, description="The version of abipy this document was built with"
    )

    pymatgen_version: str = Field(
        pmg_version, description="The version of pymatgen this document was built with"
    )

    @classmethod
    def from_mongo_oid(cls, oid: ObjectId, collection: Collection):
        """
        Return a model instance for the ObjectId oid and the MongoDB collection.
        """
        data = collection.find_one({'_id': oid})
        if not data:
            raise ValueError(f"Cannot find {self.__class__.__name__} model in collection: {collection}")

        data.pop("_id")
        data = AbipyDecoder().process_decoded(data)

        #if "scf_data" in data: print(data["scf_data"])

        return cls(**data)

    @classmethod
    def mongo_find(cls, query: dict, collection: Collection, **kwargs) -> QueryResults:
        """
        Find the models in the collection matching the mongodb query.
        """
        cursor = collection.find(filter=query, **kwargs)

        oids, models = [], []
        decoder = AbipyDecoder()
        for data in cursor:
            oids.append(data.pop("_id"))
            data = decoder.process_decoded(data)
            models.append(cls(**data))

        return QueryResults(oids, models, query)

    def mongo_insert(self, collection: Collection) -> ObjectId:
        """
        Insert the model in the collection. Return ObjectId.
        """
        # self.json build a JSON string: pydantic handles the serialization of pydantic entries
        # while arbitrary types (e.g. Strucuture) are encoded using json_encoders.
        # json loads allows us to get the final dictionary to be inserted in the MongoDB collection.
        doc = json.loads(self.json())
        return collection.insert_one(doc).inserted_id

    def mongo_full_update_oid(self, oid: ObjectId, collection: Collection) -> None:
        """
        Perform a full update of the document given the ObjectId in the collection.
        """
        # TODO: Implement Atomic transaction
        old_doc = collection.find_one({'_id': oid})
        if not old_doc:
            raise ValueError(f"Cannot find document with ObjectId: {oid}")

        new_doc = json.loads(self.json())
        old_doc.update(new_doc)

        collection.replace_one({"_id": oid}, new_doc, upsert=False)

    #def backup_oid_collection_name(self, oid: ObjectId, collection_name: str) -> str:
    #    """
    #    Write a backup file in the abipy HOME directory. Return path to the file.

    #    Useful if the MongoDB server goes down and the AbipyWorker needs
    #    to save the results somewhere on the filesystem.
    #    """
    #    bkp_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "bkp_models")
    #    if not os.path.isdir(bkp_dir): os.mkdir(bkp_dir)
    #    filename = f"{collection_name}_{str(oid)}"
    #    filepath = os.path.join(bkp_dir, filename)
    #    with open(filepath, "wt") as fp:
    #        fp.write(self.json())
    #        return filepath

#def read_backup_oid_collection_name(cls, oid, collection_name):
#    bkp_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "bkk_models")
#    filepath = os.path.join(bkp_dir, f"{collection_name}_{str(oid)}")
#    if not os.path.exists(filepath):
#    return cls.from_json_file(filepath)


def mongo_insert_models(models: List[MongoModel], collection: Collection, verbose: int = 0) -> List[ObjectId]:
    """
    Insert list of models in a collection. If verbose > 0, print elasped time.
    Return: list of ObjectId
    """
    start = time.time()
    docs = [json.loads(model.json()) for model in models]
    oids = collection.insert_many(docs).inserted_ids
    if verbose:
        print(f"Inserted {len(oids)} MongoDB documents in collection: {collection.full_name} in "
              f"{time.time() - start: .2f} [s]")
    return oids


class QueryResults:
    """
    Stores the results of a MongoDB query.
    """

    def __init__(self, oids: List[ObjectId], models: List[AbipyModel], query: dict) -> None:
        """
        Args:
            oids:
            models:
            query:
        """
        self.oids = oids
        self.models = models
        self.query = query

    @classmethod
    def empty_from_query(cls, query: dict) -> QueryResults:
        """
        Return an empty QueryResults produced by query.
        """
        return cls(oids=[], models=[], query=query)

    def __len__(self) -> int:
        return len(self.oids)

    def __bool__(self) -> bool:
        return bool(self.oids)

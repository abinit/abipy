"""
Base pydantic models used by Abipy.
"""
from __future__ import annotations

import json
import os
import inspect

from typing import List, Any, Type, TypeVar
from pprint import pprint
from pydantic import BaseModel, Field, PrivateAttr
from pydantic.main import ModelMetaclass
from bson.objectid import ObjectId
from pymongo import MongoClient
from pymongo.collection import Collection
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


class AbipyDecoder(MontyDecoder):
    """
    Extends MontyDecoder adding support for class decoding.
    """

    def process_decoded(self, d):

        if isinstance(d, dict) and "@module" in d and "@qualname" in d:
            modname = d["@module"]
            qualname = d["@qualname"]
            mod = __import__(modname, None, None, [qualname], 0)
            return getattr(mod, qualname)

        return super().process_decoded(d)


class AbipyEncoder(MontyEncoder):
    """
    Extends MontyEncoder to add support for class encoding.
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


def cls2dict(cls) -> dict:
    return {"@qualname": cls.__qualname__, "@module": cls.__module__}


M = TypeVar('M', bound='AbipyModel')


class AbipyModel(BaseModel, MSONable):
    """
    Base class providing tools to serialize Abipy/Pymatgen objects supporting the MSONable protocol.
    """

    class Config:
        # Here we register specialized encoders to convert MSONable objects to JSON strings.
        # Subclasses will inherit these json_encoders.
        json_encoders = {
            ModelMetaclass: lambda cls: cls2dict(cls),
            ObjectId: lambda oid: str(oid),
            Element: lambda o: monty_trick(o),
            Composition: lambda o: monty_trick(o),
            pmg_Structure: lambda o: monty_trick(o),
            Structure: lambda o: monty_trick(o),
            AbinitInput: lambda o: monty_trick(o),
            ElectronBands: lambda o: monty_trick(o),
        }

        # This is needed to be able to use MSONable AbiPy objects
        # such as ElectronBands that do not inherit from BaseModel
        arbitrary_types_allowed = True

        #json_loads = monty_json_loads
        #json_dumps = monty_json_dumps

    @classmethod
    def from_dict(cls: Type[M], d: dict) -> M:
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
        #return monty_trick(self)
        d = self.dict()
        #print("d in as dict:")
        #pprint(d)
        return d

    def to_json(self, **kwargs) -> str:
        """
        Returns a json string representation of the MSONable object.
        """
        # NB: we pass kwargs unlike in monty to_json
        return json.dumps(self.as_dict(), cls=AbipyEncoder, **kwargs)

    @classmethod
    def from_json_file(cls: Type[M], filepath: str) -> M:
        with open(filepath, "rt") as fp:
            return cls.from_json(fp.read())

    @classmethod
    def from_json(cls: Type[M], json_string: str) -> M:
        new = monty_json_loads(json_string)
        if isinstance(new, cls): return new
        return cls(**new)

    def json_write(self, filepath: str, **kwargs) -> None:
        with open(filepath, "wt") as fp:
            #fp.write(self.json(**kwargs))
            fp.write(self.to_json(**kwargs))

    #@abstractmethod
    #def get_panel_view(self):
    #    """Return panel object with a view of the model"""


class MongoConnector(AbipyModel):
    """
    Stores the parameters used to connect to the MongoDB server and the name of
    the (default) collection used to perform CRUD operations.
    """

    host: str = Field(..., description="Host address e.g. 'localhost'")

    port: int = Field(..., description="MongoDB server port. e.g. 27017")

    db_name: str = Field("abipy", description="Name of the MongoDB database")

    collection_name: str = Field(..., description="Name of the MongoDB collection")

    user: str = Field(None, description="User name for authentication")

    password: str = Field(None, description="Password for authentication")

    # Private attributes
    _client: Any = PrivateAttr(None)

    @classmethod
    def for_localhost(cls, collection_name: str, port: int = 27017, **kwargs) -> MongoConnector:
        """
        Build connector assuming a MongoDB server running on localhost listening on the default port
        """
        d = dict(host="localhost", port=port, collection_name=collection_name)
        d.update(kwargs)
        return cls(**d)

    def _repr_markdown_(self) -> str:
        """Markdown representation."""
        return f"""

## MongoConnector

- Host: {self.host}
- Port: {self.port}
- Database: {self.db_name}
- Collection: {self.collection_name}
- User: {self.user}

"""

    #def get_panel_view(self):
    #    import panel as pn
    #    return pn.pane.HTML(self._repr_markdown_())

    def get_client(self) -> MongoClient:
        """
        Establish a connection with the MongoDB server.

        Return: MongoClient
        """
        from pymongo.errors import ConnectionFailure

        def client_is_connected(client: MongoConnector) -> bool:
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
        self._client = MongoClient(host=self.host, port=self.port)
        return self._client

    def get_collection(self, collection_name: str = None) -> Collection:
        """
        Returns MongoDB collection from its name.
        """
        client = self.get_client()
        db = client[self.db_name]

        # Authenticate if needed
        if self.user and self.password:
            db.autenticate(self.user, password=self.password)

        return db[collection_name or self.collection_name]

    def open_mongoflow_gui(self, **serve_kwargs):
        collection = self.get_collection()
        from .flow_models import FlowModel
        flow_model = FlowModel.get_subclass_from_collection(collection)
        from abipy.panels.core import abipanel
        from abipy.panels.mongo_gui import MongoGui
        abipanel()
        app = MongoGui(self, flow_model).get_app()
        import panel as pn
        pn.serve(app, **serve_kwargs)


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

    #build_date: datetime = Field(
    #    default_factory=datetime.utcnow,
    #    description="The build date for this document",
    #)

    @classmethod
    def from_mongo_oid(cls, oid: ObjectId, collection: Collection):
        """
        Return a model instance for the ObjectId oid and the collection.
        """
        data = collection.find_one({'_id': oid})
        if data is None: return data
        data.pop("_id")
        data = AbipyDecoder().process_decoded(data)

        return cls(**data)

    @classmethod
    def mongo_find(cls, query: dict, collection: Collection, **kwargs) -> QueryResults:
        """
        Find the models in the collection matching the mongodb query.
        """
        cursor = collection.find(query, **kwargs)
        if cursor is None:
            return QueryResults.empty_from_query(query)

        oids, models = [], []
        decoder = AbipyDecoder()
        for data in cursor:
            oids.append(data.pop("_id"))
            data = decoder.process_decoded(data)
            models.append(cls(**data))

        return QueryResults(oids, models, query)

    def mongo_insert(self, collection: Collection) -> ObjectId:
        """
        Insert the model in a collection. Return ObjectId.
        """
        # Use plain json loads here as we want to insert a dictionary
        doc = json.loads(self.json())
        return collection.insert_one(doc).inserted_id

    def mongo_full_update_oid(self, oid: ObjectId, collection: Collection) -> None:
        """
        Perform a full update of the document given the ObjectId in the collection.
        """
        # TODO: Implement Atomic transaction
        old_doc = collection.find_one({'_id': oid})
        if old_doc is None:
            raise RuntimeError(f"Cannot find document with ObjectId: {oid}")

        new_doc = json.loads(self.json())
        old_doc.update(new_doc)

        collection.replace_one({"_id": oid}, new_doc, upsert=False)
        #collection.update_one({'_id': oid}, new_doc)

        #update =
        #{ $set:
        #   {
        #     quantity: 500,
        #     details: { model: "14Q3", make: "xyz" },
        #     tags: [ "coats", "outerwear", "clothing" ]
        #   }
        #}
        #self.collection.update_one(
        #    filter={"_id": oid},
        #    update={"$set": person.dict()},
        #    upsert=False
        # )

    def backup_oid_collection_name(self, oid: ObjectId, collection_name: str) -> str:
        """
        Write a backup file in the abipy HOME directory. Return path to the file.

        Useful if the MongoDB server goes down and the AbipyWorker needs
        to save the results somewhere on the filesystem.
        """
        bkp_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "bkk_models")
        if not os.path.isdir(bkp_dir): os.mkdir(bkp_dir)
        filename = f"{collection_name}_{str(oid)}"
        filepath = os.path.join(bkp_dir, filename)
        with open(filepath, "wt") as fp:
            fp.write(self.json())
            return filepath

#def read_backup_oid_collection_name(cls, oid, collection_name):
#    bkp_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "bkk_models")
#    filepath = os.path.join(bkp_dir, f"{collection_name}_{str(oid)}")
#    if not os.path.exists(filepath):
#    return cls.from_json_file(filepath)


def mongo_insert_models(models: List[MongoModel], collection: Collection) -> List[ObjectId]:
    """
    Insert list of models in a collection. Return list of objectid
    """
    docs = [json.loads(model.json()) for model in models]
    return collection.insert_many(docs).inserted_ids


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

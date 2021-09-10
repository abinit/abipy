"""
Pydantic Models.
"""
from __future__ import annotations

#import os
import json
import inspect
import traceback
import panel as pn

#from pprint import pprint
from datetime import datetime
from abc import ABC, abstractmethod
from enum import Enum #, IntEnum
from typing import List, Tuple, Dict, Optional, Any, Type, ClassVar
from pydantic import BaseModel, Field, BaseConfig, PrivateAttr #, ModelMetaClass
from pydantic.main import ModelMetaclass
from bson.objectid import ObjectId
from pymongo import MongoClient
from pymongo.collection import Collection
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.util.serialization import pmg_serialize
from pymatgen.util.typing import VectorLike, MatrixLike
from pymatgen.core import __version__ as pmg_version
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure as pmg_Structure
from abipy.core.release import __version__ as abipy_version
from abipy.core.structure import Structure
#from abipy.tools.tensors import Stress
from abipy.abio.inputs import AbinitInput
from abipy.electrons.ebands import ElectronBands
from abipy.electrons.gsr import GsrFile
from abipy.panels.core import ply
from abipy.panels.viewers import JSONViewer


class AbipyDecoder(MontyDecoder):

    def process_decoded(self, d):

        if isinstance(d, dict) and "@module" in d and "@qualname" in d:
            modname = d["@module"]
            qualname = d["@qualname"]
            mod = __import__(modname, None, None, [qualname], 0)
            #print("in my decoded with mod", mod)
            return getattr(mod, qualname)

        return super().process_decoded(d)


class AbipyEncoder(MontyEncoder):

    def default(self, o) -> dict:  # pylint: disable=E0202
        if inspect.isclass(o):
            #print("isclass", o)
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


class AbipyBaseModel(BaseModel):
    """
    Base class providing tools to serialize Abipy/Pymatgen objects
    supporting the MSONable protocol.
    """
    class Config:
        # Here we register specialized encoders to convert MSONable objects to JSON objects.
        # Subclasses will inherit these json_encoders.
        json_encoders = {
            ModelMetaclass: lambda cls: cls2dict(cls),
            Composition: lambda o: monty_trick(o),
            pmg_Structure: lambda o: monty_trick(o),
            Structure: lambda o: monty_trick(o),
            #Stress: lambda o: monty_trick(o),
            AbinitInput: lambda o: monty_trick(o),
            ElectronBands: lambda o: monty_trick(o),
        }

        # This is needed to be able to use MSONable AbiPy objects
        # such as ElectronBands that do not inherit from BaseModel
        arbitrary_types_allowed = True

        #json_loads = monty_json_loads
        #json_dumps = monty_json_dumps

    @classmethod
    def from_json_file(cls, filepath: str):
        with open(filepath, "rt") as fp:
            return cls.from_json(fp.read())

    @classmethod
    def from_json(cls, json_string: str):
        return cls(**monty_json_loads(json_string))

    def json_write(self, filepath: str, **kwargs) -> None:
        with open(filepath, "wt") as fp:
            fp.write(self.json(**kwargs))

    @classmethod
    def from_dict(cls, d: dict):
        return cls(**d)

    @pmg_serialize
    def as_dict(self) -> dict:
        return self.dict()


class MongoConnector(AbipyBaseModel):
    """
    Model with the parameters used to connect to the MongoDB server and the name of
    the (default) collection used to perform CRUD operations.
    """

    host: str = Field(..., description="Host address e.g. 'localhost'")

    port: int = Field(..., description="MongoDB server port. e.g. 27017")

    db_name: str = Field("abipy", description="Name of the MongoDB database")

    collection_name: str = Field(..., description="Name of the collection")

    user: str = Field(None, description="User name for authentication. Default: None")

    password: str = Field(None, description="Password for authentication. Default: None")

    # Private attributes
    _client: Any = PrivateAttr(None)

    @classmethod
    def for_localhost(cls, collection_name: str, port: int = 27017) -> MongoConnector:
        return cls(host="localhost", port=port, collection_name=collection_name)

    def _repr_markdown_(self) -> str:
        return f"""

## MongoConnector

- Host: {self.host}
- Port: {self.port}
- Database: {self.db_name}
- Collection: {self.collection_name}
- User: {self.user}

"""

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
        Returns MongoDB collection
        """
        client = self.get_client()
        db = client[self.db_name]

        # Authenticate if needed
        if self.user and self.password:
            db.autenticate(self.user, password=self.password)

        return db[collection_name or self.collection_name]

    def open_mongoflow_gui(self, **serve_kwargs):
        collection = self.get_collection()
        flow_model = FlowModel.get_subclass_from_collection(collection)
        from abipy.panels.core import abipanel
        from abipy.panels.mongo_gui import MongoGui
        abipanel()
        app = MongoGui(self, flow_model).get_app()
        pn.serve(app, **serve_kwargs)


class MongoModel(AbipyBaseModel):
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
        oid = ObjectId(oid)
        data = collection.find_one({'_id': oid})
        if data is None: return data
        data.pop("_id")
        data = AbipyDecoder().process_decoded(data)
        #data = monty_load(data)

        return cls(**data)
        #return cls(**dict(data, id=id))

    @classmethod
    def mongo_find(cls, query: dict, collection: Collection, **kwargs):

        cursor = collection.find(query, **kwargs)
        if cursor is None:
            return QueryResults.empty_from_query(query, collection)

        oids, models = [], []
        for data in cursor:
            oids.append(data.pop("_id"))
            #print(data)
            data = AbipyDecoder().process_decoded(data)
            #print(data)
            #data = monty_load(data)
            models.append(cls(**data))

        return QueryResults(oids, models, query, collection)

    def mongo_insert(self, collection: Collection) -> ObjectId:
        """
        Insert the model in collection. Return ObjectId.
        """
        # Use plain json loads here as we want to insert a dictionary
        doc = json.loads(self.json())
        return collection.insert_one(doc).inserted_id

    def mongo_full_update_oid(self, oid: ObjectId, collection: Collection) -> None:
        """
        Perform a full update of the model given the ObjectId in the collection.
        """
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

    #def write_backup_oid_collection_name(self, oid, collection_name):
    #    bkp_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "bkk_models")
    #    if not os.path.isdir(bkp_dir): os.mkdir(bkp_dir)
    #    filename = f"{collection_name}_{str(oid)}"
    #    with open(os.path.join(bkp_dir, filename), "wt") as fp:
    #        fp.write(self.json())

    #@classmethod
    #def read_backup_oid_collection_name(cls, oid, collection_name):
    #    bkp_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "bkk_models")
    #    filepath = os.path.join(bkp_dir, f"{collection_name}_{str(oid)}")
    #    return cls.from_json_file(filepath), filepath


def mongo_insert_models(models: List[MongoModel], collection: Collection) -> List[ObjectId]:
    docs = [json.loads(model.json()) for model in models]
    r = collection.insert_many(docs)
    return r.inserted_ids


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

# Some of these models are inspired to emmet
# or we use names and models that should be compatible.


class CrystalSystem(str, Enum):
    """
    The crystal system of the lattice
    """
    tri = "Triclinic"
    mono = "Monoclinic"
    ortho = "Orthorhombic"
    tet = "Tetragonal"
    trig = "Trigonal"
    hex_ = "Hexagonal"
    cubic = "Cubic"


class StructureData(AbipyBaseModel):
    """
    This model strores structural info set and symmetry metadata.
    """

    #mpid ? # Optional

    crystal_system: CrystalSystem = Field(
        None, title="Crystal System", description="The crystal system for this lattice"
    )

    spg_symbol: str = Field(
        None,
        title="Space Group Symbol",
        description="The spacegroup symbol for the lattice",
    )

    spg_number: int = Field(
        None,
        title="Space Group Number",
        description="The spacegroup number for the lattice",
    )

    point_group: str = Field(
        None, title="Point Group Symbol", description="The point group for the lattice"
    )

    symprec: float = Field(
        None,
        title="Symmetry Finding Precision",
        description="The precision given to spglib to determine the symmetry of this lattice",
    )

    spglib_version: str = Field(None, title="SPGLib version")

    # Structure metadata
    nsites: int = Field(None, description="Total number of sites in the structure")
    elements: List[Element] = Field(None, description="List of elements in the material")
    nelements: int = Field(None, title="Number of Elements")
    composition: Composition = Field(None, description="Full composition for the material")
    composition_reduced: Composition = Field(
        None,
        title="Reduced Composition",
        description="Simplified representation of the composition",
    )
    formula_pretty: str = Field(
        None,
        title="Pretty Formula",
        description="Cleaned representation of the formula",
    )
    formula_anonymous: str = Field(
        None,
        title="Anonymous Formula",
        description="Anonymized representation of the formula",
    )
    chemsys: str = Field(
        None,
        title="Chemical System",
        description="dash-delimited string of elements in the material",
    )
    volume: float = Field(
        None,
        title="Volume",
        description="Total volume for this structure in Angstroms^3",
    )

    density: float = Field(None, title="Density", description="Density in grams per cm^3")

    density_atomic: float = Field(
        None,
        title="Packing Density",
        description="The atomic packing density in atoms per cm^3",
    )

    structure: Structure = Field(..., description="Abipy Structure object.")

    @classmethod
    def from_structure(cls, structure: Structure) -> StructureData:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, spglib

        #symprec = SETTINGS.SYMPREC
        symprec = 0.1
        sg = SpacegroupAnalyzer(structure, symprec=symprec)
        symmetry: Dict[str, Any] = {"symprec": symprec}
        if not sg.get_symmetry_dataset():
            sg = SpacegroupAnalyzer(structure, 1e-3, 1)
            symmetry["symprec"] = 1e-3

        symmetry.update(
            {
                "spg_symbol": sg.get_space_group_symbol(),
                "spg_number": sg.get_space_group_number(),
                "point_group": sg.get_point_group_symbol(),
                "crystal_system": CrystalSystem(sg.get_crystal_system().title()),
                "hall": sg.get_hall(),
                "spglib_version": spglib.__version__,
            }
        )

        comp = structure.composition.remove_charges()
        elsyms = sorted(set([e.symbol for e in comp.elements]))
        #symmetry = SymmetryData.from_structure(structure)

        data = {
            "nsites": structure.num_sites,
            "elements": elsyms,
            "nelements": len(elsyms),
            "composition": comp,
            "composition_reduced": comp.reduced_composition,
            "formula_pretty": comp.reduced_formula,
            "formula_anonymous": comp.anonymized_formula,
            "chemsys": "-".join(elsyms),
            "volume": structure.volume,
            "density": structure.density,
            "density_atomic": structure.volume / structure.num_sites,
            #"symmetry": symmetry,
        }

        structure = Structure.as_structure(structure) #.copy()
        symmetry.update(data)

        return cls(structure=structure, **symmetry)

    def get_title(self) -> str:
        return f"Structure: {self.formula_pretty}, {self.spg_symbol} ({self.spg_number}), " + \
               f"{self.crystal_system}, natom: {self.nsites}"

    #def get_view(self):
    #    return


class ExecStatus(str, Enum):
    """
    Possible status of the execution flow.
    """
    init = "Initialized"
    errored = "Errored"
    completed = "Completed"


class FlowData(AbipyBaseModel):

    exec_status: ExecStatus = ExecStatus.init

    worker_name: str = Field(None, description="The name of th AbiPy worker running this Flow.")

    worker_hostname: str = Field(None, description="The name of the machine running this Flow.")

    workdir: str = Field(None, description="The directory where the Flow is executed.")

    created_at: datetime = Field(
        description="Timestamp for when this material document was first created",
        default_factory=datetime.utcnow,
    )

    completed_at: datetime = Field(
        None,
        description="Timestamp for the most recent calculation update for this property",
    )

    tracebacks: List[str] = Field([], description="List of exception tracebacks")


class QueryResults:
    """
    This object performs MongoDB queries in collection, convert dict to models
    and store the results
    """

    def __init__(self, oids, models, query, collection):
        self.oids = oids
        self.models = models
        self.query = query
        self.collection = collection

    @classmethod
    def empty_from_query(cls, query, collection) -> QueryResults:
        return cls(oids=(), models=(), query=query, collection=collection)

    def __len__(self) -> int:
        return len(self.oids)

    def __bool__(self) -> bool:
        return bool(self.oids)


class FlowModel(MongoModel, ABC):
    """
    Base class for models associated to a Flow calculation performed by an AbiPy Worker.
    This model implements the businness logic to get a Flow from the database,
    run it and post-process the results using the Flow API.
    The business logic is then used by the AbiPyWorker to automate multiple calculations.
    """

    flow_data: FlowData = Field(None, description="")

    input_structure_data: StructureData = Field(..., description="Input structure.")

    pseudos_specs: PseudoDojoSpecs = Field(..., description="PseudoPotential Table.")

    # Private class attributes
    _magic_key: ClassVar[str] = "__flow_model__"

    @classmethod
    def init_collection(cls, collection: Collection) -> ObjectId:
        new_dict = cls2dict(cls)

        cnt = 0
        for doc in collection.find({cls._magic_key: {"$exists": True}}):
            oid = doc.pop("_id")
            old_dict = doc.get(cls._magic_key)
            cnt += 1
            if new_dict != old_dict:
                raise ValueError()

        if cnt:
            return oid

        return collection.insert_one({cls._magic_key: new_dict}).inserted_id

    @classmethod
    def find_runnable_oid_models(cls, collection: Collection, **kwargs) -> List[Tuple[ObjectId, FlowModel]]:
        """
        Return list of (oid, model) tuple with the models that are ready to run.
        """
        # NB: None not only matches itself but also matches “does not exist.”
        # Thus, querying for a key with the value None will return all documents lacking that key
        # hence we have to check that the key is None and $exists:
        query = {"flow_data": {"$in": [None], "$exists": True}}
        cursor = collection.find(query, **kwargs)
        if cursor is None: return []

        items = []
        decoder = AbipyDecoder()
        for doc in cursor:
            oid = doc.pop("_id")
            doc = decoder.process_decoded(doc)
            #print("Found runnable doc:", doc)
            items.append((oid, cls(**doc)))

        return items

    @classmethod
    def get_subclass_from_collection(cls, collection: Collection) -> FlowModel:
        cnt = 0
        for doc in collection.find({cls._magic_key: {"$exists": True}}):
            _ = doc.pop("_id")
            data = doc.get(cls._magic_key)
            cnt += 1

        if cnt == 0:
            raise RuntimeError(f"Cannot find document with magic key: {cls._magic_key}")
        if cnt > 1:
            raise RuntimeError(f"Found {cnt} documents with magic key: {cls._magic_key}")

        sub_class = AbipyDecoder().process_decoded(data)
        return sub_class

    @abstractmethod
    def build_flow(self, workdir: str, manager):
        """
        Return Flow. Must be provided by the subclass
        """

    def build_flow_and_update_collection(self, workdir: str,
                                         oid: ObjectId,
                                         collection: Collection,
                                         abipy_worker) -> None:
        """
        API used by the AbiPy Worker to build a Flow from the model.
        Wraps build_flow implemented by the subclass.
        """
        #if self.flow_data.exec_status == ExecStatus.errored: return
        self.flow_data = flow_data = FlowData()

        # Set the workdir of the Flow in the model.
        flow_data.workdir = workdir
        flow_data.exec_status = ExecStatus.init

        flow_data.worker_name = abipy_worker.name
        from socket import gethostname
        flow_data.worker_hostname = gethostname()

        try:
            # Call the method provided by the user.
            return self.build_flow(workdir, abipy_worker.manager)

        except:
            flow_data.exec_status = ExecStatus.errored
            flow_data.tracebacks.append(traceback.format_exc())
            return None

        finally:
            self.mongo_full_update_oid(oid, collection)

    @abstractmethod
    def postprocess_flow(self, flow):
        """
        Postprocess the Flow. Must be provided by the subclass.
        """

    def postprocess_flow_and_update_collection(self, flow,
                                               oid: ObjectId,
                                               collection: Collection) -> None:
        """
        API used by the AbiPy Worker to postprocess a Flow from the model.
        Wraps postprocess_flow implented in the subclass.
        """
        flow_data = self.flow_data
        #if flow_data.exec_status == ExecStatus.errored: return
        try:
            self.postprocess_flow(flow)
            flow_data.exec_status = ExecStatus.completed
        except:
            flow_data.exec_status = ExecStatus.errored
            flow_data.tracebacks.append(traceback.format_exc())
        finally:
            flow_data.completed_at = datetime.utcnow()
            print("postprocessing_done, status:", flow_data.exec_status)
            self.mongo_full_update_oid(oid, collection)

    @classmethod
    def find(cls, query: dict, collection: Collection, **kwargs) -> QueryResults:
        cursor = collection.find(query, **kwargs)
        if cursor is None:
            return QueryResults.empty_from_query(query, collection)

        oids, models = [], []
        decoder = AbipyDecoder()
        for doc in cursor:
            oid = doc.pop("_id")
            doc = decoder.process_decoded(doc)
            models.append(cls(**doc))
            oids.append(oid)

        return QueryResults(oids, models, query, collection)

    @classmethod
    def find_by_spg_number(cls, spg_number: int, collection, **kwargs) -> QueryResults:
        query = {"input_structure_data.spg_number": int(spg_number)}
        return cls.find(query, collection, **kwargs)

    @classmethod
    def find_by_formula(cls, formula: str, collection, **kwargs) -> QueryResults:
        query = {"input_structure_data.formula_pretty": formula}
        return cls.find(query, collection, **kwargs)

    #@abstractmethod
    #def get_common_queries(self): -> list:
    #    """
    #    Return list of MongoDB queries Flow. Used by the GUI
    #    """


class GsData(AbipyBaseModel):
    """
    Ground-state results: energy, forces, stress tensor, Fermi level, band gaps.
    """
    max_force_ev_over_ang: float = Field(
            None,
            description="Max absolute Cartesian force in eV/Ang. None if forces are not available.")

    abs_pressure_gpa: float = Field(None, description="Pressure in GPa. Absolute Value.")
    pressure_gpa: float = Field(None, description="Pressure in GPa. NB: Value with sign!.")

    is_scf_run: bool = Field(None, description="True is this a Ground-state run.")

    #cart_stress_tensor_gpa: MatrixLike = Field(None, description="Cartesian stress Tensor in GPA")
    #cart_forces_ev_ang: MatrixLike = Field(None, description="Cartesian forces in eV ang^-1")

    ebands: ElectronBands = Field(..., description="Electronic bands.")

    energy: float = Field(
        None, description="The total DFT energy in eV"
    )
    energy_per_atom: float = Field(
        None, description="The DFT energy per atom in eV"
    )
    #bandgap: float = Field(None, description="The DFT bandgap for the last calculation")
    #forces: List[Vector3D] = Field(
    #    [], description="Forces on atoms from the last calculation"
    #)
    #stress: Matrix3D = Field(
    #    [], description="Stress on the unitcell from the last calculation"
    #)

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str) -> GsData:
        """
        Fill the model from the GSR filepath.
        """
        with GsrFile(gsr_filepath) as gsr:
            return cls.from_gsr(gsr)

    @classmethod
    def from_gsr(cls, gsr: GsrFile) -> GsData:
        """
        Fill the model from a |GsrFile|
        """
        # TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
        #gsr.ebands.structure.remove_site_property("cartesian_forces")
        kwargs = dict(
            ebands=gsr.ebands,
            is_scf_run=gsr.is_scf_run,
        )

        if gsr.is_scf_run:
            kwargs.update(dict(
                pressure_gpa=gsr.pressure,
                abs_pressure_gpa=abs(gsr.pressure),
                max_force_ev_over_ang=gsr.max_force,
                energy=float(gsr.energy),
                energy_per_atom=float(gsr.energy_per_atom),
                #cart_stress_tensor_gpa=
                #cart_forces_ev_over_ang=
            ))

        return cls(**kwargs)


class EbandsFlowModel(FlowModel):
    """
    This model defines the input arguments used to build a Flow for band structure calculations
    as well as the submodels used to store the final results.

    Users are supposed to use this model to initialize a MongoDB collection with all
    the input arguments that will be used to generate the flow and provide a concrete
    implementation of:

        - build_flow.
        - postprocess_flow

    The first method receives the input arguments from the MongoDB database
    and use these values to build a flow.

    The second method is invoked by the AbiPy worker when the calculation is completed.
    The function uses the Flow API to fill the ouput part of the model that
    will be then be stored in the database collection.

    NOTES: The flow should have a single structure.
    """

    ########
    # Input
    ########

    kppa: int = Field(1000, description="Defines the sampling used for the SCF run. Defaults to 1000 if not given.")

    ecut: float = Field(6, description="")

    ndivsm: int = Field(2, description="Number of divisions used to sample the smallest segment of the k-path.")

    spin_mode: str = Field("unpolarized", description="Spin polarization")

    charge: float = Field(0.0, description="Electronic charge added to the unit cell.")

    smearing: str = Field("fermi_dirac:0.1 eV", description="Smearing technique.")

    dos_kppa: int = Field(None,
        description="Scalar or List of integers with the number of k-points per atom " +
                    "to be used for the computation of the DOS (None if DOS is not wanted")

    paral_kgb: int = Field(0, description="")

    ########
    # Output
    ########

    scf_input: AbinitInput = Field(None, description="Abinit Input file generated by AbiPy.")

    scf_data: GsData = Field(None, description="Results produced by the GS SCF run.")

    nscf_kpath_data: GsData = Field(None, description="Results produced by the GS NSCF run.")

    nscf_kmesh_data: GsData = Field(None, description="Results produced by the GS NSCF run.")

    def build_flow(self, workdir: str, manager):
        """
        Build an AbiPy Flow using the input data available in the model and return it.

        Args:
            workdir: Working directory provided by the caller.
            manager: |TaskManager| object.

        Return: |Flow| object.
        """
        from abipy.abio.factories import ebands_input
        import abipy.data as abidata
        pseudos = abidata.pseudos("14si.pspnc")
        #pseudos = self.pseudo_specs.get_pseudos()  # FIXME
        #ecut = 6
        multi = ebands_input(self.input_structure_data.structure, pseudos,
                             kppa=self.kppa, nscf_nband=None, ndivsm=self.ndivsm,
                             ecut=self.ecut, pawecutdg=None, scf_nband=None, accuracy="normal",
                             spin_mode=self.spin_mode,
                             smearing=self.smearing, charge=self.charge,
                             scf_algorithm=None, dos_kppa=self.dos_kppa,
                             )

        multi.set_vars(paral_kgb=self.paral_kgb)

        self.scf_input, nscf_input = multi.split_datasets()
        from .flows import bandstructure_flow
        # TODO: Dos!
        return bandstructure_flow(workdir, self.scf_input, nscf_input, manager=manager)

    def postprocess_flow(self, flow):
        """
        Analyze the flow and fills the model with output results.
        This function is called by the worker if the flow completed succesfully.
        """
        with flow[0][0].open_gsr() as gsr:
            self.scf_data = GsData.from_gsr(gsr)

        with flow[0][1].open_gsr() as gsr:
            self.nscf_kpath_data = GsData.from_gsr(gsr)

        if self.dos_kppa is not None:
            with flow[0][2].open_gsr() as gsr:
                self.nscf_kmesh_data = GsData.from_gsr(gsr)

    def get_view(self):
        title = self.input_structure_data.get_title()
        structure = self.input_structure_data.structure
        a, b, c = structure.lattice.abc
        alpha, beta, gamma = structure.lattice.angles
        header = f"""
## {title}

- Lattice lengths: a = {a:.6f}, b = {b:.6f}, c = {c:.6f} Ang
- Lattice angles: α = {alpha:.3f}, β = {beta:.3f}, ɣ = {gamma:.3f} degrees
- Pressure: {self.scf_data.pressure_gpa:.3f} GPa
- Max |Force|: {self.scf_data.max_force_ev_over_ang:.8f} eV/Ang
- Energy: {self.scf_data.energy:.4f} eV
- Energy per atom: {self.scf_data.energy_per_atom:.4f} eV
"""
        return pn.Column(
            #self.input_structure_data.get_view(),
            header,
            pn.Row(
                ply(self.nscf_kpath_data.ebands.plotly(show=False)),
                pn.pane.HTML(self.scf_input._repr_html_()),
            ),
            "### MongoDB Document",
            JSONViewer(self.json(), depth=1),
            pn.layout.Divider(),
            sizing_mode="stretch_both",
        )

    #def get_dataframe(self, **kwargs):

    @classmethod
    def get_common_queries(cls) -> list:
        return [
            {"$and": [
                {"scf_data.abs_pressure_gpa:": {"$gt": 2}},
                {"scf_data.max_force_ev_over_ang": {"$gt": 1e-6}}]
            },
        ]


class PhononData(AbipyBaseModel):

    ddb_string: str = Field(..., description="DDB string")

    #min_phfreq_ev: float = Field(..., description="DDB string")
    #max_phfreq_ev: float = Field(..., description="DDB string")

    #phbands: PhononBands = Field(..., description="DDB string")
    #phdos: PhononDos = Field(..., description="DDB string")

    @classmethod
    def from_ddb(cls, ddb):
        kwargs = dict(
            ddb_string=ddb.get_string(),
        )

        #with ddb.anaget_phbst_and_phdos_files() as g:
        #    phbst_file, phdos_file = g[0], g[1]
        #    phbands = phbst_file.phbands
        #    phdos = phdos_file.phdos
        #    kwargs.update(dict(
        #        phbands=phbands,
        #        phdos=phdos,
        #
        #    ))

        return cls(**kwargs)


class PhononFlowModel(FlowModel):
    """
    This model defines the input arguments used to build a Flow for band structure calculations
    as well as the submodels used to store the final results.

    Users are supposed to use this model to initialize a MongoDB collection with all
    the input arguments that will be used to generate the flow and provide a concrete
    implementation of:

        - build_flow.
        - postprocess_flow

    The first method receives the input arguments from the MongoDB database
    and use these values to build a flow.

    The second method is invoked by the AbiPy worker when the calculation is completed.
    The function uses the Flow API to fill the ouput part of the model that
    will be then stored in the database collection.

    NOTES: The flow should have a single structure.
    """

    ########
    # Input
    ########

    scf_input: AbinitInput = Field(..., description="Input structure.")

    with_becs: bool = Field(..., description="Compute Born effective charges.")
    with_quad: bool = Field(..., description="Activate calculation of dynamical quadrupoles.")
    with_flexoe: bool = Field(False, description="Activate computation of flexoelectric tensor.")

    #kppa: int = Field(1000, description="")
    #ecut: float= Field(6, description="")
    #ndivsm: int = Field(2, description="")
    #spin_mode: str = Field("unpolarized", description="")
    #charge: float = Field(0.0, description="")
    #smearing: str = Field("fermi_dirac:0.1 eV", description="")
    #dos_kppa: int = Field(None, description="")
    #paral_kgb: int = Field(0, description="")

    ########
    # Output
    ########

    scf_data: GsData = Field(None, description="Results produced by the GS SCF run.")

    phonon_data: PhononData = Field(None, description="Results produced by the GS SCF run.")

    def build_flow(self, workdir, manager):
        """
        Build an AbiPy Flow using the input data available in the model and return it.

        Args:
            workdir: Working directory provided by the caller.
            manager: |TaskManager| object.

        Return: |Flow| object.
        """
        #import abipy.data as abidata
        #pseudos = abidata.pseudos("14si.pspnc")
        ##pseudos = self.pseudo_specs.get_pseudos()  # FIXME
        ##ecut = 6

        # Build input for GS calculation
        #multi.set_vars(paral_kgb=self.paral_kgb)

        #scf_input(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
        #          spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
        #          shift_mode="Monkhorst-Pack"):

        # Create flow to compute all the independent atomic perturbations
        # corresponding to a [4, 4, 4] q-mesh.
        # Electric field and Born effective charges are also computed.
        from abipy.flowtk import PhononFlow
        flow = PhononFlow.from_scf_input(workdir, self.scf_input,
                                         ph_ngqpt=(2, 2, 2),
                                         #ph_ngqpt=(4, 4, 4),
                                         with_becs=self.with_becs, with_quad=self.with_quad,
                                         with_flexoe=self.with_flexoe, manager=manager)
        return flow

    def postprocess_flow(self, flow):
        """
        Analyze the flow and fills the model with output results.
        This function is called by the worker if the flow completed succesfully.
        """
        with flow[0][0].open_gsr() as gsr:
            self.scf_data = GsData.from_gsr(gsr)

        with flow.open_final_ddb() as ddb:
            self.phonon_data = PhononData.from_ddb(ddb)

    #@classmethod
    #def get_common_queries(cls) -> list:
    #    return [
    #        {"$and": [
    #            {"scf_data.abs_pressure_gpa:": {"$gt": 2}},
    #            {"scf_data.max_force_ev_over_ang": {"$gt": 1e-6}}]
    #        },
    #    ]
"""
Pydantic Models.
"""
from __future__ import annotations

import os
import json
import inspect
import traceback
import panel as pn

from pprint import pprint
from datetime import datetime
from abc import ABC, abstractmethod
from typing import List, Optional, Any, Type
from monty import design_patterns
from pydantic import BaseModel, Field, BaseConfig, PrivateAttr #, ModelMetaClass
from pydantic.main import ModelMetaclass
from bson.objectid import ObjectId

from monty.json import MontyEncoder, MontyDecoder
from pymatgen.core import __version__ as pmg_version
from pymatgen.util.serialization import pmg_serialize
from pymatgen.util.typing import VectorLike, MatrixLike
from abipy.core.release import __version__ as abipy_version
from abipy.core.structure import Structure
from abipy.tools.tensors import Stress
from abipy.abio.inputs import AbinitInput
from abipy.electrons.ebands import ElectronBands
from abipy.panels.core import ply
from abipy.panels.viewers import JSONViewer


class AbipyDecoder(MontyDecoder):

    def process_decoded(self, d):

        if isinstance(d, dict) and "@module" in d and "@qualname" in d:
            modname = d["@module"]
            qualname = d["@qualname"]
            mod = __import__(modname, None, None, [qualname], 0)
            print("in my decoded with mod", mod)
            return getattr(mod, qualname)

        return super().process_decoded(d)



class AbipyEncoder(MontyEncoder):

    def default(self, o) -> dict:  # pylint: disable=E0202
        if inspect.isclass(o):
            print("isclass", o)
            return cls2dict(o)

        return super().default(o)


def monty_json_dumps(obj, **kwargs):
    #print("in monty json dumps with type:", type(obj))
    return json.dumps(obj, cls=AbipyEncoder, **kwargs)


def monty_json_loads(string, **kwargs):
    #print("in monty json loads with string", string)
    return json.loads(string, cls=AbipyDecoder, **kwargs)


def monty_load(obj, **kwargs):
    s = monty_json_dumps(obj)
    #print("in monty load:\n", s)
    d = monty_json_loads(s)
    print(d)
    return d


def monty_trick(obj):
    #print(f"in monty trick with {type(obj)}")
    return json.loads(monty_json_dumps(obj))


def cls2dict(cls):
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
    def from_json(cls, json_string):
        return cls(**monty_json_loads(json_string))

    def json_write(self, filepath: str, **kwargs) -> None:
        with open(filepath, "wt") as fp:
            fp.write(self.json(**kwargs))

    @classmethod
    def from_dict(cls, d):
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
    def for_localhost(cls, collection_name: str, port: int = 27017):
        return cls(host="localhost", port=port, collection_name=collection_name)

    def _repr_markdown_(self):
        return f"""

* Host: {self.host}
* Port: {self.port}
* User: {self.user}
* Database: {self.db_name}
* Collection: {self.collection_name}

"""

    def get_client(self):
        """
        Establish a connection with the MongoDB server.

        Return: MongoClient
        """
        from pymongo import MongoClient
        from pymongo.errors import ConnectionFailure

        def client_is_connected(client):
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

    def get_collection(self, collection_name=None):
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
        pn = abipanel()
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
    def from_mongo_oid(cls, oid, collection):
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
    def mongo_find_oids_models(cls, query, collection, **kwargs):

        cursor = collection.find(query, **kwargs)
        if cursor is None: return None

        oids, models = [], []
        for data in cursor:
            oids.append(data.pop("_id"))
            #print(data)
            data = AbipyDecoder().process_decoded(data)
            #print(data)
            #data = monty_load(data)
            models.append(cls(**data))

        return oids, models

    def mongo_insert(self, collection):
        """
        Insert the model in collection. Return ObjectId.
        """
        # Use plain json loads here as we want to insert a dictionary
        doc = json.loads(self.json())
        return collection.insert_one(doc).inserted_id

    def mongo_full_update_oid(self, oid, collection):
        """
        Perform a full update of the model given the ObjectId in the collection.
        """
        oid = ObjectId(oid)
        old_doc = collection.find_one({'_id': oid})
        if old_doc is None:
            raise RuntimeError(f"Cannot find document with oid: {oid}")

        new_doc = json.loads(self.json())
        old_doc.update(new_doc)

        collection.replace_one({"_id": oid}, new_doc, upsert=False)
        #collection.update_one({'_id': oid}, new_doc)

        #key = {"id": person.id}
        #self.collection.update_one(
        #    filter=key,
        #    update={"$set": person.dict()},
        #    upsert=True
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


def mongo_insert_models(models, collection):
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
    def from_table_name(cls, table_name: str):
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


class FlowData(AbipyBaseModel):

    flow_status: str = "Initialized"

    worker_name: str = Field(None, description="The name of th AbiPy worker running this Flow.")

    worker_hostname: str = Field(None, description="The name of the machine running this Flow.")

    workdir: str = Field(None, description="The directory where the Flow is executed.")

    created_at: datetime = Field(
        description="Timestamp for when this material document was first created",
        default_factory=datetime.utcnow,
    )

    completed_at: datetime = Field(None,
        description="Timestamp for the most recent calculation update for this property",
        #default_factory=datetime.utcnow,
    )

    tracebacks: List[str] = Field([], description="List of exception tracebacks")

    #abinit_version: str = Field(None, description="The version of the calculation code")


class FlowModel(MongoModel, ABC):
    """
    Base class for models associated to a Flow calculation performed by an AbiPy Worker.
    """
    flow_data: FlowData = Field(None, description="")

    pseudos_specs: PseudoDojoSpecs = Field(..., description="The input structure.")

    @classmethod
    def find_runnable_oid_models(cls, collection, limit=0):
        """
        Return list of models that are ready to run.
        """
        # NB: None not only matches itself but also matches “does not exist.”
        # Thus, querying for a key with the value None will return all documents lacking that key
        # hence we have to check that the key is None and $exists:
        query = {"flow_data": {"$in": [None], "$exists": True}}
        cursor = collection.find(query, limit=limit)
        if cursor is None: return []

        items = []
        for doc in cursor:
            oid = doc.pop("_id")
            doc = AbipyDecoder().process_decoded(doc)
            print("Found runnable doc:", doc)
            items.append((oid, cls(**doc)))

        return items

    @classmethod
    def init_collection(cls, collection):
        magic_key = "__flow_model__"
        new_dict = cls2dict(cls)

        cnt = 0
        for doc in collection.find({magic_key: {"$exists": True}}):
            oid = doc.pop("_id")
            old_dict = doc.get(magic_key)
            cnt += 1
            if new_dict != old_dict:
                raise ValueError()

        if cnt:
            return oid

        return collection.insert_one({magic_key: new_dict}).inserted_id

    @classmethod
    def get_subclass_from_collection(cls, collection):
        magic_key = "__flow_model__"

        cnt = 0
        for doc in collection.find({magic_key: {"$exists": True}}):
            oid = doc.pop("_id")
            data = doc.get(magic_key)
            cnt += 1

        if cnt == 0:
            raise RuntimeError(f"Cannot find document with magic key: {magic_key}")
        if cnt > 1:
            raise RuntimeError(f"Find {cnt} documents with magic key: {magic_key}")

        sub_class = AbipyDecoder().process_decoded(data)
        return sub_class

    @abstractmethod
    def build_flow(self, workdir, manager):
        """
        Return Flow. Must be provided by the subclass
        """

    def build_flow_and_update_collection(self, workdir, oid, collection, worker):
        """
        API used by the AbiPy Worker to build a Flow from the model.
        Wraps build_flow implemented by the subclass.
        """
        #if self.flow_data.flow_status == "errored": return
        self.flow_data = FlowData()

        # Set the workdir of the Flow in the model.
        self.flow_data.workdir = workdir
        self.flow_data.flow_status = "Initialized"

        self.flow_data.worker_name = worker.name
        from socket import gethostname
        self.flow_data.worker_hostname = gethostname()

        try:
            # Call the method provided by the user.
            return self.build_flow(workdir, worker.manager)

        except:
            self.flow_data.flow_status = "Errored"
            self.flow_data.tracebacks.append(traceback.format_exc())
            return None

        finally:
            self.mongo_full_update_oid(oid, collection)

    @abstractmethod
    def postprocess_flow(self, flow):
        """
        Postprocess the Flow. Must be provided by the subclass.
        """

    def postprocess_flow_and_update_collection(self, flow, oid, collection):
        """
        API used by the AbiPy Worker to postprocess a Flow from the model.
        Wraps postprocess_flow implented in the subclass.
        """
        #if self.flow_data.flow_status == "Errored": return
        print("in postprocess_flow")

        try:
            self.postprocess_flow(flow)
            self.flow_data.flow_status = "Completed"
        except:
            self.flow_data.flow_status = "Errored"
            self.flow_data.tracebacks.append(traceback.format_exc())
        finally:
            self.flow_data.completed_at = datetime.utcnow()
            print("postprocessing_done, status:", self.flow_data.flow_status)
            self.mongo_full_update_oid(oid, collection)

    #@abstractmethod
    #def get_common_queries(self):
    #    """
    #    Return list of MongoDB queries Flow. Used by the GUI
    #    """


class StructureData(AbipyBaseModel):

    structure: Structure = Field(..., description="Abipy Structure object.")

    @classmethod
    def from_structure(cls, structure):
        structure = Structure.as_structure(structure) #.copy()
        return cls(structure=structure)

    #def get_full_view(self):
    #    return


class GsData(AbipyBaseModel):
    """
    Ground-state results: energy, forces, stresses fermi level, gaps
    """
    max_force_ev_over_ang: float = Field(
            None,
            description="Max Cartesian force in eV/Ang. None if forces are not available.")

    pressure_gpa: float = Field(None, description="Pressure in GPa. None if forces are not available.")

    is_scf_run: float = Field(None, description="True is this a Ground-state run.")

    #cart_stress_tensor_gpa: MatrixLike = Field(None, description="Cartesian stress Tensor in GPA")
    #cart_forces_ev_ang: MatrixLike = Field(None, description="Cartesian forces in eV ang^-1")

    ebands: ElectronBands = Field(..., description="Electronic bands.")

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str):
        """
        Fill the model from the GSR filepath.
        """
        from abipy.electrons.gsr import GsrFile
        with GsrFile(gsr_filepath) as gsr:
            return cls.from_gsr(gsr)

    @classmethod
    def from_gsr(cls, gsr):
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
                max_force_ev_over_ang=gsr.max_force,
                #cart_stress_tensor_gpa=
                #cart_forces_ev_over_ang=
            ))

        return cls(**kwargs)

    #@classmethod
    #def get_common_queries(cls):
    #    return [
    #        {"$and": [{"is_scf_run:" True}, {"max_force_ev_over_ang": {"$gt": 1e-6}}]},
    #    ]


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
    will be then stored in the database collection.

    NOTES: The flow should have a single structure.
    """

    ########
    # Input
    ########

    structure_data: StructureData = Field(..., description="Input structure.")

    kppa: int = Field(1000,
        description="Defines the sampling used for the SCF run. Defaults to 1000 if not given.")

    ecut: float= Field(6, description="")

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

    scf_input: AbinitInput = Field(None, description="Input structure.")

    scf_data: GsData = Field(None, description="Results produced by the GS SCF run.")

    nscf_kpath_data: GsData = Field(None, description="Results produced by the GS NSCF run.")

    nscf_kmesh_data: GsData = Field(None, description="Results produced by the GS NSCF run.")

    def build_flow(self, workdir, manager):
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
        multi = ebands_input(self.structure_data.structure, pseudos,
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

    def get_full_view(self):
        return pn.Column(
            #self.structure_data.get_full_view(),
            ply(self.nscf_kpath_data.ebands.plotly(show=False)),
            self.scf_input._repr_html_(),
            #JSONViewer(self.json(), depth=1)
            #pn.layout.Divider(),
            sizing_mode="stretch_both",
        )


    #def get_dataframe(self, **kwargs):

    @classmethod
    def get_common_queries(cls):
        return [
            {"$and": [
                {"scf_data.pressure_gpa:": {"$gt": 2}},
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
            ddb_string = ddb.get_string(),
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

    structure_data: StructureData = Field(..., description="Input structure.")

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
        scf_input = self.scf_input
        flow = PhononFlow.from_scf_input(workdir, scf_input,
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

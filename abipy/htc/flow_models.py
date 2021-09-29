from __future__ import annotations

import traceback

import pandas as pd
import logging

from datetime import datetime
from enum import Enum
from abc import ABC, abstractmethod
from typing import List, Tuple, ClassVar, Union, TypeVar, Type, Dict  #, Optional, Any, Type,
from pydantic import Field
from bson.objectid import ObjectId
from pymongo.collection import Collection
from abipy.flowtk.tasks import Task
from abipy.flowtk.events import EventReport
from abipy.flowtk import TaskManager, Work, Flow
from abipy.abio.inputs import AbinitInput
from .base_models import AbipyModel, MongoModel, cls2dict, AbipyDecoder, QueryResults
from .structure_models import StructureData
from .pseudos_models import PseudoSpecs

logger = logging.getLogger(__name__)


class ExecStatus(str, Enum):
    """
    Possible status of the execution flow.
    """
    init = "Initialized"
    built = "Built"
    errored = "Errored"
    completed = "Completed"

    #@classmethod
    #def is_member(cls, key):
    #    return key in cls.__members__


class FlowData(AbipyModel):
    """
    Submodel storing the status the Flow and additional metadata.
    Users do not need to interact with FlowData explicitly as this task is delegated to the AbipyWorker.
    """
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


FM = TypeVar('FM', bound='FlowModel')


class FlowModel(MongoModel, ABC):
    """
    Abstract class for models associated to a Flow calculation performed by an AbiPy Worker.
    This model implements the businness logic to create a Flow from a MongoDB collection
    and post-process the results once the Flow is completed.

    Users are supposed to use this model to initialize a MongoDB collection with all
    the input arguments that will be used to generate the Flow and provide a concrete
    implementation of:

        - build_flow.
        - postprocess_flow

    The first method receives the input arguments from the MongoDB database
    and use these values to build a flow.

    The second method is invoked by the AbiPy worker when the calculation is completed.
    The function uses the Flow API to fill the ouput part of the model that
    will be then be stored in the database collection.
    The AbipyWorker uses this API to automate calculations.

    NOTES: The flow is assumed to have the same tructure and pseudopotentials.
    """

    flow_data: FlowData = Field(default_factory=FlowData,
                                description="Stores info on the Flow execution status. Set by the AbiPy Worker.")

    in_structure_data: StructureData = Field(..., description="Input structure and associated metadata. Used ")

    pseudos_specs: PseudoSpecs = Field(..., description="Pseudopotential specifications.")

    # Private class attributes
    _magic_key: ClassVar[str] = "_flow_model_"

    @classmethod
    def init_collection(cls: Type[FM], collection: Collection) -> ObjectId:
        """
        Initialize a MongoDB collection by inserting a special document containing
        the serialized FlowModel subclass.
        This function must be used when creating the collection of FlowModels for the first time.

        Return: ObiectID of the new document.
        """
        new_dict = cls2dict(cls)

        cnt, oid = 0, None
        for doc in collection.find({cls._magic_key: {"$exists": True}}):
            oid = doc.pop("_id")
            old_dict = doc.get(cls._magic_key)
            cnt += 1
            if new_dict != old_dict:
                raise TypeError(f"Cannot register new Model {new_dict}\n"
                                f"Collection is already associated to the FlowModel {old_dict}")

        if cnt:
            return oid

        return collection.insert_one({cls._magic_key: new_dict}).inserted_id

    @classmethod
    def get_subclass_from_collection(cls, collection: Collection) -> Type[FlowModel]:
        """
        Retrieve the FlowModel subclass from a MongoDB collection.
        The collections should have been initialized previously by calling `init_collection`.

        Return: FlowModel subclass.
        """
        cnt, data = 0, None
        for doc in collection.find({cls._magic_key: {"$exists": True}}):
            _ = doc.pop("_id")
            data = doc.get(cls._magic_key)
            cnt += 1

        if cnt == 0:
            raise RuntimeError(f"Cannot find document with magic key: {cls._magic_key}")
        if cnt > 1:
            raise RuntimeError(f"Found {cnt} documents with magic key: {cls._magic_key}")

        sub_class = AbipyDecoder().process_decoded(data)

        if sub_class is not cls and cls is not FlowModel:
            raise TypeError(f"Found FlowModel {sub_class} in collection\nwhile it should be {cls}")

        return sub_class

    @classmethod
    def mongo_find_completed_oids(cls, collection: Collection, **kwargs) -> List[ObjectId]:
        """
        Find all the models in collection that are completed.
        Return: list of ObjectId.
        """
        d = cls.mongo_get_status2oids(collection, **kwargs)
        return d[ExecStatus.completed]

    @classmethod
    def mongo_get_status2oids(cls, collection: Collection, **kwargs) -> Dict[ExecStatus, List[ObjectId]]:
        """
        Return dictionary mapping the Execution status to the list of list of ObjectId.
        """
        status2oids = {e: [] for e in ExecStatus}
        key = "flow_data.exec_status"
        collection.create_index(key)
        query = {key: {"$exists": True}}
        #if "add_filter" in kwargs: query.update(kwargs.get("add_filter"))
        projection = {key: 1, "_id": 1}
        cursor = collection.find(query, projection)

        for doc in cursor:
            status = ExecStatus(doc["flow_data"]["exec_status"])
            status2oids[status].append(doc["_id"])

        return status2oids

    @classmethod
    def find_runnable_oid_models(cls, collection: Collection, **kwargs: object) -> List[Tuple[ObjectId, FlowModel]]:
        """
        Return list of (oid, model) tuple with the models that are ready to run.
        """
        # NB: None not only matches itself but also matches “does not exist.”
        # Thus, querying for a key with the value None will return all documents lacking that key
        # hence we have to check that the key is None and $exists:

        # TODO: If something goes wrong when creating the Flow, flow_data is not None and flow_data.status is Init
        #query = {"flow_data": {"$in": [None], "$exists": True}}
        query = {"flow_data.exec_status": ExecStatus.init.value}
        cursor = collection.find(query, **kwargs)

        items = []
        decoder = AbipyDecoder()
        for doc in cursor:
            oid = doc.pop("_id")
            doc = decoder.process_decoded(doc)
            items.append((oid, cls(**doc)))

        return items

    @property
    def is_init(self):
        """True if the model is in `init` state"""
        return self.flow_data.exec_status == ExecStatus.init

    @property
    def is_built(self):
        """True if the model is in `built` state"""
        return self.flow_data.exec_status == ExecStatus.built

    @property
    def is_completed(self):
        """True if the model is in `completed` state"""
        return self.flow_data.exec_status == ExecStatus.completed

    @property
    def is_errored(self):
        """True if the model is in `errored` state"""
        return self.flow_data.exec_status == ExecStatus.errored

    #def create_indexes(self, collection: Collection):
    #    creates an index if an index of the same specification does not already exist.
    #    from pymongo import IndexModel, ASCENDING, DESCENDING
    #    index1 = IndexModel([("hello", DESCENDING),
    #                         ("world", ASCENDING)], name="hello_world")
    #    index2 = IndexModel([("goodbye", DESCENDING)])
    #    collection.create_indexes([index1, index2])

    @abstractmethod
    def build_flow(self, workdir: str, manager: TaskManager) -> Flow:
        """
        Use the data stored in the model to build and return a Flow to the AbipyWorker.
        Must be implemented in the concrete subclass.
        """

    def build_flow_and_update_collection(self, workdir: str, oid: ObjectId,
                                         collection: Collection, abipy_worker) -> Union[Flow, None]:
        """
        This is the API used by the AbiPyWorker to build a Flow from the model.
        Wraps build_flow implemented by the subclass to perform extra operations on the
        model, exception handling and the final insertion in the collection.
        """
        #self.flow_data = flow_data = FlowData()
        flow_data = self.flow_data

        # Set the workdir of the Flow in the model.
        flow_data.workdir = workdir
        flow_data.worker_name = abipy_worker.name
        from socket import gethostname
        flow_data.worker_hostname = gethostname()

        try:
            # Call the method provided by the concrete class to return the Flow.
            flow = self.build_flow(workdir, abipy_worker.manager)
            flow_data.exec_status = ExecStatus.built
            return flow

        except Exception:
            flow_data.exec_status = ExecStatus.errored
            exc_string = traceback.format_exc()
            logger.critical(exc_string)
            flow_data.tracebacks.append(exc_string)
            return None

        finally:
            # Update document in the MongoDB collection.
            self.mongo_full_update_oid(oid, collection)

    @abstractmethod
    def postprocess_flow(self, flow: Flow) -> None:
        """
        Postprocess the Flow and update the model with output results.
        Must be implemented in the concrete subclass.
        """

    def postprocess_flow_and_update_collection(self, flow: Flow, oid: ObjectId,
                                               collection: Collection) -> None:
        """
        API used by the AbiPy Worker to postprocess a Flow from the model.
        Wraps the postprocess_flow method implented in the subclass to add extra
        operations on the model, error handling and finally the insertion in the collection.
        """
        flow_data = self.flow_data
        try:
            self.postprocess_flow(flow)
            flow_data.exec_status = ExecStatus.completed
        except Exception:
            flow_data.exec_status = ExecStatus.errored
            flow_data.tracebacks.append(traceback.format_exc())
        finally:
            flow_data.completed_at = datetime.utcnow()
            #print("in postprocessing_flow_and_update_collection with status:", flow_data.exec_status)
            self.mongo_full_update_oid(oid, collection)

    @classmethod
    def mongo_get_crystal_systems_incoll(cls, collection: Collection, **kwargs) -> List[int]:
        """
        Return list of crystal systems in collection
        kwargs are passed to pymongo `create_index`
        """
        key = "in_structure_data.crystal_system"
        collection.create_index(key, **kwargs)
        return [k for k in collection.distinct(key) if k is not None]

    @classmethod
    def mongo_find_by_crystal_system(cls, crystal_system: str, collection: Collection, **kwargs) -> QueryResults:
        """
        Filter documents in the collection according to the crystal system.
        kwargs are passed to mongo.find.
        """
        key = "in_structure_data.crystal_system"
        collection.create_index(key)
        query = {key: crystal_system}
        return cls.mongo_find(query, collection, **kwargs)

    @classmethod
    def mongo_get_spg_numbers_incoll(cls, collection: Collection, **kwargs) -> List[int]:
        """
        Return list if all space group numbers found in the collection.
        """
        key = "in_structure_data.spg_number"
        collection.create_index(key, **kwargs)
        return [k for k in collection.distinct(key) if k is not None]

    @classmethod
    def mongo_find_by_spg_number(cls, spg_number: int, collection: Collection, **kwargs) -> QueryResults:
        """
        Filter documents in the collection according to the space group number.
        kwargs are passed to mongo.find.
        """
        key = "in_structure_data.spg_number"
        collection.create_index(key)
        query = {key: int(spg_number)}
        return cls.mongo_find(query, collection, **kwargs)

    @classmethod
    def mongo_find_by_formula(cls, formula_pretty: str, collection: Collection, **kwargs) -> QueryResults:
        """
        Filter documents in the collection according to the formula.
        kwargs are passed to pymongo `collection.find`.
        """
        key = "in_structure_data.formula_pretty"
        collection.create_index(key)
        query = {key: formula_pretty}
        return cls.mongo_find(query, collection, **kwargs)

    @classmethod
    @abstractmethod
    def get_common_queries(cls) -> List[dict]:
        """
        Return list of dictionaries with the MongoDB queries typically used to filter results.
        Empty list if no suggestion is available. Mainly used by the panel-based GUI.
        Must be implemented in the subclass.
        """

    @abstractmethod
    def get_panel_view(self):
        """Return panel object with a view of the model."""

    #@abstractmethod
    #@classmethod
    #def get_aggregate_panel_ui(cls, collection: Collection):
    #    """Return panel object with a view of the model."""

    @classmethod
    def mongo_aggregate_in_structures(cls, collection: Collection) -> pd.DataFrame:
        oids = cls.mongo_find_completed_oids(collection)
        projection = {"_id": 1, "in_structure_data": 1}
        rows = []
        for oid in oids:
            doc = collection.find_one({"_id": oid}, projection)
            structure_data = StructureData(**doc["in_structure_data"])
            d = structure_data.dict()
            d.pop("structure")
            d["_id"] = doc["_id"]
            rows.append(d)

        return pd.DataFrame(rows).set_index("_id")


class _NodeData(AbipyModel):

    node_id: str = Field(..., description="Node identifier")

    node_class: str = Field(..., description="Class of the Node")

    status: str = Field(..., description="Status of the Node")

    @classmethod
    def get_data_for_node(cls, node):
        return dict(
            node_id=node.node_id,
            node_class=node.__class__.__name__,
            status=str(node.status),
        )


class TaskData(_NodeData):
    """
    Data Model associated to an AbiPy |Task|.
    """

    input: AbinitInput = Field(..., description="Abinit input object")
    #output_str: str = Field(..., description="Abinit input")
    #log_str: str = Field(..., description="Abinit input")

    report: EventReport = Field(..., description="Number of warnings")

    num_warnings: int = Field(..., description="Number of warnings")

    num_errors: int = Field(..., description="Number of errors")

    num_comments: int = Field(..., description="Number of comments")

    #num_restarts: int = Field(..., description="Number of comments")

    #cpu_time: float = Field(..., description="CPU-time in seconds")

    #wall_time: float = Field(..., description="Wall-time in seconds")

    #mpi_nprocs: int = Field(..., description="Number of MPI processes")

    #omp_nthreads: int = Field(..., description="Number of OpenMP threads. -1 if not used.")

    @classmethod
    def from_task(cls, task: Task) -> TaskData:
        """
        Build the model from an AbiPy |Task|
        """
        data = cls.get_data_for_node(task)

        data["input"] = task.input

        # TODO: Handle None!
        report = task.get_event_report()
        data["report"] = report
        #data["report"] = report.as_dict()
        for a in ("num_errors", "num_comments", "num_warnings"):
            data[a] = getattr(report, a)

        return cls(**data)


class WorkData(_NodeData):
    """
    Data Model associated to an AbiPy |Work|
    """

    tasks: List[TaskData] = Field(..., description="List of TaskData")

    @classmethod
    def from_work(cls, work: Work) -> WorkData:
        """
        Build the model from a AbiPy |Work| instance.
        """
        data = cls.get_data_for_node(work)
        data["tasks"] = [TaskData.from_task(task) for task in work]
        return cls(**data)


class __FlowData(_NodeData):
    """
    Document associated to an AbiPy |Flow|
    """
    works: List[WorkData] = Field(..., description="")

    @classmethod
    def from_flow(cls, flow: Flow) -> __FlowData:
        """
        Build the model from an AbiPy |Flow|
        """
        data = cls.get_data_for_node(flow)
        data["works"] = [WorkData.from_work(work) for work in flow]
        return cls(**data)

    #def __getitem__(self, name):
    #    try:
    #        # Dictionary-style field of super
    #        return super().__getitem__(name)
    #    except KeyError:
    #        # Assume int or slice
    #        try:
    #            return self.works[name]
    #        except IndexError:
    #            raise

    #def delete(self):
    #    # Remove GridFs files.
    #    for work in self.works:
    #        work.outfiles.delete()
    #        #work.delete()
    #        for task in work:
    #            #task.delete()
    #            task.outfiles.delete()

    #    self.delete()

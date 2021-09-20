from __future__ import annotations

import traceback
import panel as pn
import logging

from datetime import datetime
from enum import Enum
from abc import ABC, abstractmethod
from typing import List, Tuple, ClassVar, Union, TypeVar, Type  #, Dict, Optional, Any, Type,
from pydantic import Field
from bson.objectid import ObjectId
from pymongo.collection import Collection
#from pymatgen.util.typing import VectorLike, MatrixLike
#from abipy.tools.tensors import Stress
from abipy.abio.inputs import AbinitInput
from abipy.panels.core import ply
from abipy.panels.viewers import JSONViewer
from abipy.flowtk import TaskManager, Flow, PhononFlow
from .base_models import AbipyModel, MongoModel, cls2dict, AbipyDecoder, QueryResults #, AbipyEncoder
from .structure_models import StructureData
from .pseudos_models import PseudoSpecs
from .gs_models import GsData
from .dfpt_models import PhononData
#from .worker import AbipyWorker

logger = logging.getLogger(__name__)


class ExecStatus(str, Enum):
    """
    Possible status of the execution flow.
    """
    init = "Initialized"
    errored = "Errored"
    completed = "Completed"


class FlowData(AbipyModel):

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
    Base class for models associated to a Flow calculation performed by an AbiPy Worker.
    This model implements the businness logic to create a Flow from the database,
    execute it and post-process the results.
    The AbipyWorker uses this API to automate calculations.
    """

    flow_data: FlowData = Field(None, description="")

    input_structure_data: StructureData = Field(..., description="Input structure.")

    pseudos_specs: PseudoSpecs = Field(..., description="PseudoPotential Specifications.")

    # Private class attributes
    _magic_key: ClassVar[str] = "__flow_model__"

    @classmethod
    def init_collection(cls: Type[FM], collection: Collection) -> ObjectId:
        """
        Initialize a MongoDB collection by insert a special document
        containing the serialized FlowModel class.
        Return ObiectID of the new document.
        This function must be used when crearing the collection of FlowModels for the first time.
        """
        new_dict = cls2dict(cls)

        cnt, oid = 0, None
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
    def get_subclass_from_collection(cls, collection: Collection) -> Type[FlowModel]:
        """
        Retrieve the FlowModel subclass from Collection
        The collections should have been initialized previously by calling `init_collection`.
        Return FlowModel subclass.
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
        return sub_class

    @classmethod
    def find_runnable_oid_models(cls, collection: Collection, **kwargs: object) -> List[Tuple[ObjectId, FlowModel]]:
        """
        Return list of (oid, model) tuple with the models that are ready to run.
        """
        # NB: None not only matches itself but also matches “does not exist.”
        # Thus, querying for a key with the value None will return all documents lacking that key
        # hence we have to check that the key is None and $exists:
        query = {"flow_data": {"$in": [None], "$exists": True}}
        cursor = collection.find(query, **kwargs)
        if cursor is None:
            return []

        items = []
        decoder = AbipyDecoder()
        for doc in cursor:
            oid = doc.pop("_id")
            doc = decoder.process_decoded(doc)
            #print("Found runnable doc:", doc)
            items.append((oid, cls(**doc)))

        return items

    @abstractmethod
    def build_flow(self, workdir: str, manager: TaskManager) -> Flow:
        """
        Use the data stored in the model to build and return a Flow to the AbipyWorker.
        Must be implemented in the concrete subclass.
        """

    def build_flow_and_update_collection(self, workdir: str,
                                         oid: ObjectId,
                                         collection: Collection,
                                         abipy_worker) -> Union[Flow, None]:
        """
        This is the API used by the AbiPyWorker to build a Flow from the model.
        Wraps build_flow implemented by the subclass to perform extra operations on the
        model, exception handling and the final insertion in the collection.
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
            # Call the method provided by the concrete class to return the Flow.
            return self.build_flow(workdir, abipy_worker.manager)

        except Exception:
            flow_data.exec_status = ExecStatus.errored
            exc_string = traceback.format_exc()
            logger.critical(exc_string)
            flow_data.tracebacks.append(exc_string)
            return None

        finally:
            # Insert model in the collection.
            self.mongo_full_update_oid(oid, collection)

    @abstractmethod
    def postprocess_flow(self, flow: Flow) -> None:
        """
        Postprocess the Flow.
        Must be implemented in the concrete subclass.
        """

    def postprocess_flow_and_update_collection(self, flow: Flow,
                                               oid: ObjectId,
                                               collection: Collection) -> None:
        """
        API used by the AbiPy Worker to postprocess a Flow from the model.
        Wraps the postprocess_flow method implented in the subclass to add extra
        operations on the model, error handling and finally the insertion in the collection.
        """
        flow_data = self.flow_data
        #if flow_data.exec_status == ExecStatus.errored: return
        try:
            self.postprocess_flow(flow)
            flow_data.exec_status = ExecStatus.completed
        except Exception:
            flow_data.exec_status = ExecStatus.errored
            flow_data.tracebacks.append(traceback.format_exc())
        finally:
            flow_data.completed_at = datetime.utcnow()
            print("postprocessing_done, status:", flow_data.exec_status)
            self.mongo_full_update_oid(oid, collection)

    @classmethod
    @abstractmethod
    def get_common_queries(cls) -> List[dict]:
        """
        Return list of dictionaries with the MongoDB queries typically used to filter results.
        Empty list if not suggestion is available. Mainly used by the panel-based GUI.
        Must be implemented in the subclass.
        """

    @abstractmethod
    def get_panel_view(self):
        """Return panel object with a view of the model."""

    @classmethod
    def mongo_find_by_spg_number(cls, spg_number: int, collection: Collection, **kwargs) -> QueryResults:
        """
        Filter documents in the collection according to the space group number.
        kwargs are passed to pymong collection.find.
        """
        query = {"input_structure_data.spg_number": int(spg_number)}
        return cls.mongo_find(query, collection, **kwargs)

    @classmethod
    def mongo_find_by_formula(cls, formula_pretty: str, collection: Collection, **kwargs) -> QueryResults:
        """
        Filter documents in the collection according to the formula.
        kwargs are passed to pymongo `collection.find`.
        """
        query = {"input_structure_data.formula_pretty": formula_pretty}
        return cls.mongo_find(query, collection, **kwargs)


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

    def build_flow(self, workdir: str, manager: TaskManager) -> Flow:
        """
        Build an AbiPy Flow using the input data available in the model and return it.

        Args:
            workdir: Working directory provided by the caller.
            manager: |TaskManager| object.

        Return: |Flow| object.
        """
        from abipy.abio.factories import ebands_input
        #ecut = 6
        pseudos = self.pseudos_specs.get_pseudos()

        multi = ebands_input(self.input_structure_data.structure, pseudos,
                             kppa=self.kppa, nscf_nband=None, ndivsm=self.ndivsm,
                             #ecut=self.ecut, pawecutdg=None,
                             scf_nband=None, accuracy="normal",
                             spin_mode=self.spin_mode,
                             smearing=self.smearing, charge=self.charge,
                             scf_algorithm=None, dos_kppa=self.dos_kppa,
                             )

        multi.set_vars(paral_kgb=self.paral_kgb)

        self.scf_input, nscf_input = multi.split_datasets()
        from abipy.flowtk.flows import bandstructure_flow

        # TODO: Dos!
        return bandstructure_flow(workdir, self.scf_input, nscf_input, manager=manager)

    def postprocess_flow(self, flow: Flow) -> None:
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

    def get_panel_view(self):
        """
        Return panel object with a view of the model.
        """
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
            #self.input_structure_data.get_panel_view(),
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
    def get_common_queries(cls) -> List[dict]:
        """
        Return list of dictionaries with the MongoDB queries typically used to filter results.
        Empty list if not suggestion is available. Mainly used by the panel-based GUI.
        """
        return [
            {"$and": [
                {"scf_data.abs_pressure_gpa:": {"$gt": 2}},
                {"scf_data.max_force_ev_over_ang": {"$gt": 1e-6}}]
            },
        ]


class PhononFlowModel(FlowModel):
    """
    This model defines the input arguments used to build a Flow for phonon calculations
    as well as submodels for storing the final results in the MongoDB collection.

    Users are supposed to use this model to initialize a MongoDB collection with all
    the input arguments that will be used to generate the flow and provide a concrete
    implementation of:

        - build_flow.
        - postprocess_flow
        - get_common_queries

    The first method receives the input arguments from the MongoDB database
    and use these values to build a flow.

    The second method is invoked by the AbiPy worker when the calculation is completed.
    The function uses the Flow API to fill the ouput part of the model that
    will be then stored in the database collection.
    """

    ########
    # Input
    ########

    scf_input: AbinitInput = Field(..., description="Input structure.")
    #ph_ngqpt: Tuple[int, int, int]

    with_becs: bool = Field(..., description="Compute Born effective charges.")
    with_quad: bool = Field(..., description="Activate calculation of dynamical quadrupoles.")
    with_flexoe: bool = Field(False, description="Activate computation of flexoelectric tensor.")

    ########
    # Output
    ########

    scf_data: GsData = Field(None, description="Results produced by the GS SCF run.")

    phonon_data: PhononData = Field(None, description="Results produced by the GS SCF run.")

    def build_flow(self, workdir: str, manager: TaskManager) -> PhononFlow:
        """
        Build an AbiPy Flow using the input data available in the model and return it.

        Args:
            workdir: Working directory provided by the caller.
            manager: |TaskManager| object.
        """
        # Create flow to compute all the independent atomic perturbations
        # corresponding to a [4, 4, 4] q-mesh.
        # Electric field and Born effective charges are also computed.

        flow = PhononFlow.from_scf_input(workdir, self.scf_input,
                                         ph_ngqpt=(2, 2, 2),
                                         #ph_ngqpt=(4, 4, 4),
                                         with_becs=self.with_becs, with_quad=self.with_quad,
                                         with_flexoe=self.with_flexoe, manager=manager)
        return flow

    def postprocess_flow(self, flow: PhononFlow) -> None:
        """
        Analyze the flow and fills the model with output results.
        This function is called by the worker if the flow completed succesfully.
        """
        with flow[0][0].open_gsr() as gsr:
            self.scf_data = GsData.from_gsr(gsr)

        with flow.open_final_ddb() as ddb:
            self.phonon_data = PhononData.from_ddb(ddb)

    @classmethod
    def get_common_queries(cls) -> List[dict]:
        return [
            #{"$and": [
            #    {"scf_data.abs_pressure_gpa:": {"$gt": 2}},
            #    {"scf_data.max_force_ev_over_ang": {"$gt": 1e-6}}]
            #},
        ]

    def get_panel_view(self):
        raise NotImplementedError()

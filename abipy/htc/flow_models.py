from __future__ import annotations

import traceback
import panel as pn

from datetime import datetime
from enum import Enum
from abc import ABC, abstractmethod
from typing import List, Tuple, ClassVar, Union #, Dict, Optional, Any, Type,
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
from .pseudos_models import PseudoDojoSpecs
from .gs_models import GsData
from .dfpt_models import PhononData


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

    @abstractmethod
    def build_flow(self, workdir: str, manager: TaskManager) -> Flow:
        """
        Return Flow. Must be provided by the subclass
        """

    def build_flow_and_update_collection(self, workdir: str,
                                         oid: ObjectId,
                                         collection: Collection,
                                         abipy_worker) -> Union[Flow, None]:
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
    def postprocess_flow(self, flow: Flow) -> None:
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
            return QueryResults.empty_from_query(query)

        oids, models = [], []
        decoder = AbipyDecoder()
        for doc in cursor:
            oid = doc.pop("_id")
            doc = decoder.process_decoded(doc)
            models.append(cls(**doc))
            oids.append(oid)

        return QueryResults(oids, models, query)

    @classmethod
    def find_by_spg_number(cls, spg_number: int, collection, **kwargs) -> QueryResults:
        query = {"input_structure_data.spg_number": int(spg_number)}
        return cls.find(query, collection, **kwargs)

    @classmethod
    def find_by_formula(cls, formula: str, collection, **kwargs) -> QueryResults:
        query = {"input_structure_data.formula_pretty": formula}
        return cls.find(query, collection, **kwargs)

    #@abstractmethod
    #def get_common_queries(self): -> List[dict]:
    #    """
    #    Return list of MongoDB queries Flow. Used by the GUI
    #    """


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
    def get_common_queries(cls) -> List[dict]:
        return [
            {"$and": [
                {"scf_data.abs_pressure_gpa:": {"$gt": 2}},
                {"scf_data.max_force_ev_over_ang": {"$gt": 1e-6}}]
            },
        ]


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

    def build_flow(self, workdir, manager: TaskManager) -> PhononFlow:
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

    #@classmethod
    #def get_common_queries(cls) -> List[dict]:
    #    return [
    #        {"$and": [
    #            {"scf_data.abs_pressure_gpa:": {"$gt": 2}},
    #            {"scf_data.max_force_ev_over_ang": {"$gt": 1e-6}}]
    #        },
    #    ]

from __future__ import annotations

#import logging

from typing import List  # , Tuple, ClassVar, Union, TypeVar, Type, Dict  #, Optional, Any, Type,
from pydantic import Field, validator
from abipy.abio.inputs import AbinitInput
from abipy.flowtk import TaskManager, PhononFlow
from .gs_models import GsData
from .dfpt_models import PhononData
from .flow_models import FlowModel

#logger = logging.getLogger(__name__)


class _BasePhononFlowModel(FlowModel):
    """
    This is the base class for models performing phonon band structure calculations.

    It defines ...
    as well as submodels storing the output results.
    """

    ########
    # Input
    ########

    scf_input: AbinitInput = Field(..., description="Input structure.")

    with_becs: bool = Field(..., description="Activate computation of Born effective charges.")

    with_quad: bool = Field(..., description="Activate calculation of dynamical quadrupoles.")

    with_flexoe: bool = Field(False, description="Activate computation of the flexoelectric tensor.")

    ########
    # Output
    ########

    scf_data: GsData = Field(None, description="Results produced by the GS SCF run.")

    phonon_data: PhononData = Field(None, description="Results produced by the Phonon calculation")

    #gsr_gfsd: GridFsDesc = Field(None, description="Link to a GridFS entry.")
    #ddb_gfsd: GridFsDesc = Field(None, description="Link to a GridFS entry.")
    with_dvdb: bool = Field(True, description="False if the DVDB file should not be added to GridFs")
    #dvdb_gfsd: GridFsDesc = Field(None, description="Link to a GridFS entry.")

    def build_flow(self, workdir: str, manager: TaskManager) -> PhononFlow:
        """
        Build an AbiPy Flow using the input data available in the model and return it.

        Args:
            workdir: Working directory provided by the caller.
            manager: |TaskManager| object.
        """
        # TODO: Add option to compute bands, get_dvdb filepath!
        return PhononFlow.from_scf_input(workdir, self.scf_input,
                                         ph_ngqpt=(2, 2, 2),
                                         #ph_ngqpt=(4, 4, 4),
                                         with_becs=self.with_becs, with_quad=self.with_quad,
                                         with_flexoe=self.with_flexoe, manager=manager)

    def postprocess_flow(self, flow: PhononFlow) -> None:
        """
        Analyze the flow and fills the model with output results.
        This function is called by the worker if the flow completed successfully.
        """
        with flow[0][0].open_gsr() as gsr:
            self.scf_data = GsData.from_gsr(gsr)

        with flow.open_final_ddb() as ddb:
            #self.ddb_gfsd = GridFsDesc(filepath=ddb.filepath)
            self.phonon_data = PhononData.from_ddb(ddb)

        #if self.with_dvdb
        #self.dvdb_gfsd = GridFsDesc(filepath=dvdb_filepath)

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


class PhononFlowModelWithParams(_BasePhononFlowModel):
    """
    This model generates the SCF input from meta parameters that will be stored in the database.
    """

    #kppa: float
    #qppa: float

    def __init__(self, **data):
        super().__init__(**data)

        self.make_input_from_params()
        if self.scf_input is None:
            raise ValueError("scf_input should be defined after the call to `make_input_from_params")


    def make_input_from_params(self):
        if self.scf_input is not None:
            raise ValueError(f"scf_input should be None when calling {self.__class__.__name__} "
                             "as scf_input is automatically generated from the meta-parameters.")
        pseudos = self.pseudos_specs.get_pseudos()
        structure = self.in_structure_data.structure
        # Call the factory function to build the scf_input from the meta-parameters.
        #self.scf_input =


class PhononFlowModelWithInput(_BasePhononFlowModel):
    """
    This model requires a scf AbinitInput constructed by the user.
    It is more flexible than the parameter-based version but queries related to input variables
    are more complex to perform as one should use the |AbinitInput| dictionary.
    """

    #ph_ngqpt: Tuple[int, int, int]

    @validator("scf_input")
    def check_scf_input_is_not_none(cls, v):
        if v is None:
            raise ValueError(f"Model {cls.__name__} requires an scf_input object.")
        return v

    #@root_validator
    #def check_inputs(cls, values):
    #    scf_input, nscf_input = values.get('scf_input'), values.get('nscf_input')
    #    if scf_input is None or nscf_input is None:
    #        raise ValueError(f"Constructing a {cls.__name__} model requires both scf_input and nscf_input")
    #    return values

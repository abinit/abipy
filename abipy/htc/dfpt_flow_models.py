from __future__ import annotations

from typing import List  # , Tuple, ClassVar, Union, TypeVar, Type, Dict  #, Optional, Any, Type,
from pydantic import Field, validator
from abipy.abio.inputs import AbinitInput
from abipy.flowtk import TaskManager, PhononFlow
from .base_models import MongoConnector, GfsFileDesc
from .gs_models import GsData
from .dfpt_models import PhononData
from .flow_models import FlowModel, PresetQuery


class _BasePhononFlowModel(FlowModel):
    """
    This is the base class for models performing phonon band structure calculations.

    It defines ...
    as well as submodels storing the output results.
    """

    ########
    # Input
    ########

    scf_input: AbinitInput = Field(..., description="Input for the initial SCF run")

    with_becs: bool = Field(..., description="Activate computation of Born effective charges.")

    with_quad: bool = Field(..., description="Activate calculation of dynamical quadrupoles.")

    with_flexoe: bool = Field(False, description="Activate computation of the flexoelectric tensor.")

    ########
    # Output
    ########

    scf_data: GsData = Field(None, description="Results produced by the GS SCF run.")

    phonon_data: PhononData = Field(None, description="Results produced by the Phonon calculation")

    #gsr_gfsd: GfsFileDesc = Field(None, description="Link to a GridFS entry.")
    #with_den: bool = Field(True, description="Set it to True to save the DEN file in GridFS.")
    #with_pot: bool = Field(True, description="Set it to True to save the POT file in GridFS.")

    #gsden_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the DVDB Fortran file from GridFS.")
    #gspot_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the DVDB Fortran file from GridFS.")

    #ddb_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the DDB file from GridFS.")

    with_dvdb: bool = Field(True, description="False if the DVDB file should not be added to GridFs")

    dvdb_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the DVDB Fortran file from GridFS.")

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
                                         with_becs=self.with_becs,
                                         with_quad=self.with_quad,
                                         with_flexoe=self.with_flexoe,
                                         manager=manager,
                                         )

    def postprocess_flow(self, flow: PhononFlow, mng_connector: MongoConnector) -> None:
        """
        Analyze the flow and fills the model with output results.
        MongoConnector should be used only to insert files in GridFs as the final insertion is done by the caller.
        This function is called by the worker if the flow completed successfully.
        """
        with flow[0][0].open_gsr() as gsr:
            self.scf_data = GsData.from_gsr(gsr, mng_connector, with_gsr=False)

        # Add DEN and POT files
        #self.gsden_gfsd = mng_connector.gfs_put_filepath(gsden_filepath)
        #self.gspot_gfsd = mng_connector.gfs_put_filepath(gspot_filepath)

        with flow.open_final_ddb() as ddb:
            self.phonon_data = PhononData.from_ddb(ddb, mng_connector)

        #if self.with_dvdb:
        #    dvdb_filepath = None
        #    self.dvdb_gfsd = mng_connector.gfs_put_filepath(dvdb_filepath)

    @classmethod
    def get_preset_queries(cls) -> List[PresetQuery]:
        """
        Return list of dictionaries with the MongoDB queries typically used to filter documents for this model.
        Empty list if no suggestion is available.
        """
        return [
            PresetQuery.for_large_forces_or_high_pressure("scf_data", cls),
        ]

    def get_panel_view(self, mng_connector: MongoConnector):
        """
        Return panel object with a view of the model.
        """
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

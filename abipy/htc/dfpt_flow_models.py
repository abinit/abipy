from __future__ import annotations

from typing import List  # , Tuple, ClassVar, Union, TypeVar, Type, Dict  #, Optional, Any, Type,
from pydantic import Field, validator
from abipy.abio.inputs import AbinitInput
from abipy.flowtk import TaskManager, PhononFlow
from .base_models import MongoConnector, GfsFileDesc
from .gs_models import ScfData
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

    with_becs: bool = Field(..., description="Activate computation of Born effective charges.")

    with_quad: bool = Field(..., description="Activate calculation of dynamical quadrupoles.")

    with_flexoe: bool = Field(False, description="Activate computation of the flexoelectric tensor.")

    with_gsden: bool = Field(True, description="Set it to True to save the DEN file in GridFS.")

    with_gspot: bool = Field(True, description="Set it to True to save the POT file in GridFS.")

    with_dvdb: bool = Field(True, description="False if the DVDB file should not be added to GridFs")

    ########################################
    # Input/Output depending on the subclass
    ########################################

    scf_input: AbinitInput = Field(..., description="Input for the initial SCF run")

    ########
    # Output
    ########

    scf_data: ScfData = Field(None, description="Results produced by the initial GS SCF run.")

    phonon_data: PhononData = Field(None, description="Results produced by the Phonon calculation")

    gsden_gfsd: GfsFileDesc = Field(None, description="Pointer to the GS DEN file in GridFS.")

    gspot_gfsd: GfsFileDesc = Field(None, description="Pointer to the GS KS POT file in GridFS.")

    ddb_gfsd: GfsFileDesc = Field(None, description="Pointer to the DDB file from GridFS.")

    dvdb_gfsd: GfsFileDesc = Field(None, description="Pointer to the DVDB Fortran file in GridFS.")

    def build_flow(self, workdir: str, worker: AbipyWorker) -> PhononFlow:
        """
        Build an AbiPy Flow in workdir using the input data available in the model and return it.
        """
        # TODO: Add option to compute bands, get_dvdb filepath!
        flow = PhononFlow.from_scf_input(workdir, self.scf_input,
                                         ph_ngqpt=(2, 2, 2),
                                         #ph_ngqpt=(4, 4, 4),
                                         with_becs=self.with_becs,
                                         with_quad=self.with_quad,
                                         with_flexoe=self.with_flexoe,
                                         )

        # Make sure Abinit writes the output files required by the model
        scf_input = flow[0][0].input
        if self.with_gsden: scf_input.set_vars(prtden=1)
        if self.with_gspot: scf_input.set_vars(prtpot=1)

        return flow

    def postprocess_flow(self, flow: PhononFlow, worker: AbipyWorker) -> None:
        """
        Analyze the flow and fills the model with output results.
        MongoConnector should be used only to insert files in GridFs as the final insertion is done by the caller.
        This function is called by the worker if the flow completed successfully.
        """
        mng_connector = worker.mng_connector

        scf_task = flow[0][0]
        with scf_task.open_gsr() as gsr:
            self.scf_data = ScfData.from_gsr(gsr, mng_connector, with_gsr=False)

        # Add DEN and POT files to GridFs.
        if self.with_gsden:
            self.gsden_gfsd = mng_connector.gfs_put_filepath(scf_task.outdir.need_abiext("DEN"))

        if self.with_gspot:
            self.gspot_gfsd = mng_connector.gfs_put_filepath(scf_task.outdir.need_abiext("POT"))

        with flow.open_final_ddb() as ddb:
            self.phonon_data = PhononData.from_ddb(ddb, mng_connector)
            self.ddb_gfsd = mng_connector.gfs_put_filepath(ddb.filepath)

        if self.with_dvdb:
            # Add to GridDF the DVDB required for e-ph computations.
            ph_work = flow[1]
            dvdb_filepath = ph_work.outdir.need_abiext("DVDB")
            self.dvdb_gfsd = mng_connector.gfs_put_filepath(dvdb_filepath)

    @classmethod
    def get_preset_queries(cls) -> List[PresetQuery]:
        """
        Return list of dictionaries with the MongoDB queries typically used to filter documents for this model.
        Empty list if no suggestion is available.
        """
        return [
            PresetQuery.for_large_forces_or_high_pressure("scf_data", cls),
            #PresetQuery.unstable_ph_modes("phonon_data", cls),
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

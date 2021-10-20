"""
FlowModels for ground-state computations.
"""
from __future__ import annotations

import panel as pn

from abc import ABC
from typing import List
from pydantic import Field, root_validator
from abipy.abio.inputs import AbinitInput
from abipy.panels.core import ply
#from abipy.panels.viewers import JSONViewer
from abipy.abio.factories import ebands_input
from abipy.flowtk import TaskManager, Flow
from abipy.flowtk.flows import bandstructure_flow
from abipy.htc.base_models import AbipyModel, MongoConnector, GfsFileDesc
from abipy.htc.structure_models import StructureData
from abipy.htc.flow_models import FlowModel, PresetQuery
from abipy.htc.gs_models import ScfData, NscfData, RelaxData
from abipy.htc.worker import AbipyWorker


class _EbandsFlowModel(FlowModel, ABC):
    """
    This is the base class for models performing band structure calculations.
    It defines the three input files for the SCF/NSCF/DOS calculation
    as well as sub models for storing the output results.

    The input part i.e. the set of abinit variables are customized by the concrete subclass.
    """

    with_gsr: bool = Field(False, description="True if GSR files should be saved in GridFS.")

    #with_den: bool = Field(False, description="True if DEN file should be saved in GridFS.")

    # These inputs can be generated either manually or from meta-parameters + factory functions.

    scf_input: AbinitInput = Field(None, description="Input for the GS-SCF run")

    nscf_input: AbinitInput = Field(None, description="Input for the NSCF run with k-path")

    dos_input: AbinitInput = Field(None, description="Input for the NSCF run with k-mesh (optional")

    #####################################
    # Output (filled in postprocess_flow)
    #####################################

    scf_data: ScfData = Field(None, description="Results produced by the SCF run.")

    #scf_den_gfsd: GfsFileDesc = Field(None, description="Pointer to the DEN file in GridFs. None if file is not stored")

    nscf_kpath_data: NscfData = Field(None, description="Results produced by the NSCF run with a k-path.")

    nscf_kmesh_data: NscfData = Field(None, description="Results produced by the NSCF run with a k-mesh.")

    def postprocess_flow(self, flow: Flow, worker: AbipyWorker) -> None:
        """
        Analyze the flow and fills the model with output results.
        MongoConnector should be used only to insert files in GridFs as the final insertion is done by the caller.
        This function is called by the AbiPy Worker if the flow completed successfully.
        """
        mng_connector = worker.mng_connector
        scf_task = flow[0][0]

        with scf_task.open_gsr() as gsr:
            self.scf_data = ScfData.from_gsr(gsr, mng_connector, self.with_gsr)

        with flow[0][1].open_gsr() as gsr:
            self.nscf_kpath_data = NscfData.from_gsr(gsr, mng_connector, self.with_gsr)

        if self.dos_input is not None:
            with flow[0][2].open_gsr() as gsr:
                self.nscf_kmesh_data = NscfData.from_gsr(gsr, mng_connector, self.with_gsr)

        #if self.with_den:
        #    den_filepath = scf_task.outdir.need_abiext("DEN")
        #    self.scf_den_gfsd = mng_connector.gfs_put_filepath(den_filepath)

    @classmethod
    def get_preset_queries(cls) -> List[PresetQuery]:
        """
        Return list of MongoDB queries typically used to filter documents.
        """
        return [
            PresetQuery.for_large_forces_or_high_pressure("scf_data", cls),
            PresetQuery.for_collinear_magnetic("scf_data", cls),
        ]

    def get_panel_view(self, mng_connector: MongoConnector):
        """
        Return panel object with a view of the model.
        """
        title = self.in_structure_data.get_title()
        structure = self.in_structure_data.structure
        a, b, c = structure.lattice.abc
        alpha, beta, gamma = structure.lattice.angles
        header = f"""\
## {title}

- Lattice lengths: a = {a:.6f}, b = {b:.6f}, c = {c:.6f} Ang
- Lattice angles: α = {alpha:.3f}, β = {beta:.3f}, ɣ = {gamma:.3f} degrees
"""

        if self.is_completed:
            header += f"""\
- Pressure: {self.scf_data.pressure_gpa:.3f} GPa
- Max |Force|: {self.scf_data.max_force_ev_over_ang:.8f} eV/Ang
- Energy: {self.scf_data.energy:.4f} eV
- Energy per atom: {self.scf_data.energy_per_atom:.4f} eV
"""
            ebands_kpath = self.nscf_kpath_data.ebands
            #ebands_kpath = mng_connector.gfs_get_mson_obj(self.nscf_kpath_data.ebands_gfsd)
            plotly_bands = ply(ebands_kpath.plotly(show=False))
        else:
            plotly_bands = pn.pane.Alert(f"Bands are not available because exec_status is `{self.flow_data.exec_status}`")

        return pn.Column(
            #self.in_structure_data.get_panel_view(mongo_connector),
            header,
            pn.Row(
                plotly_bands,
                pn.pane.HTML(self.scf_input._repr_html_()),
            ),
            #"### MongoDB Document",
            #JSONViewer(self.json(), depth=1),
            pn.layout.Divider(),
            sizing_mode="stretch_both",
        )

    #@classmethod
    #def mongo_aggregate_egaps(cls, collection: Collection) -> pd.DataFrame:
        #oids = cls.mng_find_completed_oids(collection)
        #    projection = [
        #        "in_structure_data
        #    }
        #
        #    fund_gap_projection = [
        #        "scf_data.fundamental_gap",
        #        "nscf_kpath_data.fundamental_gap",
        #        "nscf_kmesh_data.fundamental_gap",
        #    ]

        #    direct_gap_projection = {
        #        "scf_data.direct_gap",
        #        "nscf_kpath_data.direct_gap",
        #        "nscf_kmesh_data.direct_gap",
        #    }

        #    projection.extend(fundamental_gap_projection + direct_gap_projection)

        #    oids = cls.mng_find_completed_oids(collection)
        #    rows = []
        #    for oid in oids:
        #        doc = collection.find_one({"_id": oid}, projection)
        #        structure_data = AbipyDecoder().process_decoded(doc["in_structure_data"])
        #        #structure_data = cls.decode(doc["in_structure_data"])
        #        row = structure_data.dict4pandas()
        #        # Here I need a tool to access foo.bar instead of d["foo"]["bar"]
        #        row["fund_gap"] = min(doc[key] for key in fundamental_gap_projection)
        #        row["direct_gap"] = min(doc[key] for key in direct_gap_projection)
        #        rows.append(row)

        #    return pd.DataFrame(rows)


from abipy.htc.protocol import Protocol


class HasProtocol(AbipyModel):
    """
    Mixin class for FlowModels with a protocol object that is used to generate input files
    and pass meta parameters to the flow.
    """

    protocol: Protocol = Field(..., description="Protocol used to generate input files")

    @classmethod
    def from_structure_and_protocol(cls, structure: Structure, protocol: Protocol, **kwargs):
        """
        Build a FlowModel from a structure and the protocol.
        """
        data = dict(
                in_structure_data=StructureData.from_structure(structure),
                pseudos_specs=protocol.pseudos_specs,
                protocol=protocol
        )
        data.update(kwargs)

        return cls(**data)


class EbandsFlowModelWithProtocol(_EbandsFlowModel, HasProtocol):
    """
    This model defines the input arguments used to build a Flow for band structure calculations
    as well as the sub models used to store the final results.
    """

    def build_flow(self, workdir: str, worker: AbipyWorker) -> Flow:
        """
        Build an AbiPy Flow in `workdir` using the input data available in the model and return it.
        """
        structure = self.in_structure_data.structure
        self.scf_input, self.nscf_input = self.protocol.get_ebands_input(structure)

        #pseudos = self.pseudos_specs.get_pseudos()
        #multi = ebands_input(structure, pseudos,
        #                     kppa=self.kppa, nscf_nband=None, ndivsm=self.ndivsm,
        #                     #ecut=6, pawecutdg=None,
        #                     scf_nband=None, accuracy="normal",
        #                     spin_mode=self.spin_mode,
        #                     smearing=self.smearing, charge=self.charge,
        #                     scf_algorithm=None, dos_kppa=self.dos_kppa,
        #                     )

        #multi.set_vars(paral_kgb=self.paral_kgb)

        #if self.dos_kppa is not None:
        #    self.scf_input, self.nscf_input, self.dos_input = multi.split_datasets()
        #else:
        #    self.scf_input, self.nscf_input = multi.split_datasets()

        return bandstructure_flow(workdir, self.scf_input, self.nscf_input)
                                  #dos_inputs=self.dos_input)


class EbandsFlowModelWithParams(_EbandsFlowModel):
    """
    This model defines the input arguments used to build a Flow for band structure calculations
    as well as the sub models used to store the final results.
    """

    ##################
    # Input parameters
    ##################
    kppa: int = Field(1000, description="Defines the sampling used for the SCF run. Defaults to 1000 if not given.")

    spin_mode: str = Field("polarized", description="Spin polarization")

    charge: float = Field(0.0, description="Electronic charge added to the unit cell.")

    smearing: str = Field("fermi_dirac:0.1 eV", description="Smearing technique.")

    # TODO ndivsm < 0 is much better.
    ndivsm: int = Field(2, description="Number of divisions used to sample the smallest segment of the k-path.")

    dos_kppa: int = Field(None,
                          description="Scalar or List of integers with the number of k-points per atom " +
                                      "to be used for the computation of the DOS (None if DOS is not wanted")

    paral_kgb: int = Field(0, description="")

    def build_flow(self, workdir: str, worker: AbipyWorker) -> Flow:
        """
        Build an AbiPy Flow in `workdir` using the input data available in the model and return it.
        """
        pseudos = self.pseudos_specs.get_pseudos()
        structure = self.in_structure_data.structure

        #self.protocol.get_gs_scf_input(structure)
        #self.scf_input, self.nscf_input = self.protocol.get_ebands_input(structure)

        multi = ebands_input(structure, pseudos,
                             kppa=self.kppa, nscf_nband=None, ndivsm=self.ndivsm,
                             #ecut=6, pawecutdg=None,
                             scf_nband=None, accuracy="normal",
                             spin_mode=self.spin_mode,
                             smearing=self.smearing, charge=self.charge,
                             scf_algorithm=None, dos_kppa=self.dos_kppa,
                             )

        multi.set_vars(paral_kgb=self.paral_kgb)

        if self.dos_kppa is not None:
            self.scf_input, self.nscf_input, self.dos_input = multi.split_datasets()
        else:
            self.scf_input, self.nscf_input = multi.split_datasets()

        return bandstructure_flow(workdir, self.scf_input, self.nscf_input,
                                  dos_inputs=self.dos_input)


class EbandsFlowModelWithInputs(_EbandsFlowModel):
    """
    More flexible class for band structure calculations that requires |AbinitInput| objects
    explicitly generated by the user before calling __init__.
    """

    @root_validator
    def check_inputs(cls, values):
        """Enforce the presence of scf_input and nscf_input once the model is created."""
        scf_input, nscf_input = values.get("scf_input"), values.get("nscf_input")
        if scf_input is None or nscf_input is None:
            raise ValueError(f"Building a {cls.__name__} model requires both `scf_input` and `nscf_input`")
        return values

    def build_flow(self, workdir: str, worker: AbipyWorker) -> Flow:
        """
        Build an AbiPy Flow in `workdir` using the input data available in the model and return it.
        """
        return bandstructure_flow(workdir, self.scf_input, self.nscf_input,
                                  dos_inputs=self.dos_input)


#class RelaxFlowModelWithProtocol(FlowModel):
#
#    with_hist: bool = Field(True, description="True if HIST.nc file should be saved in GridFS.")
#
#    with_gsr: bool = Field(False, description="True if GSR.nc file should be saved in GridFS.")
#
#    relax_data: RelaxData = Field(None, description="Results of the structural relaxation.")
#
#    def build_flow(self, workdir: str, worker: AbipyWorker) -> Flow:
#        """
#        Build an AbiPy Flow in `workdir` using the input data available in the model and return it.
#        """
#        pseudos = self.pseudos_specs.get_pseudos()
#        structure = self.in_data_structure.structure
#        relax_input = self.protocol.get_relax_input(structure)
#        return Flow
#
#    def postprocess_flow(self, flow: Flow, worker: AbipyWorker) -> None:
#        """
#        Analyze the flow and fills the model with output results.
#        MongoConnector should be used only to insert files in GridFs as the final insertion is done by the caller.
#        This function is called by the AbiPy Worker if the flow completed successfully.
#        """
#        mng_connector = worker.mng_connector
#        self.relax_data = RelaxData.from_hist_and_gsr_filepaths(
#                                              hist_filepath, gsr_filepath,
#                                              mng_connector, self.with_gsr, self.with_hist)

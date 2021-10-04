"""
pydantic models storing results of ground-state calculations.
"""
from __future__ import annotations

#from abc import ABC, abstractmethod
from pydantic import Field
from typing import List
from abipy.electrons.gsr import GsrFile
from abipy.electrons.ebands import ElectronBands
from abipy.dynamics.hist import HistFile
from .base_models import AbipyModel, GfsFileDesc, GfsDesc

# TODO Add dictionary with metavariables TAGS
# mpid:

class _ModelWithGsr(AbipyModel):
    """
    Base model with an ElectronBands object and other important quantities such as
    the Fermi level and the band gaps.
    """

    ebands: ElectronBands = Field(..., description="Electronic bands.")

    ebands_gfsd: GfsDesc = Field(None, description="Metadata needed to retrieve the ebands object from GridFS.")

    #fundamental_band_gap_ev: float = Field(..., description="Fundamental band gap in eV.")

    #direct_band_gap_ev: float = Field(..., description="Direct band gap in eV.")

    #fermie_ev: float = Field(..., description="Fermi level in eV.")

    gsr_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the GSR.nc file from GridFS. "
                                                    "None if the GSR.nc file should not be stored")

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str, mng_connector: MongoConnector, with_gsr: bool):
        """
        Fill the model from the GSR filepath.
        """
        with GsrFile(gsr_filepath) as gsr:
            return cls.from_gsr(gsr, mng_connector, with_gsr)


class NscfData(_ModelWithGsr):
    """
    Model for NSCF calculations.
    """

    @classmethod
    def from_gsr(cls, gsr: GsrFile, mng_connector: MongoConnector, with_gsr: bool) -> NscfData:
        """
        Fill the model from a |GsrFile|
        """
        if gsr.is_scf_run:
            raise RuntimeError("Expecting a GSR file produced by a NSCF calculation.")

        data = dict(
            ebands=gsr.ebands,
            ebands_gfsd=mng_connector.gfs_put_mson_obj(gsr.ebands),
        )

        if with_gsr:
            data["gsr_gfsd"] = mng_connector.gfs_put_filepath(gsr.filepath)

        return cls(**data)


class _ScfBaseModel(_ModelWithGsr):
    """
    Model for SCF calculations providing energy, forces, and stress tensor besides
    the electronic band dispersion.
    """

    max_force_ev_over_ang: float = Field(...,
        description="Max absolute Cartesian force in eV/Ang. None if forces are not available.")

    abs_pressure_gpa: float = Field(..., description="Pressure in GPa. Absolute value.")

    pressure_gpa: float = Field(..., description="Pressure in GPa. NB: Value with sign!.")

    energy: float = Field(..., description="The total DFT energy in eV")

    energy_per_atom: float = Field(..., description="The DFT energy per atom in eV")

    #gap_ev: float = Field(None, description="The DFT bandgap for the last calculation")

    #forces: List[Vector3D] = Field(
    #    [], description="Forces on atoms from the last calculation"
    #)

    #cart_stress_tensor_gpa: MatrixLike = Field(None, description="Cartesian stress Tensor in GPA")

    #cart_forces_ev_over_ang: MatrixLike = Field(None, description="Cartesian forces in eV ang^-1")

    #cart_stress: Matrix3D = Field(...
    #    [], description="Stress on the unitcell from the last calculation"
    #)


class GsData(_ScfBaseModel):
    """
    Ground-state results: energy, forces, stress tensor, Fermi level, band gaps.
    """

    @classmethod
    def from_gsr(cls, gsr: GsrFile, mng_connector: MongoConnector, with_gsr: bool) -> NscfData:
        """
        Fill the model from a |GsrFile|
        """
        if not gsr.is_scf_run:
            raise RuntimeError("Expecting a GSR file produced by a SCF calculation")

        # TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
        #gsr.ebands.structure.remove_site_property("cartesian_forces")
        data = dict(
            ebands=gsr.ebands,
            ebands_gfsd=mng_connector.gfs_put_mson_obj(gsr.ebands),
        )

        data.update(dict(
            pressure_gpa=gsr.pressure,
            abs_pressure_gpa=abs(gsr.pressure),
            max_force_ev_over_ang=gsr.max_force,
            energy=float(gsr.energy),
            energy_per_atom=float(gsr.energy_per_atom),
            #cart_stress_tensor_gpa=??,
            #cart_forces_ev_over_ang=??.
        ))

        if with_gsr:
            data["gsr_gfsd"] = mng_connector.gfs_put_filepath(gsr.filepath)

        return cls(**data)


class RelaxData(_ScfBaseModel):
    """
    Model for Structural relaxation calculations
    Store energy, forces, stress tensor, Fermi level, band gaps for the final configuration
    (from _ScfBaseModel) as well as additional lists with the evolution of the most important
    physical parameters as each relaxation step.
    """

    num_iterations: int = Field(..., description="Number of relaxation steps.")

    etotal_ev_hist: List[float] = Field(..., description="Energies in eV at the different relaxation steps.")

    pressure_gpa_hist: List[float] = Field(..., description="Pressure in GPa at the different relaxation steps.")

    angles_hist: List[List[float]] = Field(...,
                                           description="Lattice angles in degrees at the different relaxation steps.")

    lenghts_ang_hist: List[List[float]] = Field(...,
                                                description="Lattice lengths in Ang at the different relaxation steps.")

    volume_ang3_hist: List[float] = Field(...,
                                          description="Unit cell volume in Ang**3 at the different relaxation steps.")

    #with_hist: bool = Field(False, description="Upload HIST.nc file to GridFS")

    gsr_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the GSR.nc file from GridFS. "
                                                    "None if the GSR.nc is not stored")

    hist_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the HIST.nc file from GridFS. "
                                                     "None if the HIST.nc is not stored")

    #@classmethod
    #def from_gsr(cls, gsr: GsrFile, mng_connector: MongoConnector, with_gsr: bool) -> NscfData:

    @classmethod
    def from_hist_and_gsr_filepaths(cls, hist_filepath: str, gsr_filepath: str,
                                    mng_connector: MongoConnector, with_gsr: bool, with_hist: bool) -> RelaxData:
        """
        Fill the model from a HIST.nc filepath.
        """
        with HistFile(hist_filepath) as hist, GsrFile(gsr_filepath) as gsr:
            return cls.from_hist_and_gsr(hist, gsr, mng_connector, with_gsr, with_hist)

    @classmethod
    def from_hist_and_gsr(cls, hist: HistFile, gsr: GsrFile, mng_connector: MongoConnector,
                          with_gsr: bool,  with_hist: bool) -> RelaxData:
        """
        Fill the model from a |HistFile|
        """
        data = dict()
        ## TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
        #gsr.ebands.structure.remove_site_property("cartesian_forces")
        #data = dict(
        #    ebands=gsr.ebands,
        #    ebands_gfsd=mng_connector.gfs_put_mson_obj(gsr.ebands),
        #    is_scf_run=gsr.is_scf_run,
        #)

        #data.update(dict(
        #    pressure_gpa=gsr.pressure,
        #    abs_pressure_gpa=abs(gsr.pressure),
        #    max_force_ev_over_ang=gsr.max_force,
        #    energy=float(gsr.energy),
        #    energy_per_atom=float(gsr.energy_per_atom),
        #    #cart_stress_tensor_gpa=
        #    #cart_forces_ev_over_ang=
        #))

        if with_gsr:
            data["gsr_gfsd"] = mng_connector.gfs_put_filepath(gsr.filepath)

        if with_hist:
            data["hist_gfsd"] = mng_connector.gfs_put_filepath(hist.filepath)

        return cls(**data)


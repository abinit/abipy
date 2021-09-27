"""pydantic models for storing results of ground-state calculations."""
from __future__ import annotations

#from abc import ABC, abstractmethod
from pydantic import Field
from typing import List
from abipy.electrons.gsr import GsrFile
from abipy.electrons.ebands import ElectronBands
from abipy.dynamics.hist import HistFile
from .base_models import AbipyModel

# TODO Add dictionary with metavariables TAGS
# mpid:


class _ModelWithEbands(AbipyModel):
    """
    Base model with an ElectronBands object and other important quantiies such as
    the Fermi level and the band gaps.
    """

    ebands: ElectronBands = Field(..., description="Electronic bands.")

    #fundamental_band_gap_ev: float = Field(..., description="Fundamental band gap in eV.")

    #direct_band_gap_ev: float = Field(..., description="Direct band gap in eV.")

    #fermie_ev: float = Field(..., description="Fermi level in eV.")


class NscfData(_ModelWithEbands):
    """
    Model for NSCF calculations.
    """

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str):
        """
        Fill the model from the GSR filepath.
        """
        with GsrFile(gsr_filepath) as gsr:
            return cls.from_gsr(gsr)

    @classmethod
    def from_gsr(cls, gsr: GsrFile) -> NscfData:
        """
        Fill the model from a |GsrFile|
        """
        if gsr.is_scf_run:
            raise RuntimeError("Expecting a GSR file produced by a NSCF calculation.")

        kwargs = dict(
            ebands=gsr.ebands,
        )

        return cls(**kwargs)


class _ScfBaseModel(_ModelWithEbands):
    """
    Model for SCF calculations providing energy, forces, and stress tensor besides
    the electronic band dispersion.
    """

    max_force_ev_over_ang: float = Field(
        None,
        description="Max absolute Cartesian force in eV/Ang. None if forces are not available.")

    abs_pressure_gpa: float = Field(None, description="Pressure in GPa. Absolute Value.")
    pressure_gpa: float = Field(None, description="Pressure in GPa. NB: Value with sign!.")

    #cart_stress_tensor_gpa: MatrixLike = Field(None, description="Cartesian stress Tensor in GPA")
    #cart_forces_ev_over_ang: MatrixLike = Field(None, description="Cartesian forces in eV ang^-1")

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


class GsData(_ScfBaseModel):
    """
    Ground-state results: energy, forces, stress tensor, Fermi level, band gaps.
    """

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str):
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
        if not gsr.is_scf_run:
            raise RuntimeError("Expecting a GSR file produced by a SCF calculation")

        # TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
        #gsr.ebands.structure.remove_site_property("cartesian_forces")
        kwargs = dict(
            ebands=gsr.ebands,
        )

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
                                                description="Lattice lenghts in Ang at the different relaxation steps.")

    volume_ang3_hist: List[float] = Field(...,
                                          description="Unit cell volume in Ang**3 at the different relaxation steps.")

    @classmethod
    def from_hist_gsr_filepaths(cls, hist_filepath: str, gsr_filepath: str) -> RelaxData:
        """
        Fill the model from a HIST.nc filepath.
        """
        with HistFile(hist_filepath) as hist_file, GsrFile(gsr_filepath) as gsr_file:
            return cls.from_hist_gsr(hist_file, gsr_file)

    @classmethod
    def from_hist_gsr(cls, hist: HistFile, gsr_file: GsrFile) -> RelaxData:
        """
        Fill the model from a |HistFile|
        """
        raise NotImplementedError()
        #    # TODO
        #    kwargs = dict()
        #    # TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
        #    gsr.ebands.structure.remove_site_property("cartesian_forces")
        #    kwargs = dict(
        #        ebands=gsr.ebands,
        #        is_scf_run=gsr.is_scf_run,
        #    )

        #    kwargs.update(dict(
        #        pressure_gpa=gsr.pressure,
        #        abs_pressure_gpa=abs(gsr.pressure),
        #        max_force_ev_over_ang=gsr.max_force,
        #        energy=float(gsr.energy),
        #        energy_per_atom=float(gsr.energy_per_atom),
        #        #cart_stress_tensor_gpa=
        #        #cart_forces_ev_over_ang=
        #    ))

        #    return cls(**kwargs)
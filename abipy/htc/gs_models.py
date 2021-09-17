"""pydantic models for storing results of ground-state calculations."""
from __future__ import annotations

from pydantic import Field
#from typing import List
from abipy.electrons.gsr import GsrFile
from abipy.electrons.ebands import ElectronBands
#from abipy.dynamics.hist import HistFile
from .base_models import AbipyModel

# TODO Add dictionary with metavariables TAGS
# mpid:


class GsData(AbipyModel):
    """
    Ground-state results: energy, forces, stress tensor, Fermi level, band gaps.
    """
    max_force_ev_over_ang: float = Field(
        None,
        description="Max absolute Cartesian force in eV/Ang. None if forces are not available.")

    abs_pressure_gpa: float = Field(None, description="Pressure in GPa. Absolute Value.")
    pressure_gpa: float = Field(None, description="Pressure in GPa. NB: Value with sign!.")

    is_scf_run: bool = Field(None, description="True is this a Ground-state run.")

    #cart_stress_tensor_gpa: MatrixLike = Field(None, description="Cartesian stress Tensor in GPA")
    #cart_forces_ev_over_ang: MatrixLike = Field(None, description="Cartesian forces in eV ang^-1")

    ebands: ElectronBands = Field(..., description="Electronic bands.")

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

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str) -> GsData:
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
        # TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
        #gsr.ebands.structure.remove_site_property("cartesian_forces")
        kwargs = dict(
            ebands=gsr.ebands,
            is_scf_run=gsr.is_scf_run,
        )

        if gsr.is_scf_run:
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


#class RelaxData(AbipyModel):
#    """
#    Ground-state results: energy, forces, stress tensor, Fermi level, band gaps.
#    """
#    max_force_ev_over_ang: float = Field(
#        None,
#        description="Max absolute Cartesian force in eV/Ang. None if forces are not available.")
#
#    abs_pressure_gpa: float = Field(None, description="Pressure in GPa. Absolute Value.")
#    pressure_gpa: float = Field(None, description="Pressure in GPa. NB: Value with sign!.")
#
#    #cart_stress_tensor_gpa: MatrixLike = Field(None, description="Cartesian stress Tensor in GPA")
#    #cart_forces_ev_ang: MatrixLike = Field(None, description="Cartesian forces in eV ang^-1")
#
#    ebands: ElectronBands = Field(..., description="Electronic bands.")
#
#    energy: float = Field(
#        None, description="The total DFT energy in eV"
#    )
#    energy_per_atom: float = Field(
#        None, description="The DFT energy per atom in eV"
#    )
#    #bandgap: float = Field(None, description="The DFT bandgap for the last calculation")
#    #forces: List[Vector3D] = Field(
#    #    [], description="Forces on atoms from the last calculation"
#    #)
#    #stress: Matrix3D = Field(
#    #    [], description="Stress on the unitcell from the last calculation"
#    #)
#
#    num_iterations: int = Field(None, description="")
#    etotal_hist: List[float]
#    pressure_gpa_hist: List[float]
#    angles_hist: List[List[float]]
#    lenghts_hist: List[List[float]]
#    volume_hist: List[float]
#
#    @classmethod
#    def from_hist_filepath(cls, hist_filepath: str) -> RelaxData:
#        """
#        Fill the model from the Hist filepath.
#        """
#        with HistFile(hist_filepath) as hist:
#            return cls.from_hist(hist)
#
#    @classmethod
#    def from_hist(cls, hist: HistFile) -> RelaxData:
#        """
#        Fill the model from a |HistFile|
#        """
#        # TODO
#        kwargs = dict()
#        # TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
#        #gsr.ebands.structure.remove_site_property("cartesian_forces")
#        #kwargs = dict(
#        #    ebands=gsr.ebands,
#        #    is_scf_run=gsr.is_scf_run,
#        #)
#        #
#        #kwargs.update(dict(
#        #    pressure_gpa=gsr.pressure,
#        #    abs_pressure_gpa=abs(gsr.pressure),
#        #    max_force_ev_over_ang=gsr.max_force,
#        #    energy=float(gsr.energy),
#        #    energy_per_atom=float(gsr.energy_per_atom),
#        #    #cart_stress_tensor_gpa=
#        #    #cart_forces_ev_over_ang=
#        #))
#
#        return cls(**kwargs)
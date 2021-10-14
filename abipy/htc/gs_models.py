"""
pydantic models storing results of ground-state calculations.
"""
from __future__ import annotations

from pydantic import Field
from typing import List
from abipy.electrons.gsr import GsrFile
from abipy.electrons.ebands import ElectronBands
from abipy.dynamics.hist import HistFile
from .base_models import AbipyModel, GfsFileDesc, GfsDesc
#from .structure_models import StructureData


class _ModelWithGsr(AbipyModel):
    """
    This is the base class that defines the dimensions and the results of an "abstract"
    GS calculations i.e. an Abinit computation that produces a GSR file.

    This class is not supposed to be instanciated directly.
    It just contains Fields that are common to specialized subclasses corresponding to Scf/NSCF/Relax runs.
    """

    ebands: ElectronBands = Field(..., description="Electronic bands.")

    nkpt: int = Field(..., description="Number of k-points in ebands")

    nsppol: int = Field(..., description="Number of independent spins in WFs")

    nspinor: int = Field(..., description="Number of spinor components in WFs")

    nspden: int = Field(..., description="Number of independent spin components in density")

    fermie_ev: float = Field(..., description="Fermi energy in eV.")

    is_metal: bool = Field(..., description="True if system is metallic")

    fundamental_gap_ev: float = Field(None, description="Fundamental band gap in eV (min over spins, if any).")

    direct_gap_ev: float = Field(None, description="Direct band gap in eV (min over spins, if any).")

    ebands_gfsd: GfsDesc = Field(None, description="Pointer to the ElectronBands object in GridFS")

    gsr_gfsd: GfsFileDesc = Field(None, description="Pointer to the GSR.nc file stored in GridFS. "
                                                    "None if the GSR.nc is not stored")

    @classmethod
    def from_gsr_filepath(cls, gsr_filepath: str, mng_connector: MongoConnector, with_gsr: bool):
        """
        Fill the model from a GSR filepath.
        """
        with GsrFile(gsr_filepath) as gsr:
            return cls.from_gsr(gsr, mng_connector, with_gsr)

    @classmethod
    def get_base_data(cls, gsr: GsrFile, mng_connector: MongoConnector, with_gsr: bool) -> dict:

        #attrs = [
        #    "ecut", "pawecutdg",
        #    "tsmear", "nkpt",
        #    "nsppol", "nspinor", "nspden",
        #]

        fundamental_gap_ev, direct_gap_ev = 0.0, 0.0
        if not gsr.ebands.is_metal:
            fundamental_gap_ev = min(gap.energy for gap in gsr.ebands.fundamental_gaps)
            direct_gap_ev = min(gap.energy for gap in gsr.ebands.direct_gaps)

        data = dict(
            nkpt=len(gsr.ebands.kpoints),
            nsppol=gsr.nsppol,
            nspinor=gsr.nspinor,
            nspden=gsr.nspden,
            ebands=gsr.ebands,
            ebands_gfsd=mng_connector.gfs_put_mson_obj(gsr.ebands),
            is_metal=gsr.ebands.is_metal,
            fermie_ev=gsr.ebands.fermie,
            fundamental_gap_ev=fundamental_gap_ev,
            direct_gap_ev=direct_gap_ev,
        )

        if with_gsr:
            # Add the GSR file to GridFS.
            data["gsr_gfsd"] = mng_connector.gfs_put_filepath(gsr.filepath)

        return data


class NscfData(_ModelWithGsr):
    """
    Model for NSCF calculations.
    """

    #dos_ef
    #dos_values

    @classmethod
    def from_gsr(cls, gsr: GsrFile, mng_connector: MongoConnector, with_gsr: bool) -> NscfData:
        """
        Fill the model from a |GsrFile|.
        """
        if gsr.is_scf_run:
            raise RuntimeError("Expecting a GSR file produced by a NSCF calculation.")

        data = cls.get_base_data(gsr, mng_connector, with_gsr)

        return cls(**data)


class ScfData(_ModelWithGsr):
    """
    Model for SCF calculations with k-points in the IBZ. This kind of sampling allows one
    to compute energies, forces, stress tensor and collinear magnetization.
    """

    max_force_ev_over_ang: float = Field(...,
        description="Max absolute Cartesian force in eV/Ang. None if forces are not available")

    abs_pressure_gpa: float = Field(..., description="Pressure in GPa. Absolute value")

    pressure_gpa: float = Field(..., description="Pressure in GPa. NB: Value with sign!")

    energy: float = Field(..., description="The total DFT energy in eV")

    energy_per_atom: float = Field(..., description="The DFT energy per atom in eV")

    #forces: List[Vector3D] = Field(
    #    [], description="Forces on atoms from the last calculation"
    #)

    #cart_stress_tensor_gpa: MatrixLike = Field(None, description="Cartesian stress Tensor in GPA")

    #cart_forces_ev_over_ang: MatrixLike = Field(None, description="Cartesian forces in eV ang^-1")

    #cart_stress: Matrix3D = Field(...
    #    [], description="Stress on the unitcell from the last calculation"
    #)

    collinear_mag: float = Field(None,
        description="Total collinear magnetization in Bohr magneton as the difference "
                    "between the spin up and spin down densities. None if undefined e.g. nspinor 2")

    @classmethod
    def from_gsr(cls, gsr: GsrFile, mng_connector: MongoConnector, with_gsr: bool) -> ScfData:
        """
        Fill the model from a |GsrFile|
        """
        if not gsr.is_scf_run:
            raise RuntimeError("Expecting a GSR file produced by a SCF calculation")

        data = cls.get_base_data(gsr, mng_connector, with_gsr)

        # TODO: Use "cartesian_forces_eV/Ang" and get rid of ArrayWithUnits
        #gsr.ebands.structure.remove_site_property("cartesian_forces")
        #data = dict(
        #    ebands=gsr.ebands,
        #    ebands_gfsd=mng_connector.gfs_put_mson_obj(gsr.ebands),
        #)

        # Add entries that are available only if SCF run.
        try:
            collinear_mag = gsr.ebands.get_collinear_mag()
        except ValueError:
            collinear_mag = None

        data.update(dict(
            pressure_gpa=gsr.pressure,
            abs_pressure_gpa=abs(gsr.pressure),
            max_force_ev_over_ang=gsr.max_force,
            energy=float(gsr.energy),
            energy_per_atom=float(gsr.energy_per_atom),
            #cart_stress_tensor_gpa=??,
            #cart_forces_ev_over_ang=??.
            collinear_mag=collinear_mag,
        ))

        return cls(**data)


class RelaxData(ScfData):
    """
    Model for structural relaxations.
    Store energy, forces, stress tensor, Fermi level, band gaps for the final configuration
    inherited from ScfData as well as additional lists with the evolution of the most important
    physical parameters at each relaxation step for plotting purposes.
    """

    #final_structure_data = StructureData = Field(..., description="Data computed from the final structure")

    num_steps: int = Field(..., description="Number of relaxation steps.")

    #etotal_ev_hist: List[float] = Field(..., description="Energies in eV at the different relaxation steps")

    #pressure_gpa_hist: List[float] = Field(..., description="Pressure in GPa at the different relaxation steps")

    #angles_hist: List[List[float]] = Field(...,
    #                                       description="Lattice angles in degrees at the different relaxation steps")

    #lenghts_ang_hist: List[List[float]] = Field(...,
    #                                            description="Lattice lengths in Ang at the different relaxation steps")

    #volume_ang3_hist: List[float] = Field(...,
    #                                      description="Unit cell volume in Ang**3 at the different relaxation steps")


    hist_gfsd: GfsFileDesc = Field(None, description="Pointer to the HIST.nc stored in GridFS. None if not stored")

    @classmethod
    def from_gsr(cls, gsr: GsrFile, mng_connector: MongoConnector, with_gsr: bool) -> NscfData:
        raise RuntimeError("For {cls.__name__}, one should call from_hist_and_gsr instead of from_gsr")

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
                          with_gsr: bool, with_hist: bool) -> RelaxData:
        """
        Fill the model from a |HistFile|
        """
        data = ScfData.from_gsr(gsr, mng_connector, with_gsr).dict()

        #an = self.get_relaxation_analyzer()
        #app("Volume change in percentage: %.2f%%" % (an.get_percentage_volume_change() * 100))
        #d = an.get_percentage_lattice_parameter_changes()
        #vals = tuple(d[k] * 100 for k in ("a", "b", "c"))
        #app("Percentage lattice parameter changes:\n\ta: %.2f%%, b: %.2f%%, c: %2.f%%" % vals)
        #an.get_percentage_bond_dist_changes(max_radius=3.0)

        data.update(dict(
            num_steps=hist.num_steps,
            #final_structure_data=StructureData.from_structure(hist.final_structure),
            #etotal_ev_hist=hist.etotals.tolist(),
            #pressure_gpa_hist=,
            #angles_hist=,
            #lenghts_ang_hist=,
            #volume_ang3_hist=,
        ))

        if with_hist:
            data["hist_gfsd"] = mng_connector.gfs_put_filepath(hist.filepath)

        return cls(**data)


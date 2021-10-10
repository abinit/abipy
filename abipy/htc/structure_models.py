"""
Pydantic models storing Structure objects and associated metadata

Some of these models are inspired to emmet or we try to use field names and models
that can be easily converted into each other if needed.
"""
from __future__ import annotations

from enum import Enum
from typing import Dict, Any, List
from pydantic import Field
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from abipy.core.structure import Structure
from abipy.htc.base_models import AbipyModel


class CrystalSystem(str, Enum):
    """
    The crystal system of the lattice
    """
    tri = "Triclinic"
    mono = "Monoclinic"
    ortho = "Orthorhombic"
    tet = "Tetragonal"
    trig = "Trigonal"
    hex_ = "Hexagonal"
    cubic = "Cubic"


class StructureData(AbipyModel):
    """
    Store structural info set and symmetry metadata.
    """

    #mpid: str = Field(None, description="MP identifier")
    # Optional or in CustomDict?

    #source_id: str  = Field(..., "Identifier e.g. mp-123")
    #source_type: str = Field(..., "source type e.g. MP for the materials project")
    #meta: Dict[str, Any] = Field(None)

    #hall:

    crystal_system: CrystalSystem = Field(
        None, title="Crystal System", description="The crystal system for this lattice"
    )

    spg_symbol: str = Field(
        None,
        title="Space Group Symbol",
        description="The spacegroup symbol for the lattice",
    )

    spg_number: int = Field(
        None,
        title="Space Group Number",
        description="The spacegroup number for the lattice",
    )

    point_group: str = Field(
        None, title="Point Group Symbol", description="The point group for the lattice"
    )

    symprec: float = Field(
        None,
        title="Symmetry Finding Precision",
        description="The precision given to spglib to determine the symmetry of this lattice",
    )

    angle_tolerance: float = Field(
        None,
        description="Tolerance on the angles",
    )

    spglib_version: str = Field(None, title="SPGLib version")

    # Structure metadata
    nsites: int = Field(None, description="Total number of sites in the structure")

    elements: List[Element] = Field(None, description="List of elements in the material")

    nelements: int = Field(None, title="Number of Elements")

    composition: Composition = Field(None, description="Full composition for the material")

    composition_reduced: Composition = Field(
        None,
        title="Reduced Composition",
        description="Simplified representation of the composition",
    )

    formula_pretty: str = Field(
        None,
        title="Pretty Formula",
        description="Cleaned representation of the formula",
    )

    formula_anonymous: str = Field(
        None,
        title="Anonymous Formula",
        description="Anonymized representation of the formula",
    )

    chemsys: str = Field(
        None,
        title="Chemical System",
        description="dash-delimited string of elements in the material",
    )

    volume: float = Field(
        None,
        title="Volume",
        description="Total volume for this structure in Angstroms^3",
    )

    density: float = Field(None, title="Density", description="Density in grams per cm^3")

    density_atomic: float = Field(
        None,
        title="Packing Density",
        description="The atomic packing density in atoms per cm^3",
    )

    structure: Structure = Field(..., description="Abipy Structure object.")

    custom_data: Dict[str, Any] = Field(None, description="Dictionary containing custom data set by the user.")

    @classmethod
    def from_mpid(cls, material_id: str) -> StructureData:
        """
        Initialize the model from a MaterialsProject id.
        """
        structure = Structure.from_mpid(material_id, final=True)
        return cls.from_structure(structure)

    @classmethod
    def from_structure(cls, structure: Structure) -> StructureData:
        """
        Initialize the model from a Structure.
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, spglib

        # 0.1 is the value used in Materials Project. 5.0 is the default value for angles.
        symprec, angle_tolerance = 0.1, 5.0
        sg = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
        #symmetry: Dict[str, Any] = {"symprec": symprec}
        if not sg.get_symmetry_dataset():
            symprec, angle_tolerance = 1e-3, 51
            sg = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)

        #symmetry.update(
        symmetry = {
                "spg_symbol": sg.get_space_group_symbol(),
                "spg_number": sg.get_space_group_number(),
                "point_group": sg.get_point_group_symbol(),
                "crystal_system": CrystalSystem(sg.get_crystal_system().title()),
                #"hall": sg.get_hall(),
                "spglib_version": spglib.__version__,
                "symprec": symprec,
                "angle_tolerance": angle_tolerance,
            }
        #)

        comp = structure.composition.remove_charges()
        elsyms = sorted(set([str(e.symbol) for e in comp.elements]))
        #symmetry = SymmetryData.from_structure(structure)

        data = {
            "nsites": structure.num_sites,
            "elements": elsyms,
            "nelements": len(elsyms),
            "composition": comp,
            "composition_reduced": comp.reduced_composition,
            "formula_pretty": comp.reduced_formula,
            "formula_anonymous": comp.anonymized_formula,
            "chemsys": "-".join(elsyms),
            "volume": structure.volume,
            "density": structure.density,
            "density_atomic": structure.volume / structure.num_sites,
            #"symmetry": symmetry,
        }

        structure = Structure.as_structure(structure)  #.copy()
        symmetry.update(data)

        return cls(structure=structure, **symmetry)

    def get_title(self) -> str:
        """Return string with metadata. Useful when creating views"""
        return f"Structure: {self.formula_pretty}, {self.spg_symbol} ({self.spg_number}), " + \
               f"{self.crystal_system}, natom: {self.nsites}"

    #def get_panel_view(self, mng_connector: MongoConnector):
    #    return

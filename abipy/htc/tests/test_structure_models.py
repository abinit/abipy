import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
from abipy.htc.structure_models import StructureData


class TestStructureData(AbipyTest):

    def test_structure_data(self):
        si = Structure.as_structure(abidata.cif_file("si.cif"))
        structure_data = StructureData.from_structure(si)
        self.assertMSONable(structure_data, test_if_subclass=True)
        assert structure_data.crystal_system == "Cubic"
        assert "Cubic" in structure_data.get_title()

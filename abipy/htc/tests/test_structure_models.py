#import mongomock
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
from abipy.htc.structure_models import StructureData


class TestStructureData(AbipyTest):

    def test_structure_data(self):
        si = Structure.as_structure(abidata.cif_file("si.cif"))
        structure_data = StructureData.from_structure(si)
        #collection = mongomock.MongoClient().db.collection
        self.assert_msonable(structure_data, test_if_subclass=True)
        assert structure_data.crystal_system == "Cubic"
        assert "Cubic" in structure_data.get_title()

        #filepath = self.get_tmpname()

        import tempfile
        import os
        with tempfile.TemporaryDirectory() as tmp:
            filepath = os.path.join(tmp, "structure_data.json")
            structure_data.json_write(filepath, indent=4)
            same_data = StructureData.from_json_file(filepath)
            assert type(structure_data) is type(same_data)
            assert same_data.structure == structure_data.structure

        # Test from_mpid.
        mp_structure_data = StructureData.from_mpid("mp-149")
        assert mp_structure_data.structure.formula == "Si2"

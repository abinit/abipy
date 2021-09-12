import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.ebands import ElectronBands
from abipy.htc.gs_models import GsData


class TestGsModel(AbipyTest):

    def test_api(self):
        filepath = abidata.ref_file("si_scf_GSR.nc")
        gs_model = GsData.from_gsr_filepath(filepath)
        assert isinstance(gs_model.ebands, ElectronBands)
        #self.assertMSONable(model)
        d = gs_model.as_dict()
        new = GsData.from_dict(d)
        gs_model.ebands.structure.remove_site_property("cartesian_forces")
        new.ebands.structure.remove_site_property("cartesian_forces")
        assert gs_model.ebands.structure == new.ebands.structure




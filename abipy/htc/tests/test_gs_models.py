import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.ebands import ElectronBands
from abipy.htc.base_models import MockedMongoConnector
from abipy.htc.gs_models import ScfData, NscfData  #, RelaxData


class TestGsModels(AbipyTest):

    def test_gs_data(self):
        """Testing ScfData model"""
        collection_name = "dummy_collection_name"
        mocked_connector = MockedMongoConnector(host="example.com", port=27017, collection_name=collection_name)

        gs_data = ScfData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"), mocked_connector, with_gsr=True)
        assert isinstance(gs_data.ebands, ElectronBands)
        #self.assertMSONable(model)
        new = ScfData.from_dict(gs_data.as_dict())
        gs_data.ebands.structure.remove_site_property("cartesian_forces")
        new.ebands.structure.remove_site_property("cartesian_forces")
        assert gs_data.ebands.structure == new.ebands.structure

        with self.assertRaises(RuntimeError):
            # Should raise as we are asking from a GS model from a NSCF calculation
            ScfData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"), mocked_connector, with_gsr=True)

        # Testing NscfData model
        with self.assertRaises(RuntimeError):
            NscfData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"), mocked_connector, with_gsr=False)

        nscf_data = NscfData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"), mocked_connector, with_gsr=False)
        assert isinstance(nscf_data.ebands, ElectronBands)

    #def test_relax_data(self):
    #    """Testing RelaxData model"""
    #    relax_data = RelaxData.from_hist_gsr_filepaths(abidata.ref_file("sic_relax_HIST.nc"))
    #    print(relax_data)


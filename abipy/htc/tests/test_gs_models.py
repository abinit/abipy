import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.ebands import ElectronBands
from abipy.htc.gs_models import GsData, NscfData  #, RelaxData


class TestGsModelsl(AbipyTest):

    def test_gs_data(self):
        """Testing GsData model"""
        gs_data = GsData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"))
        assert isinstance(gs_data.ebands, ElectronBands)
        #self.assertMSONable(model)
        new = GsData.from_dict(gs_data.as_dict())
        gs_data.ebands.structure.remove_site_property("cartesian_forces")
        new.ebands.structure.remove_site_property("cartesian_forces")
        assert gs_data.ebands.structure == new.ebands.structure

        with self.assertRaises(RuntimeError):
            GsData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"))

    def test_nscf_data(self):
        """Testing NscfData model"""
        with self.assertRaises(RuntimeError):
            NscfData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"))

        nscf_data = NscfData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"))
        assert isinstance(nscf_data.ebands, ElectronBands)

    #def test_relax_data(self):
    #    """Testing RelaxData model"""
    #    relax_data = RelaxData.from_hist_gsr_filepaths(abidata.ref_file("sic_relax_HIST.nc"))
    #    print(relax_data)


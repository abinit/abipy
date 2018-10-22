"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
#from abipy.dfpt.elastic import ElasticData


class ElasticDataFileTest(AbipyTest):

    def test_alas_elastic(self):
        """
        Testing DDB containing also third order derivatives.
        """
        self.skip_if_abinit_not_ge("8.9.3")

        with abilab.abiopen(abidata.ref_file("refs/alas_elastic_dfpt/AlAs_elastic_DDB")) as ddb:
            assert ddb.has_strain_terms(select="at_least_one")
            assert ddb.has_strain_terms(select="all")
            assert ddb.has_internalstrain_terms(select="all")
            assert ddb.has_piezoelectric_terms(select="all")
            assert ddb.has_at_least_one_atomic_perturbation()

            # Get ElasticData by calling anaddb.
            e = ddb.anaget_elastic(verbose=2)
            assert e.params["elaflag"] == 3
            assert e.params["piezoflag"] == 3
            assert e.params["instrflag"] == 1
            assert e.params["asr"] == 2 and e.params["chneut"] == 1
            # Elastic tensors.
            self.assert_almost_equal(e.elastic_relaxed[0,0,0,0], 122.23496623977118)
            assert e.elastic_clamped is not None
            assert e.elastic_stress_corr is None
            assert e.elastic_relaxed_fixed_D is None

            # Piezoelectric tensors.
            self.assert_almost_equal(e.piezo_relaxed[2,2,2], -0.041496005147475756)
            assert e.piezo_clamped is not None
            assert e.d_piezo_relaxed is None
            assert e.g_piezo_relaxed is None
            assert e.g_piezo_relaxed is None

            assert repr(e); assert str(e)
            assert e.to_string(verbose=2)
            assert e.structure.formula == "Al2 As2"
            assert e.elastic_relaxed._repr_html_()

            name_tensor_list = e.name_tensor_list(tensor_type="elastic")
            names = [nt[0] for nt in name_tensor_list]
            assert "elastic_relaxed" in names
            name_tensor_list = e.name_tensor_list(tensor_type="piezoelectric")
            names = [nt[0] for nt in name_tensor_list]
            assert "piezo_relaxed" in names
            edata_fit = e.fit_to_structure()
            assert edata_fit is not None
            edata_ieee = e.convert_to_ieee()
            assert edata_ieee is not None

            self.assertMSONable(e)

            df = e.get_elastic_tensor_dataframe(tensor_name="elastic_clamped", tol=1e-5)
            df = e.get_piezoelectric_tensor_dataframe(tensor_name="piezo_clamped", tol=1e-8)

            df = e.get_voigt_dataframe("elastic_relaxed", voigt_as_index=False, tol=1e-1)
            self.assert_almost_equal(df[(0, 0)][0], 122.23496623977118)
            df = e.get_elastic_properties_dataframe(tensor_names="elastic_relaxed", fit_to_structure=True)

"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dfpt.elastic import ElasticData


class ElaticDataFileTest(AbipyTest):

    def test_base(self):
        """Base tests for ElasticData"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")

        # Construct ElasticData by calling anaddb.
        e = ElasticData.from_file(anaddbnc_fname)
        #assert e.params["elaflag"] == 3
        #assert e.params["piezoflag"] == 3
        #assert e.params["instrflag"] == 1
        #assert e.params["asr"] == 2 and e.params["chneut"] == 1

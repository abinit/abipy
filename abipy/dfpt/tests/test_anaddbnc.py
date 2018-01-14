"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dfpt.anaddbnc import AnaddbNcFile


class AnaddbNcFileTest(AbipyTest):

    def test_base(self):
        """Base tests for AnaddbNcFile"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")

        with AnaddbNcFile(anaddbnc_fname) as anc:
            repr(anc); str(anc)
            anc.to_string(verbose=2)
            assert anc.structure.formula == "Al1 As1"
            assert anc.becs is not None
            assert anc.emacro is not None
            assert anc.emacro_rlx is not None
            assert anc.dchidt is not None
            assert anc.dchide is not None
            assert anc.oscillator_strength is not None
            assert anc.ifc is None
            assert not anc.params

            if self.has_nbformat():
                anc.write_notebook(nbpath=self.get_tmpname(text=True))

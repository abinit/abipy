from __future__ import division, print_function, unicode_literals

#import numpy as np
import abipy.data as abidata
#import abipy.core
from abipy.core.testing import *
from abipy.dynamics.hist import HistFile


class HistFileTest(AbipyTest):

    def test_hist_api(self):
        """Testing HistFile API."""
        hist = HistFile(abidata.ref_file("si_scf_GSR.nc"))
        print(hist)

        assert hist.num_steps == X
        assert len(hist.structures) == hist.num_steps
        assert hist.final_structure is hist.final_structures[-1]
        assert hist.final_structure.volume == XX
        assert len(hist.etotals) == hist.num_steps
        self.assert_equal(hist.etotals, [])

        # Test matplotlib plots.
        if self.has_matplotlib():
            hist.plot(show=False)
            hist.plot_energies(show=False)

        # Test notebook generation.
        if self.has_nbformat():
            hist.write_notebook(nbpath=self.get_tmpname(text=True))

        hist.close()

"""Tests for psps module."""
from __future__ import division, print_function, unicode_literals, absolute_import

#import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
import abipy.display.mvtk as mvtk


class MayaviTest(AbipyTest):

    def test_mayavi_toolkit(self):
        """Test mayavi toolkit."""
        if not self.has_mayavi():
            raise self.SkipTest("This test requires mayavi!")

        figure, mlab = mvtk.get_fig_mlab(figure=None)

        si_structure = self.get_abistructure_from_abiref("si_nscf_GSR.nc")

        same_fig = mvtk.plot_wigner_seitz(si_structure.lattice, figure=figure)
        assert same_fig is figure

        figure = mvtk.plot_unit_cell(si_structure.lattice)
        assert mvtk.plot_lattice_vectors(si_structure.lattice, figure=figure) is figure

        assert mvtk.plot_structure(si_structure, frac_coords=False, to_unit_cell=False, style="points+labels",
                                   unit_cell_color=(0, 0, 0), color_scheme="VESTA", figure=None, show=False)

        #mvtk.plot_labels(labels, lattice=None, coords_are_cartesian=False, figure=None, **kwargs)

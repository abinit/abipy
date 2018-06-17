# coding: utf-8
"""Tests for derivatives module."""
from __future__ import division, print_function, absolute_import, unicode_literals

import os
import numpy as np

from abipy import abilab
import abipy.data as abidata
from abipy.tools.plotting import *
from abipy.core.testing import AbipyTest


class TestPlotting(AbipyTest):
    """Test plotting module."""

    def test_set_axlims(self):
        """Testing set_axlims."""
        if not self.has_matplotlib():
            raise self.SkipTest("This test requires matplotlib")

        ax, fig, plt = get_ax_fig_plt(ax=None)

        left, right = set_axlims(ax, None, "x")
        assert left is None and right is None

        left, right = set_axlims(ax, [1], "x")
        assert left == 1 and right is None

        left, right = set_axlims(ax, 1, "x")
        assert left == 1 and right is None

        left_right = set_axlims(ax, [1, 2], "y")
        assert left_right == (1, 2)

    def test_data_from_cplx_mode(self):
        """Testing data_from_cplx_mode."""
        carr = np.empty((2, 4), dtype=np.complex)

        self.assert_equal(data_from_cplx_mode("re", carr), carr.real)
        self.assert_equal(data_from_cplx_mode("im", carr), carr.imag)
        self.assert_equal(data_from_cplx_mode("abs", carr), np.abs(carr))
        self.assert_equal(data_from_cplx_mode("angle", carr), np.angle(carr))
        with self.assertRaises(ValueError):
            data_from_cplx_mode("foo", carr)

        # Test plot_array
        if self.has_matplotlib():
            assert plot_array(carr, cplx_mode="abs", show=False)
            cvec = np.empty(10, dtype=np.complex)
            assert plot_array(cvec, cplx_mode="im", show=False)

    def test_plot_xy_with_hue(self):
        """Testing plot_xy_with_hue."""
        # Here is some sample data. The 'x2' data is slightly offset from 'x1'
        x1 = list(range(0, 100, 10))
        x2 = list(range(1, 100, 10))
        x = x1 + x2

        #The y-values I generate here mimic the general shape of my actual data
        y1 = x1[::-1]
        y2 = [i+25 for i in x1[::-1]]
        y = y1 + y2

        # Two levels of labels that will be applied to the data
        z1 = ["1"] * 10
        z2 = ["2"] * 10
        z = z1 + z2

        # A pandas data frame from the above data
        import pandas as pd
        df = pd.DataFrame({'x': x, 'y': y, 'z': z})

        from abipy.tools.plotting import plot_xy_with_hue
        if self.has_matplotlib():
            assert plot_xy_with_hue(data=df, x="x", y="y", hue="z", ax=None, show=False)
            assert plot_xy_with_hue(data=df, x="x", y="y", hue="z", ax=None, show=False,
                                    color="red", marker="v")
            assert plot_xy_with_hue(data=df, x="x", y="y", hue="z", decimals=0, ax=None, show=False,
                                    color="red", marker="v")
            with self.assertRaises(ValueError):
                plot_xy_with_hue(data=df, x="foo", y="y", hue="bar", ax=None, show=False)

    def test_array_plotter(self):
        """Testing array plotter."""
        plotter = ArrayPlotter()
        assert len(plotter) == 0
        hello = np.ones((5, 2), dtype=np.complex)
        plotter.add_array("hello", hello)
        assert "hello" in plotter.keys()
        with self.assertRaises(ValueError):
            plotter.add_array("hello", hello)

        plotter.add_arrays(["world", "foo"], [2 * np.ones((2, 2)), 4 * np.ones((4, 1))])
        assert len(plotter) == 3
        for k, v in plotter.items():
            assert k in plotter.keys()

        if self.has_matplotlib():
            plotter.plot(cplx_mode="re", show=False)

    def test_marker(self):
        """Testing Marker."""
        marker = Marker()
        assert not marker
        assert not marker.x
        x, y, s = [1, 2, 3], [4, 5, 6], [0.1, 0.2, -0.3]
        marker.extend((x, y, s))
        assert marker
        self.assert_equal(marker.x, x)
        self.assert_equal(marker.y, y)
        self.assert_equal(marker.s, s)

        pos_mark, neg_mark = marker.posneg_marker()
        assert all(s >= 0 for s in pos_mark.s)
        assert all(s < 0 for s in neg_mark.s)
        assert len(neg_mark.s) == 1 and neg_mark.s[0] == -0.3

        with self.assertRaises(TypeError):
            marker.extend((x, y))

        with self.assertRaises(TypeError):
            x.append(-1)
            marker.extend((x, y, s))

    def test_plot_cell_tools(self):
        """Testing plot_unit_cell."""
        lattice = abilab.Lattice.hexagonal(a=2, c=4)

        if self.has_matplotlib():
            fig, ax = plot_unit_cell(lattice, ax=None, linestyle="-")
            assert hasattr(fig, "show")

    def test_generic_data_file_plotter(self):
        """Testing GenericDataFilePlotter object."""
        filepath = os.path.join(abidata.dirpath, "refs", "sio2_screening", "sio2_EM1_NLF")
        plotter = GenericDataFilePlotter(filepath)
        assert plotter.to_string(verbose=2)
        assert str(plotter)
        assert len(plotter.od) == 1
        key = "# Omega [eV] Re epsilon_M IM eps_M"
        assert key in plotter.od
        assert plotter.od[key].shape == (3, 30)

        if self.has_matplotlib():
            assert plotter.plot(show=False)
            assert plotter.plot(use_index=True, show=False)

    def test_generic_data_files_plotter(self):
        """Testing GenericDataFilesPlotter object."""
        filepaths = [
                os.path.join(abidata.dirpath, "refs", "sio2_screening", "sio2_EM1_NLF"),
                os.path.join(abidata.dirpath, "refs", "sio2_screening", "sio2_EM1_NLF"),
        ]

        plotter = GenericDataFilesPlotter.from_files(filepaths)
        assert plotter.to_string(verbose=2)
        assert str(plotter)

        if self.has_matplotlib():
            assert plotter.plot(show=False)
            assert plotter.plot(use_index=True, show=False)

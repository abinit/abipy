# coding: utf-8
"""Tests for derivatives module."""
from __future__ import division, print_function, absolute_import, unicode_literals

import numpy as np

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

    def test_array_plotter(self):
        """Testing array plotter."""
        plotter = ArrayPlotter()
        assert len(plotter) == 0
        hello = np.ones((5, 2), dtype=np.complex)
        plotter.add_array("hello", hello)
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

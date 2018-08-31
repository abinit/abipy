from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np

from abipy.tools.numtools import *
from abipy.core.testing import AbipyTest


class TestNumTools(AbipyTest):
    """Test numtools."""

    def test_transpose_last3dims(self):
        """Testing transpose_last3dims"""
        arr = np.arange(120)
        arr.shape = (2, 2, 10, 3)

        same_arr = transpose_last3dims(arr)
        same_arr = transpose_last3dims(same_arr)
        assert np.all(arr == same_arr)

    def test_add_periodic_replicas(self):
        """Testing add_periodic_replicas"""
        # 1D nd 2D case
        arr = np.array([1, 2, 3, 4, 5, 6])
        new_arr = add_periodic_replicas(arr)
        assert new_arr[-1] == 1

        arr.shape = (2, 3)
        new_arr = add_periodic_replicas(arr)
        assert np.all(new_arr[-1] == [1, 2, 3, 1])
        assert np.all(new_arr[:,-1] == [1, 4, 1])

        # 4D case.
        arr = np.arange(120)
        arr.shape = (2, 2, 10, 3)

        new_arr = add_periodic_replicas(arr)
        assert np.all(new_arr[:,:-1,:-1,:-1] == arr)

        axes = [[0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2],]

        for ax in axes:
            view = np.transpose(new_arr, axes=ax)
            assert np.all(view[...,0] == view[...,-1])
            assert np.all(view[...,0,0] == view[...,-1,-1])
            assert np.all(view[...,0,0,0] == view[...,-1,-1,-1])

    def test_data_from_cplx_mode(self):
        """Testing data_from_cplx_mode."""
        carr = np.empty((2, 4), dtype=np.complex)

        self.assert_equal(data_from_cplx_mode("all", carr), carr)
        self.assert_equal(data_from_cplx_mode("re", carr), carr.real)
        self.assert_equal(data_from_cplx_mode("im", carr), carr.imag)
        self.assert_equal(data_from_cplx_mode("abs", carr), np.abs(carr))
        self.assert_equal(data_from_cplx_mode("angle", carr), np.angle(carr))

        with self.assertRaises(ValueError):
            data_from_cplx_mode("foo", carr)
        with self.assertRaises(ValueError):
            data_from_cplx_mode("angle", carr, tol=1.0)

        rarr = np.ones((2, 4), dtype=np.float)
        self.assert_equal(data_from_cplx_mode("re", rarr, tol=1.1), np.zeros_like(rarr))

    def test_special_functions(self):
        """Testing special functions."""
        assert gaussian(x=0.0, width=1.0, center=0.0, height=1.0) == 1.0

        assert lorentzian(x=0.0, width=1.0, center=0.0, height=1.0) == 1.0
        self.assert_almost_equal(lorentzian(x=0.0, width=1.0, center=0.0, height=None), 1/np.pi)

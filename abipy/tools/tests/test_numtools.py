import itertools
import numpy as np
import pytest

from pymatgen.core.lattice import Lattice
from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
from abipy.tools.numtools import *


class TestNumTools(AbipyTest):
    """Test numtools."""

    def test_print_stats_arr(self):
        arr = np.array([1, 2, 3, 4, 5])
        print_stats_arr(arr)

    def test_nparr_to_df(self):
        arr = np.array([[1, 2], [3, 4]])
        df = nparr_to_df("values", arr, ["x", "y"])
        assert len(df) == 4
        assert df.columns.tolist() == ["x", "y", "values"]

    def test_build_mesh(self):
        mesh, index = build_mesh(0.0, 5, 1.0, "centered")
        assert len(mesh) == 11
        assert index, 5

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
        carr = np.empty((2, 4), dtype=complex)

        self.assert_equal(data_from_cplx_mode("all", carr), carr)
        self.assert_equal(data_from_cplx_mode("re", carr), carr.real)
        self.assert_equal(data_from_cplx_mode("im", carr), carr.imag)
        self.assert_equal(data_from_cplx_mode("abs", carr), np.abs(carr))
        self.assert_equal(data_from_cplx_mode("angle", carr), np.angle(carr))

        with self.assertRaises(ValueError):
            data_from_cplx_mode("foo", carr)
        with self.assertRaises(ValueError):
            data_from_cplx_mode("angle", carr, tol=1.0)

        rarr = np.ones((2, 4), dtype=float)
        self.assert_equal(data_from_cplx_mode("re", rarr, tol=1.1), np.zeros_like(rarr))

    def test_is_diagonal(self):
        diag_matrix = np.diag([1, 2, 3])
        non_diag_matrix = np.array([[1, 1], [0, 1]])
        assert is_diagonal(diag_matrix)
        assert not is_diagonal(non_diag_matrix)

    def test_special_functions(self):
        """Testing special functions."""
        assert gaussian(x=0.0, width=1.0, center=0.0, height=1.0) == 1.0

        assert lorentzian(x=0.0, width=1.0, center=0.0, height=1.0) == 1.0
        self.assert_almost_equal(lorentzian(x=0.0, width=1.0, center=0.0, height=None), 1/np.pi)

    def test_iflat(self):
        nested_list = [[0], [1, 2, [3, 4]]]
        assert list(iflat(nested_list)) == [0, 1, 2, 3, 4]

    def test_grouper(self):
        assert grouper(3, "ABCDEFG", "x") == [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', 'x', 'x')]

    def test_sort_and_groupby(self):
        keys, groups = sort_and_groupby([1, 2, 1], ret_lists=True)
        assert keys, [1, 2]
        assert groups, [[1, 1], [2]]

    def test_prune_ord(self):
        assert prune_ord([1, 1, 2, 3, 3]) == [1, 2, 3]

    def test_smooth(self):
        x = np.linspace(-2, 2, 50)
        y = np.sin(x)
        smoothed = smooth(y, window_len=11, window='hanning')
        assert len(smoothed) == len(y)

    def test_find_convindex(self):
        values = [1, 1.1, 1.01, 1.001, 1.0001]
        assert find_convindex(values, tol=0.01) == 2

    def test_find_degs_sk(self):
        energies = [1, 1, 2, 3.4, 3.401]
        deg_sets = find_degs_sk(energies, atol=0.01)
        assert deg_sets == [[0, 1], [2], [3, 4]]


class TestBzRegularGridInterpolator(AbipyTest):

    def test_api(self):
        # Creates a simple cubic structure for testing.
        lattice = Lattice.cubic(1.0)  # Simple cubic lattice with a=1
        structure = Structure(lattice, species=["Si"], coords=[[0, 0, 0]])

        # Test that BzRegularGridInterpolator initializes correctly.
        ndat, nx, ny, nz = 2, 4, 5, 6
        shifts = [0, 0, 0]
        shape = (ndat, nx, ny ,nz)     # (ndat=1, nx=4, ny=4, nz=4)
        datak = np.zeros(shape)  # All zeros except one point
        datak[0, 2, 2, 2] = 1.0  # Set one known value
        datak[1, 2, 2, 2] = 2.0  # Set one known value

        # Multiple shifts should raise an error
        with pytest.raises(ValueError, match="Multiple shifts are not supported"):
            BzRegularGridInterpolator(structure, [[0, 0, 0], [0.5, 0.5, 0.5]], datak)

        # Non-zero shift should raise an error
        with pytest.raises(ValueError, match="Shift should be zero"):
            BzRegularGridInterpolator(structure, [0.1, 0.2, 0.3], datak)

        interp = BzRegularGridInterpolator(structure, shifts, datak)
        assert interp.ndat == ndat
        assert interp.dtype == datak.dtype

        # Test interpolation at known fractional coordinates."""
        values = interp.eval_kpoint([0.5, 0.5, 0.5])  # Middle of the grid

        assert isinstance(values, np.ndarray)
        assert values.shape == (ndat,)
        assert 0 <= values[0] <= 1  # Ensure interpolation is reasonable
        assert 0 <= values[1] <= 2  # Ensure interpolation is reasonable

        # Test interpolation with Cartesian coordinates.
        cart_coords = structure.reciprocal_lattice.matrix @ [0.5, 0.5, 0.5]  # Convert to Cartesian

        values = interp.eval_kpoint(cart_coords, cartesian=True)
        assert isinstance(values, np.ndarray)
        assert values.shape == (ndat,)
        assert 0 <= values[0] <= 1
        assert 0 <= values[1] <= 2

        # Test that interpolation handles periodic boundaries correctly."""
        values1 = interp.eval_kpoint([1.0, 1.0, 1.0])
        values2 = interp.eval_kpoint([0.0, 0.0, 0.0])

        np.testing.assert_allclose(values1, values2, atol=1e-6)

        # Compare interpolated and initial reference value.
        for ix, iy, iz in itertools.product(range(nx), range(ny), range(nz)):
            kpoint = [ix/nx, iy/ny, iz/nz]
            values = interp.eval_kpoint(kpoint)
            ref_values = datak[:, ix, iy, iz]
            self.assert_almost_equal(values, ref_values)

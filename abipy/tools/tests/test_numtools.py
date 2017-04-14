from __future__ import print_function, division

from abipy.tools.numtools import *
from abipy.core.testing import *


class TestTools(AbipyTest):
    """Test numtools."""

    def test_transpose_last3dims(self):
        """test transpose_last3dims"""
        arr = np.arange(120)
        arr.shape = (2,2,10,3)

        same_arr = transpose_last3dims(arr)
        same_arr = transpose_last3dims(same_arr)

        self.assertTrue( np.all(arr == same_arr) )

    def test_add_periodic_replicas(self):
        """test add_periodic_replicas"""

        # 1D nd 2D case
        arr = np.array([1,2,3,4,5,6])
        new_arr = add_periodic_replicas(arr)
        self.assertTrue(new_arr[-1] == 1)

        arr.shape = (2,3)
        new_arr = add_periodic_replicas(arr)
        self.assertTrue(np.all(new_arr[-1] == [1,2,3,1]))
        self.assertTrue(np.all(new_arr[:,-1] == [1,4,1]))

        # 4D case.
        arr = np.arange(120)
        arr.shape = (2,2,10,3)

        new_arr = add_periodic_replicas(arr)

        self.assertTrue(np.all(new_arr[:,:-1,:-1,:-1] == arr))

        axes = [[0,1,2,3], [0,2,3,1], [0,3,1,2],]

        for ax in axes:
            view = np.transpose(new_arr, axes=ax)
            self.assertTrue(np.all(view[...,0] == view[...,-1]))
            self.assertTrue(np.all(view[...,0,0] == view[...,-1,-1]))
            self.assertTrue(np.all(view[...,0,0,0] == view[...,-1,-1,-1]))


if __name__ == "__main__":
    import unittest
    unittest.main()

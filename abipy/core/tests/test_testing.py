#!/usr/bin/env python
"""Tests for core.density module"""
from __future__ import print_function, division

import os

from abipy.core.testing import *
from abipy.core.testing import input_equality_check
from abipy.abio.inputs import AbinitInput

root = os.path.dirname(__file__)


class TestTEstingTools(AbipyTest):
    """Unit tests for testing tools."""

    def test_check_input_equality(self):
        """Testing ScalarField."""

        ref_file = os.path.join(root, '..', '..', 'test_files', 'convergence_inputs_single_factory_00.json')

        input1 = AbinitInput.

        input_equality_check(ref_file, input_to_test, rtol=rtol, atol=atol, equal_nan=equal_nan)



if __name__ == "__main__":
    import unittest
    unittest.main()

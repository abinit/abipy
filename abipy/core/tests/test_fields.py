#!/usr/bin/env python
"""Tests for core.field module"""
from __future__ import print_function, division

import numpy as np
import abipy.data as abidata

from abipy.core.fields import *
from abipy.core.testing import AbipyTest


class TestScalarField(AbipyTest):
    """Unit tests for ScalarField."""

    def test_base(self):
        """Testing ScalarField."""
        structure = abidata.structure_from_ucell("Si")
        nspinor, nsppol, nspden = 1, 1, 1
        xyz_shape = (2, 3, 4)
        datar = np.zeros((nsppol,) + xyz_shape)

        field = ScalarField(nspinor, nsppol, nspden, datar, structure, iorder="c")

        print(field)
        assert field.datar.ndim ==  4
        assert len(field) == nsppol
        assert field.shape == datar.shape
        assert field.nx == 2 and field.ny == 3 and field.nz == 4
        assert field.is_collinear

        #assert field.datar_xyz.ndim == 4
        #assert field.datar_xyz.shape[-3:], xyz_shape)

        field.export(self.get_tmpname(text=True, suffix=".xsf"))

"""Tests for eph module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dfpt.eph import EphReader, EliashbergFunction


class EliashbergFunctionTest(AbipyTest):
    """Tests for EliashbergFunction."""


class EphReaderTest(AbipyTest):
    """Tests for EphReader."""

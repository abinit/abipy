# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re
import subprocess
import json
import numpy as np

import logging
logger = logging.getLogger(__file__)


class NumPyArangeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist() # or map(int, obj)

        return json.JSONEncoder.default(self, obj)


def profile(statement, global_vars, local_vars):
    """
    Run statement under profiler, supplying your own globals and locals
    """
    import pstats
    import cProfile
    import tempfile

    _, filename = tempfile.mkstemp()
    cProfile.runctx(statement, global_vars, local_vars, filename=filename)

    s = pstats.Stats(filename)
    s.strip_dirs().sort_stats("time").print_stats()

    os.remove(filename)

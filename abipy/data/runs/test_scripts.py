#!/usr/bin/env python
"""
This script runs all the python scripts located in this directory 
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os 

from subprocess import call, Popen

from abipy.core.testing import *


root = os.path.abspath(os.path.join(os.path.dirname(__file__)))


class TestScripts(AbipyTest):
    def test_all_scripts(self):
        """Testing all scripts in abipy/data/runs/"""
        retcode = call(os.path.join(root, "_runemall.py"))
        assert retcode == 0

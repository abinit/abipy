"""Test suite for abipy.parser.parser."""
import os

import unittest

from abipy.core.testing import *
from abipy.abitools import *
from abipy.abitools.parser import _parser_ok   # if not _parser_ok: Do not run the test suite.

TEST_DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')


class TestParser(AbipyTest):

    @unittest.skipIf(not _parser_ok, "Bindings for the ABINIT parser are not available")
    def test_base(self):
        filename = pjoin(TEST_DATA_PATH, "t66.in")

        dtsets = Datasets(filename)

        dt1 = dtsets[1]
        dt1.export_crystal(".xsf")

        print("dtset(acell_orig) = ", dt1["acell_orig"])

        # Test json export.
        json_name = "test-json.txt"
        dtset.export_json(json_name);

    #def test_json(self):
    #    #import ctypes
    #    #ctypes.PyDLL("./ab6_invars.so", mode=ctypes.RTLD_GLOBAL)
    #    import locale
    #    if(len(sys.argv) == 1):
    #       fname="../tests/tutorial/Input/tbs_1.in"
    #    else:
    #       fname=sys.argv[1]
    #    #print "Filename = ",fname
    #    dtset = abipy.parser.Datasets(fname)
    #    print("dtset(acell_orig) = ",dtset.dtsets[1]["acell_orig"])

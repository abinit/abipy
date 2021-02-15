"""Tests for test_testing module"""
import os
from abipy.core.testing import AbipyTest
from abipy.core.testing import input_equality_check

root = os.path.dirname(__file__)


class TestTEstingTools(AbipyTest):
    """Unit tests for testing tools."""

    def test_check_input_equality(self):
        """Testing the function to test input equality."""

        ref_file = os.path.join(root, '..', '..', 'test_files', 'convergence_inputs_single_factory_00.json')
        bad_file = os.path.join(root, '..', '..', 'test_files', 'convergence_inputs_single_factory_00-bad.json')

        input_good = self.json_read_abinit_input(ref_file)
        input_bad = self.json_read_abinit_input(bad_file)

        input_equality_check(ref_file, input_good)

        try:
            input_equality_check(ref_file, input_bad)
        except Exception as ex:
            self.assertIsInstance(ex, AssertionError)
            error_header = "Two inputs were found to be not equal:\n" \
                            "   not the same input parameters:\n" \
                            "     ['kptopt'] were found in ref but not in actual\n" \
                            "     [] were found in actual but not in ref\n\n"

            s = str(ex).replace("[u'", "['")
            assert s.startswith(error_header)
            for l in [
                "variable kptopt from the reference is not in the actual input\n\n",
                "var shiftk differs: [[0.0, 0.0, 0.0]] (reference) != [[0.1, 0.0, 0.0]] (actual)\n",
                "var ngkpt differs: [10, 10, 10] (reference) != [11, 10, 10] (actual)\n",
                "var nshiftk differs: 1 (reference) != 5 (actual)\n",
                "var charge differs: 0.0 (reference) != 0.01 (actual)\n",
                ]:
                assert l in s

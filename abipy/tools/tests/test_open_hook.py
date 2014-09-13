"""Unit tests for open_hook"""
import unittest
import tempfile

py3k = False
try:
    import __builtin__
except ImportError:
    py3k = True
    import builtins as __builtin__


@unittest.skipIf(py3k, "OpenHook not compatible with py3k")
class OpenHookTest(unittest.TestCase):
    """Test open_hook module."""
    def setUp(self):
        from abipy.tools import open_hook
        self.open_hook = open_hook
        self.open_hook.install()

    def test_open(self):
        """Test if open_hook detects open files."""
        self.assertEqual(self.open_hook.num_openfiles(), 0)

        _, fname = tempfile.mkstemp()
        fh = open(fname)

        self.open_hook.print_open_files()
        self.assertEqual(self.open_hook.num_openfiles(), 1)

    def tearDown(self):
        self.open_hook.remove()
        # Here we should have the __builtin__ version

        self.assertTrue(file is __builtin__.file)
        self.assertTrue(open is __builtin__.open)


if __name__ == "__main__":
    unittest.main()

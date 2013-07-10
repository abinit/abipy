"""Tests for htc.FilesFile."""
import warnings

from abipy.htc import *
from abipy.tests import *

# =========================================================================== #

class TestFilesFile(AbipyFileTest):
    """Unit tests for FilesFile."""

    def setUp(self):
        self.file = FilesFile('MyDir/MyName.files',
                              input='mycalc.in',
                              output='mycalc.out',
                              idat_root='i_mycalc',
                              odat_root='o_mycalc',
                              tmp_root='t_mycalc')

        self.file.pseudos = ['ps1', 'ps2']
        self.file.pseudodir = '/path/to/my/pseudodir'

    def test_check_pseudos(self):
        """Test the user is warned of pseudopotential not found."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.file.check_pseudos()
            self.assertEqual(len(w), 2)
            msg = str(w[-1].message)
            self.assertRegexpMatches(msg, "file not found")


    def test_str(self):
        """Test the FilesFile is printed correctly."""
        lookfor = """
mycalc.in
mycalc.out
i_mycalc
o_mycalc
t_mycalc
/path/to/my/pseudodir/ps1
/path/to/my/pseudodir/ps2
        """.strip()
        self.assertContains(lookfor)



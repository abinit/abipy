# coding: utf-8
"""Tests for devtools module."""
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.tools.devtools import profile, HtmlDiff


class DevtoolsTest(AbipyTest):

    def test_profile(self):
        """Testing profile function."""
        def statement():
            return 1

        stats = profile("statement()", global_vars=globals(), local_vars=locals())
        assert stats is not None

    def test_htmldiff(self):
        """Testing HtmlDiff."""
        filepaths = abidata.cif_files("si.cif", "al.cif", "gan.cif")
        diff = HtmlDiff(filepaths)

        # patch _launch_browser
        def _launch_browser(*args, **kwargs): return "patched"
        diff._launch_browser = _launch_browser
        assert diff.open_browser() == "patched"
        # This requires pygmentize package.
        retcode = diff.open_browser(diffmode="pygmentize")

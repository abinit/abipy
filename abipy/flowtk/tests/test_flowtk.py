"""Tests for flowtk __init__ module"""
from __future__ import print_function, division, unicode_literals, absolute_import

#import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest


class TestFlowtk(AbipyTest):
    """Unit tests for flowtk.__init__ module."""

    def test_flow_main(self):
        parser = flowtk.build_flow_main_parser()
        assert parser is not None

        @flowtk.flow_main
        def main(options):
            return flowtk.Flow.temporary_flow()

        mock = self.get_mock_module()
        with mock.patch('sys.argv', ["test_main.py", "prof", "--help"]):
            with self.assertRaises(SystemExit) as cm:
                main()
        assert cm.exception.code == 0

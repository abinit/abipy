"""Tests for flowtk __init__ module"""

from abipy import flowtk
from abipy.core.testing import AbipyTest


class TestFlowtk(AbipyTest):
    """Unit tests for flowtk.__init__ module."""

    def test_flow_main(self):
        """Testing flow_main decorator."""
        parser = flowtk.build_flow_main_parser()
        assert parser is not None
        # parser.parse_args("--help")

        @flowtk.flow_main
        def main(options):
            return flowtk.Flow.temporary_flow()

        mock = self.get_mock_module()
        with mock.patch("sys.argv", ["test_main.py", "prof", "--help"]), self.assertRaises(SystemExit) as cm:
            main()
        assert cm.exception.code == 0
